!> \file    mo_river_temperature.f90
!> \copydoc mo_river_temperature

!> \brief   River temperature routing.
!> \details This module contains mRM-owned thermal state and an open-water heat
!! processor for the shared river DAG traversal. The routed quantity is physical
!! thermal power relative to 0 K, in watts.
!> \version 0.1
!> \changelog
!! - 2026: v6 open-water scalar heat-balance solver.
!! - 2026: mRM-owned thermal state, processor-owned forcing caches, and scaling workspaces.
!> \authors Sebastian Mueller
!> \date    October 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_temperature
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_kind, only: dp, i4, i8
  use mo_constants, only: T0_dp, sigma_dp, Psychro_dp, cp_w_dp, nodata_dp
  use mo_pet, only: slope_satpressure
  use mo_river, only: river_t
  use mo_river_tools, only: river_muskingum_state_t, river_depth_to_flux
  use mo_river_flow, only: river_flow_state_t
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t, up_a_mean, down_nearest, down_scaling
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  use mo_message, only: error_message
  use mo_string_utils, only: n2s => num2str
  use mo_utils, only: optval
  implicit none
  private
  ! cp_w_dp in FORCES is numerically 4.19 kJ kg-1 K-1 despite its J-based comment.
  ! 1000 kg m-3 * 4.19 kJ kg-1 K-1 * 1000 J kJ-1.
  real(dp), parameter, public :: water_heat_capacity_volume = 1000.0_dp * cp_w_dp * 1000.0_dp ! J m-3 K-1
  real(dp), parameter :: min_temperature = 200.0_dp, max_temperature = 400.0_dp
  public :: river_heat_flux, solve_river_temperature, river_thermal_state_t, river_thermal_t

  !> \class river_thermal_state_t
  !> \brief mRM-owned open-water histories and published results; energy rates are W.
  type :: river_thermal_state_t
    logical :: active = .false. !< Owning container's thermal activation policy, checked during router initialization.
    real(dp) :: albedo !< Open-water shortwave albedo.
    real(dp) :: pt !< Priestley-Taylor coefficient.
    real(dp) :: emissivity !< Open-water emissivity.
    real(dp) :: kh !< Sensible heat exchange coefficient.
    real(dp) :: parameters(5) !< Width reference and four heat parameters.
    real(dp), allocatable :: area(:) !< Physical reach surface area [m2].
    real(dp), allocatable :: previous_input(:) !< Previous link energy input [W].
    real(dp), allocatable :: previous_routed(:) !< Previous routed power [W].
    real(dp), allocatable :: current_input(:) !< Current link energy input [W].
    real(dp), allocatable :: current_routed(:) !< Current routed power [W].
    real(dp), allocatable :: energy_flux(:) !< Completed mean thermal power [W].
    real(dp), allocatable :: temperature(:) !< Completed river temperature [degC].
  contains
    procedure :: allocate => river_thermal_state_allocate
    procedure :: reset => river_thermal_state_reset
    procedure :: destroy => river_thermal_state_destroy
  end type river_thermal_state_t

  !> \class river_thermal_forcing_t
  !> \brief Processor-owned endpoint map and meteorology coefficients for one routing interval.
  type, private :: river_thermal_forcing_t
    integer(i8), allocatable :: cell_in(:) !< Upstream endpoint cell on the meteo grid.
    integer(i8), allocatable :: cell_out(:) !< Downstream endpoint cell on the meteo grid.
    real(dp), allocatable :: ta(:) !< Mean endpoint air temperature [K].
    real(dp), allocatable :: r0(:) !< Incoming radiation after albedo [W m-2].
    real(dp), allocatable :: eta(:) !< Priestley-Taylor multiplier.
    real(dp), allocatable :: tswitch(:) !< Radiation branch boundary [K].
    real(dp) :: albedo !< Parameter used to prepare incoming radiation.
    real(dp) :: pt !< Parameter used to prepare latent heat exchange.
    real(dp) :: emissivity !< Parameter used to prepare the radiation boundary.
  contains
    procedure :: init => river_thermal_forcing_init
    procedure :: update => river_thermal_forcing_update
    procedure :: destroy => river_thermal_forcing_destroy
  end type river_thermal_forcing_t

  !> \class river_thermal_t
  !> \brief Open-water processor with owned forcing and energy scratch.
  type :: river_thermal_t
    type(river_t), pointer :: river => null() !< Borrowed network; the mRM-owned target outlives this processor.
    type(river_thermal_forcing_t) :: forcing !< Cached per-reach meteorology.
    real(dp), allocatable :: meteo_work(:, :) !< Reusable unpacked L1 forcing workspace.
    real(dp), allocatable :: scale_grid(:, :) !< Reusable unpacked L1 runoff workspace.
    type(scaler_t) :: meteo_scaler !< L1-to-L3 meteorology remapping cache.
    real(dp), allocatable :: acc_depth(:) !< Accumulated runoff temperature-depth [K mm].
    real(dp), allocatable :: acc_air(:) !< Accumulated air temperature [degC].
    real(dp), allocatable :: acc_ssrd(:) !< Accumulated shortwave radiation [W m-2].
    real(dp), allocatable :: acc_strd(:) !< Accumulated longwave radiation [W m-2].
    real(dp), allocatable :: air(:) !< Mapped air temperature [degC].
    real(dp), allocatable :: ssrd(:) !< Mapped shortwave radiation [W m-2].
    real(dp), allocatable :: strd(:) !< Mapped longwave radiation [W m-2].
    real(dp), allocatable :: local_energy(:) !< Local lateral thermal power [W].
    real(dp), allocatable :: scaled_flux(:) !< L3 temperature-weighted runoff rate [K m3 s-1].
    real(dp), allocatable :: scale_work(:) !< L1 temperature-weighted volume work buffer.
    real(dp), allocatable :: sum_energy(:) !< Sum of substep input power [W].
    real(dp), allocatable :: sum_discharge(:) !< Sum of substep discharge [m3 s-1].
  contains
    procedure :: init => river_thermal_init
    procedure :: to_restart_dataset => river_thermal_to_restart_dataset
    procedure :: destroy => river_thermal_destroy
    procedure :: accumulate => river_thermal_accumulate
    procedure :: prepare => river_thermal_prepare
    procedure, private :: map_meteo => river_thermal_map_meteo
    procedure :: advance => river_thermal_advance
    procedure :: complete => river_thermal_complete
    procedure :: route_nodes => river_thermal_route_nodes
  end type river_thermal_t
contains
  !> \brief Initialize the endpoint map and parameters used to prepare forcing.
  subroutine river_thermal_forcing_init(self, river, albedo, pt, emissivity)
    class(river_thermal_forcing_t), intent(inout), target :: self !< Forcing cache to initialize.
    type(river_t), intent(in) :: river !< Routing network providing endpoint cells.
    real(dp), intent(in) :: albedo !< Open-water albedo for this initialization.
    real(dp), intent(in) :: pt !< Priestley-Taylor parameter for this initialization.
    real(dp), intent(in) :: emissivity !< Water emissivity for this initialization.
    integer(i8) :: n, downstream
    call self%destroy()
    allocate(self%cell_in(river%n_nodes), self%cell_out(river%n_nodes), self%ta(river%n_nodes), &
      self%r0(river%n_nodes), self%eta(river%n_nodes), self%tswitch(river%n_nodes))
    self%albedo = albedo
    self%pt = pt
    self%emissivity = emissivity
    do n = 1_i8, river%n_nodes
      downstream = n
      if (.not.river%is_sink(n)) downstream = river%down(n)
      self%cell_in(n) = n
      self%cell_out(n) = downstream
      if (allocated(river%node_cell)) then
        self%cell_in(n) = river%node_cell(n)
        self%cell_out(n) = river%node_cell(downstream)
      end if
    end do
  end subroutine river_thermal_forcing_init

  !> \brief Prepare reach meteorology once, before all numerical routing substeps.
  subroutine river_thermal_forcing_update(self, air, shortwave, longwave)
    class(river_thermal_forcing_t), intent(inout), target :: self !< Per-reach cache to refresh.
    real(dp), intent(in) :: air(:) !< Mapped air temperature [degC].
    real(dp), intent(in) :: shortwave(:) !< Mapped shortwave [W m-2].
    real(dp), intent(in) :: longwave(:) !< Mapped longwave [W m-2].
    integer(i8) :: n, cell_in, cell_out
    real(dp) :: sw, lw, delta
    !$omp parallel do default(shared) private(cell_in,cell_out,sw,lw,delta) schedule(static)
    do n = 1_i8, size(self%ta, kind=i8)
      cell_in = self%cell_in(n)
      cell_out = self%cell_out(n)
      self%ta(n) = 0.5_dp*(air(cell_in)+air(cell_out)) + T0_dp
      sw = 0.5_dp*(shortwave(cell_in)+shortwave(cell_out))
      lw = 0.5_dp*(longwave(cell_in)+longwave(cell_out))
      delta = slope_satpressure(self%ta(n)-T0_dp)
      self%eta(n) = self%pt*delta / (Psychro_dp+delta)
      self%r0(n) = (1.0_dp-self%albedo)*sw + lw
      self%tswitch(n) = min_temperature
      if (self%r0(n) > 0.0_dp) self%tswitch(n) = sqrt(sqrt(self%r0(n)/(self%emissivity*sigma_dp)))
    end do
    !$omp end parallel do
  end subroutine river_thermal_forcing_update

  !> \brief Release the processor-owned forcing cache.
  subroutine river_thermal_forcing_destroy(self)
    class(river_thermal_forcing_t), intent(inout), target :: self !< Forcing cache to release.
    if (allocated(self%cell_in)) deallocate(self%cell_in)
    if (allocated(self%cell_out)) deallocate(self%cell_out)
    if (allocated(self%ta)) deallocate(self%ta)
    if (allocated(self%r0)) deallocate(self%r0)
    if (allocated(self%eta)) deallocate(self%eta)
    if (allocated(self%tswitch)) deallocate(self%tswitch)
  end subroutine river_thermal_forcing_destroy

  !> \brief Remap forcing using packed downscaling or a reusable unpacked workspace.
  subroutine river_thermal_map_meteo(self, input, output)
    class(river_thermal_t), intent(inout), target :: self !< Processor owning meteo support and workspace.
    real(dp), intent(in) :: input(:) !< Packed L1 meteorology.
    real(dp), intent(out) :: output(:) !< Packed L3 meteorology.
    if (self%meteo_scaler%scaling_mode == down_scaling) then
      call self%meteo_scaler%downscale_nearest(input, output)
    else
      call self%meteo_scaler%source_grid%unpack_into(input, self%meteo_work)
      call self%meteo_scaler%execute(self%meteo_work, output)
    end if
  end subroutine river_thermal_map_meteo

  !> \brief Allocate mRM-owned histories and stable published arrays.
  subroutine river_thermal_state_allocate(self, n_nodes)
    class(river_thermal_state_t), intent(inout), target :: self !< Thermal state to initialize.
    integer(i8), intent(in) :: n_nodes !< Number of routing nodes.
    call self%destroy()
    allocate(self%area(n_nodes), self%previous_input(n_nodes), self%previous_routed(n_nodes), &
      self%current_input(n_nodes), self%current_routed(n_nodes), self%energy_flux(n_nodes), &
      self%temperature(n_nodes))
    call self%reset()
  end subroutine river_thermal_state_allocate

  !> \brief Reset histories and completed outputs without replacing exchange targets.
  subroutine river_thermal_state_reset(self)
    class(river_thermal_state_t), intent(inout), target :: self !< Thermal state to reset in place.
    self%active = .false.
    if (allocated(self%area)) self%area = 0.0_dp
    if (allocated(self%previous_input)) self%previous_input = 0.0_dp
    if (allocated(self%previous_routed)) self%previous_routed = 0.0_dp
    if (allocated(self%current_input)) self%current_input = 0.0_dp
    if (allocated(self%current_routed)) self%current_routed = 0.0_dp
    if (allocated(self%energy_flux)) self%energy_flux = 0.0_dp
    if (allocated(self%temperature)) self%temperature = nodata_dp
  end subroutine river_thermal_state_reset

  !> \brief Release mRM-owned thermal histories and published arrays.
  subroutine river_thermal_state_destroy(self)
    class(river_thermal_state_t), intent(inout), target :: self !< Thermal state to release.
    if (allocated(self%area)) deallocate(self%area)
    if (allocated(self%previous_input)) deallocate(self%previous_input)
    if (allocated(self%previous_routed)) deallocate(self%previous_routed)
    if (allocated(self%current_input)) deallocate(self%current_input)
    if (allocated(self%current_routed)) deallocate(self%current_routed)
    if (allocated(self%energy_flux)) deallocate(self%energy_flux)
    if (allocated(self%temperature)) deallocate(self%temperature)
    self%active = .false.
  end subroutine river_thermal_state_destroy

  !> \brief Initialize scratch from prepared state parameters and area, optionally restoring restart histories.
  subroutine river_thermal_init(self, state, river, input_grid, runoff_cells, needs_scale_work, nc)
    class(river_thermal_t), intent(inout), target :: self !< Thermal processor.
    type(river_thermal_state_t), intent(inout), target :: state !< mRM-owned thermal state.
    type(river_t), intent(in), target :: river !< Routing network.
    type(grid_t), pointer, intent(in) :: input_grid !< L1 forcing grid.
    integer(i8), intent(in) :: runoff_cells !< Cells on L3 runoff support.
    logical, intent(in) :: needs_scale_work !< Whether runoff upscaling needs an L1 buffer.
    type(NcDataset), intent(in), optional :: nc !< Open restart dataset replacing fresh state histories.
    type(NcVariable) :: var
    real(dp) :: saved_parameters(5)
    integer(i8) :: n

    call self%destroy()
    self%river => river
    n = river%n_nodes
    if (.not.allocated(state%area)) call error_message("river_temperature: mRM-owned thermal state is not allocated")
    if (size(state%area, kind=i8) /= n) &
      call error_message("river_temperature: mRM-owned thermal state is not allocated")
    allocate(self%acc_depth(input_grid%ncells), self%acc_air(input_grid%ncells), &
      self%acc_ssrd(input_grid%ncells), self%acc_strd(input_grid%ncells))
    allocate(self%air(river%grid%ncells), self%ssrd(river%grid%ncells), self%strd(river%grid%ncells))
    allocate(self%local_energy(n), self%sum_energy(n), self%sum_discharge(n))
    if (river%scc) allocate(self%scaled_flux(runoff_cells))
    if (needs_scale_work) allocate(self%scale_work(input_grid%ncells), self%scale_grid(input_grid%nx,input_grid%ny))
    call self%meteo_scaler%init(input_grid, river%grid, &
      upscaling_operator=up_a_mean, downscaling_operator=down_nearest)
    if (self%meteo_scaler%scaling_mode /= down_scaling) allocate(self%meteo_work(input_grid%nx,input_grid%ny))
    self%acc_depth = 0.0_dp
    self%acc_air = 0.0_dp
    self%acc_ssrd = 0.0_dp
    self%acc_strd = 0.0_dp
    self%local_energy = 0.0_dp
    self%sum_energy = 0.0_dp
    self%sum_discharge = 0.0_dp
    state%albedo = state%parameters(2)
    state%pt = state%parameters(3)
    state%emissivity = state%parameters(4)
    state%kh = state%parameters(5)
    if (present(nc)) then
      if (.not.nc%hasVariable("thermal_parameters") .or. .not.nc%hasVariable("thermal_surface_area") .or. &
          .not.nc%hasVariable("thermal_previous_input") .or. &
          .not.nc%hasVariable("thermal_previous_routed") .or. .not.nc%hasVariable("thermal_energy_flux") .or. &
          .not.nc%hasVariable("thermal_temperature")) &
        call error_message("river_router: temperature restart state is incomplete")
      var = nc%getVariable("thermal_parameters"); call var%readInto(saved_parameters)
      if (any(abs(saved_parameters - state%parameters) > 1.0e-12_dp)) &
        call error_message("river_router: temperature parameters differ from restart")
      var = nc%getVariable("thermal_surface_area"); call var%readInto(state%area)
      var = nc%getVariable("thermal_previous_input"); call var%readInto(state%previous_input)
      var = nc%getVariable("thermal_previous_routed"); call var%readInto(state%previous_routed)
      var = nc%getVariable("thermal_energy_flux"); call var%readInto(state%energy_flux)
      var = nc%getVariable("thermal_temperature"); call var%readInto(state%temperature)
    end if
    if (any(state%area < 0.0_dp)) call error_message("river_router: negative river surface area")
    call self%forcing%init(river, state%albedo, state%pt, state%emissivity)
  end subroutine river_thermal_init

  !> \brief Persist thermal state using the existing restart variable names.
  subroutine river_thermal_to_restart_dataset(self, state, nc, node_dim)
    class(river_thermal_t), intent(in), target :: self !< Thermal processor.
    type(river_thermal_state_t), intent(in), target :: state !< mRM-owned thermal state.
    type(NcDataset), intent(inout) :: nc !< Open restart dataset.
    type(NcDimension), intent(in) :: node_dim !< River node dimension.
    type(NcDimension) :: parameter_dim
    type(NcVariable) :: nc_var
    parameter_dim = nc%setDimension("thermal_parameter", 5_i4)
    nc_var = nc%setVariable("thermal_parameters", "f64", [parameter_dim])
    call nc_var%setData(state%parameters)
    nc_var = nc%setVariable("thermal_surface_area", "f64", [node_dim])
    call nc_var%setAttribute("units", "m2")
    call nc_var%setData(state%area)
    nc_var = nc%setVariable("thermal_previous_input", "f64", [node_dim])
    call nc_var%setAttribute("units", "W")
    call nc_var%setData(state%previous_input)
    nc_var = nc%setVariable("thermal_previous_routed", "f64", [node_dim])
    call nc_var%setAttribute("units", "W")
    call nc_var%setData(state%previous_routed)
    nc_var = nc%setVariable("thermal_energy_flux", "f64", [node_dim])
    call nc_var%setAttribute("units", "W")
    call nc_var%setData(state%energy_flux)
    nc_var = nc%setVariable("thermal_temperature", "f64", [node_dim])
    call nc_var%setAttribute("units", "degC")
    call nc_var%setData(state%temperature)
  end subroutine river_thermal_to_restart_dataset

  !> \brief Release processor-owned forcing and energy scratch.
  subroutine river_thermal_destroy(self)
    class(river_thermal_t), intent(inout), target :: self !< Thermal processor to release.
    if (allocated(self%acc_depth)) deallocate(self%acc_depth)
    if (allocated(self%acc_air)) deallocate(self%acc_air)
    if (allocated(self%acc_ssrd)) deallocate(self%acc_ssrd)
    if (allocated(self%acc_strd)) deallocate(self%acc_strd)
    if (allocated(self%air)) deallocate(self%air)
    if (allocated(self%ssrd)) deallocate(self%ssrd)
    if (allocated(self%strd)) deallocate(self%strd)
    if (allocated(self%local_energy)) deallocate(self%local_energy)
    if (allocated(self%scaled_flux)) deallocate(self%scaled_flux)
    if (allocated(self%scale_work)) deallocate(self%scale_work)
    if (allocated(self%sum_energy)) deallocate(self%sum_energy)
    if (allocated(self%sum_discharge)) deallocate(self%sum_discharge)
    if (allocated(self%meteo_work)) deallocate(self%meteo_work)
    if (allocated(self%scale_grid)) deallocate(self%scale_grid)
    call self%forcing%destroy()
    call self%meteo_scaler%destroy()
    nullify(self%river)
  end subroutine river_thermal_destroy

  !> \brief Accumulate runoff energy and meteorology until the next routing result.
  subroutine river_thermal_accumulate(self, temperature_depth, air_temperature, shortwave, longwave)
    class(river_thermal_t), intent(inout), target :: self !< Thermal processor.
    real(dp), optional, intent(in) :: temperature_depth(:) !< L1 runoff temperature-depth [K mm].
    real(dp), optional, intent(in) :: air_temperature(:) !< L1 air temperature [degC].
    real(dp), optional, intent(in) :: shortwave(:) !< L1 incoming shortwave [W m-2].
    real(dp), optional, intent(in) :: longwave(:) !< L1 incoming longwave [W m-2].
    if (.not.present(temperature_depth) .or. .not.present(air_temperature) .or. &
        .not.present(shortwave) .or. .not.present(longwave)) &
      call error_message("river_router: active temperature routing requires runoff energy and all meteorology")
    if (size(temperature_depth) /= size(self%acc_depth) .or. size(air_temperature) /= size(self%acc_depth) .or. &
        size(shortwave) /= size(self%acc_depth) .or. size(longwave) /= size(self%acc_depth)) &
      call error_message("river_router: temperature input size differs from runoff input grid")
    block
      integer(i8) :: c
      !$omp simd
      do c = 1_i8, size(self%acc_depth, kind=i8)
        self%acc_depth(c) = self%acc_depth(c) + temperature_depth(c)
        self%acc_air(c) = self%acc_air(c) + air_temperature(c)
        self%acc_ssrd(c) = self%acc_ssrd(c) + shortwave(c)
        self%acc_strd(c) = self%acc_strd(c) + longwave(c)
      end do
    end block
  end subroutine river_thermal_accumulate

  !> \brief Map runoff energy and mean meteorology, then reset completed-interval accumulators.
  subroutine river_thermal_prepare(self, scaler, muskingum, input_step, node_runoff_cell)
    class(river_thermal_t), intent(inout), target :: self !< Thermal processor.
    type(scaler_t), intent(inout) :: scaler !< Router-owned L1-to-L3 depth-to-flux map.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Shared coefficients and completed routing cadence.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    integer(i8), intent(in), optional :: node_runoff_cell(:) !< Router-owned mapping from nodes to land runoff cells.
    integer(i8) :: n, c

    if (muskingum%accumulations > 1_i4) &
      self%acc_depth = self%acc_depth / muskingum%accumulations
    if (self%river%scc) then
      ! Q and Q*T use the same map with processor-owned output buffers.
      call river_depth_to_flux(scaler, self%acc_depth, self%scaled_flux, input_step, self%scale_work, self%scale_grid)
      !$omp parallel do default(shared) private(n,c) schedule(static)
      do n = 1_i8, self%river%n_nodes
        c = self%river%node_cell(n)
        self%local_energy(n) = 0.0_dp
        if (present(node_runoff_cell)) c = node_runoff_cell(n)
        if (c > 0_i8) self%local_energy(n) = &
          water_heat_capacity_volume * self%scaled_flux(c) * self%river%area_fraction(n)
      end do
      !$omp end parallel do
    else
      call river_depth_to_flux(scaler, self%acc_depth, self%local_energy, input_step, self%scale_work, self%scale_grid)
      self%local_energy = water_heat_capacity_volume * self%local_energy
    end if
    if (muskingum%accumulations > 1_i4) then
      !$omp simd
      do c = 1_i8, size(self%acc_air, kind=i8)
        self%acc_air(c) = self%acc_air(c) / muskingum%accumulations
        self%acc_ssrd(c) = self%acc_ssrd(c) / muskingum%accumulations
        self%acc_strd(c) = self%acc_strd(c) / muskingum%accumulations
      end do
    end if
    call self%map_meteo(self%acc_air, self%air)
    call self%map_meteo(self%acc_ssrd, self%ssrd)
    call self%map_meteo(self%acc_strd, self%strd)
    call self%forcing%update(self%air, self%ssrd, self%strd)
    !$omp simd
    do c = 1_i8, size(self%acc_air, kind=i8)
      self%acc_depth(c) = 0.0_dp
      self%acc_air(c) = 0.0_dp
      self%acc_ssrd(c) = 0.0_dp
      self%acc_strd(c) = 0.0_dp
    end do
    !$omp simd
    do c = 1_i8, size(self%sum_energy, kind=i8)
      self%sum_energy(c) = 0.0_dp
      self%sum_discharge(c) = 0.0_dp
    end do
  end subroutine river_thermal_prepare

  !> \brief Swap thermal history buffers after a routing substep without copying them.
  subroutine river_thermal_advance(self, state)
    class(river_thermal_t), intent(in), target :: self !< Thermal processor.
    type(river_thermal_state_t), intent(inout), target :: state !< Thermal histories to advance.
    real(dp), allocatable :: tmp(:)
    call move_alloc(state%previous_input, tmp)
    call move_alloc(state%current_input, state%previous_input)
    call move_alloc(tmp, state%current_input)
    call move_alloc(state%previous_routed, tmp)
    call move_alloc(state%current_routed, state%previous_routed)
    call move_alloc(tmp, state%current_routed)
  end subroutine river_thermal_advance

  !> \brief Publish mean thermal power and its discharge-weighted temperature.
  subroutine river_thermal_complete(self, state, muskingum)
    class(river_thermal_t), intent(in), target :: self !< Thermal processor.
    type(river_thermal_state_t), intent(inout), target :: state !< Thermal output state.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Number of substeps in the completed interval.
    integer(i8) :: n
    !$omp simd
    do n = 1_i8, size(self%sum_energy, kind=i8)
      state%energy_flux(n) = self%sum_energy(n) / muskingum%iterations
      if (self%sum_discharge(n) > 0.0_dp) then
        state%temperature(n) = self%sum_energy(n) / (water_heat_capacity_volume*self%sum_discharge(n)) - T0_dp
      else
        state%temperature(n) = nodata_dp
      end if
    end do
  end subroutine river_thermal_complete

  !> \brief Route a complete level or consecutive serial levels in topological order.
  !> \details Parallel callers share a team and bounds; the loop barrier completes the range before returning.
  recursive subroutine river_thermal_route_nodes(self, first, last, flow, state, muskingum, parallel)
    class(river_thermal_t), intent(inout), target :: self !< Thermal processor with mapped forcing.
    integer(i8), intent(in) :: first !< First topological level to route.
    integer(i8), intent(in) :: last !< Last topological level to route.
    type(river_flow_state_t), intent(in), target :: flow !< Discharge already updated for this range.
    type(river_thermal_state_t), intent(inout), target :: state !< Thermal histories to update.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Read-only routing coefficients and cadence.
    logical, optional, intent(in) :: parallel !< Collective static worksharing; default is a serial topological traversal.
    integer(i8) :: j, first_node, last_node

    first_node = self%river%order%level_start(first)
    last_node = self%river%order%level_end(last)

    if (optval(parallel,.false.)) then
      ! Complete energy throughout this range before downstream levels read it.
      !$omp do schedule(static)
      do j = first_node, last_node
        call river_thermal_process_node(j)
      end do
      !$omp end do
    else
      do j = first_node, last_node
        call river_thermal_process_node(j)
      end do
    end if
  contains
    !> \brief Apply the node balance using the enclosing range's shared state and coefficient arrays.
    subroutine river_thermal_process_node(position)
      integer(i8), intent(in) :: position !< Position in the node order.
      integer(i8) :: n
      real(dp) :: qin, qout, previous_qin, previous_qout
      integer(i8) :: m
      integer(i8), pointer, contiguous :: edges(:)
      real(dp) :: adv, tup, h0, surface, solved, nu1, nu2
      integer(i4) :: status

      n = self%river%order%id(position)
      qin = flow%current_input(n)
      qout = flow%current_routed(n)
      previous_qin = flow%previous_input(n)
      previous_qout = flow%previous_routed(n)
      adv = self%local_energy(n)
      call self%river%src_view(n, edges)
      do m = 1_i8, size(edges, kind=i8)
        adv = adv + state%current_routed(edges(m))
      end do
      self%sum_energy(n) = self%sum_energy(n) + adv
      self%sum_discharge(n) = self%sum_discharge(n) + qin
      if (abs(qin) <= 1.0e-12_dp) then
        if (abs(previous_qin) > 1.0e-12_dp .or. abs(previous_qout) > 1.0e-12_dp .or. &
            abs(qout) > 1.0e-12_dp .or. abs(adv) > 1.0e-6_dp) &
          call error_message("river_router: zero-flow wet reach has no thermal storage model, node ", n2s(n))
        state%current_input(n) = 0.0_dp
        state%current_routed(n) = 0.0_dp
        return
      end if
      if (qin < 0.0_dp) call error_message("river_router: negative node discharge for thermal route")
      tup = adv / (water_heat_capacity_volume * qin)
      if (self%river%is_sink(n)) then
        state%current_input(n) = adv
        state%current_routed(n) = 0.0_dp
        return
      end if
      if (qout <= 0.0_dp) call error_message("river_router: wet reach has nonpositive routed discharge, node ", n2s(n))
      nu1 = muskingum%nu1(n)
      nu2 = muskingum%nu2(n)
      h0 = (1.0_dp - nu1) * state%previous_routed(n) + &
        (nu1 - nu2) * state%previous_input(n) + nu2 * adv
      call solve_heat_balance(qout, h0, nu2 * state%area(n) * 0.5_dp, &
        state%area(n), tup, self%forcing%ta(n), self%forcing%r0(n), self%forcing%eta(n), &
        self%forcing%tswitch(n), state%emissivity, state%kh, &
        solved, surface, status)
      if (status /= 0_i4) call error_message("river_router: temperature root failed at node ", n2s(n), &
        ", status ", n2s(status))
      state%current_input(n) = adv + surface
      state%current_routed(n) = water_heat_capacity_volume * qout * solved
    end subroutine river_thermal_process_node
  end subroutine river_thermal_route_nodes

  !> \brief Evaluate atmospheric heat flux into open water in W m-2.
  real(dp) function river_heat_flux(t, ta, shortwave, longwave, albedo, pt, emissivity, kh) result(flux)
    real(dp), intent(in) :: t !< Water temperature [K].
    real(dp), intent(in) :: ta !< Air temperature [K].
    real(dp), intent(in) :: shortwave !< Incoming shortwave [W m-2].
    real(dp), intent(in) :: longwave !< Incoming longwave [W m-2].
    real(dp), intent(in) :: albedo !< Open-water shortwave albedo.
    real(dp), intent(in) :: pt !< Priestley-Taylor coefficient.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp) :: r0, delta, eta
    delta = slope_satpressure(ta - T0_dp)
    eta = pt * delta / (Psychro_dp + delta)
    r0 = (1.0_dp - albedo) * shortwave + longwave
    flux = heat_flux_precomputed(t, ta, r0, eta, emissivity, kh)
  end function river_heat_flux

  !> \brief Evaluate heat flux with meteorology-dependent terms computed once for a reach.
  pure real(dp) function heat_flux_precomputed(t, ta, r0, eta, emissivity, kh) result(flux)
    real(dp), intent(in) :: t !< Water temperature [K].
    real(dp), intent(in) :: ta !< Air temperature [K].
    real(dp), intent(in) :: r0 !< Incoming radiation after albedo [W m-2].
    real(dp), intent(in) :: eta !< Priestley-Taylor multiplier.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp) :: netrad
    netrad = r0 - emissivity * sigma_dp * t**4
    flux = netrad - eta * max(netrad, 0.0_dp) - kh * (t - ta)
  end function heat_flux_precomputed

  !> \brief Solve the branch-aware heat balance for downstream water temperature.
  !> \details Q is the link outflow in m3 s-1; H0 is the known routed power in W;
  !! D is nu2 times half the surface area in m2. Status 0 is success, 1 no root,
  !! and 2 multiple admissible roots. Scalar root updates are safeguarded by a bracket.
  subroutine solve_river_temperature(q, h0, d, area, tup, ta, shortwave, longwave, albedo, pt, emissivity, kh, &
                                     temperature, surface_source, status, used_fallback)
    real(dp), intent(in) :: q !< Routed link discharge [m3 s-1].
    real(dp), intent(in) :: h0 !< Known Muskingum thermal term [W].
    real(dp), intent(in) :: d !< nu2 times half the physical reach area [m2].
    real(dp), intent(in) :: area !< Physical reach surface area [m2].
    real(dp), intent(in) :: tup !< Incoming water temperature [K].
    real(dp), intent(in) :: ta !< Mean endpoint air temperature [K].
    real(dp), intent(in) :: shortwave !< Mean endpoint shortwave [W m-2].
    real(dp), intent(in) :: longwave !< Mean endpoint longwave [W m-2].
    real(dp), intent(in) :: albedo !< Open-water shortwave albedo.
    real(dp), intent(in) :: pt !< Priestley-Taylor coefficient.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp), intent(out) :: temperature !< Solved outgoing water temperature [K].
    real(dp), intent(out) :: surface_source !< Atmospheric source added to link input [W].
    integer(i4), intent(out) :: status !< Zero for one valid root.
    logical, optional, intent(out) :: used_fallback !< Whether bracketed solving was used.
    real(dp) :: r0, eta, delta, tswitch
    if (q <= 0.0_dp .or. area < 0.0_dp) then
      temperature = 0.0_dp
      surface_source = 0.0_dp
      status = 1_i4
      if (present(used_fallback)) used_fallback = .false.
      return
    end if
    delta = slope_satpressure(ta - T0_dp)
    eta = pt * delta / (Psychro_dp + delta)
    r0 = (1.0_dp - albedo) * shortwave + longwave
    tswitch = min_temperature
    if (r0 > 0.0_dp) tswitch = sqrt(sqrt(r0 / (emissivity * sigma_dp)))
    call solve_heat_balance(q, h0, d, area, tup, ta, r0, eta, tswitch, emissivity, kh, &
      temperature, surface_source, status, used_fallback)
  end subroutine solve_river_temperature

  !> \brief Solve the physical balance with forcing coefficients prepared by the caller.
  subroutine solve_heat_balance(q, h0, d, area, tup, ta, r0, eta, tswitch, emissivity, kh, &
                                temperature, surface_source, status, used_fallback)
    real(dp), intent(in) :: q !< Routed discharge [m3 s-1].
    real(dp), intent(in) :: h0 !< Known routed power [W].
    real(dp), intent(in) :: d !< Muskingum coefficient times half the reach area [m2].
    real(dp), intent(in) :: area !< Physical reach surface area [m2].
    real(dp), intent(in) :: tup !< Incoming water temperature [K].
    real(dp), intent(in) :: ta !< Mean endpoint air temperature [K].
    real(dp), intent(in) :: r0 !< Incoming radiation after albedo [W m-2].
    real(dp), intent(in) :: eta !< Priestley-Taylor multiplier.
    real(dp), intent(in) :: tswitch !< Radiation branch boundary [K].
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat coefficient [W m-2 K-1].
    real(dp), intent(out) :: temperature !< Solved water temperature [K].
    real(dp), intent(out) :: surface_source !< Atmospheric source [W].
    integer(i4), intent(out) :: status !< Zero for one valid root, one for none, two for multiple.
    logical, optional, intent(out) :: used_fallback !< Whether any bisection was needed.
    real(dp) :: lo, hi, a, a4, a1, a0, b, fup, cvq
    real(dp) :: cuts(3), gx, gl, gu, mid, root, tol, tol_upper, stationary
    real(dp) :: accepted(4)
    integer(i4) :: branch, nc, nroot, j, k
    logical :: valid, bisected

    ! Default to failure and zero outputs so every early return has defined results.
    temperature = 0.0_dp
    surface_source = 0.0_dp
    status = 1_i4
    if (present(used_fallback)) used_fallback = .false.
    if (q <= 0.0_dp .or. area < 0.0_dp) return

    ! These terms are shared by both radiation branches.  The known upstream
    ! surface flux contributes to the constant part of the implicit heat balance.
    cvq = water_heat_capacity_volume * q
    fup = heat_flux_precomputed(tup, ta, r0, eta, emissivity, kh)
    b = h0 + d * fup
    accepted = 0.0_dp
    nroot = 0_i4

    ! The max(netrad,0) latent-heat term gives two quartics: one while
    ! net radiation is positive, and one after it reaches zero.
    do branch = 1_i4, 2_i4
      if (branch == 1_i4) then
        lo = min_temperature
        hi = min(max_temperature, tswitch)
        a = 1.0_dp - eta
      else
        lo = max(min_temperature, tswitch)
        hi = max_temperature
        a = 1.0_dp
      end if
      ! The radiative coefficient is 1-eta while net radiation drives
      ! evaporation; it is 1 when negative radiation switches latent loss off.
      if (hi < lo) cycle

      ! Expand g(T) = Cv*Q*T - H0 - D*(F(Tup)+F(T)) into
      ! a4*T**4 + a1*T + a0 on this branch; only the radiation
      ! coefficient a differs between the sunlit and dark branches.
      a4 = d * a * emissivity * sigma_dp
      a1 = cvq + d * kh
      a0 = -(b + d * (a * r0 + kh * ta))
      ! A quartic can turn once on the branch. Split at its stationary
      ! temperature, when present, so each interval is monotonic and can
      ! be checked for a root by its endpoint residuals.
      cuts = [lo, hi, hi]
      nc = 2_i4
      ! A stationary point solves g'(T) = 4*a4*T**3 + a1 = 0.
      if (abs(a4) > tiny(a4)) then
        stationary = -a1 / (4.0_dp * a4)
        if (stationary > 0.0_dp) then
          stationary = stationary**(1.0_dp / 3.0_dp)
          if (stationary > lo .and. stationary < hi) then
            cuts = [lo, stationary, hi]
            nc = 3_i4
          end if
        end if
      end if
      ! Check each monotonic interval independently. Record roots at cuts
      ! (including a tangent root at a stationary point) before testing for
      ! a sign change in the open interval.
      do j = 1_i4, nc - 1_i4
        call evaluate_balance(cuts(j), gl, tol)
        call evaluate_balance(cuts(j+1), gu, tol_upper)
        if (abs(gl) <= tol) call add_root(cuts(j))
        if (abs(gu) <= tol_upper) call add_root(cuts(j+1))
        if ((gl >= 0.0_dp .and. gu >= 0.0_dp) .or. (gl <= 0.0_dp .and. gu <= 0.0_dp)) cycle
        valid = .false.
        ! Avoid iteration when the quartic term is negligible at every
        ! permitted temperature; solve the remaining linear balance directly.
        if (abs(a4 * max_temperature**4) <= 1.0e-12_dp * &
            max(abs(a1 * max_temperature), abs(a0), 1.0_dp)) then
          ! The quartic contribution is below the residual tolerance: solve the linear limit directly.
          if (abs(a1) > tiny(a1)) then
            root = -a0 / a1
            valid = ieee_is_finite(root)
          end if
        else
          call sparse_quartic_root(a4, a1, a0, tup, cuts(j), cuts(j+1), root, valid, bisected)
          if (present(used_fallback)) used_fallback = used_fallback .or. bisected
        end if
        ! An iterative candidate is provisional: require it to stay on the
        ! interval that produced the bracket and satisfy the original W balance.
        if (valid) valid = root >= cuts(j) .and. root <= cuts(j+1)
        if (valid) then
          call evaluate_balance(root, gx, tol)
          valid = abs(gx) <= tol
        end if
        if (.not.valid) then
          ! Keep the sign-changing interval if the fast solve fails validation.
          ! Bisection stops early when both residual and width tolerances pass;
          ! the final candidate is checked against the physical residual below.
          if (present(used_fallback)) used_fallback = .true.
          lo = cuts(j)
          hi = cuts(j+1)
          do k = 1_i4, 100_i4
            mid = 0.5_dp * (lo + hi)
            call evaluate_balance(mid, gx, tol)
            if (abs(gx) <= tol .and. hi - lo <= 1.0e-6_dp) exit
            if ((gl <= 0.0_dp .and. gx >= 0.0_dp) .or. (gl >= 0.0_dp .and. gx <= 0.0_dp)) then
              hi = mid
            else
              lo = mid
              gl = gx
            end if
          end do
          root = mid
        end if
        call evaluate_balance(root, gx, tol)
        if (abs(gx) <= tol) call add_root(root)
      end do
    end do
    ! Publish only an unambiguous physical root. The atmospheric source is
    ! the trapezoidal mean of the endpoint heat fluxes times full reach area.
    if (nroot == 1_i4) then
      status = 0_i4
      temperature = accepted(1)
      surface_source = area * 0.5_dp * (fup + &
        heat_flux_precomputed(temperature, ta, r0, eta, emissivity, kh))
    else if (nroot > 1_i4) then
      status = 2_i4
    end if
  contains
    !> \brief Evaluate the original heat balance and its relative energy tolerance together.
    subroutine evaluate_balance(t, value, tolerance)
      real(dp), intent(in) :: t !< Trial temperature [K].
      real(dp), intent(out) :: value !< Physical heat-balance residual [W].
      real(dp), intent(out) :: tolerance !< Allowed energy residual [W].
      real(dp) :: power, source
      ! This residual is in watts and matches the routed physical balance;
      ! scale its tolerance to the magnitudes of the competing energy terms.
      power = cvq * t
      source = d * (fup + heat_flux_precomputed(t, ta, r0, eta, emissivity, kh))
      value = power - h0 - source
      tolerance = 1.0e-10_dp * max(abs(power), abs(h0), abs(source), 1.0_dp)
    end subroutine evaluate_balance

    !> \brief Add a unique, physically valid temperature root to the accepted candidates.
    subroutine add_root(t)
      real(dp), intent(in) :: t !< Trial water temperature [K].
      integer(i4) :: m
      if (.not.ieee_is_finite(t)) return
      if (t < min_temperature .or. t > max_temperature) return
      ! Adjacent intervals share cut temperatures; do not count the same root twice.
      do m = 1_i4, nroot
        if (abs(t - accepted(m)) < 1.0e-6_dp) return
      end do
      if (nroot < size(accepted)) then
        nroot = nroot + 1_i4
        accepted(nroot) = t
      end if
    end subroutine add_root
  end subroutine solve_heat_balance

  !> \brief Find a root of a4*T**4+a1*T+a0 on one monotonic branch without allocating arrays.
  subroutine sparse_quartic_root(a4, a1, a0, seed, lower, upper, root, valid, bisected)
    real(dp), intent(in) :: a4 !< Quartic coefficient.
    real(dp), intent(in) :: a1 !< Linear coefficient.
    real(dp), intent(in) :: a0 !< Constant coefficient.
    real(dp), intent(in) :: seed !< Incoming water-temperature seed [K].
    real(dp), intent(in) :: lower !< Lower branch endpoint [K].
    real(dp), intent(in) :: upper !< Upper branch endpoint [K].
    real(dp), intent(out) :: root !< Candidate root [K].
    logical, intent(out) :: valid !< Whether the normalized polynomial converged.
    logical, intent(out) :: bisected !< Whether a safeguarded midpoint was needed.
    real(dp) :: c4, c1, c0, scale, x, x2, x3, f, df, step, candidate, lo, hi, fl, checkpoint
    integer(i4) :: iteration

    ! Solve in x = T/300 K so quartic and linear terms have comparable scales.
    c4 = a4 * 300.0_dp**4
    c1 = a1 * 300.0_dp
    scale = max(abs(c4), abs(c1), abs(a0))
    ! Keep the seed as the failure result; the caller can reject it via valid.
    root = seed
    valid = .false.
    bisected = .false.
    if (.not.ieee_is_finite(scale) .or. scale <= 0.0_dp) return
    c4 = c4 / scale
    c1 = c1 / scale
    c0 = a0 / scale
    lo = lower / 300.0_dp
    hi = upper / 300.0_dp
    ! Clamp the upstream-temperature seed to this monotonic branch's bracket.
    x = min(max(seed / 300.0_dp, lo), hi)
    ! The branch endpoints bracket at most one root; fl tracks the lower endpoint.
    fl = (c4 * lo**3 + c1) * lo + c0
    checkpoint = hi - lo
    do iteration = 1_i4, 50_i4
      ! Reuse powers for both polynomial value and its Newton derivative.
      x2 = x * x
      x3 = x2 * x
      f = (c4 * x3 + c1) * x + c0
      df = 4.0_dp * c4 * x3 + c1
      if (abs(f) < 1.0e-12_dp) then
        root = x * 300.0_dp
        valid = .true.
        return
      end if
      if ((fl <= 0.0_dp .and. f <= 0.0_dp) .or. (fl >= 0.0_dp .and. f >= 0.0_dp)) then
        lo = x
        fl = f
      else
        hi = x
      end if
      ! Newton is fast near a root; the bracket midpoint remains the safe default.
      candidate = 0.5_dp * (lo + hi)
      if (abs(df) > epsilon(df) * max(abs(4.0_dp*c4*x3), abs(c1), tiny(df))) then
        step = f / df
        candidate = x - step
      else
        bisected = .true.
      end if
      ! Periodically enforce bracket progress even if all proposed root updates stay inside it.
      if (.not.ieee_is_finite(candidate) .or. candidate <= lo .or. candidate >= hi .or. candidate == x) then
        candidate = 0.5_dp * (lo + hi)
        bisected = .true.
      else if (mod(iteration, 4_i4) == 0_i4) then
        if (hi - lo > 0.5_dp * checkpoint) then
          candidate = 0.5_dp * (lo + hi)
          bisected = .true.
        end if
        checkpoint = hi - lo
      end if
      x = candidate
    end do
    ! Return the last estimate for diagnostics, but leave valid false on exhaustion.
    root = x * 300.0_dp
  end subroutine sparse_quartic_root
end module mo_river_temperature
