!> \file    mo_river_router.f90
!> \copydoc mo_river_router

!> \brief   River router.
!> \details This module contains a router for river networks.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_router

  use mo_kind, only: i4, i8, dp
  use mo_constants, only: nodata_dp
  use mo_utils, only: optval
  use mo_string_utils, only: n2s => num2str
  use mo_river, only: river_t
  use mo_river_tools, only: river_muskingum_state_t, routing_steps, routing_space_weight
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t, up_sum, down_nearest, down_scaling
  use mo_message, only: error_message, warn_message, message
  use mo_river_temperature, only: river_thermal_t, river_thermal_state_t
  use mo_river_flow, only: river_flow_t, river_flow_state_t
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension

  implicit none
  private

  public :: routing_steps, routing_space_weight
  ! time slice selector in muskingum schema for link_inflow and link_routed
  integer(i4), public, parameter :: current = 1_i4 !< selector for current time step
  integer(i4), public, parameter :: previous = 2_i4 !< selector for previous time step
  real(dp), public, parameter :: liter_to_m3 = 0.001_dp !< conversion factor from liter to m3

  !> \class river_router_t
  !> \brief Coordinate discharge and thermal worksharing over each DAG level.
  type, public :: river_router_t
    type(river_t), pointer :: river => null() !< Borrowed routing network.
    type(grid_t), pointer :: input_grid => null() !< Borrowed L1 runoff support.
    type(grid_t), pointer :: runoff_grid => null() !< Borrowed L3 runoff support.
    type(scaler_t) :: scaler !< Shared L1-to-L3 depth-to-flux map.
    integer(i4) :: input_step = 0_i4 !< Represented runoff support [h] shared by all processors.
    integer(i4) :: input_count = 0_i4 !< Model updates accumulated for the next completed routing interval.
    integer(i8), allocatable :: node_runoff_cell(:) !< Shared mapping from river nodes to land runoff cells.
    type(river_muskingum_state_t) :: muskingum !< Router-owned coefficients and cadence.
    type(river_flow_t) :: flow !< Discharge routing and owned runoff scratch.
    type(river_thermal_t) :: thermal !< Thermal routing and owned forcing scratch.
    integer(i8) :: omp_level_thresh = 0_i8 !< Minimum level size for parallel routing.
    logical, allocatable :: parallel_level(:) !< Eligible network levels for the configured thread count.
    integer(i8), allocatable :: serial_level_end(:) !< Last level in each consecutive serial run.
  contains
    procedure, public :: init => river_router_init
    procedure, public :: update => river_router_update
    procedure, public :: destroy => river_router_destroy
    procedure, private :: setup_grids => river_router_setup_grids
    procedure, private :: setup_parallelization => river_router_setup_parallelization
    procedure, private :: route => river_router_route
    ! restart
    procedure, public :: to_restart_dataset => river_router_to_restart_dataset, to_restart_file => river_router_to_restart_file
    procedure, public :: from_restart_dataset => river_router_from_restart_dataset, from_restart_file => river_router_from_restart_file
    generic, public :: to_restart => to_restart_dataset, to_restart_file
    generic, public :: from_restart => from_restart_dataset, from_restart_file
  end type river_router_t

contains

  !> \brief Initialize the shared runoff remapping support.
  subroutine river_router_setup_grids(this, river, input_grid, input_step, runoff_grid)
    class(river_router_t), intent(inout), target :: this !< Router to configure.
    type(river_t), pointer, intent(in) :: river !< Routing network.
    type(grid_t), pointer, intent(in) :: input_grid !< L1 runoff support.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    type(grid_t), pointer, optional, intent(in) :: runoff_grid !< Optional land-only L3 runoff support.
    integer(i8), allocatable :: full_to_land(:)
    integer(i8) :: n

    this%input_step = input_step
    this%input_count = 0_i4
    this%river => river
    this%input_grid => input_grid
    if (present(runoff_grid)) then
      this%runoff_grid => runoff_grid
    else
      this%runoff_grid => river%grid
    end if
    call this%scaler%init(this%input_grid, this%runoff_grid, &
      upscaling_operator=up_sum, downscaling_operator=down_nearest)
    if (allocated(this%node_runoff_cell)) deallocate(this%node_runoff_cell)
    if (.not.associated(this%river%grid, this%runoff_grid)) then
      allocate(full_to_land(this%river%grid%ncells), this%node_runoff_cell(this%river%n_nodes))
      call this%river%grid%gen_id_map(this%runoff_grid, full_to_land)
      !$omp parallel do default(shared) private(n) schedule(static)
      do n = 1_i8, this%river%n_nodes
        this%node_runoff_cell(n) = full_to_land(this%river%node_cell(n))
      end do
      !$omp end parallel do
    end if
  end subroutine river_router_setup_grids

  !> \brief Release router-owned support and detach the borrowed network.
  subroutine river_router_destroy(this)
    class(river_router_t), intent(inout), target :: this !< Router to destroy.
    call this%flow%destroy()
    call this%thermal%destroy()
    call this%muskingum%destroy()
    call this%scaler%destroy()
    nullify(this%river, this%input_grid, this%runoff_grid)
    if (allocated(this%parallel_level)) deallocate(this%parallel_level)
    if (allocated(this%serial_level_end)) deallocate(this%serial_level_end)
    if (allocated(this%node_runoff_cell)) deallocate(this%node_runoff_cell)
    this%input_step = 0_i4
    this%input_count = 0_i4
    this%omp_level_thresh = 0_i8
  end subroutine river_router_destroy

  !> \brief Initialize cadence and processors using domain-owned discharge and thermal state.
  subroutine river_router_init(this, river, flow_state, thermal_state, celerity, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, model_step, &
    runoff_grid, lake_nodes)
    implicit none
    class(river_router_t), intent(inout), target :: this !< Router coordinating routing work.
    type(river_t), pointer, intent(in) :: river !< river definition
    type(river_flow_state_t), intent(inout), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(inout), target :: thermal_state !< mRM-owned thermal state and activation policy.
    real(dp), intent(in) :: celerity(:) !< celerity of the link starting at each river node
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: -1 - threads*8, specials: 0 - all serial, 1 - all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh !< Minimum river-level size for parallel routing.
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)
    type(grid_t), pointer, optional, intent(in) :: runoff_grid !< optional land-only L3 runoff grid
    integer(i8), optional, intent(in) :: lake_nodes(:) !< canonical lake nodes in point-set order

    integer(i4) :: model_step_, input_step_

    input_step_ = optval(input_step, 1_i4)
    model_step_ = optval(model_step, input_step_)
    call this%setup_grids(river, input_grid, input_step_, runoff_grid)
    call this%flow%init(flow_state, this%river, this%input_grid%ncells, &
      this%runoff_grid%ncells, this%river%scc, this%scaler%scaling_mode /= down_scaling, &
      lake_nodes, input_grid=this%input_grid)

    ! setup muskingum parameters
    if (.not.allocated(this%river%link_length)) call error_message("river_router: link_length not available")
    call this%muskingum%init(this%river%link_length, this%river%is_sink, celerity, &
      this%input_step, model_step_, max_route_step)
    ! setup parallelization
    call this%setup_parallelization(root_levels, omp_level_thresh)
    ! The owning container supplies activation, parameters, and fresh surface area.
    if (thermal_state%active) then
      if (this%flow%lake_routing) call error_message("river_router: temperature routing needs a lake outflow water temperature")
      call this%thermal%init(thermal_state, this%river, this%input_grid, this%runoff_grid%ncells, &
        this%scaler%scaling_mode /= down_scaling)
    else
      call this%thermal%destroy()
    end if
  end subroutine river_router_init

  !> \brief Restore flow and any active thermal processor from a restart file.
  subroutine river_router_from_restart_file(this, path, river, flow_state, thermal_state, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, read_fluxes, &
    model_step, runoff_grid, lake_nodes)
    use mo_netcdf, only : NcDataset
    implicit none
    class(river_router_t), intent(inout), target :: this !< Router coordinating routing work.
    character(*), intent(in) :: path !< NetCDF file path
    type(river_t), pointer, intent(in) :: river !< river definition
    type(river_flow_state_t), intent(inout), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(inout), target :: thermal_state !< mRM-owned thermal state and activation policy.
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: threads * 8, 0 to run all in serial, 1 to run all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh !< Minimum river-level size for parallel routing.
    logical, optional, intent(in) :: read_fluxes !< whether to read fluxes from restart file (default: .true.)
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)
    type(grid_t), pointer, optional, intent(in) :: runoff_grid !< optional land-only L3 runoff grid
    integer(i8), optional, intent(in) :: lake_nodes(:) !< canonical lake nodes in point-set order
    type(NcDataset) :: nc
    nc = NcDataset(path, "r")
    call this%from_restart_dataset( &
      nc=nc, river=river, flow_state=flow_state, thermal_state=thermal_state, input_grid=input_grid, input_step=input_step, max_route_step=max_route_step, &
      root_levels=root_levels, omp_level_thresh=omp_level_thresh, read_fluxes=read_fluxes, model_step=model_step, &
      runoff_grid=runoff_grid, lake_nodes=lake_nodes)
    call nc%close()
  end subroutine river_router_from_restart_file

  !> \brief Restore cadence and initialize processors using an open restart dataset.
  subroutine river_router_from_restart_dataset(this, nc, river, flow_state, thermal_state, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, read_fluxes, &
    model_step, runoff_grid, lake_nodes)
    implicit none
    class(river_router_t), intent(inout), target :: this !< Router coordinating routing work.
    type(NcDataset), intent(in) :: nc !< restart dataset
    type(river_t), pointer, intent(in) :: river !< river definition
    type(river_flow_state_t), intent(inout), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(inout), target :: thermal_state !< mRM-owned thermal state and activation policy.
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: threads * 8, 0 to run all in serial, 1 to run all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh !< Minimum river-level size for parallel routing.
    logical, optional, intent(in) :: read_fluxes !< whether to read fluxes from restart file (default: .true.)
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)
    type(grid_t), pointer, optional, intent(in) :: runoff_grid !< optional land-only L3 runoff grid
    integer(i8), optional, intent(in) :: lake_nodes(:) !< canonical lake nodes in point-set order

    type(NcVariable) :: var
    integer(i4) :: model_step_, input_step_

    input_step_ = optval(input_step, 1_i4)
    model_step_ = optval(model_step, input_step_)
    call this%setup_grids(river, input_grid, input_step_, runoff_grid)
    call this%flow%init(flow_state, this%river, this%input_grid%ncells, &
      this%runoff_grid%ncells, this%river%scc, this%scaler%scaling_mode /= down_scaling, &
      lake_nodes, input_grid=this%input_grid, nc=nc, read_fluxes=read_fluxes)
    call this%muskingum%allocate(this%river%n_nodes)

    ! get muskingum parameter nu1
    if (nc%hasVariable("nu1")) then
      var = nc%getVariable("nu1")
      call var%readInto(this%muskingum%nu1)
    else
      call warn_message("nu1 not found in restart file, cannot initialize river router from restart without nu1.")
      error stop 1
    end if

    ! get muskingum parameter nu2
    if (nc%hasVariable("nu2")) then
      var = nc%getVariable("nu2")
      call var%readInto(this%muskingum%nu2)
    else
      call warn_message("nu2 not found in restart file, cannot initialize river router from restart without nu2.")
      error stop 1
    end if

    ! get routing step
    if (nc%hasVariable("route_step")) then
      var = nc%getVariable("route_step")
      call var%readInto(this%muskingum%routing_substep)
    else
      call warn_message("route_step not found in restart file, cannot initialize river router from restart without route_step.")
      error stop 1
    end if

    if (present(max_route_step)) then
      if (this%muskingum%routing_substep > max_route_step) then
        call warn_message("Routing substep in restart file (", n2s(this%muskingum%routing_substep), &
          "s) is larger than specified max_route_step (", n2s(max_route_step), "s).")
        error stop 1
      end if
    end if

    call this%muskingum%setup_timing(this%input_step, model_step_)

    ! setup parallelization
    call this%setup_parallelization(root_levels, omp_level_thresh)

    ! The owning container supplies activation, parameters, and fresh surface area.
    if (thermal_state%active) then
      if (this%flow%lake_routing) call error_message("river_router: temperature routing needs a lake outflow water temperature")
      call this%thermal%init(thermal_state, this%river, this%input_grid, this%runoff_grid%ncells, &
        this%scaler%scaling_mode /= down_scaling, nc, read_fluxes)
    else
      call this%thermal%destroy()
    end if
  end subroutine river_router_from_restart_dataset

  !> \brief Write river router state to restart file
  subroutine river_router_to_restart_file(this, path, flow_state, thermal_state, append)
    use mo_netcdf, only: NcDataset
    implicit none
    class(river_router_t), intent(in), target :: this !< Initialized router.
    character(*), intent(in) :: path !< NetCDF file path
    type(river_flow_state_t), intent(in), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(in), target :: thermal_state !< mRM-owned thermal histories.
    logical, optional, intent(in) :: append !< whether NetCDF file should be opened in append mode (default .false.)
    type(NcDataset) :: nc
    character(1) :: fmode
    fmode = "w"
    if ( optval(append, .false.) ) fmode = "a"
    nc = NcDataset(path, fmode)
    call this%to_restart_dataset(nc, flow_state, thermal_state)
    call nc%close()
  end subroutine river_router_to_restart_file

  !> \brief Write river router state to restart file
  subroutine river_router_to_restart_dataset(this, nc, flow_state, thermal_state)
    use mo_netcdf, only: NcDataset, NcVariable
    implicit none
    class(river_router_t), intent(in), target :: this !< Initialized router.
    type(NcDataset), intent(inout) :: nc !< restart dataset
    type(river_flow_state_t), intent(in), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(in), target :: thermal_state !< mRM-owned thermal histories.
    type(NcVariable) :: nc_var
    type(NcDimension) :: node_dim, dims(0)

    if (this%input_count /= 0_i4) call error_message( &
      "Cannot write routing restart inside an unfinished routing step; choose a routing boundary.")

    if ( .not.nc%hasDimension("node") ) then
      if (this%river%n_nodes > int(huge(1_i4),i8)) call error_message("river_router%to_restart: river network too large to write to restart.")
      node_dim = nc%setDimension("node", int(this%river%n_nodes, i4))  ! only works if network is not to huge for i4
    else
      node_dim = nc%getDimension("node")
    end if

    call this%flow%to_restart_dataset(flow_state, nc, node_dim)

    ! write muskingum parameter nu1
    if ( allocated(this%muskingum%nu1) ) then
      call message("writing nu1 to restart_file")
      nc_var = nc%setVariable("nu1", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "muskingum parameter nu1")
      call nc_var%setAttribute("units", "1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%muskingum%nu1)
    end if

    ! write muskingum parameter nu2
    if ( allocated(this%muskingum%nu2) ) then
      call message("writing nu2 to restart_file")
      nc_var = nc%setVariable("nu2", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "muskingum parameter nu2")
      call nc_var%setAttribute("units", "1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%muskingum%nu2)
    end if

    ! numerical routing substep (used to determine nu1 and nu2)
    nc_var = nc%setVariable("route_step", "f64", dims(:0)) ! scalar
    call nc_var%setAttribute("long_name", "numerical routing substep")
    call nc_var%setAttribute("units", "s")
    call nc_var%setData(this%muskingum%routing_substep)
    if (thermal_state%active) call this%thermal%to_restart_dataset(thermal_state, nc, node_dim)

  end subroutine river_router_to_restart_dataset

  !> \brief Setup parallelization of river routing based on river level sizes and user settings.
  subroutine river_router_setup_parallelization(this, root_levels, omp_level_thresh)
    !$ use omp_lib, only: omp_get_num_threads
    implicit none
    class(river_router_t), intent(inout), target :: this !< Router coordinating routing work.
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: -1 - threads*8, specials: 0 - all serial, 1 - all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh !< Minimum river-level size for parallel routing.
    integer(i8) :: i, last, n_levels
    integer(i4) :: threads

    if (.not.allocated(this%river%order%id)) call this%river%calc_order(root=root_levels, sort=.true.)
    threads = 1_i4
    !$omp parallel
    !$omp single
    !$ threads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    this%omp_level_thresh = optval(omp_level_thresh, -1_i8)
    if (this%omp_level_thresh < 0_i8) this%omp_level_thresh = int(threads, i8)*8_i8
    if (allocated(this%parallel_level)) deallocate(this%parallel_level)
    if (allocated(this%serial_level_end)) deallocate(this%serial_level_end)
    n_levels = this%river%order%n_levels
    allocate(this%parallel_level(n_levels), source=.false.)
    allocate(this%serial_level_end(n_levels))
    if (threads > 1_i4 .and. this%omp_level_thresh > 0_i8) then
      do i = 1_i8, n_levels
        this%parallel_level(i) = this%river%order%level_size(i) >= this%omp_level_thresh
      end do
    end if

    ! One single section can process several small levels without an intervening barrier.
    i = 1_i8
    do while (i <= n_levels)
      last = i
      if (.not.this%parallel_level(i)) then
        do while (last < n_levels)
          if (this%parallel_level(last + 1_i8)) exit
          last = last + 1_i8
        end do
      end if
      this%serial_level_end(i:last) = last
      i = last + 1_i8
    end do
  end subroutine river_router_setup_parallelization

  !> \brief Update routing for one time step
  subroutine river_router_update(this, input_runoff, flow_state, thermal_state, completed, lake_outflow, &
    temperature_depth, air_temperature, shortwave, longwave)
    implicit none
    class(river_router_t), intent(inout), target :: this !< Router coordinating routing work.
    real(dp), dimension(this%input_grid%ncells), intent(in) :: input_runoff !< L1 runoff depth [mm].
    type(river_flow_state_t), intent(inout), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(inout), target :: thermal_state !< mRM-owned thermal histories.
    logical, intent(out) :: completed !< whether a new routing-step mean was produced
    real(dp), optional, intent(in) :: lake_outflow(:) !< prescribed lake outflow in point-set order [m3 s-1]
    real(dp), optional, intent(in) :: temperature_depth(:) !< temperature-weighted runoff [K mm]
    real(dp), optional, intent(in) :: air_temperature(:) !< L1 air temperature [degC].
    real(dp), optional, intent(in) :: shortwave(:) !< L1 shortwave radiation [W m-2].
    real(dp), optional, intent(in) :: longwave(:) !< L1 longwave radiation [W m-2].
    integer(i4) :: i
    completed = .false.
    if (thermal_state%active) call this%thermal%accumulate(temperature_depth, air_temperature, shortwave, longwave)
    if (this%flow%lake_routing) then
      if (.not.present(lake_outflow) .or. .not.allocated(flow_state%lake_inflow)) then
        call error_message("river_router: lake topology requires lake inflow and outflow arrays")
      end if
      if (size(lake_outflow) /= this%flow%n_lakes .or. size(flow_state%lake_inflow) /= this%flow%n_lakes) then
        call error_message("river_router: lake flux array size does not match lake topology")
      end if
    end if
    call this%flow%accumulate(input_runoff)
    this%input_count = this%input_count + 1_i4
    if (this%input_count < this%muskingum%accumulations) return ! not yet ready to route
    this%input_count = 0_i4
    call this%flow%prepare(this%scaler, this%muskingum, this%input_step, this%node_runoff_cell)
    if (thermal_state%active) call this%thermal%prepare(this%scaler, this%muskingum, this%input_step, this%node_runoff_cell)
    ! Route every substep before publishing either completed mean.
    do i = 1_i4, this%muskingum%iterations
      call this%route(flow_state, thermal_state, lake_outflow)
      call this%flow%advance(flow_state)
      if (thermal_state%active) call this%thermal%advance(thermal_state)
    end do
    call this%flow%complete(flow_state, this%muskingum)
    if (thermal_state%active) call this%thermal%complete(thermal_state, this%muskingum)
    completed = .true.
  end subroutine river_router_update

  !> \brief Coordinate processor-owned static loops and synchronize completed DAG levels.
  subroutine river_router_route(this, flow_state, thermal_state, lake_outflow)
    class(river_router_t), intent(inout), target :: this !< Router scheduling each DAG level.
    type(river_flow_state_t), intent(inout), target :: flow_state !< mRM-owned discharge histories.
    type(river_thermal_state_t), intent(inout), target :: thermal_state !< mRM-owned thermal histories.
    real(dp), optional, intent(in) :: lake_outflow(:) !< Prescribed lake release [m3 s-1].
    integer(i8) :: i, last, n_levels
    logical :: parallel

    n_levels = this%river%order%n_levels

    if (any(this%parallel_level)) then
      !$omp parallel default(shared) private(i,last,parallel)
      i = 1_i8
      do while (i <= n_levels)
        parallel = this%parallel_level(i)
        last = this%serial_level_end(i)
        ! Finish the team before the final serial group; a wholly serial network never starts it.
        if (.not.parallel .and. last == n_levels) exit

        if (parallel) then
          ! Each processor completes the level before the next dependent processor starts.
          call this%flow%route_nodes(i, last, flow_state, this%muskingum, lake_outflow, parallel=.true.)
          if (thermal_state%active) call this%thermal%route_nodes(i, last, flow_state, thermal_state, this%muskingum, parallel=.true.)
        else
          ! One thread processes the complete serial group in topological order.
          !$omp single
          call this%flow%route_nodes(i, last, flow_state, this%muskingum, lake_outflow)
          if (thermal_state%active) call this%thermal%route_nodes(i, last, flow_state, thermal_state, this%muskingum)
          !$omp end single
        end if
        i = last + 1_i8
      end do
      !$omp end parallel
    end if

    if (.not.this%parallel_level(n_levels)) then
      ! Walk cached groups to locate the serial tail; entirely serial networks need no walk.
      i = 1_i8
      do while (this%serial_level_end(i) < n_levels)
        i = this%serial_level_end(i) + 1_i8
      end do
      call this%flow%route_nodes(i, n_levels, flow_state, this%muskingum, lake_outflow)
      if (thermal_state%active) call this%thermal%route_nodes(i, n_levels, flow_state, thermal_state, this%muskingum)
    end if
  end subroutine river_router_route

end module mo_river_router
