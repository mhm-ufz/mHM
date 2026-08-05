!> \file    mo_mrm_container.f90
!> \copydoc mo_mrm_container

!> \brief   Module for a mHM process container.
!> \version 0.1
!> \changelog
!! - Stephan Thober Sep 2026
!!   - initial version using river dag
!> \authors Sebastian Mueller, Stephan Thober
!> \date    Aug 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_mrm_container
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_kind, only: i2, i4, i8, dp
  use mo_nml, only: position_nml ! , close_nml
  use mo_exchange_type, only: exchange_t
  use mo_message, only: message, error_message
  use mo_river, only: river_t
  use mo_river_upscaler, only: river_upscaler_t
  use mo_river_router, only: river_router_t
  use mo_river_output, only: river_output_dataset
  use mo_grid, only: grid_t
  use mo_grid_io, only: output_dataset
  use mo_utils, only: is_close
  use mo_string_utils, only: n2s => num2str
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  use nml_helper, only: NML_OK
  use nml_config_mrm, only: nml_config_mrm_t
  use nml_output_mrm, only: nml_output_mrm_t

  character(len=*), parameter :: s = "mrm" !< logging scope
  public :: derive_mrm_output_timing

  !> \class   mrm_t
  !> \brief   Class for a single mRM process container.
  type, public :: mrm_t
    type(nml_config_mrm_t)     :: config              !< configuration of the mRM process container
    type(nml_output_mrm_t)     :: output_config       !< output configuration of the mRM process container
    type(exchange_t), pointer  :: exchange => null()  !< exchange container of the domain
    type(grid_t)               :: level3              !< mrm grid
    type(river_t)              :: river_l0            !< level-0 river network (for upscaling)
    type(river_t)              :: river               !< upscaled river network
    type(river_router_t)       :: router              !< river router
    type(river_upscaler_t)     :: upscaler            !< river upscaler for upscaling from level-0 to level-3 river network
    type(output_dataset)       :: ds_out              !< output dataset for gridded outputs
    type(river_output_dataset) :: ds_node_out         !< output dataset for river node based outputs
    real(dp), allocatable      :: discharge(:)        !< discharge array for all river nodes
    logical                    :: active = .false.    !< whether mRM participates in the configured domain
    logical                    :: scc_active = .false. !< whether scc based based upscaling is active
    logical                    :: read_restart = .false. !< whether to read restart file
    character(:), allocatable  :: restart_input_path  !< path to restart file to read
    logical                    :: write_restart = .false. !< whether to write restart file
    character(:), allocatable  :: restart_output_path !< path to restart file to write
    logical                    :: output_active = .false. !< whether output is enabled
    logical                    :: output_node_active = .false. !< whether node based output is enabled
    character(:), allocatable  :: output_path         !< path to output file
    character(:), allocatable  :: output_node_path    !< path to node output file
  contains
    procedure :: set_dims => mrm_set_dims
    procedure :: configure => mrm_configure
    procedure :: connect => mrm_connect
    procedure :: initialize => mrm_initialize
    procedure :: update => mrm_update
    procedure :: finalize => mrm_finalize
    procedure :: create_restart => mrm_create_restart
    procedure :: create_output => mrm_create_output
    procedure, private :: configure_parameters => mrm_configure_parameters
    procedure, private :: get_routing_steps => mrm_get_routing_steps
    procedure, private :: read_public_discharge => mrm_read_public_discharge
    procedure, private :: validate_timing => mrm_validate_timing
    procedure, private :: at_output_boundary => mrm_at_output_boundary
    procedure, private :: add_output => mrm_add_output
    procedure, private :: update_output => mrm_update_output
  end type mrm_t

contains

  !> \brief Derive the relation between a fixed-hour mRM output cadence and completed routing results.
  subroutine derive_mrm_output_timing(routing_step, output_frequency, routing_results, output_records, status)
    implicit none
    integer(i4), intent(in) :: routing_step !< [h] interval represented by one completed routing result
    integer(i4), intent(in) :: output_frequency !< [h] fixed output interval
    integer(i4), intent(out) :: routing_results !< completed routing results contributing to one output record
    integer(i4), intent(out) :: output_records !< output records filled by one completed routing result
    integer(i4), intent(out) :: status !< 0 success; positive values identify an invalid timing input

    routing_results = 0_i4
    output_records = 0_i4
    status = 0_i4
    if (routing_step < 1_i4) then
      status = 1_i4
      return
    end if
    if (output_frequency < 1_i4) then
      status = 2_i4
      return
    end if
    if (output_frequency >= routing_step) then
      if (mod(output_frequency, routing_step) /= 0_i4) then
        status = 3_i4
        return
      end if
      routing_results = output_frequency / routing_step
      output_records = 1_i4
    else
      if (mod(routing_step, output_frequency) /= 0_i4) then
        status = 3_i4
        return
      end if
      routing_results = 1_i4
      output_records = routing_step / output_frequency
    end if
  end subroutine derive_mrm_output_timing

  !> \brief Set runtime dimensions for generated mRM namelists.
  subroutine mrm_set_dims(self)
    class(mrm_t), intent(inout), target :: self
    character(1024) :: errmsg
    integer :: status

    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Error setting mRM config dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mrm_set_dims

  !> \brief Create a restart file for the mRM process container.
  subroutine mrm_create_restart(self)
    class(mrm_t), intent(inout), target :: self
    type(NcDataset) :: nc
    type(NcVariable) :: nc_var
    type(NcDimension) :: dims(0)

    if (self%router%input_count /= 0_i4) then
      log_fatal(*) "mRM restart can only be written at a completed routing boundary."
      error stop 1
    end if
    if ((self%output_active .or. self%output_node_active) .and. &
      .not.self%at_output_boundary(self%exchange%time)) then
      log_fatal(*) "mRM restart can only be written at an active mRM output boundary."
      error stop 1
    end if
    log_info(*) "Write mRM restart to file: ", self%restart_output_path
    nc = NcDataset(self%restart_output_path, "w")
    call self%level3%to_restart(nc)
    call self%river%to_restart(nc)
    call self%router%to_restart(nc)
    nc_var = nc%setVariable("mrm_discharge", "f64", [nc%getDimension("node")])
    call nc_var%setAttribute("long_name", "most recently completed mean routing discharge")
    call nc_var%setAttribute("units", "m3 s-1")
    call nc_var%setData(self%discharge)
    nc_var = nc%setVariable("mrm_meta", "i8", dims(:0)) ! scalar integer to indicate scc river
    call nc_var%setAttribute("routing_case", self%exchange%config%processes%routing)
    call nc_var%setAttribute("routing_gamma", self%exchange%parameters%get_process("routing"))
    call nc_var%setAttribute("time_stamp", self%exchange%time%str())
    call nc%close()
  end subroutine mrm_create_restart

  !> \brief Configure the mRM process container.
  subroutine mrm_configure(self, file, out_file)
    class(mrm_t), intent(inout), target :: self
    character(*), intent(in), optional :: file !< file containing the config_mrm namelists
    character(*), intent(in), optional :: out_file !< file containing the output_mrm namelists
    integer(i4) :: id(1) ! domain id
    integer(i4) :: case
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer :: status

    log_info(*) "Configure mRM"
    self%active = any([ &
      self%exchange%config%processes%routing, &
      self%exchange%config%processes%temperature_routing &
    ] /= 0_i4)
    if (.not.self%active) return
    ! get domain id
    id(1) = self%exchange%nml_domain_id

    ! read and check config
    if (present(file)) then
      path = self%exchange%get_path(file) ! get absolute path relative to cwd
      log_info(*) "Read mRM config: ", path
      status = self%config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading mRM config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%is_configured) then
      log_fatal(*) "mRM configuration not set."
      error stop 1
    end if
    status = self%config%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mRM config not valid: ", trim(errmsg)
      error stop 1
    end if

    ! output
    self%output_active = .true.
    if (present(out_file)) then
      path = self%exchange%get_path(out_file, root=.true.)
      log_info(*) "Read mRM output config: ", path
      status = self%output_config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        self%output_active = .false.
        log_warn(*) "mRM output disabled, config not found: ", trim(errmsg)
      end if
    end if
    if (self%output_config%is_configured) then
      status = self%output_config%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mRM output config invalid: ", trim(errmsg)
        error stop 1
      end if
    else
      self%output_active = .false.
      log_warn(*) "mRM output disabled, config not set."
    end if
    self%output_node_active = self%output_active ! further controlled by file presence in connect subroutine

    ! set output paths
    if (self%output_active) then
      status = self%config%is_set("output_path", idx=id, errmsg=errmsg)
      self%output_active = self%output_active .and. (status == NML_OK)
      if (status /= NML_OK) then
        log_warn(*) "mRM output disabled, path not set for domain ", n2s(id(1)), ": ", trim(errmsg)
      else
        self%output_path = self%exchange%get_path(self%config%output_path(id(1))) ! resolve relative path
      end if
    end if
    if (self%output_node_active) then
      status = self%config%is_set("output_node_path", idx=id, errmsg=errmsg)
      self%output_node_active = self%output_node_active .and. (status == NML_OK)
      if (status /= NML_OK) then
        log_warn(*) "mRM node output disabled, path not set for domain ", n2s(id(1)), ": ", trim(errmsg)
      else
        self%output_node_path = self%exchange%get_path(self%config%output_node_path(id(1))) ! resolve relative path
      end if
    end if

    ! check routing case
    call self%configure_parameters()

    case = self%exchange%config%processes%routing
    select case (case)
      ! case (1_i4)
      !   log_info(*) "mRM routing case 1"
      case (2_i4)
        log_info(*) "mRM routing case 2: constant celerity"
      case (3_i4)
        log_info(*) "mRM routing case 3: variable celerity based on slope"
      case default
        log_fatal(*) "mRM routing case ", n2s(case), " not implemented."
        error stop 1
    end select
  end subroutine mrm_configure

  !> \brief Read selected mRM parameter namelists and register them in execution order.
  subroutine mrm_configure_parameters(self)
    class(mrm_t), intent(inout), target :: self
    character(1024) :: errmsg
    integer :: status

    associate(parameter_config => self%exchange%config%parameters)
      select case (self%exchange%config%processes%routing)
        case (1_i4)
          if (.not.parameter_config%routing_1%is_configured .and. allocated(self%exchange%parameter_file)) then
            status = parameter_config%routing_1%from_file(file=self%exchange%parameter_file, errmsg=errmsg)
            call check_parameter_status(status, "routing_1", "read", errmsg)
          end if
          status = parameter_config%routing_1%is_valid(errmsg=errmsg)
          call check_parameter_status(status, "routing_1", "validate", errmsg)
          call self%exchange%parameters%add_process("routing", [ &
            parameter_config%routing_1%travel_time_constant, &
            parameter_config%routing_1%travel_time_river_length, &
            parameter_config%routing_1%travel_time_river_slope, &
            parameter_config%routing_1%travel_time_impervious, &
            parameter_config%routing_1%attenuation_river_slope], [character(64) :: &
            "travel_time_constant", "travel_time_river_length", "travel_time_river_slope", &
            "travel_time_impervious", "attenuation_river_slope"], group="routing_1")
        case (2_i4)
          if (.not.parameter_config%routing_2%is_configured .and. allocated(self%exchange%parameter_file)) then
            status = parameter_config%routing_2%from_file(file=self%exchange%parameter_file, errmsg=errmsg)
            call check_parameter_status(status, "routing_2", "read", errmsg)
          end if
          status = parameter_config%routing_2%is_valid(errmsg=errmsg)
          call check_parameter_status(status, "routing_2", "validate", errmsg)
          call self%exchange%parameters%add_process("routing", [parameter_config%routing_2%streamflow_celerity], &
            [character(64) :: "streamflow_celerity"], group="routing_2")
        case (3_i4)
          if (.not.parameter_config%routing_3%is_configured .and. allocated(self%exchange%parameter_file)) then
            status = parameter_config%routing_3%from_file(file=self%exchange%parameter_file, errmsg=errmsg)
            call check_parameter_status(status, "routing_3", "read", errmsg)
          end if
          status = parameter_config%routing_3%is_valid(errmsg=errmsg)
          call check_parameter_status(status, "routing_3", "validate", errmsg)
          call self%exchange%parameters%add_process("routing", [parameter_config%routing_3%slope_factor], &
            [character(64) :: "slope_factor"], group="routing_3")
      end select

      if (self%exchange%config%processes%temperature_routing == 1_i4) then
        if (.not.parameter_config%river_temperature_1%is_configured .and. allocated(self%exchange%parameter_file)) then
          status = parameter_config%river_temperature_1%from_file(file=self%exchange%parameter_file, errmsg=errmsg)
          call check_parameter_status(status, "river_temperature_1", "read", errmsg)
        end if
        status = parameter_config%river_temperature_1%is_valid(errmsg=errmsg)
        call check_parameter_status(status, "river_temperature_1", "validate", errmsg)
        call self%exchange%parameters%add_process("temperature_routing", [ &
          parameter_config%river_temperature_1%albedo_water, &
          parameter_config%river_temperature_1%pt_a_water, &
          parameter_config%river_temperature_1%emissivity_water, &
          parameter_config%river_temperature_1%turbulent_heat_exchange_coefficient], [character(64) :: &
          "albedo_water", "pt_a_water", "emissivity_water", "turbulent_heat_exchange_coefficient"], &
          group="river_temperature_1")
      end if
    end associate
  end subroutine mrm_configure_parameters

  !> \brief Fail with context when a generated routing parameter operation fails.
  subroutine check_parameter_status(status, block, action, errmsg)
    integer, intent(in) :: status
    character(*), intent(in) :: block
    character(*), intent(in) :: action
    character(*), intent(in) :: errmsg

    if (status /= NML_OK) then
      log_fatal(*) "mRM: failed to ", trim(action), " parameter block '", trim(block), "': ", trim(errmsg)
      error stop 1
    end if
  end subroutine check_parameter_status

  ! read initial values and populate exchange
  subroutine mrm_connect(self)
    use mo_datetime, only: datetime, timedelta, one_hour, one_day
    use mo_grid, only: grid_t
    use mo_river, only: river_t
    use mo_river_tools, only: read_scc_gauges
    use mo_os, only: path_ext, path_isfile
    use mo_string_utils, only: n2s => num2str

    implicit none

    class(mrm_t), target, intent(inout) :: self
    logical, allocatable        :: scc_latlon ! allocatable to be able to make it "not present" if not allocated
    character(:), allocatable   :: file
    real(dp), allocatable       :: scc_gauges(:,:)
    integer(i4)                 :: id(1)
    logical                     :: const_celerity
    integer(i4)                 :: model_step

    integer :: status
    character(1024) :: errmsg

    log_info(*) "Connect mRM"

    model_step = int(self%exchange%step / one_hour(), i4)

    ! get domain id
    id(1) = self%exchange%nml_domain_id
    ! check if scc_gauges_path is given
    self%scc_active = self%config%is_set("scc_gauges_path", idx=id, errmsg=errmsg) == NML_OK
    ! check routing case
    const_celerity = (self%exchange%config%processes%routing == 2_i4)
    ! get restart setting
    self%read_restart = self%config%read_restart(id(1))
    self%write_restart = self%config%write_restart(id(1))
    if (self%read_restart) then
      status = self%config%is_set("restart_input_path", idx=id, errmsg=errmsg)
      if (status /= NML_OK) call error_message("mRM restart input path not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg))
      self%restart_input_path = self%exchange%get_path(self%config%restart_input_path(id(1)))
    end if
    if (self%write_restart) then
      status = self%config%is_set("restart_output_path", idx=id, errmsg=errmsg)
      if (status /= NML_OK) call error_message("mRM restart output path not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg))
      self%restart_output_path = self%exchange%get_path(self%config%restart_output_path(id(1)))
    end if

    ! Runoff may be provided by dynamic input and connected during update.
    call self%exchange%runoff_total%require("mRM", .true., check_data=.false.)
    call self%exchange%fdir%require("mRM", .not.self%read_restart)
    call self%exchange%slope%require("mRM", .not.const_celerity .and. .not.self%read_restart)

    ! derive level-3 grid
    if (self%read_restart) then
      scope_info(s,*) "Read mRM grid from restart file: ", self%restart_input_path
      call self%level3%from_restart(self%restart_input_path)
    else
      if (.not.ieee_is_finite(self%exchange%level3_resolution) .or. self%exchange%level3_resolution <= 0.0_dp) then
        log_fatal(*) "mRM: level3 resolution not configured (expected from config_resolution/route)."
        error stop 1
      end if
      scope_info(s,*) "Derive mRM grid from level-0 grid with resolution: ", self%exchange%level3_resolution
      call self%exchange%level0%gen_grid(self%level3, target_resolution=self%exchange%level3_resolution)
    end if
    ! if (self%level3%has_aux_coords()) call self%level3%estimate_aux_vertices()
    scope_debug(s,*) "level0 ncells", n2s(self%exchange%level0%ncells)
    scope_debug(s,*) "level0 cellsize", n2s(self%exchange%level0%cellsize)
    scope_debug(s,*) "level1 ncells", n2s(self%exchange%level1%ncells)
    scope_debug(s,*) "level1 cellsize", n2s(self%exchange%level1%cellsize)
    scope_debug(s,*) "level3 ncells", n2s(self%level3%ncells)
    scope_debug(s,*) "level3 cellsize", n2s(self%level3%cellsize)

    ! create rivers
    if (self%read_restart) then
      call self%river%from_restart_file(self%restart_input_path, self%level3)
    else if (is_close(self%level3%cellsize, self%exchange%level0%cellsize)) then
      ! TODO: the upscaler should handle also the case of no upscaling (level0 == level11)
      scope_info(s,*) "level-0 and level-3 river network are equal of size:", n2s(self%exchange%level3%ncells)
      call self%river%from_fdir(int(self%exchange%fdir%data, i2), self%level3)
    else
      scope_info(s,*) "Create level-0 river network of size:", n2s(self%exchange%level0%ncells)
      ! TODO: make fdir i2
      call self%river_l0%from_fdir(int(self%exchange%fdir%data, i2), self%exchange%level0)
      scope_info(s,*) "Calculate facc on level-0"
      call self%river_l0%calc_order()
      call self%river_l0%calc_facc()
      ! check SCC config
      if (self%scc_active) then
        file = self%exchange%get_path(self%config%scc_gauges_path(id(1)))
        scope_info(s,*) "Read SCC gauges from file: ", file
        allocate(scc_latlon)  ! if not allocated, it is not present as optional argument
        call read_scc_gauges(file, scc_gauges, scc_latlon)
      end if
      ! scc_gauges/scc_latlon not present if not allocated
      scope_info(s,*) "Initialize upscaler and upscale river network to level-3"
      call self%upscaler%init(self%river_l0, self%river, self%level3, scc_gauges, scc_latlon)
      if (self%config%is_set("diagnostics_path", idx=id) == NML_OK) then
        file = self%exchange%get_path(self%config%diagnostics_path(id(1)))
        log_info(*) "Write mRM upscaling diagnostics to file: ", file
        call self%river_l0%export( &
          path        = file, &
          sub_map     = self%upscaler%scc_map, &
          leaving     = self%upscaler%leaving_cells, &
          stream_mask = self%upscaler%stream_mask, &
          stream_sub  = self%upscaler%stream_sub, &
          highlight   = self%upscaler%is_link_start.or.self%river_l0%is_sink, &
          factor      = self%upscaler%upscaler%factor &
        )
      end if
    end if

    ! populate exchange type
    allocate(self%discharge(self%river%n_nodes))
    call self%exchange%discharge%publish_local("mRM", self%discharge, model_step)
    self%exchange%level3 => self%level3
  end subroutine mrm_connect

  ! set initial values like timestep 0
  subroutine mrm_initialize(self)
    use mo_datetime, only: datetime, timedelta, one_hour, one_day
    class(mrm_t), target, intent(inout) :: self

    integer(i4)           :: id(1)
    integer(i4)           :: input_step, model_step
    logical               :: const_celerity
    real(dp), allocatable :: gamma(:)

    log_info(*) "Initialize mRM"

    ! get domain id
    id(1) = self%exchange%nml_domain_id
    ! calculate celerity
    gamma = self%exchange%parameters%get_process("routing")
    const_celerity = (self%exchange%config%processes%routing == 2_i4)
    call self%get_routing_steps(input_step, model_step)

    if (self%read_restart) then
      if (.not.self%exchange%parameters%is_default("routing")) then
        log_fatal(*) "mRM: routing parameters cannot be overridden when routing state is read from restart."
        error stop 1
      end if
      scope_info(s,*) "Read routing state from restart file: ", self%restart_input_path
      ! TODO: warn about gamma mismatch between restart and config
      call self%router%from_restart_file( &
        path              = self%restart_input_path, &
        river             = self%river, &
        input_grid        = self%exchange%level1, &
        input_step        = input_step, &
        model_step        = model_step, &
        max_route_step    = real(self%config%max_route_step(id(1)), dp), &
        root_levels       = self%config%river_net_order_root_based(id(1)), &
        omp_level_thresh  = int(self%config%river_net_omp_level_min(id(1)), i8), &
        read_fluxes       = self%config%read_restart_fluxes(id(1)))
    else
      if (allocated(self%river_l0%celerity)) deallocate(self%river_l0%celerity)
      if (allocated(self%river%celerity)) deallocate(self%river%celerity)
      ! NOTE: if slope data pointer is null (i.e. slope not provided), optional slope will be seen as "not present"
      if (is_close(self%level3%cellsize, self%exchange%level0%cellsize)) then
        call self%river%calc_celerity(gamma=gamma(1), slope=self%exchange%slope%data, constant_celerity=const_celerity)
      else
        call self%upscaler%calc_celerity(gamma=gamma(1), slope=self%exchange%slope%data, constant_celerity=const_celerity)
      end if
      scope_info(s,*) "Initialize router"
      call self%router%init( &
        river            = self%river, &
        input_grid       = self%exchange%level1, &
        input_step       = input_step, &
        model_step       = model_step, &
        max_route_step   = real(self%config%max_route_step(id(1)), dp), &
        root_levels      = self%config%river_net_order_root_based(id(1)), &
        omp_level_thresh = int(self%config%river_net_omp_level_min(id(1)), i8))
    end if

    scope_debug(s,*) "router%routing_substep: ", self%router%routing_substep
    scope_debug(s,*) "router%routing_step: ", self%router%routing_step
    scope_debug(s,*) "last level in parallel: ", self%router%last_parallel_level, "/", self%router%river%order%n_levels
    call self%exchange%discharge%set_stepping("mRM", self%router%routing_step)
    if (self%read_restart) then
      call self%read_public_discharge()
    else
      self%discharge = 0.0_dp
    end if

    call self%validate_timing()

    call self%create_output()

  end subroutine mrm_initialize

  !> \brief Convert runoff support and the global model cadence to fixed hours for routing.
  subroutine mrm_get_routing_steps(self, input_step, model_step)
    use mo_datetime, only: one_hour
    use mo_grid_io, only: daily, monthly, yearly, no_time, varying
    class(mrm_t), intent(in), target :: self
    integer(i4), intent(out) :: input_step, model_step

    model_step = int(self%exchange%step / one_hour(), i4)
    if (model_step < 1_i4 .or. self%exchange%step /= model_step * one_hour()) then
      log_fatal(*) "mRM requires the global model step to be a positive whole number of hours."
      error stop 1
    end if

    select case (self%exchange%runoff_total%stepping)
      case (daily)
        input_step = 24_i4
      case (1_i4:)
        input_step = self%exchange%runoff_total%stepping
      case (no_time)
        log_fatal(*) "mRM cannot route static runoff; runoff stepping metadata must describe a fixed interval."
        error stop 1
      case (monthly)
        log_fatal(*) "mRM cannot route monthly runoff; runoff support must be a fixed number of hours."
        error stop 1
      case (yearly)
        log_fatal(*) "mRM cannot route yearly runoff; runoff support must be a fixed number of hours."
        error stop 1
      case (varying)
        log_fatal(*) "mRM cannot route irregular runoff; runoff support must be a fixed number of hours."
        error stop 1
      case default
        log_fatal(*) "mRM runoff stepping metadata is missing or unsupported: ", self%exchange%runoff_total%stepping
        error stop 1
    end select

    if (input_step < model_step) then
      log_fatal(*) "mRM runoff support (", input_step, "h) is finer than the model step (", model_step, &
        "h); input aggregation is not implemented."
      error stop 1
    end if
    if (mod(input_step, model_step) /= 0_i4) then
      log_fatal(*) "mRM runoff support (", input_step, "h) must be divisible by the model step (", model_step, "h)."
      error stop 1
    end if
  end subroutine mrm_get_routing_steps

  !> \brief Restore the mRM-owned published discharge from a routing-boundary restart.
  subroutine mrm_read_public_discharge(self)
    class(mrm_t), intent(inout), target :: self
    type(NcDataset) :: nc
    type(NcVariable) :: nc_var
    integer(i4) :: id(1)

    id(1) = self%exchange%nml_domain_id
    if (.not.self%config%read_restart_fluxes(id(1))) then
      self%discharge = 0.0_dp
      return
    end if

    nc = NcDataset(self%restart_input_path, "r")
    if (nc%hasVariable("mrm_discharge")) then
      nc_var = nc%getVariable("mrm_discharge")
      call nc_var%readInto(self%discharge)
    else
      self%discharge = self%router%previous_discharge
      log_warn(*) "mRM restart has no published discharge; using the final internal routing state."
    end if
    call nc%close()
  end subroutine mrm_read_public_discharge

  !> \brief Validate routing/output compatibility and boundary-only restart rules.
  subroutine mrm_validate_timing(self)
    use mo_datetime, only: one_hour
    use mo_grid_io, only: daily, monthly, yearly, no_time
    class(mrm_t), intent(in), target :: self
    integer(i4) :: frequency
    integer(i4) :: output_records
    integer(i4) :: routing_results
    integer(i4) :: status
    integer(i4) :: total_hours
    logical :: output_enabled
    logical :: output_boundary

    total_hours = nint((self%exchange%end_time - self%exchange%start_time) / one_hour(), i4)
    if (mod(total_hours, self%router%routing_step) /= 0_i4) then
      log_fatal(*) "mRM simulation end must be a routing boundary: duration ", total_hours, &
        "h is not divisible by routing_step=", self%router%routing_step, "h."
      error stop 1
    end if

    output_enabled = self%output_active .or. self%output_node_active
    if (.not.output_enabled) return
    frequency = self%output_config%output_frequency
    output_boundary = self%at_output_boundary(self%exchange%end_time)
    select case (frequency)
      case (daily)
        if (mod(24_i4, self%router%routing_step) /= 0_i4) then
          log_fatal(*) "Daily mRM output requires routing_step to divide 24h."
          error stop 1
        end if
        if (.not.self%exchange%start_time%is_new_day()) then
          log_fatal(*) "Daily mRM output requires the simulation to start at a day boundary."
          error stop 1
        end if
      case (monthly)
        if (mod(24_i4, self%router%routing_step) /= 0_i4) then
          log_fatal(*) "Monthly mRM output requires routing_step to divide 24h."
          error stop 1
        end if
        if (.not.self%exchange%start_time%is_new_month()) then
          log_fatal(*) "Monthly mRM output requires the simulation to start at a month boundary."
          error stop 1
        end if
      case (yearly)
        if (mod(24_i4, self%router%routing_step) /= 0_i4) then
          log_fatal(*) "Yearly mRM output requires routing_step to divide 24h."
          error stop 1
        end if
        if (.not.self%exchange%start_time%is_new_year()) then
          log_fatal(*) "Yearly mRM output requires the simulation to start at a year boundary."
          error stop 1
        end if
      case (no_time)
        continue
      case default
        call derive_mrm_output_timing( &
          self%router%routing_step, frequency, routing_results, output_records, status)
        if (status /= 0_i4) then
          log_fatal(*) "mRM output_frequency=", frequency, "h and routing_step=", &
            self%router%routing_step, "h must be whole multiples of one another."
          error stop 1
        end if
    end select

    if (self%write_restart .and. .not.output_boundary) then
      log_fatal(*) "mRM restart output time must coincide with an mRM output boundary."
      error stop 1
    end if
  end subroutine mrm_validate_timing

  !> \brief Report whether a timestamp is a configured mRM output boundary.
  logical function mrm_at_output_boundary(self, time) result(boundary)
    use mo_datetime, only: datetime, one_hour
    use mo_grid_io, only: daily, monthly, yearly, no_time
    class(mrm_t), intent(in), target :: self
    type(datetime), intent(in) :: time
    integer(i4) :: elapsed_hours

    if (.not.self%output_active .and. .not.self%output_node_active) then
      boundary = .true.
      return
    end if
    select case (self%output_config%output_frequency)
      case (daily)
        boundary = time%is_new_day()
      case (monthly)
        boundary = time%is_new_month()
      case (yearly)
        boundary = time%is_new_year()
      case (no_time)
        boundary = time == self%exchange%end_time
      case default
        elapsed_hours = nint((time - self%exchange%start_time) / one_hour(), i4)
        boundary = elapsed_hours >= 0_i4 .and. &
          mod(elapsed_hours, self%output_config%output_frequency) == 0_i4
    end select
  end function mrm_at_output_boundary

  ! perform routing within time loop
  subroutine mrm_update(self)
    class(mrm_t), target, intent(inout) :: self
    logical :: completed
    log_trace(*) "Update mRM"

    ! route runoff
    call self%router%update(self%exchange%runoff_total%data, self%discharge, completed)
    if (completed) call self%update_output()
  end subroutine mrm_update

  !> \brief Add one newly completed routing-step mean to all active mRM outputs.
  subroutine mrm_add_output(self)
    class(mrm_t), intent(inout), target :: self

    if (self%output_node_active) call self%ds_node_out%update("discharge", self%discharge)
    if (self%output_active) then
      if (self%scc_active) then
        call self%ds_out%update("discharge", self%river%select_cell_values(self%discharge))
      else
        call self%ds_out%update("discharge", self%discharge)
      end if
    end if
  end subroutine mrm_add_output

  !> \brief Backfill or aggregate completed routing results at the configured mRM output cadence.
  subroutine mrm_update_output(self)
    use mo_datetime, only: datetime, timedelta, one_hour
    use mo_grid_io, only: daily, monthly, yearly, no_time
    class(mrm_t), intent(inout), target :: self
    type(datetime) :: write_time
    integer(i4) :: elapsed_hours
    integer(i4) :: frequency
    integer(i4) :: i
    integer(i4) :: records
    integer(i4) :: routing_results
    integer(i4) :: status
    logical :: write_stamp

    if (.not.self%output_active .and. .not.self%output_node_active) return
    frequency = self%output_config%output_frequency

    ! A completed routing mean can be repeated over finer fixed-hour output records.
    if (frequency > 0_i4 .and. frequency < self%router%routing_step) then
      call derive_mrm_output_timing( &
        self%router%routing_step, frequency, routing_results, records, status)
      if (status /= 0_i4) then
        log_fatal(*) "Invalid fixed-hour mRM output cadence after initialization."
        error stop 1
      end if
      do i = 1_i4, records
        write_time = self%exchange%time
        call write_time%add(timedelta(hours=-(self%router%routing_step - i * frequency)))
        call self%add_output()
        if (self%output_node_active) call self%ds_node_out%write(write_time)
        if (self%output_active) call self%ds_out%write(write_time)
      end do
      return
    end if

    ! Equal or coarser outputs receive exactly one sample per completed routing step.
    call self%add_output()
    write_stamp = .false.
    select case(frequency)
      case(daily)
        if (self%exchange%time%is_new_day()) write_stamp = .true.
      case(monthly)
        if (self%exchange%time%is_new_month()) write_stamp = .true.
      case(yearly)
        if (self%exchange%time%is_new_year()) write_stamp = .true.
      case(no_time) ! once
        if (self%exchange%time == self%exchange%end_time) write_stamp = .true.
      case default ! fixed hours
        elapsed_hours = nint((self%exchange%time - self%exchange%start_time) / one_hour(), i4)
        if (mod(elapsed_hours, frequency) == 0_i4) write_stamp = .true.
    end select

    if (write_stamp) then
      if (self%output_node_active) call self%ds_node_out%write(self%exchange%time)
      if (self%output_active) call self%ds_out%write(self%exchange%time)
    end if
  end subroutine mrm_update_output

  subroutine mrm_finalize(self)
    class(mrm_t), intent(inout), target :: self
    log_info(*) "Finalize mRM"
    if (self%write_restart) call self%create_restart()
    if (self%output_active) then
      call self%ds_out%close()
      log_info(*) "Close mRM grid output file: ", self%output_path
    else
      log_info(*) "No mRM grid output file will be written"
    end if
    if (self%output_node_active) then
      call self%ds_node_out%close()
      log_info(*) "Close mRM node output file: ", self%output_node_path
    else
      log_info(*) "No mRM node output file will be written"
    end if
  end subroutine mrm_finalize

  subroutine mrm_cleanup(self)
    class(mrm_t), intent(inout), target :: self
    log_info(*) "Cleanup mRM"
    ! deallocate arrays, close files, ...
    call self%upscaler%destroy()
    call self%river_l0%clean()
  end subroutine mrm_cleanup

  subroutine mrm_create_output(self)
    use mo_grid_io, only: var, time_units_delta
    class(mrm_t), intent(inout), target :: self

    integer(i4) :: timestamp
    character(:), allocatable :: delta, dtype
    type(var), allocatable :: vars(:)

    ! shortcut
    if (.not.self%output_active .and. .not.self%output_node_active) return

    ! general config
    timestamp = self%output_config%output_time_reference
    delta = time_units_delta(self%output_config%output_frequency, timestamp)

    ! create output variables
    dtype = "f64"
    if (.not.self%output_config%output_double_precision) dtype = "f32"
    allocate(vars(0))
    if (self%output_config%out_Qrouted) vars = [vars, var(name="discharge", units="m3 s-1", dtype=dtype, avg=.true.)]

    ! create grid based output
    if (self%output_active) then
      log_info(*) "Create mRM grid based output file: ", self%output_path
      ! create output dataset
      call self%ds_out%init( &
        path        = self%output_path, &
        grid        = self%level3, &
        vars        = vars, &
        start_time  = self%exchange%start_time, &
        delta       = delta, &
        timestamp   = timestamp, &
        deflate_level = self%output_config%output_deflate_level)
    end if

    ! create node based output
    if (self%output_node_active) then
      log_info(*) "Create mRM node based output file: ", self%output_node_path
      call self%ds_node_out%init( &
        path        = self%output_node_path, &
        river       = self%river, &
        vars        = vars, &
        start_time  = self%exchange%start_time, &
        delta       = delta, &
        timestamp   = timestamp, &
        deflate_level = self%output_config%output_deflate_level)
    end if

  end subroutine mrm_create_output
end module mo_mrm_container
