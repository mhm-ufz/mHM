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
    logical                    :: scc_active          !< whether scc based based upscaling is active
    logical                    :: read_restart        !< whether to read restart file
    character(:), allocatable  :: restart_input_path  !< path to restart file to read
    logical                    :: write_restart       !< whether to write restart file
    character(:), allocatable  :: restart_output_path !< path to restart file to write
    logical                    :: output_active       !< whether output is enabled
    logical                    :: output_node_active  !< whether node based output is enabled
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
  end type mrm_t

contains

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
    class(mrm_t), intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: nc_var
    type(NcDimension) :: dims(0)
    log_info(*) "Write mRM restart to file: ", self%restart_output_path
    nc = NcDataset(self%restart_output_path, "w")
    call self%level3%to_restart(nc)
    call self%river%to_restart(nc)
    call self%router%to_restart(nc)
    nc_var = nc%setVariable("mrm_meta", "i8", dims(:0)) ! scalar integer to indicate scc river
    call nc_var%setAttribute("routing_case", self%exchange%parameters%config%processes%routing)
    call nc_var%setAttribute("routing_gamma", self%exchange%parameters%get_process(8_i4))
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
      log_info(*) "Read mRM output config: ", out_file ! out file is already absolute path
      status = self%output_config%from_file(file=out_file, errmsg=errmsg)
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
    case = self%exchange%parameters%config%processes%routing
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

    integer :: status
    character(1024) :: errmsg

    log_info(*) "Connect mRM"

    ! get domain id
    id(1) = self%exchange%nml_domain_id
    ! check if scc_gauges_path is given
    self%scc_active = self%config%is_set("scc_gauges_path", idx=id, errmsg=errmsg) == NML_OK
    ! check routing case
    const_celerity = (self%exchange%parameters%config%processes%routing == 2_i4)
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
    self%exchange%discharge%provided = .true.
    self%exchange%discharge%data => self%discharge
    self%exchange%level3 => self%level3
  end subroutine mrm_connect

  ! set initial values like timestep 0
  subroutine mrm_initialize(self)
    use mo_datetime, only: datetime, timedelta, one_hour, one_day
    class(mrm_t), target, intent(inout) :: self

    integer(i4)           :: id(1)
    logical               :: const_celerity
    real(dp), allocatable :: gamma(:)

    log_info(*) "Initialize mRM"

    ! get domain id
    id(1) = self%exchange%nml_domain_id
    ! calculate celerity
    gamma = self%exchange%parameters%get_process(8_i4)  ! routing still process 8
    const_celerity = (self%exchange%parameters%config%processes%routing == 2_i4)

    if (self%read_restart) then
      scope_info(s,*) "Read routing state from restart file: ", self%restart_input_path
      ! TODO: warn about gamma mismatch between restart and config
      call self%router%from_restart_file( &
        path              = self%restart_input_path, &
        river             = self%river, &
        input_grid        = self%exchange%level1, &
        input_step        = int(self%exchange%step/one_hour(), i4), &
        max_route_step    = real(self%config%max_route_step(id(1)), dp), &
        root_levels       = self%config%river_net_order_root_based(id(1)), &
        omp_level_thresh  = int(self%config%river_net_omp_level_min(id(1)), i8), &
        read_fluxes       = self%config%read_restart_fluxes(id(1)))
    else
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
        input_step       = int(self%exchange%step/one_hour(), i4), &
        max_route_step   = real(self%config%max_route_step(id(1)), dp), &
        root_levels      = self%config%river_net_order_root_based(id(1)), &
        omp_level_thresh = int(self%config%river_net_omp_level_min(id(1)), i8))
    end if

    scope_debug(s,*) "router%step: ", self%router%step
    scope_debug(s,*) "last level in parallel: ", self%router%last_parallel_level, "/", self%router%river%order%n_levels
    if (self%router%step > 3600.0_dp) then
      log_error(*) "mRM routing time step is larger than 1 hour. This is not yet supported."
      stop 1
    end if

    call self%create_output()

  end subroutine mrm_initialize

  ! perform routing within time loop
  subroutine mrm_update(self)
    use mo_grid_io, only: daily, monthly, yearly, no_time
    class(mrm_t), target, intent(inout) :: self
    logical :: write_stamp
    log_trace(*) "Update mRM"

    ! route runoff
    call self%router%update(self%exchange%runoff_total%data, self%discharge)

    ! write time-stamp depending on config
    write_stamp = .false.
    select case(self%output_config%output_frequency)
      case(daily)
        if (self%exchange%time%is_new_day()) write_stamp = .true.
      case(monthly)
        if (self%exchange%time%is_new_month()) write_stamp = .true.
      case(yearly)
        if (self%exchange%time%is_new_year()) write_stamp = .true.
      case(no_time) ! once
        if (self%exchange%time == self%exchange%end_time) write_stamp = .true.
      case default ! every n time steps
        if (mod(self%exchange%step_count, self%output_config%output_frequency) == 0_i4) write_stamp = .true.
    end select

    ! update output
    if (self%output_node_active) then
      call self%ds_node_out%update("discharge", self%discharge)
      if (write_stamp) call self%ds_node_out%write(self%exchange%time)
    end if
    if (self%output_active) then
      if (self%scc_active) then
        call self%ds_out%update("discharge", self%river%select_cell_values(self%discharge))
      else
        call self%ds_out%update("discharge", self%discharge)
      end if
      if (write_stamp) call self%ds_out%write(self%exchange%time)
    end if
  end subroutine mrm_update

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
