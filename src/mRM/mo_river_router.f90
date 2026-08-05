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
  use mo_utils, only: equal, optval
  use mo_string_utils, only: n2s => num2str
  use mo_river, only: river_t
  use mo_grid, only: grid_t, bottom_up, cartesian
  use mo_grid_scaler, only: scaler_t, up_sum, down_nearest, down_scaling
  use mo_message, only: error_message, warn_message, message
  use mo_datetime, only: HOUR_SECONDS
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension

  implicit none
  private

  !> \name Routing Constants
  !> \brief Constants controlling the routing algorithms.
  !!@{
  !> valid routing time steps [s] (integer minute fractions of hour or integer hour fractions of day)
  real(dp), public, dimension(19), parameter :: routing_steps = &
    [   60.0_dp,   120.0_dp,   180.0_dp,   240.0_dp,   300.0_dp,   360.0_dp,            &
       600.0_dp,   720.0_dp,   900.0_dp,  1200.0_dp,  1800.0_dp,  3600.0_dp,            &
      7200.0_dp, 10800.0_dp, 14400.0_dp, 21600.0_dp, 28800.0_dp, 43200.0_dp, 86400.0_dp ]
  !> space weighting of routing is set to 0, since this parameter has no effect on the routing results, see Thober et al. 2017
  real(dp), public, parameter :: routing_space_weight = 0.0_dp
  ! time slice selector in muskingum schema for link_inflow and link_routed
  integer(i4), public, parameter :: current = 1_i4 !< selector for current time step
  integer(i4), public, parameter :: previous = 2_i4 !< selector for previous time step
  real(dp), public, parameter :: liter_to_m3 = 0.001_dp !< conversion factor from liter to m3
  public :: derive_routing_timing
  !!@}

  !> \class inflow_t
  !> \brief Inflow handler
  type, public :: inflow_t
    integer(i8) :: count = 0_i8 !< number of inflow nodes
    integer(i8), allocatable :: nodes(:) !< list of inflow nodes [id]
    logical, allocatable :: headwater(:) !< whether to consider headwater cells of inflow gauge
    real(dp), allocatable :: rate(:) !< inflow rate [m3 s-1]
  ! contains
    ! procedure, public :: init => inflow_init
    ! procedure, public :: add => inflow_add
  end type inflow_t

  !> \class river_router_t
  !> \brief River network router
  type, public :: river_router_t
    type(river_t), pointer :: river => null() !< river definition to route on
    type(grid_t), pointer :: input_grid => null() !< grid the input (e.g. runoff) is defined on
    type(scaler_t) :: scaler !< input rescaler
    integer(i4) :: input_step !< [h] time step size of the input in hours
    integer(i4) :: input_count !< [-] number of accumulated model updates for one routing step
    integer(i4) :: routing_step !< [h] interval represented by one completed routing result
    real(dp) :: routing_substep !< [s] numerical routing step for meeting the CFL condition
    integer(i4) :: iterations !< number of numerical substeps for each routing step
    integer(i4) :: accumulations !< number of model updates for each routing step
    type(inflow_t) :: inflow_handler !< inflow handler
    !> [mm] accumulated runoff size(input_grid\%ncells)
    real(dp), allocatable :: acc_runoff(:)
    !> [m3 s-1] runoff flux for current cell as intermediate result for SCC, size(river\%grid\%ncells)
    real(dp), allocatable :: scaled_runoff(:)
    !> [m3 s-1] runoff flux to route for current river node, size(river\%n_nodes)
    real(dp), allocatable :: runoff(:)
    !> [m3 s-1] inflow flux to link starting at current node from current time step, size(river\%n_nodes)
    real(dp), allocatable :: discharge(:)
    !> [m3 s-1] inflow flux to link starting at current node from previous time step, size(river\%n_nodes)
    real(dp), allocatable :: previous_discharge(:)
    !> [m3 s-1] routed flux from link starting at current node from current time step, size(river\%n_nodes)
    real(dp), allocatable :: tributary(:)
    !> [m3 s-1] routed flux from link starting at current node from previous time step, size(river\%n_nodes)
    real(dp), allocatable :: previous_tributary(:)
    real(dp), allocatable :: nu1(:) !< Muskingum parameter nu1, size(river\%n_nodes) -> C1 = nu2, C2 = nu1-nu2, C3=1-nu1
    real(dp), allocatable :: nu2(:) !< Muskingum parameter nu2, size(river\%n_nodes) -> C1 = nu2, C2 = nu1-nu2, C3=1-nu1
    !> minimum size of river-levels to route in parallel (default: threads * 8, 0 to run all in serial, 1 to run all in parallel)
    integer(i8) :: omp_level_thresh = 0_i8
    integer(i8) :: last_parallel_level = 0_i8 !< last level to run in parallel
  contains
    procedure, public :: init => river_router_init
    procedure, public :: update => river_router_update
    procedure, public :: allocate => river_router_allocate
    procedure, public :: deallocate => river_router_deallocate
    procedure, private :: setup_grids => river_router_setup_grids
    procedure, private :: setup_muskingum => river_router_setup_muskingum
    procedure, private :: setup_timing => river_router_setup_timing
    procedure, private :: setup_parallelization => river_router_setup_parallelization
    procedure, private :: scale_runoff => river_router_scale_runoff
    procedure, private :: route => river_router_route
    procedure, private :: time_shift => river_router_time_shift
    ! restart
    procedure, public :: to_restart_dataset => river_router_to_restart_dataset, to_restart_file => river_router_to_restart_file
    procedure, public :: from_restart_dataset => river_router_from_restart_dataset, from_restart_file => river_router_from_restart_file
    generic, public :: to_restart => to_restart_dataset, to_restart_file
    generic, public :: from_restart => from_restart_dataset, from_restart_file
  end type river_router_t

contains

  !> \brief Derive the completed routing step and its update/iteration counts.
  subroutine derive_routing_timing(input_step, model_step, routing_substep, routing_step, iterations, accumulations, status)
    implicit none
    integer(i4), intent(in) :: input_step !< [h] runoff support interval
    integer(i4), intent(in) :: model_step !< [h] model update interval
    real(dp), intent(in) :: routing_substep !< [s] numerical routing substep
    integer(i4), intent(out) :: routing_step !< [h] completed routing interval
    integer(i4), intent(out) :: iterations !< numerical routing iterations per routing step
    integer(i4), intent(out) :: accumulations !< model updates per routing step
    integer(i4), intent(out) :: status !< 0 success; positive values identify an invalid timing input
    integer(i8) :: a, b, divisor
    integer(i8) :: model_seconds, routing_seconds, substep_seconds

    routing_step = 0_i4
    iterations = 0_i4
    accumulations = 0_i4
    status = 0_i4
    if (input_step < 1_i4) then
      status = 1_i4
      return
    end if
    if (model_step < 1_i4) then
      status = 2_i4
      return
    end if
    if (input_step < model_step) then
      status = 3_i4
      return
    end if
    if (mod(input_step, model_step) /= 0_i4) then
      status = 4_i4
      return
    end if

    substep_seconds = nint(routing_substep, i8)
    if (substep_seconds < 1_i8 .or. .not.equal(routing_substep, real(substep_seconds, dp))) then
      status = 5_i4
      return
    end if
    model_seconds = int(model_step, i8) * int(HOUR_SECONDS, i8)
    a = model_seconds
    b = substep_seconds
    do while (b /= 0_i8)
      divisor = mod(a, b)
      a = b
      b = divisor
    end do
    routing_seconds = (model_seconds / a) * substep_seconds
    if (mod(routing_seconds, int(HOUR_SECONDS, i8)) /= 0_i8 .or. &
      routing_seconds / int(HOUR_SECONDS, i8) > int(huge(routing_step), i8)) then
      status = 6_i4
      return
    end if

    routing_step = int(routing_seconds / int(HOUR_SECONDS, i8), i4)
    iterations = int(routing_seconds / substep_seconds, i4)
    accumulations = routing_step / model_step
  end subroutine derive_routing_timing

  !> \brief Set river and input grid for router and initialize scaler if not yet initialized or input grid has changed.
  subroutine river_router_setup_grids(this, river, input_grid)
    implicit none
    class(river_router_t), intent(inout) :: this
    type(river_t), pointer, intent(in) :: river !< river definition
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    this%river => river
    this%input_grid => input_grid
    if (.not.associated(this%scaler%source_grid,this%input_grid)) then
      ! if scaler is not yet initialized or input grid has changed, initialize scaler to only do it once
      ! TODO: this is dangerous, since calling init twice is not possible
      ! TODO: scaler should also safely re-allocate attributes when being initialized
      call this%scaler%init( &
        source_grid=this%input_grid, &
        target_grid=this%river%grid, &
        upscaling_operator=up_sum, &       ! if L3 coarser than L1: sum runoff
        downscaling_operator=down_nearest) ! if L3 finer than L1: distribute same value on fine cells
    end if
  end subroutine river_router_setup_grids

  !> \brief Deallocate all arrays in river router.
  subroutine river_router_deallocate(this)
    implicit none
    class(river_router_t), intent(inout) :: this
    if (allocated(this%scaled_runoff)) deallocate(this%scaled_runoff)
    if (allocated(this%runoff)) deallocate(this%runoff)
    if (allocated(this%acc_runoff)) deallocate(this%acc_runoff)
    if (allocated(this%discharge)) deallocate(this%discharge)
    if (allocated(this%previous_discharge)) deallocate(this%previous_discharge)
    if (allocated(this%tributary)) deallocate(this%tributary)
    if (allocated(this%previous_tributary)) deallocate(this%previous_tributary)
    if (allocated(this%nu1)) deallocate(this%nu1)
    if (allocated(this%nu2)) deallocate(this%nu2)
  end subroutine river_router_deallocate

  !> \brief Allocate arrays in river router based on river and input grid size.
  !> \details Deallocate first to avoid memory leaks if allocate is called multiple times.
  subroutine river_router_allocate(this)
    implicit none
    class(river_router_t), intent(inout) :: this
    ! only need scaled runoff as intermediate result in case of SCC
    call this%deallocate() ! deallocate first to avoid memory leaks if allocate is called multiple times
    if (this%river%scc) allocate(this%scaled_runoff(this%river%grid%ncells))
    allocate(this%runoff(this%river%n_nodes))
    allocate(this%acc_runoff(this%input_grid%ncells))
    allocate(this%discharge(this%river%n_nodes))
    allocate(this%previous_discharge(this%river%n_nodes))
    allocate(this%tributary(this%river%n_nodes))
    allocate(this%previous_tributary(this%river%n_nodes))
    allocate(this%nu1(this%river%n_nodes))
    allocate(this%nu2(this%river%n_nodes))
  end subroutine river_router_allocate

  !> \brief Setup river upscaler from fine river and coarse target grid.
  subroutine river_router_init(this, river, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, model_step)
    implicit none
    class(river_router_t), intent(inout) :: this
    type(river_t), pointer, intent(in) :: river !< river definition
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    ! type(inflow_t), intent(in), optional :: inflow_handler !< inflow specifications
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: -1 - threads*8, specials: 0 - all serial, 1 - all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)

    integer(i4) :: model_step_

    this%input_step = optval(input_step, 1_i4)
    model_step_ = optval(model_step, this%input_step)
    call this%setup_grids(river, input_grid)
    call this%allocate() ! allocate arrays based on river and input grid size

    ! TODO: if (present(inflow_handler)) this%inflow_handler = inflow_handler

    ! initial states
    this%input_count = 0_i4
    this%discharge(:) = 0.0_dp
    this%previous_discharge(:) = 0.0_dp
    this%tributary(:) = 0.0_dp
    this%previous_tributary(:) = 0.0_dp

    ! setup muskingum parameters
    call this%setup_muskingum(max_route_step, model_step_)
    ! setup parallelization
    call this%setup_parallelization(root_levels, omp_level_thresh)
  end subroutine river_router_init

  subroutine river_router_from_restart_file(this, path, river, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, read_fluxes, model_step)
    use mo_netcdf, only : NcDataset
    implicit none
    class(river_router_t), intent(inout) :: this
    character(*), intent(in) :: path !< NetCDF file path
    type(river_t), pointer, intent(in) :: river !< river definition
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    ! type(inflow_t), intent(in), optional :: inflow_handler !< inflow specifications
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: threads * 8, 0 to run all in serial, 1 to run all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh
    logical, optional, intent(in) :: read_fluxes !< whether to read fluxes from restart file (default: .true.)
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)
    type(NcDataset) :: nc
    nc = NcDataset(path, "r")
    call this%from_restart_dataset( &
      nc=nc, river=river, input_grid=input_grid, input_step=input_step, max_route_step=max_route_step, &
      root_levels=root_levels, omp_level_thresh=omp_level_thresh, read_fluxes=read_fluxes, model_step=model_step)
    call nc%close()
  end subroutine river_router_from_restart_file

  !> \brief Setup river router from restart file
  subroutine river_router_from_restart_dataset(this, nc, river, input_grid, input_step, max_route_step, root_levels, omp_level_thresh, read_fluxes, model_step)
    use mo_utils, only: locate
    !$ use omp_lib, only: omp_get_num_threads
    implicit none
    class(river_router_t), intent(inout) :: this
    type(NcDataset), intent(in) :: nc !< restart dataset
    type(river_t), pointer, intent(in) :: river !< river definition
    type(grid_t), pointer, intent(in) :: input_grid !< input grid
    integer(i4), intent(in), optional :: input_step !< [h] input time step size (1 by default)
    ! type(inflow_t), intent(in), optional :: inflow_handler !< inflow specifications
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: threads * 8, 0 to run all in serial, 1 to run all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh
    logical, optional, intent(in) :: read_fluxes !< whether to read fluxes from restart file (default: .true.)
    integer(i4), intent(in), optional :: model_step !< [h] time between update calls (input_step by default)

    type(NcVariable) :: var
    logical :: do_read_fluxes
    integer(i4) :: model_step_

    do_read_fluxes = optval(read_fluxes, .true.)

    this%input_step = optval(input_step, 1_i4)
    model_step_ = optval(model_step, this%input_step)
    call this%setup_grids(river, input_grid)
    call this%allocate() ! allocate arrays based on river and input grid size

    ! TODO: if (present(inflow_handler)) this%inflow_handler = inflow_handler

    ! initial states
    this%input_count = 0_i4
    this%acc_runoff(:) = 0.0_dp

    ! Current arrays are scratch buffers overwritten by the next routing calculation.
    this%discharge(:) = 0.0_dp

    ! get previous discharge
    if (do_read_fluxes .and. nc%hasVariable("previous_discharge")) then
      var = nc%getVariable("previous_discharge")
      call var%readInto(this%previous_discharge)
    else
      ! TODO: warn about this
      this%previous_discharge(:) = 0.0_dp
    end if

    this%tributary(:) = 0.0_dp

    ! get previous tributary
    if (do_read_fluxes .and. nc%hasVariable("previous_tributary")) then
      var = nc%getVariable("previous_tributary")
      call var%readInto(this%previous_tributary)
    else
      ! TODO: warn about this
      this%previous_tributary(:) = 0.0_dp
    end if

    ! get muskingum parameter nu1
    if (nc%hasVariable("nu1")) then
      var = nc%getVariable("nu1")
      call var%readInto(this%nu1)
    else
      call warn_message("nu1 not found in restart file, cannot initialize river router from restart without nu1.")
      error stop 1
    end if

    ! get muskingum parameter nu2
    if (nc%hasVariable("nu2")) then
      var = nc%getVariable("nu2")
      call var%readInto(this%nu2)
    else
      call warn_message("nu2 not found in restart file, cannot initialize river router from restart without nu2.")
      error stop 1
    end if

    ! get routing step
    if (nc%hasVariable("route_step")) then
      var = nc%getVariable("route_step")
      call var%readInto(this%routing_substep)
    else
      call warn_message("route_step not found in restart file, cannot initialize river router from restart without route_step.")
      error stop 1
    end if

    if (present(max_route_step)) then
      if (this%routing_substep > max_route_step) then
        call warn_message("Routing substep in restart file (", n2s(this%routing_substep), &
          "s) is larger than specified max_route_step (", n2s(max_route_step), "s).")
        error stop 1
      end if
    end if

    call this%setup_timing(model_step_)

    ! setup parallelization
    call this%setup_parallelization(root_levels, omp_level_thresh)

  end subroutine river_router_from_restart_dataset

  !> \brief Write river router state to restart file
  subroutine river_router_to_restart_file(this, path, append)
    use mo_netcdf, only: NcDataset
    implicit none
    class(river_router_t), intent(in) :: this
    character(*), intent(in) :: path !< NetCDF file path
    logical, optional, intent(in) :: append !< whether NetCDF file should be opened in append mode (default .false.)
    type(NcDataset) :: nc
    character(1) :: fmode
    fmode = "w"
    if ( optval(append, .false.) ) fmode = "a"
    nc = NcDataset(path, fmode)
    call this%to_restart_dataset(nc)
    call nc%close()
  end subroutine river_router_to_restart_file

  !> \brief Write river router state to restart file
  subroutine river_router_to_restart_dataset(this, nc)
    use mo_netcdf, only: NcDataset, NcVariable
    implicit none
    class(river_router_t), intent(in) :: this
    type(NcDataset), intent(inout) :: nc !< restart dataset
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

    ! Write only the previous arrays: the current arrays are scratch buffers that
    ! are overwritten before they are read by the next routing calculation.
    if ( allocated(this%previous_discharge) ) then
      call message("writing previous_discharge to restart_file")
      nc_var = nc%setVariable("previous_discharge", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "previous_discharge")
      call nc_var%setAttribute("units", "m3 s-1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%previous_discharge)
    end if

    ! write previous_tributary
    if ( allocated(this%previous_tributary) ) then
      call message("writing previous_tributary to restart_file")
      nc_var = nc%setVariable("previous_tributary", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "previous_tributary")
      call nc_var%setAttribute("units", "m3 s-1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%previous_tributary)
    end if

    ! write muskingum parameter nu1
    if ( allocated(this%nu1) ) then
      call message("writing nu1 to restart_file")
      nc_var = nc%setVariable("nu1", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "muskingum parameter nu1")
      call nc_var%setAttribute("units", "1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%nu1)
    end if

    ! write muskingum parameter nu2
    if ( allocated(this%nu2) ) then
      call message("writing nu2 to restart_file")
      nc_var = nc%setVariable("nu2", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "muskingum parameter nu2")
      call nc_var%setAttribute("units", "1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%nu2)
    end if

    ! numerical routing substep (used to determine nu1 and nu2)
    nc_var = nc%setVariable("route_step", "f64", dims(:0)) ! scalar
    call nc_var%setAttribute("long_name", "numerical routing substep")
    call nc_var%setAttribute("units", "s")
    call nc_var%setData(this%routing_substep)

  end subroutine river_router_to_restart_dataset

  !> \brief calculate the muskingum parameters nu1 and nu2
  subroutine river_router_setup_muskingum(this, max_route_step, model_step)
    use mo_utils, only: locate
    implicit none
    class(river_router_t), intent(inout) :: this
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    integer(i4), intent(in) :: model_step !< [h] time between router update calls
    real(dp), allocatable :: k(:)
    real(dp) :: xi
    integer(i4) :: step_id

    xi = routing_space_weight ! NOTE: fixed for now

    if (.not.allocated(this%river%link_length)) call error_message("river_router%setup_muskingum: link_length not available")
    if (.not.allocated(this%river%celerity)) call error_message("river_router%setup_muskingum: celerity not available")

    ! wave travel time parameter [s]
    allocate(k(this%river%n_nodes), source=(this%river%link_length / this%river%celerity))

    ! set min wave travel time to min routing step
    step_id = max(1_i4, locate(routing_steps, minval(k, mask=.not.this%river%is_sink)))
    this%routing_substep = routing_steps(step_id)
    if (present(max_route_step)) this%routing_substep = min(this%routing_substep, max_route_step)
    if (.not.any(equal(this%routing_substep, routing_steps))) &
      call error_message("routing substep is invalid: ", n2s(this%routing_substep))

    ! muskingum parameters
    this%nu1 = this%routing_substep / ( k * (1.0_dp - xi) + this%routing_substep / 2.0_dp )
    this%nu2 = 1.0_dp - this%nu1 * k / this%routing_substep

    call this%setup_timing(model_step)
  end subroutine river_router_setup_muskingum

  !> \brief Derive completed routing cadence from the model and numerical routing steps.
  subroutine river_router_setup_timing(this, model_step)
    implicit none
    class(river_router_t), intent(inout) :: this
    integer(i4), intent(in) :: model_step !< [h] time between router update calls
    integer(i4) :: status

    call derive_routing_timing( &
      this%input_step, model_step, this%routing_substep, this%routing_step, this%iterations, this%accumulations, status)
    select case (status)
      case (1_i4)
        call error_message("Runoff input step must be positive.")
      case (2_i4)
        call error_message("Model step must be positive.")
      case (3_i4)
        call error_message("Runoff input step must not be finer than the model step.")
      case (4_i4)
        call error_message("Runoff input step must be divisible by the model step.")
      case (5_i4)
        call error_message("Routing substep must be a positive whole number of seconds.")
      case (6_i4)
        call error_message("Could not derive a whole-hour routing step.")
    end select
    if (mod(this%routing_step, model_step) /= 0_i4) call error_message( &
      "Routing step must be divisible by the model step.")
    if (mod(this%routing_step * HOUR_SECONDS, nint(this%routing_substep, i4)) /= 0_i4) call error_message( &
      "Routing substep must divide the completed routing step.")
  end subroutine river_router_setup_timing

  !> \brief Setup parallelization of river routing based on river level sizes and user settings.
  subroutine river_router_setup_parallelization(this, root_levels, omp_level_thresh)
    !$ use omp_lib, only: omp_get_num_threads
    implicit none
    class(river_router_t), intent(inout) :: this
    logical, intent(in), optional :: root_levels !< order levels as distance from graph roots (default: .false.)
    !> minimum size of river-levels to route in parallel (default: -1 - threads*8, specials: 0 - all serial, 1 - all in parallel)
    integer(i8), optional, intent(in) :: omp_level_thresh
    integer(i8) :: omp_lvl_thr
    integer(i8) :: i
    integer(i8), pointer :: level_size(:)

    if (.not.allocated(this%river%order%id)) call this%river%calc_order(root_levels)

    this%last_parallel_level = 0_i8 ! default is to run all levels in serial
    omp_lvl_thr = optval(omp_level_thresh, -1_i8) ! -1 means default

    !$omp parallel
    !$ this%omp_level_thresh = int(omp_get_num_threads() * 8, kind=i8)
    !$ if (omp_lvl_thr >= 0_i8) this%omp_level_thresh = omp_lvl_thr
    !$omp end parallel

    ! determine last level to run in parallel
    if (this%river%order%n_levels == 1_i8) then
      if (this%omp_level_thresh > 0_i8) this%last_parallel_level = 1_i8
      return
    end if
    level_size => this%river%order%level_size
    ! check if levels are sorted in size (only true for leave based ordering)
    if (all(level_size(:this%river%order%n_levels-1_i8) >= level_size(2_i8:))) then
      if (optval(root_levels, .false.)) then
        call warn_message("River levels were requested as distance from graph roots, but level sizes are sorted in descending order.")
        call warn_message("This is unexpected, maybe your restart file has inconsistent ordering compared to namelist settings.")
      end if
      if (this%omp_level_thresh > 0_i8) then
        do i = 1_i8, this%river%order%n_levels
          this%last_parallel_level = i
          if (this%river%order%level_size(i) < this%omp_level_thresh) exit
        end do
      end if
    else
      ! root based levels are not sorted in size, so run all in parallel if wanted
      if (this%omp_level_thresh > 0_i8) this%last_parallel_level = this%river%order%n_levels
    end if
  end subroutine river_router_setup_parallelization

  !> \brief Update routing for one time step
  subroutine river_router_update(this, input_runoff, discharge, completed)
    implicit none
    class(river_router_t), intent(inout) :: this
    real(dp), dimension(this%input_grid%ncells), intent(in) :: input_runoff
    real(dp), dimension(this%river%n_nodes), intent(inout) :: discharge
    logical, intent(out) :: completed !< whether a new routing-step mean was produced
    integer(i4) :: i
    integer(i8) :: n, c
    completed = .false.
    ! accumulate input runoff
    this%input_count = this%input_count + 1_i4
    if (this%input_count == 1_i4) then
      this%acc_runoff = input_runoff ! first accumulation as runoff reset
    else
      this%acc_runoff = this%acc_runoff + input_runoff ! accumulate
    end if
    if (this%input_count < this%accumulations) return ! not yet ready to route
    ! reset input counter
    this%input_count = 0_i4
    ! Average repeated runoff amounts across model updates before converting the represented input support to a rate.
    if (this%accumulations > 1_i4) this%acc_runoff = this%acc_runoff / this%accumulations
    ! distribute runoff in case of SCC
    if (this%river%scc) then
      this%scaled_runoff = this%scale_runoff(this%acc_runoff)
      ! Keep SIMD but avoid the combined OpenMP 4.x form that NAG rejects here.
      !$omp parallel default(none) shared(this) private(c, n)
      !$omp do simd schedule(static)
      do n = 1_i8, this%river%n_nodes
        c = this%river%node_cell(n)
        this%runoff(n) = this%scaled_runoff(c) * this%river%area_fraction(n)
      end do
      !$omp end do simd
      !$omp end parallel
    else
      this%runoff = this%scale_runoff(this%acc_runoff)
    end if
    ! route
    do i = 1_i4, this%iterations
      call this%route(this%discharge, this%tributary)
      if (i == 1_i4) then
        discharge = this%discharge
      else
        discharge = discharge + this%discharge
      end if
      ! store discharge and tributary for next time-step
      call this%time_shift()
    end do
    ! average discharge flux over output time step iterations
    if (this%iterations > 1_i4) discharge = discharge / this%iterations
    completed = .true.
  end subroutine river_router_update

  !> \brief Setup river upscaler from fine river and coarse target grid.
  function river_router_scale_runoff(this, input_runoff) result(scaled_runoff_flux)
    implicit none
    class(river_router_t), intent(inout) :: this
    real(dp), dimension(this%input_grid%ncells), intent(in) :: input_runoff
    real(dp), dimension(this%river%grid%ncells) :: scaled_runoff_flux
    ! rescaled result in [liter] = [mm m2]
    if (this%scaler%scaling_mode == down_scaling) then
      call this%scaler%downscale_nearest(input_runoff, scaled_runoff_flux)
      scaled_runoff_flux = scaled_runoff_flux * this%river%grid%cell_area
    else
      call this%scaler%execute(input_runoff * this%input_grid%cell_area, scaled_runoff_flux)
    end if
    ! from [liter TS-1] to [m3 s-1]
    scaled_runoff_flux = scaled_runoff_flux * (liter_to_m3 / (this%input_step * HOUR_SECONDS))
  end function river_router_scale_runoff

  !> \brief Execute routing for a single routing time step.
  subroutine river_router_route(this, discharge, tributary)
    use mo_dag, only: node
    implicit none
    class(river_router_t), intent(in) :: this
    real(dp), dimension(this%river%n_nodes), intent(out) :: discharge, tributary
    integer(i8) :: i, j, n_levels
    logical, pointer, contiguous :: is_sink(:)
    integer(i8), pointer, contiguous :: id(:), level_start(:), level_end(:)

    ! pointers for speeding up dereferencing attributes (river is a pointer, so this works)
    id => this%river%order%id
    level_start => this%river%order%level_start
    level_end => this%river%order%level_end
    is_sink => this%river%is_sink
    n_levels = this%river%order%n_levels

    ! parallel routing on levels with size above threshold
    if (this%last_parallel_level > 0_i8) then
      !$omp parallel default(shared) private(i)
      do i = 1_i8, this%last_parallel_level
        !$omp do simd schedule(static)
        do j = level_start(i), level_end(i)
          call process_node(j)
        end do
        !$omp end do simd
      end do
      !$omp end parallel
    end if

    ! serial routing of nodes from remaining levels
    if (this%last_parallel_level == n_levels) return ! no node left
    do j = level_start(this%last_parallel_level + 1_i8), level_end(n_levels)
      call process_node(j)
    end do

  contains
    subroutine process_node(ni)
      integer(i8), intent(in) :: ni
      integer(i8) :: n
      integer(i4) :: n_edges, m
      integer(i8), pointer, contiguous :: edges(:)
      real(dp) :: inflow, prev_inflow, prev_routed
      n = id(ni)
      ! add runoff to inflow
      inflow = this%runoff(n)
      ! add upstream routed flux to inflow
      call this%river%src_view(n, edges)
      n_edges = size(edges)
      !$omp simd reduction(+:inflow)
      do m = 1_i4, n_edges
        inflow = inflow + tributary(edges(m))
      end do
      ! TODO: add/replace inflow with inflow handler at inflow gauges
      ! discharge is the accumulated inflow at current node
      discharge(n) = inflow
      if (is_sink(n)) then
        ! no routing at sinks
        tributary(n) = 0.0_dp
      else
        ! muskingum schema for routed flux
        prev_routed = this%previous_tributary(n)
        prev_inflow = this%previous_discharge(n)
        tributary(n) = prev_routed + this%nu1(n) * (prev_inflow - prev_routed) + this%nu2(n) * (inflow - prev_inflow)
      end if
    end subroutine process_node
  end subroutine river_router_route

  !> \brief Move current discharge and tributary to previous arrays.
  subroutine river_router_time_shift(this)
    implicit none
    class(river_router_t), intent(inout) :: this
    real(dp), allocatable :: tmp(:)
    ! swap discharge and previous_discharge
    call move_alloc(this%previous_discharge, tmp)
    call move_alloc(this%discharge, this%previous_discharge)
    call move_alloc(tmp, this%discharge)
    ! swap tributary and previous_tributary
    call move_alloc(this%previous_tributary, tmp)
    call move_alloc(this%tributary, this%previous_tributary)
    call move_alloc(tmp, this%tributary)
  end subroutine river_router_time_shift
end module mo_river_router
