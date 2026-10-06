!> \file    mo_river_tools.f90
!> \copydoc mo_river_tools

!> \brief   River related tools.
!> \details This module contains river data tools and shared routing state types.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_tools

  use mo_kind, only: i1, i2, i4, i8, dp
  use mo_message, only: error_message
  use mo_utils, only: equal, locate
  use mo_string_utils, only: n2s => num2str
  use mo_datetime, only: HOUR_SECONDS
  use mo_grid_scaler, only: scaler_t, down_scaling
  use mo_grid, only: grid_t, check_factor, id_bounds, bottom_up

  implicit none
  private
  public :: fdir_upscaling
  public :: unique_ids, river_depth_to_flux

  !> valid routing time steps [s] (integer minute fractions of hour or integer hour fractions of day)
  real(dp), public, dimension(19), parameter :: routing_steps = &
    [   60.0_dp,   120.0_dp,   180.0_dp,   240.0_dp,   300.0_dp,   360.0_dp,            &
       600.0_dp,   720.0_dp,   900.0_dp,  1200.0_dp,  1800.0_dp,  3600.0_dp,            &
      7200.0_dp, 10800.0_dp, 14400.0_dp, 21600.0_dp, 28800.0_dp, 43200.0_dp, 86400.0_dp ]
  !> space weighting of routing is set to 0, since this parameter has no effect on the routing results, see Thober et al. 2017
  real(dp), public, parameter :: routing_space_weight = 0.0_dp

  !> \class river_muskingum_state_t
  !> \brief Router-owned coefficients and completed routing cadence.
  type, public :: river_muskingum_state_t
    real(dp), allocatable :: nu1(:) !< First Muskingum coefficient per node.
    real(dp), allocatable :: nu2(:) !< Second Muskingum coefficient per node.
    integer(i4) :: routing_step = 0_i4 !< Completed routing interval [h].
    integer(i4) :: iterations = 0_i4 !< Numerical substeps per completed interval.
    integer(i4) :: accumulations = 0_i4 !< Model updates per completed interval.
    real(dp) :: routing_substep = 0.0_dp !< Numerical substep [s].
  contains
    procedure :: init => river_muskingum_state_init
    procedure :: derive_timing => river_muskingum_state_derive_timing
    procedure :: setup_timing => river_muskingum_state_setup_timing
    procedure :: allocate => river_muskingum_state_allocate
    procedure :: destroy => river_muskingum_state_destroy
  end type river_muskingum_state_t

  !> \name D8 direction values
  !> \brief Constants describing the D8 routing direction of a cell.
  !!@{
  integer(i2), public, parameter :: sink = 0_i2 !< sink
  integer(i2), public, parameter :: d8_E = 1_i2 !< east
  integer(i2), public, parameter :: d8_SE = 2_i2 !< south-east
  integer(i2), public, parameter :: d8_S = 4_i2 !< south
  integer(i2), public, parameter :: d8_SW = 8_i2 !< south-west
  integer(i2), public, parameter :: d8_W = 16_i2 !< west
  integer(i2), public, parameter :: d8_NW = 32_i2 !< north-west
  integer(i2), public, parameter :: d8_N = 64_i2 !< north
  integer(i2), public, parameter :: d8_NE = 128_i2 !< north-east
  integer(i2), public, dimension(8), parameter :: d8_all = [d8_E, d8_SE, d8_S, d8_SW, d8_W, d8_NW, d8_N, d8_NE]
  integer(i2), public, dimension(8), parameter :: d8_back = [d8_W, d8_NW, d8_N, d8_NE, d8_E, d8_SE, d8_S, d8_SW]
  integer(i2), public, dimension(3), parameter :: d8_leave_E = [d8_NE, d8_E, d8_SE]
  integer(i2), public, dimension(3), parameter :: d8_leave_S = [d8_SE, d8_S, d8_SW]
  integer(i2), public, dimension(3), parameter :: d8_leave_W = [d8_SW, d8_W, d8_NW]
  integer(i2), public, dimension(3), parameter :: d8_leave_N = [d8_NW, d8_N, d8_NE]
  integer(i2), public, dimension(5), parameter :: d8_leave_SE = [d8_NE, d8_E, d8_SE, d8_S, d8_SW]
  integer(i2), public, dimension(5), parameter :: d8_leave_SW = [d8_SE, d8_S, d8_SW, d8_W, d8_NW]
  integer(i2), public, dimension(5), parameter :: d8_leave_NE = [d8_NW, d8_N, d8_NE, d8_E, d8_SE]
  integer(i2), public, dimension(5), parameter :: d8_leave_NW = [d8_SW, d8_W, d8_NW, d8_N, d8_NE]
  !!@}

  !> \name LDD direction values
  !> \brief Constants describing the local drain direction routing direction of a cell.
  !!@{
  integer(i1), public, parameter :: o = 0_i1 !< no-data shortcut
  integer(i1), public, parameter :: ldd_sink = 5_i1 !< sink
  integer(i1), public, parameter :: ldd_E = 6_i1 !< east
  integer(i1), public, parameter :: ldd_SE = 3_i1 !< south-east
  integer(i1), public, parameter :: ldd_S = 2_i1 !< south
  integer(i1), public, parameter :: ldd_SW = 1_i1 !< south-west
  integer(i1), public, parameter :: ldd_W = 4_i1 !< west
  integer(i1), public, parameter :: ldd_NW = 7_i1 !< north-west
  integer(i1), public, parameter :: ldd_N = 8_i1 !< north
  integer(i1), public, parameter :: ldd_NE = 9_i1 !< north-east
  integer(i1), public, dimension(8), parameter :: ldd_all = [ldd_E, ldd_SE, ldd_S, ldd_SW, ldd_W, ldd_NW, ldd_N, ldd_NE]
  integer(i1), public, dimension(8), parameter :: ldd_back = [ldd_W, ldd_NW, ldd_N, ldd_NE, ldd_E, ldd_SE, ldd_S, ldd_SW]
  integer(i4), public, dimension(9), parameter :: ldd_dx = [-1_i4, 0_i4, 1_i4, -1_i4, 0_i4, 1_i4, -1_i4, 0_i4, 1_i4]
  integer(i4), public, dimension(9), parameter :: ldd_dy = [-1_i4, -1_i4, -1_i4, 0_i4, 0_i4, 0_i4, 1_i4, 1_i4, 1_i4]
  integer(i2), public, dimension(9), parameter :: map_ldd_to_d8 = [d8_SW, d8_S, d8_SE, d8_W, sink, d8_E, d8_NW, d8_N, d8_NE]
  integer(i1), public, dimension(128), parameter :: map_d8_to_ldd = [ &
    ldd_E, &
    ldd_SE, o, &
    ldd_S,  o,o,o, &
    ldd_SW, o,o,o,o,o,o,o, &
    ldd_W,  o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_NW, o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_N,  o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
            o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_NE]
  !!@}

  interface unique_ids
    module procedure unique_ids_i4
    module procedure unique_ids_i8
  end interface unique_ids

contains

  !> \brief Allocate Muskingum coefficients and derive cadence from link travel times.
  subroutine river_muskingum_state_init(self, link_length, is_sink, celerity, input_step, model_step, max_route_step)
    class(river_muskingum_state_t), intent(inout), target :: self !< Coefficients and routing cadence to initialize.
    real(dp), intent(in) :: link_length(:) !< Routing link lengths [m].
    logical, intent(in) :: is_sink(:) !< Nodes without a downstream link.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    real(dp), intent(in) :: celerity(:) !< celerity of the link starting at each river node
    real(dp), optional, intent(in) :: max_route_step !< [s] maximum routing time step (default: 86400.0)
    integer(i4), intent(in) :: model_step !< [h] time between router update calls
    real(dp), allocatable :: k(:)
    real(dp) :: xi
    integer(i4) :: step_id
    integer(i8) :: n_nodes, n

    xi = routing_space_weight ! NOTE: fixed for now

    n_nodes = size(link_length, kind=i8)
    if (size(celerity, kind=i8) /= n_nodes .or. size(is_sink, kind=i8) /= n_nodes) &
      call error_message("river_muskingum_state: node array sizes differ")
    call self%allocate(n_nodes)

    ! wave travel time parameter [s]
    allocate(k(n_nodes))
    !$omp parallel do default(shared) schedule(static)
    do n = 1_i8, n_nodes
      k(n) = link_length(n) / celerity(n)
    end do
    !$omp end parallel do

    ! set min wave travel time to min routing step
    step_id = max(1_i4, locate(routing_steps, minval(k, mask=.not.is_sink)))
    self%routing_substep = routing_steps(step_id)
    if (present(max_route_step)) self%routing_substep = min(self%routing_substep, max_route_step)
    if (.not.any(equal(self%routing_substep, routing_steps))) &
      call error_message("routing substep is invalid: ", n2s(self%routing_substep))

    ! muskingum parameters
    !$omp parallel do default(shared) schedule(static)
    do n = 1_i8, n_nodes
      self%nu1(n) = self%routing_substep / (k(n) * (1.0_dp - xi) + self%routing_substep / 2.0_dp)
      self%nu2(n) = 1.0_dp - self%nu1(n) * k(n) / self%routing_substep
    end do
    !$omp end parallel do

    call self%setup_timing(input_step, model_step)
  end subroutine river_muskingum_state_init

  !> \brief Derive completed routing cadence from the model and numerical routing steps.
  subroutine river_muskingum_state_setup_timing(self, input_step, model_step)
    class(river_muskingum_state_t), intent(inout), target :: self !< Routing cadence to derive and validate.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    integer(i4), intent(in) :: model_step !< [h] time between router update calls
    integer(i4) :: status

    call self%derive_timing(input_step, model_step, status)
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
    if (mod(self%routing_step, model_step) /= 0_i4) call error_message( &
      "Routing step must be divisible by the model step.")
    if (mod(self%routing_step * HOUR_SECONDS, nint(self%routing_substep, i4)) /= 0_i4) call error_message( &
      "Routing substep must divide the completed routing step.")
  end subroutine river_muskingum_state_setup_timing

  !> \brief Derive the completed routing step and its update/iteration counts.
  subroutine river_muskingum_state_derive_timing(self, input_step, model_step, status)
    class(river_muskingum_state_t), intent(inout), target :: self !< Numerical substep and derived cadence.
    integer(i4), intent(in) :: input_step !< [h] runoff support interval
    integer(i4), intent(in) :: model_step !< [h] model update interval
    integer(i4), intent(out) :: status !< 0 success; positive values identify an invalid timing input
    integer(i8) :: a, b, divisor
    integer(i8) :: model_seconds, routing_seconds, substep_seconds

    self%routing_step = 0_i4
    self%iterations = 0_i4
    self%accumulations = 0_i4
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

    substep_seconds = nint(self%routing_substep, i8)
    if (substep_seconds < 1_i8 .or. .not.equal(self%routing_substep, real(substep_seconds, dp))) then
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
      routing_seconds / int(HOUR_SECONDS, i8) > int(huge(self%routing_step), i8)) then
      status = 6_i4
      return
    end if

    self%routing_step = int(routing_seconds / int(HOUR_SECONDS, i8), i4)
    self%iterations = int(routing_seconds / substep_seconds, i4)
    self%accumulations = self%routing_step / model_step
  end subroutine river_muskingum_state_derive_timing

  !> \brief Map L1 runoff depth or temperature-depth to L3 volume rate.
  subroutine river_depth_to_flux(scaler, input_depth, scaled_flux, input_step, scale_work, scale_grid)
    type(scaler_t), intent(inout) :: scaler !< Shared L1-to-L3 remapping cache.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    real(dp), intent(in) :: input_depth(:) !< L1 runoff depth or temperature-depth [mm or K mm].
    real(dp), intent(out) :: scaled_flux(:) !< L3 rate on runoff support [m3 s-1 or K m3 s-1].
    real(dp), optional, intent(inout) :: scale_work(:) !< L1 volume work buffer when upscaling.
    real(dp), optional, intent(inout) :: scale_grid(:, :) !< Reusable unpacked L1 volume workspace.
    integer(i8) :: c
    real(dp) :: rate_factor
    if (size(input_depth) /= scaler%source_grid%ncells .or. &
        size(scaled_flux) /= scaler%target_grid%ncells) &
      call error_message("river_router: depth-to-flux support size mismatch")
    rate_factor = 0.001_dp / (input_step * HOUR_SECONDS)
    if (scaler%scaling_mode == down_scaling) then
      call scaler%downscale_nearest(input_depth, scaled_flux)
      !$omp parallel do default(shared) schedule(static)
      do c = 1_i8, size(scaled_flux, kind=i8)
        scaled_flux(c) = (scaled_flux(c) * scaler%target_grid%cell_area(c)) * rate_factor
      end do
      !$omp end parallel do
    else
      if (.not.present(scale_work) .or. .not.present(scale_grid)) &
        call error_message("river_router: missing depth-to-flux work buffer")
      !$omp parallel do default(shared) schedule(static)
      do c = 1_i8, size(input_depth, kind=i8)
        scale_work(c) = input_depth(c) * scaler%source_grid%cell_area(c)
      end do
      !$omp end parallel do
      call scaler%source_grid%unpack_into(scale_work, scale_grid)
      call scaler%execute(scale_grid, scaled_flux)
      !$omp parallel do default(shared) schedule(static)
      do c = 1_i8, size(scaled_flux, kind=i8)
        scaled_flux(c) = scaled_flux(c) * rate_factor
      end do
      !$omp end parallel do
    end if
  end subroutine river_depth_to_flux

  !> \brief Allocate router-owned Muskingum coefficients.
  subroutine river_muskingum_state_allocate(self, n_nodes)
    class(river_muskingum_state_t), intent(inout), target :: self !< Coefficients to initialize.
    integer(i8), intent(in) :: n_nodes !< Number of routing nodes.
    integer(i8) :: n
    call self%destroy()
    allocate(self%nu1(n_nodes), self%nu2(n_nodes))
    !$omp parallel do default(shared) schedule(static)
    do n = 1_i8, n_nodes
      self%nu1(n) = 0.0_dp
      self%nu2(n) = 0.0_dp
    end do
    !$omp end parallel do
  end subroutine river_muskingum_state_allocate

  !> \brief Release router-owned coefficients and cadence.
  subroutine river_muskingum_state_destroy(self)
    class(river_muskingum_state_t), intent(inout), target :: self !< Coefficients to release.
    if (allocated(self%nu1)) deallocate(self%nu1)
    if (allocated(self%nu2)) deallocate(self%nu2)
    self%routing_step = 0_i4
    self%routing_substep = 0.0_dp
    self%iterations = 0_i4
    self%accumulations = 0_i4
  end subroutine river_muskingum_state_destroy


  !> \brief Check whether all 32-bit integer IDs are unique.
  logical function unique_ids_i4(ids) result(unique)
    integer(i4), intent(in) :: ids(:)
    integer(i8) :: i, j

    unique = .true.
    !$omp parallel do default(shared) private(j) reduction(.and.:unique) schedule(guided)
    do i = 1_i8, size(ids, kind=i8) - 1_i8
      do j = i + 1_i8, size(ids, kind=i8)
        unique = unique .and. ids(i) /= ids(j)
      end do
    end do
    !$omp end parallel do
  end function unique_ids_i4

  !> \brief Check whether all 64-bit integer IDs are unique.
  logical function unique_ids_i8(ids) result(unique)
    integer(i8), intent(in) :: ids(:)
    integer(i8) :: i, j

    unique = .true.
    !$omp parallel do default(shared) private(j) reduction(.and.:unique) schedule(guided)
    do i = 1_i8, size(ids, kind=i8) - 1_i8
      do j = i + 1_i8, size(ids, kind=i8)
        unique = unique .and. ids(i) /= ids(j)
      end do
    end do
    !$omp end parallel do
  end function unique_ids_i8

  !> \brief Upscale flow direction from fine to coarse grid
  subroutine fdir_upscaling(fine_grid, fine_fdir, fine_facc, coarse_grid, fdir)
    type(grid_t), intent(in) :: fine_grid
    integer(i2), intent(in) :: fine_fdir(:)
    integer(i4), intent(in) :: fine_facc(:)
    type(grid_t), intent(in) :: coarse_grid
    integer(i2), intent(out) :: fdir(:)

    logical :: is_bottom_up
    integer(i2), allocatable :: in_fdir(:,:)
    integer(i4), allocatable :: in_facc(:,:)
    integer(i8) :: i
    integer(i4) :: factor, loc(2), ix, iy, yl, yu, xl, xu, yn, ys
    integer(i2) :: dir

    call check_factor(fine_grid%cellsize, coarse_grid%cellsize, factor=factor)

    allocate(in_fdir(fine_grid%nx, fine_grid%ny))
    allocate(in_facc(fine_grid%nx, fine_grid%ny))

    call fine_grid%unpack_into(fine_fdir, in_fdir)
    call fine_grid%unpack_into(fine_facc, in_facc)
    is_bottom_up = (fine_grid%y_direction == bottom_up)

    ! determine coarse fdir by finding the fine cell with max facc in each coarse cell
    !$omp parallel do default(shared) private(yl,yu,xl,xu,yn,ys,ix,iy,loc,dir) schedule(static)
    do i = 1_i8, coarse_grid%ncells
      call id_bounds(factor, coarse_grid%cell_ij(i,1_i8), coarse_grid%cell_ij(i,2_i8), &
        coarse_grid%y_direction, coarse_grid%ny, fine_grid%y_direction, fine_grid%nx, fine_grid%ny, &
        xl, xu, yl, yu)
      if (is_bottom_up) then
        yn = yu ! north y-bound is up
        ys = yl ! south y-bound is down
      else ! top_down
        yn = yl ! north y-bound is down
        ys = yu ! south y-bound is up
      end if
      loc = maxloc(in_facc(xl:xu,yl:yu), mask=fine_grid%mask(xl:xu,yl:yu))
      ix = xl + loc(1) - 1_i4
      iy = yl + loc(2) - 1_i4
      dir = in_fdir(ix,iy)
      if (dir == 0_i2) then ! coarse sink
        fdir(i) = 0_i2
      else if (ix==xu.and.iy==yn) then ! NE
        select case(dir)
          case(d8_NW, d8_N)
            fdir(i) = d8_N
          case(d8_NE)
            fdir(i) = d8_NE
          case(d8_E, d8_SE)
            fdir(i) = d8_E
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (ix==xu.and.iy==ys) then ! SE
        select case(dir)
          case(d8_NE, d8_E)
            fdir(i) = d8_E
          case(d8_SE)
            fdir(i) = d8_SE
          case(d8_S, d8_SW)
            fdir(i) = d8_S
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (ix==xl.and.iy==ys) then ! SW
        select case(dir)
          case(d8_SE, d8_S)
            fdir(i) = d8_S
          case(d8_SW)
            fdir(i) = d8_SW
          case(d8_W, d8_NW)
            fdir(i) = d8_W
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (ix==xl.and.iy==yn) then ! NW
        select case(dir)
          case(d8_SW, d8_W)
            fdir(i) = d8_W
          case(d8_NW)
            fdir(i) = d8_NW
          case(d8_N, d8_NE)
            fdir(i) = d8_N
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (ix==xu) then ! E
        select case (dir)
          case(d8_NE, d8_E, d8_SE)
            fdir(i) = d8_E
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (ix==xl) then ! W
        select case (dir)
          case(d8_SW, d8_W, d8_NW)
            fdir(i) = d8_W
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (iy==yn) then ! N
        select case (dir)
          case(d8_NW, d8_N, d8_NE)
            fdir(i) = d8_N
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else if (iy==ys) then ! S
        select case (dir)
          case(d8_SW, d8_S, d8_SE)
            fdir(i) = d8_S
          case default
            call error_message("fdir_upscaling: invalid fine flow direction for cell ")
        end select
      else
        call error_message("fdir_upscaling: could not determine coarse flow direction for cell ")
      end if
    end do
    !$omp end parallel do

    deallocate(in_fdir, in_facc)

  end subroutine fdir_upscaling

end module mo_river_tools
