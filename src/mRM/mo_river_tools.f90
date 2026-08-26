!> \file    mo_river_tools.f90
!> \copydoc mo_river_tools

!> \brief   River related tools.
!> \details This module contains tools to deal with river related data.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_tools

  use mo_kind, only: i1, i2, i4, i8
  use mo_message, only: error_message
  use mo_grid, only: grid_t, check_factor, id_bounds, bottom_up

  implicit none
  private
  public :: fdir_upscaling
  public :: unique_ids

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
