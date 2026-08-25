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

  use mo_kind, only: dp, i2, i4, i8
  use mo_message, only: error_message
  use mo_grid, only: grid_t, check_factor, coarse_ij, id_bounds, bottom_up
  use mo_river, only: river_t, d8_E, d8_SE, d8_S, d8_SW, d8_W, d8_NW, d8_N, d8_NE

  implicit none
  private
  public :: fdir_upscaling

contains

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
