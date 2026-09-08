!> \file    mo_lake_topology.f90
!> \copydoc mo_lake_topology

!> \brief   Labeled lake-footprint and restart helpers.
!> \details Shared topology validation, area derivation, and restart persistence for lake processes.
!> \version 0.1
!> \changelog
!! - Pallav Shrestha (2018-2023): lake-aware topology concepts.
!! - Sebastian Mueller (2026): exchange-side restart topology helpers.
!> \authors Pallav Shrestha, Sebastian Mueller
!> \date    September 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
#include "logging.h"
module mo_lake_topology
  use mo_logging
  use mo_kind, only: i4, i8, dp
  use mo_grid, only: grid_t
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  implicit none
  private
  public :: lake_derive_area, lake_read_topology, lake_write_topology
contains
  !> \brief Derive lake areas from the labeled packed level-0 footprint.
  subroutine lake_derive_area(grid, lake_map, lake_ids, lake_area)
    type(grid_t), intent(in) :: grid
    integer(i8), intent(in) :: lake_map(:), lake_ids(:)
    real(dp), allocatable, intent(out) :: lake_area(:)
    integer(i8) :: i, j
    if (size(lake_map, kind=i8) /= grid%ncells) then
      log_fatal(*) "Lake topology: map does not match grid."
      error stop 1
    end if
    allocate(lake_area(size(lake_ids)), source=0.0_dp)
    do i = 1_i8, size(lake_map, kind=i8)
      j = findloc(lake_ids, lake_map(i), dim=1, kind=i8)
      if (j < 1_i8) then
        log_fatal(*) "Lake topology: map contains an unknown stable lake ID."
        error stop 1
      end if
      lake_area(j) = lake_area(j) + grid%cell_area(i)
    end do
    if (any(lake_area <= 0.0_dp)) then
      log_fatal(*) "Lake topology: every lake must have positive area."
      error stop 1
    end if
  end subroutine lake_derive_area

  !> \brief Write the packed lake grid and stable lake IDs to a restart dataset.
  subroutine lake_write_topology(nc, grid, lake_map)
    type(NcDataset), intent(inout) :: nc
    type(grid_t), intent(in) :: grid
    integer(i8), intent(in) :: lake_map(:)
    type(NcDimension) :: lake_cell_dim
    type(NcVariable) :: var
    call grid%to_restart(nc)
    lake_cell_dim = nc%setDimension("lake_cell", int(size(lake_map), i4))
    var = nc%setVariable("lake_map", "i64", [lake_cell_dim])
    call var%setData(lake_map)
  end subroutine lake_write_topology

  !> \brief Restore and validate the packed lake topology from a restart dataset.
  subroutine lake_read_topology(nc, grid, lake_map)
    type(NcDataset), intent(inout) :: nc
    type(grid_t), intent(out) :: grid
    integer(i8), allocatable, intent(out) :: lake_map(:)
    type(NcVariable) :: var
    if (.not.nc%hasVariable("lake_map")) then
      log_fatal(*) "Lake restart is missing level-0 topology."
      error stop 1
    end if
    call grid%from_restart(nc)
    var = nc%getVariable("lake_map")
    call var%getData(lake_map)
    if (size(lake_map, kind=i8) /= grid%ncells .or. any(lake_map <= 0_i8)) then
      log_fatal(*) "Lake restart contains invalid level-0 topology."
      error stop 1
    end if
  end subroutine lake_read_topology
end module mo_lake_topology
