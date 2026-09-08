!> \file    mo_lake_forcing.f90
!> \copydoc mo_lake_forcing

!> \brief   Sparse level-2 support for lake-point meteorological forcings.
!> \details Builds count-weighted lake supports from labeled level-0 lake cells.
!> \version 0.1
!> \changelog
!! - Pallav Shrestha (2018-2023): lake-aware routing and coupling concepts.
!! - Sebastian Mueller (2026): sparse lake meteorological forcing support.
!> \authors Pallav Shrestha, Sebastian Mueller
!> \date    September 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
#include "logging.h"
module mo_lake_forcing
  use mo_logging
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t, down_nearest
  use mo_kind, only: i4, i8, dp
  use mo_orderpack, only: sort
  implicit none
  private
  public :: lake_support_t, lake_support_build, lake_unique_counts
  type :: lake_support_t
    integer(i8), allocatable :: l2_ids(:)
    real(dp), allocatable :: weights(:)
  contains
    procedure :: aggregate => lake_support_aggregate
  end type lake_support_t
contains

  !> \brief Return sorted unique values and occurrence counts under an optional mask.
  subroutine lake_unique_counts(values, unique, counts, mask)
    integer(i4), intent(in) :: values(:)
    integer(i8), allocatable, intent(out) :: unique(:), counts(:)
    logical, optional, intent(in) :: mask(:)
    integer(i4), allocatable :: selected(:)
    integer(i8) :: i, n_unique

    if (present(mask)) then
      if (size(mask, kind=i8) /= size(values, kind=i8)) then
        log_fatal(*) "Lake forcing: unique-value mask has an incompatible size."
        error stop 1
      end if
      selected = pack(values, mask)
    else
      selected = values
    end if
    if (size(selected, kind=i8) == 0_i8) then
      allocate(unique(0), counts(0))
      return
    end if
    call sort(selected)
    n_unique = 1_i8 + count(selected(2:) /= selected(:size(selected) - 1), kind=i8)
    allocate(unique(n_unique), counts(n_unique))
    n_unique = 1_i8
    unique(1) = int(selected(1), i8)
    counts(1) = 1_i8
    do i = 2_i8, size(selected, kind=i8)
      if (selected(i) == selected(i - 1_i8)) then
        counts(n_unique) = counts(n_unique) + 1_i8
      else
        n_unique = n_unique + 1_i8
        unique(n_unique) = int(selected(i), i8)
        counts(n_unique) = 1_i8
      end if
    end do
  end subroutine lake_unique_counts

  !> \brief Build count-weighted sparse level-2 support for every configured lake.
  subroutine lake_support_build(level2, level0_lake, lake_map, lake_ids, support)
    type(grid_t), target, intent(in) :: level2, level0_lake
    integer(i8), intent(in) :: lake_map(:), lake_ids(:)
    type(lake_support_t), allocatable, intent(out) :: support(:)
    type(scaler_t) :: scaler
    integer(i4), allocatable :: l2_ids(:)
    integer(i8), allocatable :: counts(:)
    integer(i8) :: lake, n, total
    if (size(lake_map, kind=i8) /= level0_lake%ncells) then
      log_fatal(*) "Lake forcing: map does not match level-0 lake grid."
      error stop 1
    end if
    call scaler%init(level2, level0_lake, downscaling_operator=down_nearest)
    if (any(scaler%id_map < 1_i8) .or. any(scaler%id_map > int(huge(0_i4), i8))) then
      log_fatal(*) "Lake forcing: lake cells are not covered by valid level-2 cells."
      error stop 1
    end if
    l2_ids = int(scaler%id_map, i4)
    allocate(support(size(lake_ids)))
    total = 0_i8
    !$omp parallel do default(shared) private(lake, counts, n) reduction(+:total) schedule(dynamic)
    do lake = 1_i8, size(lake_ids, kind=i8)
      call lake_unique_counts(l2_ids, support(lake)%l2_ids, counts, mask=lake_map == lake_ids(lake))
      n = sum(counts)
      total = total + n
      if (n == 0_i8) cycle
      allocate(support(lake)%weights(size(counts)))
      support(lake)%weights = real(counts, dp) / real(n, dp)
    end do
    !$omp end parallel do
    do lake = 1_i8, size(support, kind=i8)
      if (.not.allocated(support(lake)%l2_ids)) then
        log_fatal(*) "Lake forcing: lake map has unknown, duplicate, or empty lake IDs."
        error stop 1
      end if
    end do
    if (total /= level0_lake%ncells) then
      log_fatal(*) "Lake forcing: lake map has unknown, duplicate, or empty lake IDs."
      error stop 1
    end if
  end subroutine lake_support_build

  !> \brief Aggregate one packed level-2 field over this lake support.
  function lake_support_aggregate(self, field) result(result)
    class(lake_support_t), intent(in) :: self
    real(dp), intent(in) :: field(:)
    real(dp) :: result
    result = sum(field(self%l2_ids) * self%weights)
  end function lake_support_aggregate
end module mo_lake_forcing
