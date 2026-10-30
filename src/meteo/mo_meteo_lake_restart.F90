!> \file mo_meteo_lake_restart.F90
!> \copydoc mo_meteo_lake_restart

!> \brief Versioned, owned restart definition for sparse lake meteorological support.
!> \version 0.1
!> \changelog
!! - 2026: compact lake support persistence and compatibility validation.
!> \date October 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_meteo
#include "logging.h"
module mo_meteo_lake_restart
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_kind, only: dp, i1, i4, i8
  use mo_grid, only: grid_t
  use mo_lake_forcing, only: lake_support_t, lake_grids_match
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  implicit none
  private
  public :: lake_forcing_restart_t, lake_grids_match, lake_supports_match

  !> \class lake_forcing_restart_t
  !> \brief Owned spatial definition and active daily-to-hourly disaggregation settings.
  type :: lake_forcing_restart_t
    type(grid_t) :: grid !< native forcing geometry, mask and packing
    integer(i8), allocatable :: ids(:) !< stable lake IDs
    type(lake_support_t), allocatable :: support(:) !< count-weighted sparse mapping
    real(dp), allocatable :: pre_weights(:,:,:) !< aggregated monthly/hourly precipitation weights
    real(dp), allocatable :: pet_weights(:,:,:) !< aggregated monthly/hourly PET weights
    real(dp) :: frac_pre(12) = 0.0_dp !< active monthly night precipitation fractions
    real(dp) :: frac_pet(12) = 0.0_dp !< active monthly night PET fractions
    logical :: disaggregate_pre = .false. !< daily precipitation is disaggregated
    logical :: disaggregate_pet = .false. !< daily PET is disaggregated
    logical :: weighted = .false. !< use hourly weight cubes instead of day/night fractions
    logical :: loaded = .false. !< definition was read successfully
  contains
    procedure :: read => lake_forcing_restart_read
    procedure :: write => lake_forcing_restart_write
    procedure :: validate => lake_forcing_restart_validate
    procedure :: destroy => lake_forcing_restart_destroy
  end type lake_forcing_restart_t
contains
  !> \brief Compare mappings in the same stable-ID order.
  logical function lake_supports_match(a, b) result(matches)
    type(lake_support_t), intent(in) :: a(:) !< first set of lake supports
    type(lake_support_t), intent(in) :: b(:) !< second set of lake supports
    integer :: i
    matches = .false.
    if (size(a) /= size(b)) return
    do i = 1, size(a)
      if (size(a(i)%l2_ids) /= size(b(i)%l2_ids)) return
      if (any(a(i)%l2_ids /= b(i)%l2_ids)) return
      if (any(abs(a(i)%weights-b(i)%weights) > 1.0e-12_dp)) return
    end do
    matches = .true.
  end function lake_supports_match

  !> \brief Reject malformed or empty sparse mappings and temporal weights.
  subroutine lake_forcing_restart_validate(self)
    class(lake_forcing_restart_t), target, intent(in) :: self
    integer :: i, j, n
    if (.not.ieee_is_finite(self%grid%cellsize)) call lake_restart_invalid()
    if (.not.ieee_is_finite(self%grid%xllcorner) .or. .not.ieee_is_finite(self%grid%yllcorner)) call lake_restart_invalid()
    if (self%grid%cellsize <= 0.0_dp .or. self%grid%ncells < 1_i8) call lake_restart_invalid()
    if (.not.allocated(self%ids) .or. .not.allocated(self%support)) call lake_restart_invalid()
    n = size(self%ids)
    if (n < 1 .or. size(self%support) /= n) call lake_restart_invalid()
    if (any(self%ids <= 0_i8)) call lake_restart_invalid()
    do i = 1, n
      if (count(self%ids == self%ids(i)) /= 1) call lake_restart_invalid()
      if (.not.allocated(self%support(i)%l2_ids)) call lake_restart_invalid()
      if (.not.allocated(self%support(i)%weights)) call lake_restart_invalid()
      if (size(self%support(i)%l2_ids) < 1) call lake_restart_invalid()
      if (size(self%support(i)%l2_ids) /= size(self%support(i)%weights)) call lake_restart_invalid()
      if (any(self%support(i)%l2_ids < 1_i8) .or. &
          any(self%support(i)%l2_ids > self%grid%ncells)) call lake_restart_invalid()
      do j = 2, size(self%support(i)%l2_ids)
        if (self%support(i)%l2_ids(j) <= self%support(i)%l2_ids(j-1)) call lake_restart_invalid()
      end do
      if (.not.all(ieee_is_finite(self%support(i)%weights))) call lake_restart_invalid()
      if (any(self%support(i)%weights < 0.0_dp)) call lake_restart_invalid()
      if (abs(sum(self%support(i)%weights)-1.0_dp) > 1.0e-12_dp) call lake_restart_invalid()
    end do
    if (self%disaggregate_pre) then
      if (self%weighted) then
        if (.not.allocated(self%pre_weights)) call lake_restart_invalid()
        call lake_validate_temporal(self%pre_weights, n)
      else
        if (.not.all(ieee_is_finite(self%frac_pre))) call lake_restart_invalid()
        if (any(self%frac_pre < 0.0_dp) .or. any(self%frac_pre > 1.0_dp)) call lake_restart_invalid()
      end if
    end if
    if (self%disaggregate_pet) then
      if (self%weighted) then
        if (.not.allocated(self%pet_weights)) call lake_restart_invalid()
        call lake_validate_temporal(self%pet_weights, n)
      else
        if (.not.all(ieee_is_finite(self%frac_pet))) call lake_restart_invalid()
        if (any(self%frac_pet < 0.0_dp) .or. any(self%frac_pet > 1.0_dp)) call lake_restart_invalid()
      end if
    end if
  end subroutine lake_forcing_restart_validate

  !> \brief Check a monthly/hourly flux weight cube independently of forcing values.
  subroutine lake_validate_temporal(weights, n)
    real(dp), intent(in) :: weights(:,:,:) !< normalized hourly flux weights
    integer, intent(in) :: n !< number of lakes
    if (any(shape(weights) /= [n,12,24])) call lake_restart_invalid()
    if (.not.all(ieee_is_finite(weights))) call lake_restart_invalid()
    if (any(weights < 0.0_dp)) call lake_restart_invalid()
    ! Source weight files may be single precision; retain their values and allow rounding.
    if (any(abs(sum(weights, dim=3)-1.0_dp) > 1.0e-6_dp)) call lake_restart_invalid()
  end subroutine lake_validate_temporal

  !> \brief Write static support, with no forcing values or simulation-time constraint.
  subroutine lake_forcing_restart_write(self, path, creation_time)
    class(lake_forcing_restart_t), target, intent(in) :: self
    character(*), intent(in) :: path !< explicit destination file
    character(*), intent(in) :: creation_time !< informational model time
    type(NcDataset) :: nc
    type(NcVariable) :: var
    type(NcDimension) :: lake_dim, entry_dim, month_dim, hour_dim, scalar_dims(0)
    integer(i8), allocatable :: starts(:), counts(:), indices(:)
    real(dp), allocatable :: weights(:)
    integer(i8) :: first, last
    integer :: i
    call self%validate()
    counts = [(size(self%support(i)%l2_ids, kind=i8), i=1,size(self%ids))]
    allocate(starts(size(counts)), indices(sum(counts)), weights(sum(counts)))
    first = 1_i8
    do i = 1, size(counts)
      starts(i) = first
      last = first + counts(i) - 1_i8
      indices(first:last) = self%support(i)%l2_ids
      weights(first:last) = self%support(i)%weights
      first = last + 1_i8
    end do
    nc = NcDataset(path, "w")
    call self%grid%to_restart(nc, prefix="L2_", mask=.true., area=.false.)
    lake_dim = nc%setDimension("lake", size(self%ids))
    entry_dim = nc%setDimension("support_entry", size(indices))
    month_dim = nc%setDimension("month", 12)
    hour_dim = nc%setDimension("hour", 24)
    var = nc%setVariable("lake_id", "i64", [lake_dim])
    call var%setData(self%ids)
    var = nc%setVariable("support_start", "i64", [lake_dim])
    call var%setData(starts)
    call var%setAttribute("start_index", 1_i4)
    var = nc%setVariable("support_count", "i64", [lake_dim])
    call var%setData(counts)
    var = nc%setVariable("support_l2_index", "i64", [entry_dim])
    call var%setData(indices)
    call var%setAttribute("start_index", 1_i4)
    var = nc%setVariable("support_weight", "f64", [entry_dim])
    call var%setData(weights)
    var = nc%setVariable("frac_night_pre", "f64", [month_dim])
    call var%setData(self%frac_pre)
    var = nc%setVariable("frac_night_pet", "f64", [month_dim])
    call var%setData(self%frac_pet)
    if (allocated(self%pre_weights)) then
      var = nc%setVariable("pre_weights", "f64", [lake_dim, month_dim, hour_dim])
      call var%setData(self%pre_weights)
    end if
    if (allocated(self%pet_weights)) then
      var = nc%setVariable("pet_weights", "f64", [lake_dim, month_dim, hour_dim])
      call var%setData(self%pet_weights)
    end if
    var = nc%setVariable("meteo_lake_meta", "i32", scalar_dims)
    call var%setData(1_i4)
    call var%setAttribute("format_version", 1_i4)
    call var%setAttribute("creation_time", creation_time)
    call var%setAttribute("weighted", merge(1_i4,0_i4,self%weighted))
    call var%setAttribute("disaggregate_pre", merge(1_i4,0_i4,self%disaggregate_pre))
    call var%setAttribute("disaggregate_pet", merge(1_i4,0_i4,self%disaggregate_pet))
    call nc%close()
  end subroutine lake_forcing_restart_write

  !> \brief Read owned definition and validate ragged offsets before using indices.
  subroutine lake_forcing_restart_read(self, path)
    class(lake_forcing_restart_t), target, intent(inout) :: self
    character(*), intent(in) :: path !< explicit source file
    type(NcDataset) :: nc
    type(NcVariable) :: var
    integer(i8), allocatable :: starts(:), counts(:), indices(:), column_counts(:), cumulative(:)
    integer(i1), allocatable :: mask(:,:)
    real(dp), allocatable :: weights(:), fractions(:), x_axis(:), y_axis(:)
    character(:), allocatable :: x_name, y_name
    integer(i4) :: flag, version
    integer(i8) :: first, last, expected
    integer :: i
    character(32), parameter :: required(*) = [character(32) :: "meteo_lake_meta", "lake_id", &
      "support_start", "support_count", "support_l2_index", "support_weight", "L2_mask", &
      "L2_mask_col_cnt", "L2_mask_cum_col_cnt", "frac_night_pre", "frac_night_pet"]
    call self%destroy()
    nc = NcDataset(path, "r")
    do i = 1, size(required)
      if (.not.nc%hasVariable(trim(required(i)))) call lake_restart_invalid()
    end do
    var = nc%getVariable("meteo_lake_meta")
    if (.not.var%hasAttribute("format_version")) call lake_restart_invalid()
    if (.not.var%hasAttribute("weighted") .or. .not.var%hasAttribute("disaggregate_pre") .or. &
        .not.var%hasAttribute("disaggregate_pet")) call lake_restart_invalid()
    call var%getAttribute("format_version", version)
    if (version /= 1_i4) call lake_restart_invalid()
    call var%getAttribute("weighted", flag)
    if (flag /= 0_i4 .and. flag /= 1_i4) call lake_restart_invalid()
    self%weighted = flag == 1_i4
    call var%getAttribute("disaggregate_pre", flag)
    if (flag /= 0_i4 .and. flag /= 1_i4) call lake_restart_invalid()
    self%disaggregate_pre = flag == 1_i4
    call var%getAttribute("disaggregate_pet", flag)
    if (flag /= 0_i4 .and. flag /= 1_i4) call lake_restart_invalid()
    self%disaggregate_pet = flag == 1_i4
    x_name = "L2_x"
    y_name = "L2_y"
    if (.not.nc%hasVariable(x_name)) then
      x_name = "L2_lon"
      y_name = "L2_lat"
    end if
    if (.not.nc%hasVariable(x_name) .or. .not.nc%hasVariable(y_name)) call lake_restart_invalid()
    var = nc%getVariable(x_name)
    call var%getData(x_axis)
    var = nc%getVariable(y_name)
    call var%getData(y_axis)
    if (.not.all(ieee_is_finite(x_axis)) .or. .not.all(ieee_is_finite(y_axis))) call lake_restart_invalid()
    ! Validate packing counts before the grid helper constructs cell_ij from them.
    var = nc%getVariable("L2_mask")
    call var%getData(mask)
    if (any(shape(mask) /= [size(x_axis),size(y_axis)])) call lake_restart_invalid()
    var = nc%getVariable("L2_mask_col_cnt")
    call var%getData(column_counts)
    var = nc%getVariable("L2_mask_cum_col_cnt")
    call var%getData(cumulative)
    if (size(column_counts) /= size(mask,2) .or. size(cumulative) /= size(mask,2)) call lake_restart_invalid()
    if (any(mask /= 0_i1 .and. mask /= 1_i1)) call lake_restart_invalid()
    expected = 0_i8
    do i = 1, size(column_counts)
      if (column_counts(i) /= count(mask(:,i) == 1_i1, kind=i8)) call lake_restart_invalid()
      if (cumulative(i) /= expected) call lake_restart_invalid()
      expected = expected + column_counts(i)
    end do
    call self%grid%from_restart(nc, prefix="L2_")
    if (any(abs(x_axis-self%grid%x_axis()) > 1.0e-10_dp*max(1.0_dp,self%grid%cellsize))) call lake_restart_invalid()
    if (any(abs(y_axis-self%grid%y_axis()) > 1.0e-10_dp*max(1.0_dp,self%grid%cellsize))) call lake_restart_invalid()
    var = nc%getVariable("lake_id")
    call var%getData(self%ids)
    var = nc%getVariable("support_start")
    call var%getData(starts)
    var = nc%getVariable("support_count")
    call var%getData(counts)
    var = nc%getVariable("support_l2_index")
    call var%getData(indices)
    var = nc%getVariable("support_weight")
    call var%getData(weights)
    if (size(starts) /= size(self%ids) .or. size(counts) /= size(self%ids)) call lake_restart_invalid()
    if (size(weights) /= size(indices)) call lake_restart_invalid()
    allocate(self%support(size(self%ids)))
    expected = 1_i8
    do i = 1, size(self%ids)
      first = starts(i)
      if (first /= expected .or. counts(i) < 1_i8) call lake_restart_invalid()
      if (counts(i) > size(indices,kind=i8)-first+1_i8) call lake_restart_invalid()
      last = first + counts(i) - 1_i8
      self%support(i)%l2_ids = indices(first:last)
      self%support(i)%weights = weights(first:last)
      expected = last + 1_i8
    end do
    if (expected /= size(indices,kind=i8)+1_i8) call lake_restart_invalid()
    var = nc%getVariable("frac_night_pre")
    call var%getData(fractions)
    if (size(fractions) /= 12) call lake_restart_invalid()
    self%frac_pre = fractions
    var = nc%getVariable("frac_night_pet")
    call var%getData(fractions)
    if (size(fractions) /= 12) call lake_restart_invalid()
    self%frac_pet = fractions
    if (nc%hasVariable("pre_weights")) then
      var = nc%getVariable("pre_weights")
      call var%getData(self%pre_weights)
    end if
    if (nc%hasVariable("pet_weights")) then
      var = nc%getVariable("pet_weights")
      call var%getData(self%pet_weights)
    end if
    call nc%close()
    call self%validate()
    self%loaded = .true.
  end subroutine lake_forcing_restart_read

  !> \brief Release definition after consumers have detached.
  subroutine lake_forcing_restart_destroy(self)
    class(lake_forcing_restart_t), target, intent(inout) :: self
    call self%grid%destroy()
    if (allocated(self%ids)) deallocate(self%ids)
    if (allocated(self%support)) deallocate(self%support)
    if (allocated(self%pre_weights)) deallocate(self%pre_weights)
    if (allocated(self%pet_weights)) deallocate(self%pet_weights)
    self%loaded = .false.
    self%weighted = .false.
    self%disaggregate_pre = .false.
    self%disaggregate_pet = .false.
    self%frac_pre = 0.0_dp
    self%frac_pet = 0.0_dp
  end subroutine lake_forcing_restart_destroy

  !> \brief Report a corrupt or unsupported lake forcing restart.
  subroutine lake_restart_invalid()
    log_fatal(*) "Meteo: malformed lake support restart (IDs, packing, indices, weights or settings)."
    error stop 1
  end subroutine lake_restart_invalid
end module mo_meteo_lake_restart
