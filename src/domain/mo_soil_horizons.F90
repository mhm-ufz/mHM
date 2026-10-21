!> \file mo_soil_horizons.F90
!> \copydoc mo_soil_horizons

!> \brief Soil-horizon coordinate metadata and restart compatibility.
!> \details Shares the v6 mHM restart-bound validation with MPR and layered outputs.
!> \authors Sebastian Mueller
!> \version 0.1
!> \changelog
!! - 2026: shared validation of fixed interfaces and the unspecified-bottom convention.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_domain
#include "logging.h"
module mo_soil_horizons
  use mo_logging
  use mo_kind, only: dp, i4
  use mo_netcdf, only: NcDataset, NcVariable
  use mo_utils, only: is_close, optval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: write_horizon_metadata, read_horizon_bounds, validate_horizon_bounds, compare_horizon_bounds
contains
  !> \brief Add soil-specific coordinate metadata and optional general axis attributes.
  subroutine write_horizon_metadata(nc, name, variable_bottom, general)
    type(NcDataset), intent(in) :: nc !< Open dataset.
    character(*), intent(in) :: name !< Coordinate variable name.
    logical, intent(in) :: variable_bottom !< Final layer has a soil-specific bottom depth.
    logical, intent(in), optional :: general !< Write axis, positive, and standard_name attributes (default true).
    type(NcVariable) :: axis
    axis = nc%getVariable(name)
    call axis%setAttribute('units', 'mm')
    if (optval(general, .true.)) then
      call axis%setAttribute('positive', 'down')
      call axis%setAttribute('axis', 'Z')
      call axis%setAttribute('standard_name', 'depth')
    end if
    call axis%setAttribute('long_name', 'soil layer top depth')
    if (variable_bottom) then
      call axis%setAttribute('comment', &
        'mHM convention: coincident bounds on the final layer indicate an unspecified, spatially variable bottom depth; ' // &
        'they do not imply zero physical thickness.')
    end if
  end subroutine write_horizon_metadata

  !> \brief Read bounds independently of legacy bottom or current top coordinate values.
  subroutine read_horizon_bounds(nc, name, n, variable_bottom, bounds)
    type(NcDataset), intent(in) :: nc !< Open restart dataset.
    character(*), intent(in) :: name !< Coordinate name.
    integer(i4), intent(in) :: n !< Required number of layers.
    logical, intent(in) :: variable_bottom !< Configured soil-specific bottom-depth convention.
    real(dp), allocatable, intent(out) :: bounds(:) !< Restored boundaries [mm], with a repeated variable endpoint.
    type(NcVariable) :: axis, bvar
    character(256) :: bname
    character(32) :: old_bottom_type
    real(dp), allocatable :: raw(:, :), pairs(:, :)
    integer(i4), allocatable :: shp(:)
    integer(i4) :: k
    if (.not.nc%hasVariable(name)) then
      log_fatal(*) 'Soil horizons: required coordinate missing: ', name
      error stop 1
    end if
    axis = nc%getVariable(name)
    if (.not.axis%hasAttribute('bounds')) then
      log_fatal(*) 'Soil horizons: coordinate is missing bounds metadata.'
      error stop 1
    end if
    call axis%getAttribute('bounds', bname)
    if (.not.nc%hasVariable(trim(bname))) then
      log_fatal(*) 'Soil horizons: required bounds variable missing: ', trim(bname)
      error stop 1
    end if
    bvar = nc%getVariable(trim(bname))
    shp = bvar%getShape()
    if (size(shp) /= 2) then
      log_fatal(*) 'Soil horizons: bounds must have rank two.'
      error stop 1
    end if
    call bvar%getData(raw)
    if (size(raw, 1) == 2 .and. size(raw, 2) == n) then
      pairs = raw
    else if (size(raw, 2) == 2 .and. size(raw, 1) == n) then
      pairs = transpose(raw)
    else
      log_fatal(*) 'Soil horizons: bounds shape does not match the configured layer count.'
      error stop 1
    end if
    if (any(.not.ieee_is_finite(pairs))) then
      log_fatal(*) 'Soil horizons: bounds must be finite.'
      error stop 1
    end if
    allocate(bounds(n + 1))
    bounds(1) = pairs(1, 1)
    bounds(2:) = pairs(2, :)
    do k = 2, n
      if (.not.is_close(pairs(1, k), bounds(k))) then
        log_fatal(*) 'Soil horizons: non-contiguous bounds.'
        error stop 1
      end if
    end do
    ! Earlier v6 files may declare a fixed or nominal endpoint explicitly.
    if (axis%hasAttribute('bottom_boundary_type')) then
      call axis%getAttribute('bottom_boundary_type', old_bottom_type)
      select case (old_bottom_type)
      case ('fixed', 'domain_maximum', 'legacy_nominal')
      case default
        log_fatal(*) 'Soil horizons: unknown legacy bottom_boundary_type: ', trim(old_bottom_type)
        error stop 1
      end select
      if (variable_bottom .eqv. (old_bottom_type == 'fixed')) then
        log_fatal(*) 'Soil horizons: incompatible fixed/variable bottom-boundary conventions.'
        error stop 1
      end if
    end if
    if (.not.variable_bottom .and. bounds(n + 1) == bounds(n)) then
      log_fatal(*) 'Soil horizons: incompatible fixed/variable bottom-boundary conventions.'
      error stop 1
    end if
    call validate_horizon_bounds(bounds, variable_bottom)
    ! Legacy nominal endpoints convey no physical L1 depth; normalize before publication.
    if (variable_bottom) bounds(n + 1) = bounds(n)
  end subroutine read_horizon_bounds

  !> \brief Validate fixed interfaces, permitting an unspecified variable terminal depth.
  subroutine validate_horizon_bounds(bounds, variable_bottom)
    real(dp), intent(in) :: bounds(:) !< Complete boundary vector [mm].
    logical, intent(in) :: variable_bottom !< Permit coincident final bounds or an older nominal endpoint.
    integer :: last
    if (size(bounds) < 2) then
      log_fatal(*) 'Soil horizons: at least one layer is required.'
      error stop 1
    end if
    if (any(.not.ieee_is_finite(bounds))) then
      log_fatal(*) 'Soil horizons: bounds must be finite.'
      error stop 1
    end if
    last = size(bounds)
    if (variable_bottom) last = last - 1
    if (.not.is_close(bounds(1), 0.0_dp) .or. any(bounds(2:last) <= bounds(:last-1))) then
      log_fatal(*) 'Soil horizons: bounds must start at zero and increase strictly.'
      error stop 1
    end if
    if (bounds(size(bounds)) < bounds(size(bounds) - 1)) then
      log_fatal(*) 'Soil horizons: terminal bound must not lie above the final layer top.'
      error stop 1
    end if
  end subroutine validate_horizon_bounds

  !> \brief Compare fixed interfaces, ignoring the terminal bound when soil depth varies.
  subroutine compare_horizon_bounds(actual, expected, variable_bottom)
    real(dp), intent(in) :: actual(:) !< Restored boundaries [mm].
    real(dp), intent(in) :: expected(:) !< Configured or connected boundaries [mm].
    logical, intent(in) :: variable_bottom !< Final layer has a soil-specific bottom depth.
    integer :: last
    if (size(actual) /= size(expected)) then
      log_fatal(*) 'Soil horizons: layer counts do not match.'
      error stop 1
    end if
    last = size(actual) - 1
    if (.not.variable_bottom) last = size(actual)
    if (any(.not.is_close(actual(:last), expected(:last)))) then
      log_fatal(*) 'Soil horizons: fixed internal interfaces do not match the current soil axis.'
      error stop 1
    end if
  end subroutine compare_horizon_bounds
end module mo_soil_horizons
