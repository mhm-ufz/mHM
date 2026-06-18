!> \file nml_parameter_routing1.f90
!> \copydoc nml_routing1

!> \brief Routing - Case 1
!> \details Parameters for routing (case 1 - Muskingum).
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_routing1
  use nml_helper, only: &
    nml_file_t, &
    nml_line_buffer, &
    NML_OK, &
    NML_ERR_FILE_NOT_FOUND, &
    NML_ERR_OPEN, &
    NML_ERR_NOT_OPEN, &
    NML_ERR_NML_NOT_FOUND, &
    NML_ERR_READ, &
    NML_ERR_CLOSE, &
    NML_ERR_REQUIRED, &
    NML_ERR_ENUM, &
    NML_ERR_BOUNDS, &
    NML_ERR_NOT_SET, &
    NML_ERR_INVALID_NAME, &
    NML_ERR_INVALID_INDEX, &
    idx_check, &
    to_lower, &
    NML_ERR_PARTLY_SET
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  !> \class nml_routing1_t
  !> \brief Routing - Case 1
  !> \details Parameters for routing (case 1 - Muskingum).
  type, public :: nml_routing1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: muskingumtraveltime_constant !< Muskingum travel time constant
    real(dp), dimension(5) :: muskingumtraveltime_riverlength !< Muskingum travel time river length
    real(dp), dimension(5) :: muskingumtraveltime_riverslope !< Muskingum travel time river slope
    real(dp), dimension(5) :: muskingumtraveltime_impervious !< Muskingum travel time impervious
    real(dp), dimension(5) :: muskingumattenuation_riverslope !< Muskingum attenuation river slope
  contains
    procedure :: init => nml_routing1_init
    procedure :: from_file => nml_routing1_from_file
    procedure :: set => nml_routing1_set
    procedure :: is_set => nml_routing1_is_set
    procedure :: is_valid => nml_routing1_is_valid
  end type nml_routing1_t

contains

  !> \brief Initialize defaults and sentinels for routing1
  integer function nml_routing1_init(this, errmsg) result(status)
    class(nml_routing1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%muskingumtraveltime_constant = ieee_value(this%muskingumtraveltime_constant, ieee_quiet_nan) ! sentinel for required real array
    this%muskingumtraveltime_riverlength = ieee_value(this%muskingumtraveltime_riverlength, ieee_quiet_nan) ! sentinel for required real array
    this%muskingumtraveltime_riverslope = ieee_value(this%muskingumtraveltime_riverslope, ieee_quiet_nan) ! sentinel for required real array
    this%muskingumtraveltime_impervious = ieee_value(this%muskingumtraveltime_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%muskingumattenuation_riverslope = ieee_value(this%muskingumattenuation_riverslope, ieee_quiet_nan) ! sentinel for required real array
  end function nml_routing1_init


  !> \brief Read routing1 namelist from file
  integer function nml_routing1_from_file(this, file, errmsg) result(status)
    class(nml_routing1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: muskingumtraveltime_constant
    real(dp), dimension(5) :: muskingumtraveltime_riverlength
    real(dp), dimension(5) :: muskingumtraveltime_riverslope
    real(dp), dimension(5) :: muskingumtraveltime_impervious
    real(dp), dimension(5) :: muskingumattenuation_riverslope
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /routing1/ &
      muskingumtraveltime_constant, &
      muskingumtraveltime_riverlength, &
      muskingumtraveltime_riverslope, &
      muskingumtraveltime_impervious, &
      muskingumattenuation_riverslope

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    muskingumtraveltime_constant = this%muskingumtraveltime_constant
    muskingumtraveltime_riverlength = this%muskingumtraveltime_riverlength
    muskingumtraveltime_riverslope = this%muskingumtraveltime_riverslope
    muskingumtraveltime_impervious = this%muskingumtraveltime_impervious
    muskingumattenuation_riverslope = this%muskingumattenuation_riverslope

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("routing1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=routing1, iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) then
      status = NML_ERR_READ
      if (present(errmsg)) errmsg = trim(iomsg)
      close_status = nml%close()
      return
    end if
    close_status = nml%close(errmsg=errmsg)
    if (close_status /= NML_OK) then
      status = close_status
      return
    end if

    ! assign values
    this%muskingumtraveltime_constant = muskingumtraveltime_constant
    this%muskingumtraveltime_riverlength = muskingumtraveltime_riverlength
    this%muskingumtraveltime_riverslope = muskingumtraveltime_riverslope
    this%muskingumtraveltime_impervious = muskingumtraveltime_impervious
    this%muskingumattenuation_riverslope = muskingumattenuation_riverslope

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing1_from_file

  !> \brief Set routing1 values
  integer function nml_routing1_set(this, &
    muskingumtraveltime_constant, &
    muskingumtraveltime_riverlength, &
    muskingumtraveltime_riverslope, &
    muskingumtraveltime_impervious, &
    muskingumattenuation_riverslope, &
    errmsg) result(status)

    class(nml_routing1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: muskingumtraveltime_constant !< Muskingum travel time constant
    real(dp), dimension(:), intent(in) :: muskingumtraveltime_riverlength !< Muskingum travel time river length
    real(dp), dimension(:), intent(in) :: muskingumtraveltime_riverslope !< Muskingum travel time river slope
    real(dp), dimension(:), intent(in) :: muskingumtraveltime_impervious !< Muskingum travel time impervious
    real(dp), dimension(:), intent(in) :: muskingumattenuation_riverslope !< Muskingum attenuation river slope
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(muskingumtraveltime_constant, 1) > size(this%muskingumtraveltime_constant, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'muskingumtraveltime_constant'"
      return
    end if
    lb__1 = lbound(this%muskingumtraveltime_constant, 1)
    ub__1 = lb__1 + size(muskingumtraveltime_constant, 1) - 1
    this%muskingumtraveltime_constant(lb__1:ub__1) = muskingumtraveltime_constant
    if (size(muskingumtraveltime_riverlength, 1) > size(this%muskingumtraveltime_riverlength, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'muskingumtraveltime_riverlength'"
      return
    end if
    lb__1 = lbound(this%muskingumtraveltime_riverlength, 1)
    ub__1 = lb__1 + size(muskingumtraveltime_riverlength, 1) - 1
    this%muskingumtraveltime_riverlength(lb__1:ub__1) = muskingumtraveltime_riverlength
    if (size(muskingumtraveltime_riverslope, 1) > size(this%muskingumtraveltime_riverslope, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'muskingumtraveltime_riverslope'"
      return
    end if
    lb__1 = lbound(this%muskingumtraveltime_riverslope, 1)
    ub__1 = lb__1 + size(muskingumtraveltime_riverslope, 1) - 1
    this%muskingumtraveltime_riverslope(lb__1:ub__1) = muskingumtraveltime_riverslope
    if (size(muskingumtraveltime_impervious, 1) > size(this%muskingumtraveltime_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'muskingumtraveltime_impervious'"
      return
    end if
    lb__1 = lbound(this%muskingumtraveltime_impervious, 1)
    ub__1 = lb__1 + size(muskingumtraveltime_impervious, 1) - 1
    this%muskingumtraveltime_impervious(lb__1:ub__1) = muskingumtraveltime_impervious
    if (size(muskingumattenuation_riverslope, 1) > size(this%muskingumattenuation_riverslope, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'muskingumattenuation_riverslope'"
      return
    end if
    lb__1 = lbound(this%muskingumattenuation_riverslope, 1)
    ub__1 = lb__1 + size(muskingumattenuation_riverslope, 1) - 1
    this%muskingumattenuation_riverslope(lb__1:ub__1) = muskingumattenuation_riverslope

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing1_set

  !> \brief Check whether a namelist value was set
  integer function nml_routing1_is_set(this, name, idx, errmsg) result(status)
    class(nml_routing1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field name
    integer, intent(in), optional :: idx(:) !< optional field index values
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if
    select case (to_lower(trim(name)))
    case ("muskingumtraveltime_constant")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%muskingumtraveltime_constant), ubound(this%muskingumtraveltime_constant), &
          "muskingumTravelTime_constant", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%muskingumtraveltime_constant(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%muskingumtraveltime_constant))) status = NML_ERR_NOT_SET
      end if
    case ("muskingumtraveltime_riverlength")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%muskingumtraveltime_riverlength), ubound(this%muskingumtraveltime_riverlength), &
          "muskingumTravelTime_riverLength", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%muskingumtraveltime_riverlength(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%muskingumtraveltime_riverlength))) status = NML_ERR_NOT_SET
      end if
    case ("muskingumtraveltime_riverslope")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%muskingumtraveltime_riverslope), ubound(this%muskingumtraveltime_riverslope), &
          "muskingumTravelTime_riverSlope", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%muskingumtraveltime_riverslope(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%muskingumtraveltime_riverslope))) status = NML_ERR_NOT_SET
      end if
    case ("muskingumtraveltime_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%muskingumtraveltime_impervious), ubound(this%muskingumtraveltime_impervious), &
          "muskingumTravelTime_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%muskingumtraveltime_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%muskingumtraveltime_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("muskingumattenuation_riverslope")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%muskingumattenuation_riverslope), ubound(this%muskingumattenuation_riverslope), &
          "muskingumAttenuation_riverSlope", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%muskingumattenuation_riverslope(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%muskingumattenuation_riverslope))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_routing1_is_set

  !> \brief Validate required values and constraints
  integer function nml_routing1_is_valid(this, errmsg) result(status)
    class(nml_routing1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required arrays
    if (all(ieee_is_nan(this%muskingumtraveltime_constant(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: muskingumTravelTime_constant"
      return
    end if
    if (any(ieee_is_nan(this%muskingumtraveltime_constant(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: muskingumTravelTime_constant"
      return
    end if
    if (all(ieee_is_nan(this%muskingumtraveltime_riverlength(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: muskingumTravelTime_riverLength"
      return
    end if
    if (any(ieee_is_nan(this%muskingumtraveltime_riverlength(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: muskingumTravelTime_riverLength"
      return
    end if
    if (all(ieee_is_nan(this%muskingumtraveltime_riverslope(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: muskingumTravelTime_riverSlope"
      return
    end if
    if (any(ieee_is_nan(this%muskingumtraveltime_riverslope(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: muskingumTravelTime_riverSlope"
      return
    end if
    if (all(ieee_is_nan(this%muskingumtraveltime_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: muskingumTravelTime_impervious"
      return
    end if
    if (any(ieee_is_nan(this%muskingumtraveltime_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: muskingumTravelTime_impervious"
      return
    end if
    if (all(ieee_is_nan(this%muskingumattenuation_riverslope(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: muskingumAttenuation_riverSlope"
      return
    end if
    if (any(ieee_is_nan(this%muskingumattenuation_riverslope(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: muskingumAttenuation_riverSlope"
      return
    end if
  end function nml_routing1_is_valid

end module nml_routing1
