!> \file nml_parameter_snow1.f90
!> \copydoc nml_snow1

!> \brief Snow - Case 1
!> \details Parameters for Snow module.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_snow1
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

  !> \class nml_snow1_t
  !> \brief Snow - Case 1
  !> \details Parameters for Snow module.
  type, public :: nml_snow1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: snowtreshholdtemperature !< Threshold for rain/snow partitioning [degC].
    real(dp), dimension(5) :: degreedayfactor_forest !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(5) :: degreedayfactor_impervious !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(5) :: degreedayfactor_pervious !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(5) :: increasedegreedayfactorbyprecip !< Increase of degree day factor if there is precipitation [degC-1].
    real(dp), dimension(5) :: maxdegreedayfactor_forest !< Maximum values for degree day factor [m degC-1].
    real(dp), dimension(5) :: maxdegreedayfactor_impervious !< Maximum values for degree day factor [m degC-1].
    real(dp), dimension(5) :: maxdegreedayfactor_pervious !< Maximum values for degree day factor [m degC-1].
  contains
    procedure :: init => nml_snow1_init
    procedure :: from_file => nml_snow1_from_file
    procedure :: set => nml_snow1_set
    procedure :: is_set => nml_snow1_is_set
    procedure :: is_valid => nml_snow1_is_valid
  end type nml_snow1_t

contains

  !> \brief Initialize defaults and sentinels for snow1
  integer function nml_snow1_init(this, errmsg) result(status)
    class(nml_snow1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%snowtreshholdtemperature = ieee_value(this%snowtreshholdtemperature, ieee_quiet_nan) ! sentinel for required real array
    this%degreedayfactor_forest = ieee_value(this%degreedayfactor_forest, ieee_quiet_nan) ! sentinel for required real array
    this%degreedayfactor_impervious = ieee_value(this%degreedayfactor_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%degreedayfactor_pervious = ieee_value(this%degreedayfactor_pervious, ieee_quiet_nan) ! sentinel for required real array
    this%increasedegreedayfactorbyprecip = ieee_value(this%increasedegreedayfactorbyprecip, ieee_quiet_nan) ! sentinel for required real array
    this%maxdegreedayfactor_forest = ieee_value(this%maxdegreedayfactor_forest, ieee_quiet_nan) ! sentinel for required real array
    this%maxdegreedayfactor_impervious = ieee_value(this%maxdegreedayfactor_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%maxdegreedayfactor_pervious = ieee_value(this%maxdegreedayfactor_pervious, ieee_quiet_nan) ! sentinel for required real array
  end function nml_snow1_init


  !> \brief Read snow1 namelist from file
  integer function nml_snow1_from_file(this, file, errmsg) result(status)
    class(nml_snow1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: snowtreshholdtemperature
    real(dp), dimension(5) :: degreedayfactor_forest
    real(dp), dimension(5) :: degreedayfactor_impervious
    real(dp), dimension(5) :: degreedayfactor_pervious
    real(dp), dimension(5) :: increasedegreedayfactorbyprecip
    real(dp), dimension(5) :: maxdegreedayfactor_forest
    real(dp), dimension(5) :: maxdegreedayfactor_impervious
    real(dp), dimension(5) :: maxdegreedayfactor_pervious
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /snow1/ &
      snowtreshholdtemperature, &
      degreedayfactor_forest, &
      degreedayfactor_impervious, &
      degreedayfactor_pervious, &
      increasedegreedayfactorbyprecip, &
      maxdegreedayfactor_forest, &
      maxdegreedayfactor_impervious, &
      maxdegreedayfactor_pervious

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    snowtreshholdtemperature = this%snowtreshholdtemperature
    degreedayfactor_forest = this%degreedayfactor_forest
    degreedayfactor_impervious = this%degreedayfactor_impervious
    degreedayfactor_pervious = this%degreedayfactor_pervious
    increasedegreedayfactorbyprecip = this%increasedegreedayfactorbyprecip
    maxdegreedayfactor_forest = this%maxdegreedayfactor_forest
    maxdegreedayfactor_impervious = this%maxdegreedayfactor_impervious
    maxdegreedayfactor_pervious = this%maxdegreedayfactor_pervious

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("snow1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=snow1, iostat=iostat, iomsg=iomsg)
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
    this%snowtreshholdtemperature = snowtreshholdtemperature
    this%degreedayfactor_forest = degreedayfactor_forest
    this%degreedayfactor_impervious = degreedayfactor_impervious
    this%degreedayfactor_pervious = degreedayfactor_pervious
    this%increasedegreedayfactorbyprecip = increasedegreedayfactorbyprecip
    this%maxdegreedayfactor_forest = maxdegreedayfactor_forest
    this%maxdegreedayfactor_impervious = maxdegreedayfactor_impervious
    this%maxdegreedayfactor_pervious = maxdegreedayfactor_pervious

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_snow1_from_file

  !> \brief Set snow1 values
  integer function nml_snow1_set(this, &
    snowtreshholdtemperature, &
    degreedayfactor_forest, &
    degreedayfactor_impervious, &
    degreedayfactor_pervious, &
    increasedegreedayfactorbyprecip, &
    maxdegreedayfactor_forest, &
    maxdegreedayfactor_impervious, &
    maxdegreedayfactor_pervious, &
    errmsg) result(status)

    class(nml_snow1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: snowtreshholdtemperature !< Threshold for rain/snow partitioning [degC].
    real(dp), dimension(:), intent(in) :: degreedayfactor_forest !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(:), intent(in) :: degreedayfactor_impervious !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(:), intent(in) :: degreedayfactor_pervious !< Degree day factors to determine melting flux [m degC-1].
    real(dp), dimension(:), intent(in) :: increasedegreedayfactorbyprecip !< Increase of degree day factor if there is precipitation [degC-1].
    real(dp), dimension(:), intent(in) :: maxdegreedayfactor_forest !< Maximum values for degree day factor [m degC-1].
    real(dp), dimension(:), intent(in) :: maxdegreedayfactor_impervious !< Maximum values for degree day factor [m degC-1].
    real(dp), dimension(:), intent(in) :: maxdegreedayfactor_pervious !< Maximum values for degree day factor [m degC-1].
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(snowtreshholdtemperature, 1) > size(this%snowtreshholdtemperature, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'snowtreshholdtemperature'"
      return
    end if
    lb__1 = lbound(this%snowtreshholdtemperature, 1)
    ub__1 = lb__1 + size(snowtreshholdtemperature, 1) - 1
    this%snowtreshholdtemperature(lb__1:ub__1) = snowtreshholdtemperature
    if (size(degreedayfactor_forest, 1) > size(this%degreedayfactor_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'degreedayfactor_forest'"
      return
    end if
    lb__1 = lbound(this%degreedayfactor_forest, 1)
    ub__1 = lb__1 + size(degreedayfactor_forest, 1) - 1
    this%degreedayfactor_forest(lb__1:ub__1) = degreedayfactor_forest
    if (size(degreedayfactor_impervious, 1) > size(this%degreedayfactor_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'degreedayfactor_impervious'"
      return
    end if
    lb__1 = lbound(this%degreedayfactor_impervious, 1)
    ub__1 = lb__1 + size(degreedayfactor_impervious, 1) - 1
    this%degreedayfactor_impervious(lb__1:ub__1) = degreedayfactor_impervious
    if (size(degreedayfactor_pervious, 1) > size(this%degreedayfactor_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'degreedayfactor_pervious'"
      return
    end if
    lb__1 = lbound(this%degreedayfactor_pervious, 1)
    ub__1 = lb__1 + size(degreedayfactor_pervious, 1) - 1
    this%degreedayfactor_pervious(lb__1:ub__1) = degreedayfactor_pervious
    if (size(increasedegreedayfactorbyprecip, 1) > size(this%increasedegreedayfactorbyprecip, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'increasedegreedayfactorbyprecip'"
      return
    end if
    lb__1 = lbound(this%increasedegreedayfactorbyprecip, 1)
    ub__1 = lb__1 + size(increasedegreedayfactorbyprecip, 1) - 1
    this%increasedegreedayfactorbyprecip(lb__1:ub__1) = increasedegreedayfactorbyprecip
    if (size(maxdegreedayfactor_forest, 1) > size(this%maxdegreedayfactor_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'maxdegreedayfactor_forest'"
      return
    end if
    lb__1 = lbound(this%maxdegreedayfactor_forest, 1)
    ub__1 = lb__1 + size(maxdegreedayfactor_forest, 1) - 1
    this%maxdegreedayfactor_forest(lb__1:ub__1) = maxdegreedayfactor_forest
    if (size(maxdegreedayfactor_impervious, 1) > size(this%maxdegreedayfactor_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'maxdegreedayfactor_impervious'"
      return
    end if
    lb__1 = lbound(this%maxdegreedayfactor_impervious, 1)
    ub__1 = lb__1 + size(maxdegreedayfactor_impervious, 1) - 1
    this%maxdegreedayfactor_impervious(lb__1:ub__1) = maxdegreedayfactor_impervious
    if (size(maxdegreedayfactor_pervious, 1) > size(this%maxdegreedayfactor_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'maxdegreedayfactor_pervious'"
      return
    end if
    lb__1 = lbound(this%maxdegreedayfactor_pervious, 1)
    ub__1 = lb__1 + size(maxdegreedayfactor_pervious, 1) - 1
    this%maxdegreedayfactor_pervious(lb__1:ub__1) = maxdegreedayfactor_pervious

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_snow1_set

  !> \brief Check whether a namelist value was set
  integer function nml_snow1_is_set(this, name, idx, errmsg) result(status)
    class(nml_snow1_t), intent(in) :: this !< namelist instance
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
    case ("snowtreshholdtemperature")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%snowtreshholdtemperature), ubound(this%snowtreshholdtemperature), &
          "snowTreshholdTemperature", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%snowtreshholdtemperature(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%snowtreshholdtemperature))) status = NML_ERR_NOT_SET
      end if
    case ("degreedayfactor_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%degreedayfactor_forest), ubound(this%degreedayfactor_forest), &
          "degreeDayFactor_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%degreedayfactor_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%degreedayfactor_forest))) status = NML_ERR_NOT_SET
      end if
    case ("degreedayfactor_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%degreedayfactor_impervious), ubound(this%degreedayfactor_impervious), &
          "degreeDayFactor_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%degreedayfactor_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%degreedayfactor_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("degreedayfactor_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%degreedayfactor_pervious), ubound(this%degreedayfactor_pervious), &
          "degreeDayFactor_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%degreedayfactor_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%degreedayfactor_pervious))) status = NML_ERR_NOT_SET
      end if
    case ("increasedegreedayfactorbyprecip")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%increasedegreedayfactorbyprecip), ubound(this%increasedegreedayfactorbyprecip), &
          "increaseDegreeDayFactorByPrecip", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%increasedegreedayfactorbyprecip(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%increasedegreedayfactorbyprecip))) status = NML_ERR_NOT_SET
      end if
    case ("maxdegreedayfactor_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%maxdegreedayfactor_forest), ubound(this%maxdegreedayfactor_forest), &
          "maxDegreeDayFactor_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%maxdegreedayfactor_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%maxdegreedayfactor_forest))) status = NML_ERR_NOT_SET
      end if
    case ("maxdegreedayfactor_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%maxdegreedayfactor_impervious), ubound(this%maxdegreedayfactor_impervious), &
          "maxDegreeDayFactor_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%maxdegreedayfactor_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%maxdegreedayfactor_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("maxdegreedayfactor_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%maxdegreedayfactor_pervious), ubound(this%maxdegreedayfactor_pervious), &
          "maxDegreeDayFactor_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%maxdegreedayfactor_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%maxdegreedayfactor_pervious))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_snow1_is_set

  !> \brief Validate required values and constraints
  integer function nml_snow1_is_valid(this, errmsg) result(status)
    class(nml_snow1_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%snowtreshholdtemperature(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: snowTreshholdTemperature"
      return
    end if
    if (any(ieee_is_nan(this%snowtreshholdtemperature(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: snowTreshholdTemperature"
      return
    end if
    if (all(ieee_is_nan(this%degreedayfactor_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: degreeDayFactor_forest"
      return
    end if
    if (any(ieee_is_nan(this%degreedayfactor_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: degreeDayFactor_forest"
      return
    end if
    if (all(ieee_is_nan(this%degreedayfactor_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: degreeDayFactor_impervious"
      return
    end if
    if (any(ieee_is_nan(this%degreedayfactor_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: degreeDayFactor_impervious"
      return
    end if
    if (all(ieee_is_nan(this%degreedayfactor_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: degreeDayFactor_pervious"
      return
    end if
    if (any(ieee_is_nan(this%degreedayfactor_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: degreeDayFactor_pervious"
      return
    end if
    if (all(ieee_is_nan(this%increasedegreedayfactorbyprecip(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: increaseDegreeDayFactorByPrecip"
      return
    end if
    if (any(ieee_is_nan(this%increasedegreedayfactorbyprecip(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: increaseDegreeDayFactorByPrecip"
      return
    end if
    if (all(ieee_is_nan(this%maxdegreedayfactor_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: maxDegreeDayFactor_forest"
      return
    end if
    if (any(ieee_is_nan(this%maxdegreedayfactor_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: maxDegreeDayFactor_forest"
      return
    end if
    if (all(ieee_is_nan(this%maxdegreedayfactor_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: maxDegreeDayFactor_impervious"
      return
    end if
    if (any(ieee_is_nan(this%maxdegreedayfactor_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: maxDegreeDayFactor_impervious"
      return
    end if
    if (all(ieee_is_nan(this%maxdegreedayfactor_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: maxDegreeDayFactor_pervious"
      return
    end if
    if (any(ieee_is_nan(this%maxdegreedayfactor_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: maxDegreeDayFactor_pervious"
      return
    end if
  end function nml_snow1_is_valid

end module nml_snow1
