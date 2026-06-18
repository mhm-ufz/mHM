!> \file nml_parameter_rivertemp1.f90
!> \copydoc nml_rivertemp1

!> \brief River Temperature - Case 1
!> \details Parameters for River Temperature case 1.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_rivertemp1
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

  !> \class nml_rivertemp1_t
  !> \brief River Temperature - Case 1
  !> \details Parameters for River Temperature case 1.
  type, public :: nml_rivertemp1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: albedo_water !< Albedo of open water [-]
    real(dp), dimension(5) :: pt_a_water !< Priestley-Taylor coefficient for open water [-]
    real(dp), dimension(5) :: emissivity_water !< Emissivity of open water [-]
    real(dp), dimension(5) :: turb_heat_ex_coeff !< Turbulent heat exchange coefficient for open water [W m-2 K-1]
  contains
    procedure :: init => nml_rivertemp1_init
    procedure :: from_file => nml_rivertemp1_from_file
    procedure :: set => nml_rivertemp1_set
    procedure :: is_set => nml_rivertemp1_is_set
    procedure :: is_valid => nml_rivertemp1_is_valid
  end type nml_rivertemp1_t

contains

  !> \brief Initialize defaults and sentinels for rivertemp1
  integer function nml_rivertemp1_init(this, errmsg) result(status)
    class(nml_rivertemp1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%albedo_water = ieee_value(this%albedo_water, ieee_quiet_nan) ! sentinel for required real array
    this%pt_a_water = ieee_value(this%pt_a_water, ieee_quiet_nan) ! sentinel for required real array
    this%emissivity_water = ieee_value(this%emissivity_water, ieee_quiet_nan) ! sentinel for required real array
    this%turb_heat_ex_coeff = ieee_value(this%turb_heat_ex_coeff, ieee_quiet_nan) ! sentinel for required real array
  end function nml_rivertemp1_init


  !> \brief Read rivertemp1 namelist from file
  integer function nml_rivertemp1_from_file(this, file, errmsg) result(status)
    class(nml_rivertemp1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: albedo_water
    real(dp), dimension(5) :: pt_a_water
    real(dp), dimension(5) :: emissivity_water
    real(dp), dimension(5) :: turb_heat_ex_coeff
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /rivertemp1/ &
      albedo_water, &
      pt_a_water, &
      emissivity_water, &
      turb_heat_ex_coeff

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    albedo_water = this%albedo_water
    pt_a_water = this%pt_a_water
    emissivity_water = this%emissivity_water
    turb_heat_ex_coeff = this%turb_heat_ex_coeff

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("rivertemp1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=rivertemp1, iostat=iostat, iomsg=iomsg)
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
    this%albedo_water = albedo_water
    this%pt_a_water = pt_a_water
    this%emissivity_water = emissivity_water
    this%turb_heat_ex_coeff = turb_heat_ex_coeff

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_rivertemp1_from_file

  !> \brief Set rivertemp1 values
  integer function nml_rivertemp1_set(this, &
    albedo_water, &
    pt_a_water, &
    emissivity_water, &
    turb_heat_ex_coeff, &
    errmsg) result(status)

    class(nml_rivertemp1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: albedo_water !< Albedo of open water [-]
    real(dp), dimension(:), intent(in) :: pt_a_water !< Priestley-Taylor coefficient for open water [-]
    real(dp), dimension(:), intent(in) :: emissivity_water !< Emissivity of open water [-]
    real(dp), dimension(:), intent(in) :: turb_heat_ex_coeff !< Turbulent heat exchange coefficient for open water [W m-2 K-1]
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(albedo_water, 1) > size(this%albedo_water, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'albedo_water'"
      return
    end if
    lb__1 = lbound(this%albedo_water, 1)
    ub__1 = lb__1 + size(albedo_water, 1) - 1
    this%albedo_water(lb__1:ub__1) = albedo_water
    if (size(pt_a_water, 1) > size(this%pt_a_water, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pt_a_water'"
      return
    end if
    lb__1 = lbound(this%pt_a_water, 1)
    ub__1 = lb__1 + size(pt_a_water, 1) - 1
    this%pt_a_water(lb__1:ub__1) = pt_a_water
    if (size(emissivity_water, 1) > size(this%emissivity_water, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'emissivity_water'"
      return
    end if
    lb__1 = lbound(this%emissivity_water, 1)
    ub__1 = lb__1 + size(emissivity_water, 1) - 1
    this%emissivity_water(lb__1:ub__1) = emissivity_water
    if (size(turb_heat_ex_coeff, 1) > size(this%turb_heat_ex_coeff, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'turb_heat_ex_coeff'"
      return
    end if
    lb__1 = lbound(this%turb_heat_ex_coeff, 1)
    ub__1 = lb__1 + size(turb_heat_ex_coeff, 1) - 1
    this%turb_heat_ex_coeff(lb__1:ub__1) = turb_heat_ex_coeff

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_rivertemp1_set

  !> \brief Check whether a namelist value was set
  integer function nml_rivertemp1_is_set(this, name, idx, errmsg) result(status)
    class(nml_rivertemp1_t), intent(in) :: this !< namelist instance
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
    case ("albedo_water")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%albedo_water), ubound(this%albedo_water), &
          "albedo_water", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%albedo_water(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%albedo_water))) status = NML_ERR_NOT_SET
      end if
    case ("pt_a_water")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pt_a_water), ubound(this%pt_a_water), &
          "pt_a_water", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pt_a_water(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pt_a_water))) status = NML_ERR_NOT_SET
      end if
    case ("emissivity_water")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%emissivity_water), ubound(this%emissivity_water), &
          "emissivity_water", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%emissivity_water(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%emissivity_water))) status = NML_ERR_NOT_SET
      end if
    case ("turb_heat_ex_coeff")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%turb_heat_ex_coeff), ubound(this%turb_heat_ex_coeff), &
          "turb_heat_ex_coeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%turb_heat_ex_coeff(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%turb_heat_ex_coeff))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_rivertemp1_is_set

  !> \brief Validate required values and constraints
  integer function nml_rivertemp1_is_valid(this, errmsg) result(status)
    class(nml_rivertemp1_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%albedo_water(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: albedo_water"
      return
    end if
    if (any(ieee_is_nan(this%albedo_water(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: albedo_water"
      return
    end if
    if (all(ieee_is_nan(this%pt_a_water(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: pt_a_water"
      return
    end if
    if (any(ieee_is_nan(this%pt_a_water(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: pt_a_water"
      return
    end if
    if (all(ieee_is_nan(this%emissivity_water(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: emissivity_water"
      return
    end if
    if (any(ieee_is_nan(this%emissivity_water(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: emissivity_water"
      return
    end if
    if (all(ieee_is_nan(this%turb_heat_ex_coeff(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: turb_heat_ex_coeff"
      return
    end if
    if (any(ieee_is_nan(this%turb_heat_ex_coeff(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: turb_heat_ex_coeff"
      return
    end if
  end function nml_rivertemp1_is_valid

end module nml_rivertemp1
