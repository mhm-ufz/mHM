!> \file nml_parameter_pet2.f90
!> \copydoc nml_pet2

!> \brief PET - Case 2
!> \details Parameters for PET (case 2 - Priestley-Taylor).
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet2
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

  !> \class nml_pet2_t
  !> \brief PET - Case 2
  !> \details Parameters for PET (case 2 - Priestley-Taylor).
  type, public :: nml_pet2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: priestleytaylorcoeff !< Priestley-Taylor coefficient
    real(dp), dimension(5) :: priestleytaylorlaicorr !< Priestley-Taylor LAI correction factor
  contains
    procedure :: init => nml_pet2_init
    procedure :: from_file => nml_pet2_from_file
    procedure :: set => nml_pet2_set
    procedure :: is_set => nml_pet2_is_set
    procedure :: is_valid => nml_pet2_is_valid
  end type nml_pet2_t

contains

  !> \brief Initialize defaults and sentinels for pet2
  integer function nml_pet2_init(this, errmsg) result(status)
    class(nml_pet2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%priestleytaylorcoeff = ieee_value(this%priestleytaylorcoeff, ieee_quiet_nan) ! sentinel for required real array
    this%priestleytaylorlaicorr = ieee_value(this%priestleytaylorlaicorr, ieee_quiet_nan) ! sentinel for required real array
  end function nml_pet2_init


  !> \brief Read pet2 namelist from file
  integer function nml_pet2_from_file(this, file, errmsg) result(status)
    class(nml_pet2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: priestleytaylorcoeff
    real(dp), dimension(5) :: priestleytaylorlaicorr
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet2/ &
      priestleytaylorcoeff, &
      priestleytaylorlaicorr

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    priestleytaylorcoeff = this%priestleytaylorcoeff
    priestleytaylorlaicorr = this%priestleytaylorlaicorr

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("pet2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet2, iostat=iostat, iomsg=iomsg)
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
    this%priestleytaylorcoeff = priestleytaylorcoeff
    this%priestleytaylorlaicorr = priestleytaylorlaicorr

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet2_from_file

  !> \brief Set pet2 values
  integer function nml_pet2_set(this, &
    priestleytaylorcoeff, &
    priestleytaylorlaicorr, &
    errmsg) result(status)

    class(nml_pet2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: priestleytaylorcoeff !< Priestley-Taylor coefficient
    real(dp), dimension(:), intent(in) :: priestleytaylorlaicorr !< Priestley-Taylor LAI correction factor
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(priestleytaylorcoeff, 1) > size(this%priestleytaylorcoeff, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'priestleytaylorcoeff'"
      return
    end if
    lb__1 = lbound(this%priestleytaylorcoeff, 1)
    ub__1 = lb__1 + size(priestleytaylorcoeff, 1) - 1
    this%priestleytaylorcoeff(lb__1:ub__1) = priestleytaylorcoeff
    if (size(priestleytaylorlaicorr, 1) > size(this%priestleytaylorlaicorr, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'priestleytaylorlaicorr'"
      return
    end if
    lb__1 = lbound(this%priestleytaylorlaicorr, 1)
    ub__1 = lb__1 + size(priestleytaylorlaicorr, 1) - 1
    this%priestleytaylorlaicorr(lb__1:ub__1) = priestleytaylorlaicorr

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet2_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet2_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet2_t), intent(in) :: this !< namelist instance
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
    case ("priestleytaylorcoeff")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%priestleytaylorcoeff), ubound(this%priestleytaylorcoeff), &
          "PriestleyTaylorCoeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%priestleytaylorcoeff(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%priestleytaylorcoeff))) status = NML_ERR_NOT_SET
      end if
    case ("priestleytaylorlaicorr")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%priestleytaylorlaicorr), ubound(this%priestleytaylorlaicorr), &
          "PriestleyTaylorLAIcorr", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%priestleytaylorlaicorr(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%priestleytaylorlaicorr))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet2_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet2_is_valid(this, errmsg) result(status)
    class(nml_pet2_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%priestleytaylorcoeff(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PriestleyTaylorCoeff"
      return
    end if
    if (any(ieee_is_nan(this%priestleytaylorcoeff(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PriestleyTaylorCoeff"
      return
    end if
    if (all(ieee_is_nan(this%priestleytaylorlaicorr(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PriestleyTaylorLAIcorr"
      return
    end if
    if (any(ieee_is_nan(this%priestleytaylorlaicorr(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PriestleyTaylorLAIcorr"
      return
    end if
  end function nml_pet2_is_valid

end module nml_pet2
