!> \file nml_parameter_interception1.f90
!> \copydoc nml_interception1

!> \brief Interception - Case 1
!> \details Parameters for interception.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_interception1
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

  !> \class nml_interception1_t
  !> \brief Interception - Case 1
  !> \details Parameters for interception.
  type, public :: nml_interception1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: canopyinterceptionfactor !< Multiplier to relate LAI to interception storage
  contains
    procedure :: init => nml_interception1_init
    procedure :: from_file => nml_interception1_from_file
    procedure :: set => nml_interception1_set
    procedure :: is_set => nml_interception1_is_set
    procedure :: is_valid => nml_interception1_is_valid
  end type nml_interception1_t

contains

  !> \brief Initialize defaults and sentinels for interception1
  integer function nml_interception1_init(this, errmsg) result(status)
    class(nml_interception1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%canopyinterceptionfactor = ieee_value(this%canopyinterceptionfactor, ieee_quiet_nan) ! sentinel for required real array
  end function nml_interception1_init


  !> \brief Read interception1 namelist from file
  integer function nml_interception1_from_file(this, file, errmsg) result(status)
    class(nml_interception1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: canopyinterceptionfactor
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /interception1/ &
      canopyinterceptionfactor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    canopyinterceptionfactor = this%canopyinterceptionfactor

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("interception1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=interception1, iostat=iostat, iomsg=iomsg)
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
    this%canopyinterceptionfactor = canopyinterceptionfactor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_interception1_from_file

  !> \brief Set interception1 values
  integer function nml_interception1_set(this, &
    canopyinterceptionfactor, &
    errmsg) result(status)

    class(nml_interception1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: canopyinterceptionfactor !< Multiplier to relate LAI to interception storage
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(canopyinterceptionfactor, 1) > size(this%canopyinterceptionfactor, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'canopyinterceptionfactor'"
      return
    end if
    lb__1 = lbound(this%canopyinterceptionfactor, 1)
    ub__1 = lb__1 + size(canopyinterceptionfactor, 1) - 1
    this%canopyinterceptionfactor(lb__1:ub__1) = canopyinterceptionfactor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_interception1_set

  !> \brief Check whether a namelist value was set
  integer function nml_interception1_is_set(this, name, idx, errmsg) result(status)
    class(nml_interception1_t), intent(in) :: this !< namelist instance
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
    case ("canopyinterceptionfactor")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%canopyinterceptionfactor), ubound(this%canopyinterceptionfactor), &
          "canopyInterceptionFactor", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%canopyinterceptionfactor(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%canopyinterceptionfactor))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_interception1_is_set

  !> \brief Validate required values and constraints
  integer function nml_interception1_is_valid(this, errmsg) result(status)
    class(nml_interception1_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%canopyinterceptionfactor(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: canopyInterceptionFactor"
      return
    end if
    if (any(ieee_is_nan(this%canopyinterceptionfactor(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: canopyInterceptionFactor"
      return
    end if
  end function nml_interception1_is_valid

end module nml_interception1
