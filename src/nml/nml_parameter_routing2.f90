!> \file nml_parameter_routing2.f90
!> \copydoc nml_routing2

!> \brief Routing - Case 2
!> \details Parameters for routing (case 2 - constant celerity).
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_routing2
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

  !> \class nml_routing2_t
  !> \brief Routing - Case 2
  !> \details Parameters for routing (case 2 - constant celerity).
  type, public :: nml_routing2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: streamflow_celerity !< Streamflow celerity
  contains
    procedure :: init => nml_routing2_init
    procedure :: from_file => nml_routing2_from_file
    procedure :: set => nml_routing2_set
    procedure :: is_set => nml_routing2_is_set
    procedure :: is_valid => nml_routing2_is_valid
  end type nml_routing2_t

contains

  !> \brief Initialize defaults and sentinels for routing2
  integer function nml_routing2_init(this, errmsg) result(status)
    class(nml_routing2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%streamflow_celerity = ieee_value(this%streamflow_celerity, ieee_quiet_nan) ! sentinel for required real array
  end function nml_routing2_init


  !> \brief Read routing2 namelist from file
  integer function nml_routing2_from_file(this, file, errmsg) result(status)
    class(nml_routing2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: streamflow_celerity
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /routing2/ &
      streamflow_celerity

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    streamflow_celerity = this%streamflow_celerity

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("routing2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=routing2, iostat=iostat, iomsg=iomsg)
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
    this%streamflow_celerity = streamflow_celerity

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing2_from_file

  !> \brief Set routing2 values
  integer function nml_routing2_set(this, &
    streamflow_celerity, &
    errmsg) result(status)

    class(nml_routing2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: streamflow_celerity !< Streamflow celerity
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(streamflow_celerity, 1) > size(this%streamflow_celerity, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'streamflow_celerity'"
      return
    end if
    lb__1 = lbound(this%streamflow_celerity, 1)
    ub__1 = lb__1 + size(streamflow_celerity, 1) - 1
    this%streamflow_celerity(lb__1:ub__1) = streamflow_celerity

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing2_set

  !> \brief Check whether a namelist value was set
  integer function nml_routing2_is_set(this, name, idx, errmsg) result(status)
    class(nml_routing2_t), intent(in) :: this !< namelist instance
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
    case ("streamflow_celerity")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%streamflow_celerity), ubound(this%streamflow_celerity), &
          "streamflow_celerity", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%streamflow_celerity(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%streamflow_celerity))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_routing2_is_set

  !> \brief Validate required values and constraints
  integer function nml_routing2_is_valid(this, errmsg) result(status)
    class(nml_routing2_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%streamflow_celerity(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: streamflow_celerity"
      return
    end if
    if (any(ieee_is_nan(this%streamflow_celerity(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: streamflow_celerity"
      return
    end if
  end function nml_routing2_is_valid

end module nml_routing2
