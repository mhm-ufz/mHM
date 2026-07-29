!> \file nml_parameter_routing_2.f90
!> \copydoc nml_routing_2

!> \brief Routing - Case 2
!> \details Parameters for constant-celerity routing.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_routing_2
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
    parameter_t
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  !> \class nml_routing_2_t
  !> \brief Routing - Case 2
  !> \details Parameters for constant-celerity routing.
  type, public :: nml_routing_2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: streamflow_celerity !< Streamflow celerity
  contains
    procedure :: init => nml_routing_2_init
    procedure :: init_type => nml_routing_2_init_type
    procedure :: from_file => nml_routing_2_from_file
    procedure :: set => nml_routing_2_set
    procedure :: is_set => nml_routing_2_is_set
    procedure :: is_valid => nml_routing_2_is_valid
  end type nml_routing_2_t

contains

  !> \brief Initialize defaults and sentinels for routing_2
  integer function nml_routing_2_init(this, errmsg) result(status)
    class(nml_routing_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      streamflow_celerity=this%streamflow_celerity, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_routing_2_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_routing_2_init_type(this, &
    streamflow_celerity, &
    errmsg) result(status)
    class(nml_routing_2_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: streamflow_celerity !< Streamflow celerity
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(streamflow_celerity)) then
      streamflow_celerity%value = ieee_value(streamflow_celerity%value, ieee_quiet_nan) ! sentinel for derived component value
      streamflow_celerity%optimize = .false.
      streamflow_celerity%lower_bound = ieee_value(streamflow_celerity%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      streamflow_celerity%upper_bound = ieee_value(streamflow_celerity%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      streamflow_celerity%lower_bound = 0.1_dp
      streamflow_celerity%upper_bound = 15.0_dp
    end if
  end function nml_routing_2_init_type


  !> \brief Read routing_2 namelist from file
  integer function nml_routing_2_from_file(this, file, errmsg) result(status)
    class(nml_routing_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: streamflow_celerity
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /routing_2/ &
      streamflow_celerity

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    streamflow_celerity = this%streamflow_celerity

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("routing_2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=routing_2, iostat=iostat, iomsg=iomsg)
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
  end function nml_routing_2_from_file

  !> \brief Set routing_2 values
  integer function nml_routing_2_set(this, &
    streamflow_celerity, &
    errmsg) result(status)

    class(nml_routing_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: streamflow_celerity !< Streamflow celerity

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%streamflow_celerity = streamflow_celerity

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing_2_set

  !> \brief Check whether a namelist value was set
  integer function nml_routing_2_is_set(this, name, idx, errmsg) result(status)
    class(nml_routing_2_t), intent(in) :: this !< namelist instance
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
    case ("streamflow_celerity%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'streamflow_celerity'"
        return
      end if
      if (ieee_is_nan(this%streamflow_celerity%value)) status = NML_ERR_NOT_SET
    case ("streamflow_celerity%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'streamflow_celerity'"
        return
      end if
    case ("streamflow_celerity%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'streamflow_celerity'"
        return
      end if
    case ("streamflow_celerity%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'streamflow_celerity'"
        return
      end if
    case ("streamflow_celerity")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'streamflow_celerity'"
        return
      end if
      if (ieee_is_nan(this%streamflow_celerity%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_routing_2_is_set

  !> \brief Validate required values and constraints
  integer function nml_routing_2_is_valid(this, errmsg) result(status)
    class(nml_routing_2_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required parameters
    istat = this%is_set("streamflow_celerity", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: streamflow_celerity"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_routing_2_is_valid

end module nml_routing_2
