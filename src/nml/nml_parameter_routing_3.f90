!> \file nml_parameter_routing_3.f90
!> \copydoc nml_routing_3

!> \brief Routing - Case 3
!> \details Parameters for varying-celerity routing.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_routing_3
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
    to_lower
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp
  use mo_parameter_types, only: parameter_t

  implicit none

  !> \class nml_routing_3_t
  !> \brief Routing - Case 3
  !> \details Parameters for varying-celerity routing.
  type, public :: nml_routing_3_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: slope_factor !< Slope factor
  contains
    procedure :: init => nml_routing_3_init
    procedure :: init_type => nml_routing_3_init_type
    procedure :: from_file => nml_routing_3_from_file
    procedure :: set => nml_routing_3_set
    procedure :: is_set => nml_routing_3_is_set
    procedure :: is_valid => nml_routing_3_is_valid
  end type nml_routing_3_t

contains

  !> \brief Initialize defaults and sentinels for routing_3
  integer function nml_routing_3_init(this, errmsg) result(status)
    class(nml_routing_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      slope_factor=this%slope_factor, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_routing_3_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_routing_3_init_type(this, &
    slope_factor, &
    errmsg) result(status)
    class(nml_routing_3_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: slope_factor !< Slope factor
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(slope_factor)) then
      slope_factor%value = ieee_value(slope_factor%value, ieee_quiet_nan) ! sentinel for derived component value
      slope_factor%optimize = .false.
      slope_factor%min = ieee_value(slope_factor%min, ieee_quiet_nan) ! sentinel for derived component min
      slope_factor%max = ieee_value(slope_factor%max, ieee_quiet_nan) ! sentinel for derived component max
      slope_factor%min = 0.1_dp
      slope_factor%max = 100.0_dp
    end if
  end function nml_routing_3_init_type


  !> \brief Read routing_3 namelist from file
  integer function nml_routing_3_from_file(this, file, errmsg) result(status)
    class(nml_routing_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: slope_factor
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /routing_3/ &
      slope_factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    slope_factor = this%slope_factor

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("routing_3", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=routing_3, iostat=iostat, iomsg=iomsg)
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
    this%slope_factor = slope_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing_3_from_file

  !> \brief Set routing_3 values
  integer function nml_routing_3_set(this, &
    slope_factor, &
    errmsg) result(status)

    class(nml_routing_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: slope_factor !< Slope factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%slope_factor = slope_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing_3_set

  !> \brief Check whether a namelist value was set
  integer function nml_routing_3_is_set(this, name, idx, errmsg) result(status)
    class(nml_routing_3_t), intent(in) :: this !< namelist instance
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
    case ("slope_factor%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slope_factor'"
        return
      end if
      if (ieee_is_nan(this%slope_factor%value)) status = NML_ERR_NOT_SET
    case ("slope_factor%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slope_factor'"
        return
      end if
    case ("slope_factor%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slope_factor'"
        return
      end if
    case ("slope_factor%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slope_factor'"
        return
      end if
    case ("slope_factor")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slope_factor'"
        return
      end if
      if (ieee_is_nan(this%slope_factor%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_routing_3_is_set

  !> \brief Validate required values and constraints
  integer function nml_routing_3_is_valid(this, errmsg) result(status)
    class(nml_routing_3_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("slope_factor", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: slope_factor"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_routing_3_is_valid

end module nml_routing_3
