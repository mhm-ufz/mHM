!> \file nml_parameter_routing_1.f90
!> \copydoc nml_routing_1

!> \brief Routing - Case 1
!> \details Parameters for Muskingum routing.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_routing_1
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

  !> \class nml_routing_1_t
  !> \brief Routing - Case 1
  !> \details Parameters for Muskingum routing.
  type, public :: nml_routing_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: travel_time_constant !< Muskingum travel-time constant
    type(parameter_t) :: travel_time_river_length !< River-length contribution to Muskingum travel time
    type(parameter_t) :: travel_time_river_slope !< River-slope contribution to Muskingum travel time
    type(parameter_t) :: travel_time_impervious !< Impervious-area contribution to Muskingum travel time
    type(parameter_t) :: attenuation_river_slope !< River-slope contribution to Muskingum attenuation
  contains
    procedure :: init => nml_routing_1_init
    procedure :: init_type => nml_routing_1_init_type
    procedure :: from_file => nml_routing_1_from_file
    procedure :: set => nml_routing_1_set
    procedure :: is_set => nml_routing_1_is_set
    procedure :: is_valid => nml_routing_1_is_valid
  end type nml_routing_1_t

contains

  !> \brief Initialize defaults and sentinels for routing_1
  integer function nml_routing_1_init(this, errmsg) result(status)
    class(nml_routing_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      travel_time_constant=this%travel_time_constant, &
      travel_time_river_length=this%travel_time_river_length, &
      travel_time_river_slope=this%travel_time_river_slope, &
      travel_time_impervious=this%travel_time_impervious, &
      attenuation_river_slope=this%attenuation_river_slope, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_routing_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_routing_1_init_type(this, &
    travel_time_constant, &
    travel_time_river_length, &
    travel_time_river_slope, &
    travel_time_impervious, &
    attenuation_river_slope, &
    errmsg) result(status)
    class(nml_routing_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: travel_time_constant !< Muskingum travel-time constant
    type(parameter_t), intent(inout), optional :: travel_time_river_length !< River-length contribution to Muskingum travel time
    type(parameter_t), intent(inout), optional :: travel_time_river_slope !< River-slope contribution to Muskingum travel time
    type(parameter_t), intent(inout), optional :: travel_time_impervious !< Impervious-area contribution to Muskingum travel time
    type(parameter_t), intent(inout), optional :: attenuation_river_slope !< River-slope contribution to Muskingum attenuation
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(travel_time_constant)) then
      travel_time_constant%value = ieee_value(travel_time_constant%value, ieee_quiet_nan) ! sentinel for derived component value
      travel_time_constant%optimize = .false.
      travel_time_constant%lower_bound = ieee_value(travel_time_constant%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      travel_time_constant%upper_bound = ieee_value(travel_time_constant%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      travel_time_constant%lower_bound = 0.31_dp
      travel_time_constant%upper_bound = 0.35_dp
    end if
    if (present(travel_time_river_length)) then
      travel_time_river_length%value = ieee_value(travel_time_river_length%value, ieee_quiet_nan) ! sentinel for derived component value
      travel_time_river_length%optimize = .false.
      travel_time_river_length%lower_bound = ieee_value(travel_time_river_length%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      travel_time_river_length%upper_bound = ieee_value(travel_time_river_length%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      travel_time_river_length%lower_bound = 0.07_dp
      travel_time_river_length%upper_bound = 0.08_dp
    end if
    if (present(travel_time_river_slope)) then
      travel_time_river_slope%value = ieee_value(travel_time_river_slope%value, ieee_quiet_nan) ! sentinel for derived component value
      travel_time_river_slope%optimize = .false.
      travel_time_river_slope%lower_bound = ieee_value(travel_time_river_slope%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      travel_time_river_slope%upper_bound = ieee_value(travel_time_river_slope%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      travel_time_river_slope%lower_bound = 1.95_dp
      travel_time_river_slope%upper_bound = 2.1_dp
    end if
    if (present(travel_time_impervious)) then
      travel_time_impervious%value = ieee_value(travel_time_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      travel_time_impervious%optimize = .false.
      travel_time_impervious%lower_bound = ieee_value(travel_time_impervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      travel_time_impervious%upper_bound = ieee_value(travel_time_impervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      travel_time_impervious%lower_bound = 0.09_dp
      travel_time_impervious%upper_bound = 0.11_dp
    end if
    if (present(attenuation_river_slope)) then
      attenuation_river_slope%value = ieee_value(attenuation_river_slope%value, ieee_quiet_nan) ! sentinel for derived component value
      attenuation_river_slope%optimize = .false.
      attenuation_river_slope%lower_bound = ieee_value(attenuation_river_slope%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      attenuation_river_slope%upper_bound = ieee_value(attenuation_river_slope%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      attenuation_river_slope%lower_bound = 0.01_dp
      attenuation_river_slope%upper_bound = 0.5_dp
    end if
  end function nml_routing_1_init_type


  !> \brief Read routing_1 namelist from file
  integer function nml_routing_1_from_file(this, file, errmsg) result(status)
    class(nml_routing_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: travel_time_constant
    type(parameter_t) :: travel_time_river_length
    type(parameter_t) :: travel_time_river_slope
    type(parameter_t) :: travel_time_impervious
    type(parameter_t) :: attenuation_river_slope
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /routing_1/ &
      travel_time_constant, &
      travel_time_river_length, &
      travel_time_river_slope, &
      travel_time_impervious, &
      attenuation_river_slope

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    travel_time_constant = this%travel_time_constant
    travel_time_river_length = this%travel_time_river_length
    travel_time_river_slope = this%travel_time_river_slope
    travel_time_impervious = this%travel_time_impervious
    attenuation_river_slope = this%attenuation_river_slope

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("routing_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=routing_1, iostat=iostat, iomsg=iomsg)
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
    this%travel_time_constant = travel_time_constant
    this%travel_time_river_length = travel_time_river_length
    this%travel_time_river_slope = travel_time_river_slope
    this%travel_time_impervious = travel_time_impervious
    this%attenuation_river_slope = attenuation_river_slope

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing_1_from_file

  !> \brief Set routing_1 values
  integer function nml_routing_1_set(this, &
    travel_time_constant, &
    travel_time_river_length, &
    travel_time_river_slope, &
    travel_time_impervious, &
    attenuation_river_slope, &
    errmsg) result(status)

    class(nml_routing_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: travel_time_constant !< Muskingum travel-time constant
    type(parameter_t), intent(in) :: travel_time_river_length !< River-length contribution to Muskingum travel time
    type(parameter_t), intent(in) :: travel_time_river_slope !< River-slope contribution to Muskingum travel time
    type(parameter_t), intent(in) :: travel_time_impervious !< Impervious-area contribution to Muskingum travel time
    type(parameter_t), intent(in) :: attenuation_river_slope !< River-slope contribution to Muskingum attenuation

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%travel_time_constant = travel_time_constant
    this%travel_time_river_length = travel_time_river_length
    this%travel_time_river_slope = travel_time_river_slope
    this%travel_time_impervious = travel_time_impervious
    this%attenuation_river_slope = attenuation_river_slope

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_routing_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_routing_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_routing_1_t), intent(in) :: this !< namelist instance
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
    case ("travel_time_constant%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_constant'"
        return
      end if
      if (ieee_is_nan(this%travel_time_constant%value)) status = NML_ERR_NOT_SET
    case ("travel_time_constant%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_constant'"
        return
      end if
    case ("travel_time_constant%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_constant'"
        return
      end if
    case ("travel_time_constant%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_constant'"
        return
      end if
    case ("travel_time_constant")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_constant'"
        return
      end if
      if (ieee_is_nan(this%travel_time_constant%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("travel_time_river_length%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_length'"
        return
      end if
      if (ieee_is_nan(this%travel_time_river_length%value)) status = NML_ERR_NOT_SET
    case ("travel_time_river_length%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_length'"
        return
      end if
    case ("travel_time_river_length%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_length'"
        return
      end if
    case ("travel_time_river_length%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_length'"
        return
      end if
    case ("travel_time_river_length")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_length'"
        return
      end if
      if (ieee_is_nan(this%travel_time_river_length%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("travel_time_river_slope%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_slope'"
        return
      end if
      if (ieee_is_nan(this%travel_time_river_slope%value)) status = NML_ERR_NOT_SET
    case ("travel_time_river_slope%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_slope'"
        return
      end if
    case ("travel_time_river_slope%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_slope'"
        return
      end if
    case ("travel_time_river_slope%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_slope'"
        return
      end if
    case ("travel_time_river_slope")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_river_slope'"
        return
      end if
      if (ieee_is_nan(this%travel_time_river_slope%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("travel_time_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_impervious'"
        return
      end if
      if (ieee_is_nan(this%travel_time_impervious%value)) status = NML_ERR_NOT_SET
    case ("travel_time_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_impervious'"
        return
      end if
    case ("travel_time_impervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_impervious'"
        return
      end if
    case ("travel_time_impervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_impervious'"
        return
      end if
    case ("travel_time_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'travel_time_impervious'"
        return
      end if
      if (ieee_is_nan(this%travel_time_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("attenuation_river_slope%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'attenuation_river_slope'"
        return
      end if
      if (ieee_is_nan(this%attenuation_river_slope%value)) status = NML_ERR_NOT_SET
    case ("attenuation_river_slope%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'attenuation_river_slope'"
        return
      end if
    case ("attenuation_river_slope%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'attenuation_river_slope'"
        return
      end if
    case ("attenuation_river_slope%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'attenuation_river_slope'"
        return
      end if
    case ("attenuation_river_slope")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'attenuation_river_slope'"
        return
      end if
      if (ieee_is_nan(this%attenuation_river_slope%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_routing_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_routing_1_is_valid(this, errmsg) result(status)
    class(nml_routing_1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("travel_time_constant", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: travel_time_constant"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("travel_time_river_length", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: travel_time_river_length"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("travel_time_river_slope", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: travel_time_river_slope"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("travel_time_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: travel_time_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("attenuation_river_slope", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: attenuation_river_slope"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_routing_1_is_valid

end module nml_routing_1
