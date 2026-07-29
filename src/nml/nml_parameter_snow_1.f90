!> \file nml_parameter_snow_1.f90
!> \copydoc nml_snow_1

!> \brief Snow - Case 1
!> \details Parameters for the degree-day snow module.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_snow_1
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

  !> \class nml_snow_1_t
  !> \brief Snow - Case 1
  !> \details Parameters for the degree-day snow module.
  type, public :: nml_snow_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: snow_threshold_temperature !< Threshold for rain and snow partitioning [degC]
    type(parameter_t) :: degree_day_factor_forest !< Degree-day factor for forest [m degC-1]
    type(parameter_t) :: degree_day_factor_impervious !< Degree-day factor for impervious areas [m degC-1]
    type(parameter_t) :: degree_day_factor_pervious !< Degree-day factor for pervious areas [m degC-1]
    type(parameter_t) :: degree_day_factor_precipitation !< Precipitation-dependent degree-day factor increase [degC-1]
    type(parameter_t) :: max_degree_day_factor_forest !< Maximum degree-day factor for forest [m degC-1]
    type(parameter_t) :: max_degree_day_factor_impervious !< Maximum degree-day factor for impervious areas [m degC-1]
    type(parameter_t) :: max_degree_day_factor_pervious !< Maximum degree-day factor for pervious areas [m degC-1]
  contains
    procedure :: init => nml_snow_1_init
    procedure :: init_type => nml_snow_1_init_type
    procedure :: from_file => nml_snow_1_from_file
    procedure :: set => nml_snow_1_set
    procedure :: is_set => nml_snow_1_is_set
    procedure :: is_valid => nml_snow_1_is_valid
  end type nml_snow_1_t

contains

  !> \brief Initialize defaults and sentinels for snow_1
  integer function nml_snow_1_init(this, errmsg) result(status)
    class(nml_snow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      snow_threshold_temperature=this%snow_threshold_temperature, &
      degree_day_factor_forest=this%degree_day_factor_forest, &
      degree_day_factor_impervious=this%degree_day_factor_impervious, &
      degree_day_factor_pervious=this%degree_day_factor_pervious, &
      degree_day_factor_precipitation=this%degree_day_factor_precipitation, &
      max_degree_day_factor_forest=this%max_degree_day_factor_forest, &
      max_degree_day_factor_impervious=this%max_degree_day_factor_impervious, &
      max_degree_day_factor_pervious=this%max_degree_day_factor_pervious, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_snow_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_snow_1_init_type(this, &
    snow_threshold_temperature, &
    degree_day_factor_forest, &
    degree_day_factor_impervious, &
    degree_day_factor_pervious, &
    degree_day_factor_precipitation, &
    max_degree_day_factor_forest, &
    max_degree_day_factor_impervious, &
    max_degree_day_factor_pervious, &
    errmsg) result(status)
    class(nml_snow_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: snow_threshold_temperature !< Threshold for rain and snow partitioning [degC]
    type(parameter_t), intent(inout), optional :: degree_day_factor_forest !< Degree-day factor for forest [m degC-1]
    type(parameter_t), intent(inout), optional :: degree_day_factor_impervious !< Degree-day factor for impervious areas [m degC-1]
    type(parameter_t), intent(inout), optional :: degree_day_factor_pervious !< Degree-day factor for pervious areas [m degC-1]
    type(parameter_t), intent(inout), optional :: degree_day_factor_precipitation !< Precipitation-dependent degree-day factor increase [degC-1]
    type(parameter_t), intent(inout), optional :: max_degree_day_factor_forest !< Maximum degree-day factor for forest [m degC-1]
    type(parameter_t), intent(inout), optional :: max_degree_day_factor_impervious !< Maximum degree-day factor for impervious areas [m degC-1]
    type(parameter_t), intent(inout), optional :: max_degree_day_factor_pervious !< Maximum degree-day factor for pervious areas [m degC-1]
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(snow_threshold_temperature)) then
      snow_threshold_temperature%value = ieee_value(snow_threshold_temperature%value, ieee_quiet_nan) ! sentinel for derived component value
      snow_threshold_temperature%optimize = .false.
      snow_threshold_temperature%lower_bound = ieee_value(snow_threshold_temperature%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      snow_threshold_temperature%upper_bound = ieee_value(snow_threshold_temperature%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      snow_threshold_temperature%lower_bound = -2.0_dp
      snow_threshold_temperature%upper_bound = 2.0_dp
    end if
    if (present(degree_day_factor_forest)) then
      degree_day_factor_forest%value = ieee_value(degree_day_factor_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      degree_day_factor_forest%optimize = .false.
      degree_day_factor_forest%lower_bound = ieee_value(degree_day_factor_forest%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      degree_day_factor_forest%upper_bound = ieee_value(degree_day_factor_forest%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      degree_day_factor_forest%lower_bound = 0.0001_dp
      degree_day_factor_forest%upper_bound = 4.0_dp
    end if
    if (present(degree_day_factor_impervious)) then
      degree_day_factor_impervious%value = ieee_value(degree_day_factor_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      degree_day_factor_impervious%optimize = .false.
      degree_day_factor_impervious%lower_bound = ieee_value(degree_day_factor_impervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      degree_day_factor_impervious%upper_bound = ieee_value(degree_day_factor_impervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      degree_day_factor_impervious%lower_bound = 0.0_dp
      degree_day_factor_impervious%upper_bound = 1.0_dp
    end if
    if (present(degree_day_factor_pervious)) then
      degree_day_factor_pervious%value = ieee_value(degree_day_factor_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      degree_day_factor_pervious%optimize = .false.
      degree_day_factor_pervious%lower_bound = ieee_value(degree_day_factor_pervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      degree_day_factor_pervious%upper_bound = ieee_value(degree_day_factor_pervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      degree_day_factor_pervious%lower_bound = 0.0_dp
      degree_day_factor_pervious%upper_bound = 2.0_dp
    end if
    if (present(degree_day_factor_precipitation)) then
      degree_day_factor_precipitation%value = ieee_value(degree_day_factor_precipitation%value, ieee_quiet_nan) ! sentinel for derived component value
      degree_day_factor_precipitation%optimize = .false.
      degree_day_factor_precipitation%lower_bound = ieee_value(degree_day_factor_precipitation%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      degree_day_factor_precipitation%upper_bound = ieee_value(degree_day_factor_precipitation%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      degree_day_factor_precipitation%lower_bound = 0.1_dp
      degree_day_factor_precipitation%upper_bound = 0.9_dp
    end if
    if (present(max_degree_day_factor_forest)) then
      max_degree_day_factor_forest%value = ieee_value(max_degree_day_factor_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      max_degree_day_factor_forest%optimize = .false.
      max_degree_day_factor_forest%lower_bound = ieee_value(max_degree_day_factor_forest%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      max_degree_day_factor_forest%upper_bound = ieee_value(max_degree_day_factor_forest%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      max_degree_day_factor_forest%lower_bound = 0.0_dp
      max_degree_day_factor_forest%upper_bound = 8.0_dp
    end if
    if (present(max_degree_day_factor_impervious)) then
      max_degree_day_factor_impervious%value = ieee_value(max_degree_day_factor_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      max_degree_day_factor_impervious%optimize = .false.
      max_degree_day_factor_impervious%lower_bound = ieee_value(max_degree_day_factor_impervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      max_degree_day_factor_impervious%upper_bound = ieee_value(max_degree_day_factor_impervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      max_degree_day_factor_impervious%lower_bound = 0.0_dp
      max_degree_day_factor_impervious%upper_bound = 8.0_dp
    end if
    if (present(max_degree_day_factor_pervious)) then
      max_degree_day_factor_pervious%value = ieee_value(max_degree_day_factor_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      max_degree_day_factor_pervious%optimize = .false.
      max_degree_day_factor_pervious%lower_bound = ieee_value(max_degree_day_factor_pervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      max_degree_day_factor_pervious%upper_bound = ieee_value(max_degree_day_factor_pervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      max_degree_day_factor_pervious%lower_bound = 0.0_dp
      max_degree_day_factor_pervious%upper_bound = 8.0_dp
    end if
  end function nml_snow_1_init_type


  !> \brief Read snow_1 namelist from file
  integer function nml_snow_1_from_file(this, file, errmsg) result(status)
    class(nml_snow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: snow_threshold_temperature
    type(parameter_t) :: degree_day_factor_forest
    type(parameter_t) :: degree_day_factor_impervious
    type(parameter_t) :: degree_day_factor_pervious
    type(parameter_t) :: degree_day_factor_precipitation
    type(parameter_t) :: max_degree_day_factor_forest
    type(parameter_t) :: max_degree_day_factor_impervious
    type(parameter_t) :: max_degree_day_factor_pervious
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /snow_1/ &
      snow_threshold_temperature, &
      degree_day_factor_forest, &
      degree_day_factor_impervious, &
      degree_day_factor_pervious, &
      degree_day_factor_precipitation, &
      max_degree_day_factor_forest, &
      max_degree_day_factor_impervious, &
      max_degree_day_factor_pervious

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    snow_threshold_temperature = this%snow_threshold_temperature
    degree_day_factor_forest = this%degree_day_factor_forest
    degree_day_factor_impervious = this%degree_day_factor_impervious
    degree_day_factor_pervious = this%degree_day_factor_pervious
    degree_day_factor_precipitation = this%degree_day_factor_precipitation
    max_degree_day_factor_forest = this%max_degree_day_factor_forest
    max_degree_day_factor_impervious = this%max_degree_day_factor_impervious
    max_degree_day_factor_pervious = this%max_degree_day_factor_pervious

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("snow_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=snow_1, iostat=iostat, iomsg=iomsg)
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
    this%snow_threshold_temperature = snow_threshold_temperature
    this%degree_day_factor_forest = degree_day_factor_forest
    this%degree_day_factor_impervious = degree_day_factor_impervious
    this%degree_day_factor_pervious = degree_day_factor_pervious
    this%degree_day_factor_precipitation = degree_day_factor_precipitation
    this%max_degree_day_factor_forest = max_degree_day_factor_forest
    this%max_degree_day_factor_impervious = max_degree_day_factor_impervious
    this%max_degree_day_factor_pervious = max_degree_day_factor_pervious

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_snow_1_from_file

  !> \brief Set snow_1 values
  integer function nml_snow_1_set(this, &
    snow_threshold_temperature, &
    degree_day_factor_forest, &
    degree_day_factor_impervious, &
    degree_day_factor_pervious, &
    degree_day_factor_precipitation, &
    max_degree_day_factor_forest, &
    max_degree_day_factor_impervious, &
    max_degree_day_factor_pervious, &
    errmsg) result(status)

    class(nml_snow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: snow_threshold_temperature !< Threshold for rain and snow partitioning [degC]
    type(parameter_t), intent(in) :: degree_day_factor_forest !< Degree-day factor for forest [m degC-1]
    type(parameter_t), intent(in) :: degree_day_factor_impervious !< Degree-day factor for impervious areas [m degC-1]
    type(parameter_t), intent(in) :: degree_day_factor_pervious !< Degree-day factor for pervious areas [m degC-1]
    type(parameter_t), intent(in) :: degree_day_factor_precipitation !< Precipitation-dependent degree-day factor increase [degC-1]
    type(parameter_t), intent(in) :: max_degree_day_factor_forest !< Maximum degree-day factor for forest [m degC-1]
    type(parameter_t), intent(in) :: max_degree_day_factor_impervious !< Maximum degree-day factor for impervious areas [m degC-1]
    type(parameter_t), intent(in) :: max_degree_day_factor_pervious !< Maximum degree-day factor for pervious areas [m degC-1]

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%snow_threshold_temperature = snow_threshold_temperature
    this%degree_day_factor_forest = degree_day_factor_forest
    this%degree_day_factor_impervious = degree_day_factor_impervious
    this%degree_day_factor_pervious = degree_day_factor_pervious
    this%degree_day_factor_precipitation = degree_day_factor_precipitation
    this%max_degree_day_factor_forest = max_degree_day_factor_forest
    this%max_degree_day_factor_impervious = max_degree_day_factor_impervious
    this%max_degree_day_factor_pervious = max_degree_day_factor_pervious

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_snow_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_snow_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_snow_1_t), intent(in) :: this !< namelist instance
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
    case ("snow_threshold_temperature%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow_threshold_temperature'"
        return
      end if
      if (ieee_is_nan(this%snow_threshold_temperature%value)) status = NML_ERR_NOT_SET
    case ("snow_threshold_temperature%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow_threshold_temperature'"
        return
      end if
    case ("snow_threshold_temperature%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow_threshold_temperature'"
        return
      end if
    case ("snow_threshold_temperature%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow_threshold_temperature'"
        return
      end if
    case ("snow_threshold_temperature")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow_threshold_temperature'"
        return
      end if
      if (ieee_is_nan(this%snow_threshold_temperature%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("degree_day_factor_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_forest'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_forest%value)) status = NML_ERR_NOT_SET
    case ("degree_day_factor_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_forest'"
        return
      end if
    case ("degree_day_factor_forest%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_forest'"
        return
      end if
    case ("degree_day_factor_forest%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_forest'"
        return
      end if
    case ("degree_day_factor_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_forest'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("degree_day_factor_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_impervious'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_impervious%value)) status = NML_ERR_NOT_SET
    case ("degree_day_factor_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_impervious'"
        return
      end if
    case ("degree_day_factor_impervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_impervious'"
        return
      end if
    case ("degree_day_factor_impervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_impervious'"
        return
      end if
    case ("degree_day_factor_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_impervious'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("degree_day_factor_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_pervious'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_pervious%value)) status = NML_ERR_NOT_SET
    case ("degree_day_factor_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_pervious'"
        return
      end if
    case ("degree_day_factor_pervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_pervious'"
        return
      end if
    case ("degree_day_factor_pervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_pervious'"
        return
      end if
    case ("degree_day_factor_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_pervious'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("degree_day_factor_precipitation%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_precipitation'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_precipitation%value)) status = NML_ERR_NOT_SET
    case ("degree_day_factor_precipitation%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_precipitation'"
        return
      end if
    case ("degree_day_factor_precipitation%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_precipitation'"
        return
      end if
    case ("degree_day_factor_precipitation%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_precipitation'"
        return
      end if
    case ("degree_day_factor_precipitation")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'degree_day_factor_precipitation'"
        return
      end if
      if (ieee_is_nan(this%degree_day_factor_precipitation%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("max_degree_day_factor_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_forest'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_forest%value)) status = NML_ERR_NOT_SET
    case ("max_degree_day_factor_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_forest'"
        return
      end if
    case ("max_degree_day_factor_forest%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_forest'"
        return
      end if
    case ("max_degree_day_factor_forest%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_forest'"
        return
      end if
    case ("max_degree_day_factor_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_forest'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("max_degree_day_factor_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_impervious'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_impervious%value)) status = NML_ERR_NOT_SET
    case ("max_degree_day_factor_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_impervious'"
        return
      end if
    case ("max_degree_day_factor_impervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_impervious'"
        return
      end if
    case ("max_degree_day_factor_impervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_impervious'"
        return
      end if
    case ("max_degree_day_factor_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_impervious'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("max_degree_day_factor_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_pervious'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_pervious%value)) status = NML_ERR_NOT_SET
    case ("max_degree_day_factor_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_pervious'"
        return
      end if
    case ("max_degree_day_factor_pervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_pervious'"
        return
      end if
    case ("max_degree_day_factor_pervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_pervious'"
        return
      end if
    case ("max_degree_day_factor_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_degree_day_factor_pervious'"
        return
      end if
      if (ieee_is_nan(this%max_degree_day_factor_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_snow_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_snow_1_is_valid(this, errmsg) result(status)
    class(nml_snow_1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("snow_threshold_temperature", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: snow_threshold_temperature"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("degree_day_factor_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: degree_day_factor_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("degree_day_factor_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: degree_day_factor_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("degree_day_factor_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: degree_day_factor_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("degree_day_factor_precipitation", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: degree_day_factor_precipitation"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("max_degree_day_factor_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: max_degree_day_factor_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("max_degree_day_factor_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: max_degree_day_factor_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("max_degree_day_factor_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: max_degree_day_factor_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_snow_1_is_valid

end module nml_snow_1
