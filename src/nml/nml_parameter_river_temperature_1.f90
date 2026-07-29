!> \file nml_parameter_river_temperature_1.f90
!> \copydoc nml_river_temperature_1

!> \brief River temperature - Case 1
!> \details Parameters for river-temperature routing.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_river_temperature_1
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

  !> \class nml_river_temperature_1_t
  !> \brief River temperature - Case 1
  !> \details Parameters for river-temperature routing.
  type, public :: nml_river_temperature_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: albedo_water !< Albedo of open water [-]
    type(parameter_t) :: pt_a_water !< Priestley-Taylor coefficient for open water [-]
    type(parameter_t) :: emissivity_water !< Emissivity of open water [-]
    type(parameter_t) :: turbulent_heat_exchange_coefficient !< Turbulent heat-exchange coefficient [W m-2 K-1]
  contains
    procedure :: init => nml_river_temperature_1_init
    procedure :: init_type => nml_river_temperature_1_init_type
    procedure :: from_file => nml_river_temperature_1_from_file
    procedure :: set => nml_river_temperature_1_set
    procedure :: is_set => nml_river_temperature_1_is_set
    procedure :: is_valid => nml_river_temperature_1_is_valid
  end type nml_river_temperature_1_t

contains

  !> \brief Initialize defaults and sentinels for river_temperature_1
  integer function nml_river_temperature_1_init(this, errmsg) result(status)
    class(nml_river_temperature_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      albedo_water=this%albedo_water, &
      pt_a_water=this%pt_a_water, &
      emissivity_water=this%emissivity_water, &
      turbulent_heat_exchange_coefficient=this%turbulent_heat_exchange_coefficient, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_river_temperature_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_river_temperature_1_init_type(this, &
    albedo_water, &
    pt_a_water, &
    emissivity_water, &
    turbulent_heat_exchange_coefficient, &
    errmsg) result(status)
    class(nml_river_temperature_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: albedo_water !< Albedo of open water [-]
    type(parameter_t), intent(inout), optional :: pt_a_water !< Priestley-Taylor coefficient for open water [-]
    type(parameter_t), intent(inout), optional :: emissivity_water !< Emissivity of open water [-]
    type(parameter_t), intent(inout), optional :: turbulent_heat_exchange_coefficient !< Turbulent heat-exchange coefficient [W m-2 K-1]
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(albedo_water)) then
      albedo_water%value = ieee_value(albedo_water%value, ieee_quiet_nan) ! sentinel for derived component value
      albedo_water%optimize = .false.
      albedo_water%lower_bound = ieee_value(albedo_water%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      albedo_water%upper_bound = ieee_value(albedo_water%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      albedo_water%lower_bound = 0.03_dp
      albedo_water%upper_bound = 0.2_dp
    end if
    if (present(pt_a_water)) then
      pt_a_water%value = ieee_value(pt_a_water%value, ieee_quiet_nan) ! sentinel for derived component value
      pt_a_water%optimize = .false.
      pt_a_water%lower_bound = ieee_value(pt_a_water%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pt_a_water%upper_bound = ieee_value(pt_a_water%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pt_a_water%lower_bound = 1.2_dp
      pt_a_water%upper_bound = 2.0_dp
    end if
    if (present(emissivity_water)) then
      emissivity_water%value = ieee_value(emissivity_water%value, ieee_quiet_nan) ! sentinel for derived component value
      emissivity_water%optimize = .false.
      emissivity_water%lower_bound = ieee_value(emissivity_water%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      emissivity_water%upper_bound = ieee_value(emissivity_water%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      emissivity_water%lower_bound = 0.95_dp
      emissivity_water%upper_bound = 0.99_dp
    end if
    if (present(turbulent_heat_exchange_coefficient)) then
      turbulent_heat_exchange_coefficient%value = ieee_value(turbulent_heat_exchange_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      turbulent_heat_exchange_coefficient%optimize = .false.
      turbulent_heat_exchange_coefficient%lower_bound = ieee_value(turbulent_heat_exchange_coefficient%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      turbulent_heat_exchange_coefficient%upper_bound = ieee_value(turbulent_heat_exchange_coefficient%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      turbulent_heat_exchange_coefficient%lower_bound = 10.0_dp
      turbulent_heat_exchange_coefficient%upper_bound = 50.0_dp
    end if
  end function nml_river_temperature_1_init_type


  !> \brief Read river_temperature_1 namelist from file
  integer function nml_river_temperature_1_from_file(this, file, errmsg) result(status)
    class(nml_river_temperature_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: albedo_water
    type(parameter_t) :: pt_a_water
    type(parameter_t) :: emissivity_water
    type(parameter_t) :: turbulent_heat_exchange_coefficient
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /river_temperature_1/ &
      albedo_water, &
      pt_a_water, &
      emissivity_water, &
      turbulent_heat_exchange_coefficient

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    albedo_water = this%albedo_water
    pt_a_water = this%pt_a_water
    emissivity_water = this%emissivity_water
    turbulent_heat_exchange_coefficient = this%turbulent_heat_exchange_coefficient

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("river_temperature_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=river_temperature_1, iostat=iostat, iomsg=iomsg)
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
    this%turbulent_heat_exchange_coefficient = turbulent_heat_exchange_coefficient

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_river_temperature_1_from_file

  !> \brief Set river_temperature_1 values
  integer function nml_river_temperature_1_set(this, &
    albedo_water, &
    pt_a_water, &
    emissivity_water, &
    turbulent_heat_exchange_coefficient, &
    errmsg) result(status)

    class(nml_river_temperature_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: albedo_water !< Albedo of open water [-]
    type(parameter_t), intent(in) :: pt_a_water !< Priestley-Taylor coefficient for open water [-]
    type(parameter_t), intent(in) :: emissivity_water !< Emissivity of open water [-]
    type(parameter_t), intent(in) :: turbulent_heat_exchange_coefficient !< Turbulent heat-exchange coefficient [W m-2 K-1]

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%albedo_water = albedo_water
    this%pt_a_water = pt_a_water
    this%emissivity_water = emissivity_water
    this%turbulent_heat_exchange_coefficient = turbulent_heat_exchange_coefficient

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_river_temperature_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_river_temperature_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_river_temperature_1_t), intent(in) :: this !< namelist instance
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
    case ("albedo_water%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_water'"
        return
      end if
      if (ieee_is_nan(this%albedo_water%value)) status = NML_ERR_NOT_SET
    case ("albedo_water%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_water'"
        return
      end if
    case ("albedo_water%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_water'"
        return
      end if
    case ("albedo_water%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_water'"
        return
      end if
    case ("albedo_water")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_water'"
        return
      end if
      if (ieee_is_nan(this%albedo_water%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("pt_a_water%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pt_a_water'"
        return
      end if
      if (ieee_is_nan(this%pt_a_water%value)) status = NML_ERR_NOT_SET
    case ("pt_a_water%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pt_a_water'"
        return
      end if
    case ("pt_a_water%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pt_a_water'"
        return
      end if
    case ("pt_a_water%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pt_a_water'"
        return
      end if
    case ("pt_a_water")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pt_a_water'"
        return
      end if
      if (ieee_is_nan(this%pt_a_water%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("emissivity_water%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'emissivity_water'"
        return
      end if
      if (ieee_is_nan(this%emissivity_water%value)) status = NML_ERR_NOT_SET
    case ("emissivity_water%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'emissivity_water'"
        return
      end if
    case ("emissivity_water%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'emissivity_water'"
        return
      end if
    case ("emissivity_water%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'emissivity_water'"
        return
      end if
    case ("emissivity_water")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'emissivity_water'"
        return
      end if
      if (ieee_is_nan(this%emissivity_water%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("turbulent_heat_exchange_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'turbulent_heat_exchange_coefficient'"
        return
      end if
      if (ieee_is_nan(this%turbulent_heat_exchange_coefficient%value)) status = NML_ERR_NOT_SET
    case ("turbulent_heat_exchange_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'turbulent_heat_exchange_coefficient'"
        return
      end if
    case ("turbulent_heat_exchange_coefficient%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'turbulent_heat_exchange_coefficient'"
        return
      end if
    case ("turbulent_heat_exchange_coefficient%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'turbulent_heat_exchange_coefficient'"
        return
      end if
    case ("turbulent_heat_exchange_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'turbulent_heat_exchange_coefficient'"
        return
      end if
      if (ieee_is_nan(this%turbulent_heat_exchange_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_river_temperature_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_river_temperature_1_is_valid(this, errmsg) result(status)
    class(nml_river_temperature_1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("albedo_water", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: albedo_water"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("pt_a_water", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pt_a_water"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("emissivity_water", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: emissivity_water"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("turbulent_heat_exchange_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: turbulent_heat_exchange_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_river_temperature_1_is_valid

end module nml_river_temperature_1
