!> \file nml_parameter_pet_3.f90
!> \copydoc nml_pet_3

!> \brief PET - Case 3
!> \details Parameters for the Penman-Monteith PET method.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet_3
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

  !> \class nml_pet_3_t
  !> \brief PET - Case 3
  !> \details Parameters for the Penman-Monteith PET method.
  type, public :: nml_pet_3_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: canopy_height_forest !< Canopy height for forest
    type(parameter_t) :: canopy_height_impervious !< Canopy height for impervious areas
    type(parameter_t) :: canopy_height_pervious !< Canopy height for pervious areas
    type(parameter_t) :: displacement_height_coefficient !< Displacement-height coefficient
    type(parameter_t) :: momentum_roughness_length_coefficient !< Momentum roughness-length coefficient
    type(parameter_t) :: heat_roughness_length_coefficient !< Heat roughness-length coefficient
    type(parameter_t) :: stomatal_resistance !< Stomatal resistance
  contains
    procedure :: init => nml_pet_3_init
    procedure :: init_type => nml_pet_3_init_type
    procedure :: from_file => nml_pet_3_from_file
    procedure :: set => nml_pet_3_set
    procedure :: is_set => nml_pet_3_is_set
    procedure :: is_valid => nml_pet_3_is_valid
  end type nml_pet_3_t

contains

  !> \brief Initialize defaults and sentinels for pet_3
  integer function nml_pet_3_init(this, errmsg) result(status)
    class(nml_pet_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      canopy_height_forest=this%canopy_height_forest, &
      canopy_height_impervious=this%canopy_height_impervious, &
      canopy_height_pervious=this%canopy_height_pervious, &
      displacement_height_coefficient=this%displacement_height_coefficient, &
      momentum_roughness_length_coefficient=this%momentum_roughness_length_coefficient, &
      heat_roughness_length_coefficient=this%heat_roughness_length_coefficient, &
      stomatal_resistance=this%stomatal_resistance, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_pet_3_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_pet_3_init_type(this, &
    canopy_height_forest, &
    canopy_height_impervious, &
    canopy_height_pervious, &
    displacement_height_coefficient, &
    momentum_roughness_length_coefficient, &
    heat_roughness_length_coefficient, &
    stomatal_resistance, &
    errmsg) result(status)
    class(nml_pet_3_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: canopy_height_forest !< Canopy height for forest
    type(parameter_t), intent(inout), optional :: canopy_height_impervious !< Canopy height for impervious areas
    type(parameter_t), intent(inout), optional :: canopy_height_pervious !< Canopy height for pervious areas
    type(parameter_t), intent(inout), optional :: displacement_height_coefficient !< Displacement-height coefficient
    type(parameter_t), intent(inout), optional :: momentum_roughness_length_coefficient !< Momentum roughness-length coefficient
    type(parameter_t), intent(inout), optional :: heat_roughness_length_coefficient !< Heat roughness-length coefficient
    type(parameter_t), intent(inout), optional :: stomatal_resistance !< Stomatal resistance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(canopy_height_forest)) then
      canopy_height_forest%value = ieee_value(canopy_height_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      canopy_height_forest%optimize = .false.
      canopy_height_forest%lower_bound = ieee_value(canopy_height_forest%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      canopy_height_forest%upper_bound = ieee_value(canopy_height_forest%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      canopy_height_forest%lower_bound = 15.0_dp
      canopy_height_forest%upper_bound = 40.0_dp
    end if
    if (present(canopy_height_impervious)) then
      canopy_height_impervious%value = ieee_value(canopy_height_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      canopy_height_impervious%optimize = .false.
      canopy_height_impervious%lower_bound = ieee_value(canopy_height_impervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      canopy_height_impervious%upper_bound = ieee_value(canopy_height_impervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      canopy_height_impervious%lower_bound = 0.01_dp
      canopy_height_impervious%upper_bound = 0.5_dp
    end if
    if (present(canopy_height_pervious)) then
      canopy_height_pervious%value = ieee_value(canopy_height_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      canopy_height_pervious%optimize = .false.
      canopy_height_pervious%lower_bound = ieee_value(canopy_height_pervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      canopy_height_pervious%upper_bound = ieee_value(canopy_height_pervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      canopy_height_pervious%lower_bound = 0.1_dp
      canopy_height_pervious%upper_bound = 5.0_dp
    end if
    if (present(displacement_height_coefficient)) then
      displacement_height_coefficient%value = ieee_value(displacement_height_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      displacement_height_coefficient%optimize = .false.
      displacement_height_coefficient%lower_bound = ieee_value(displacement_height_coefficient%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      displacement_height_coefficient%upper_bound = ieee_value(displacement_height_coefficient%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      displacement_height_coefficient%lower_bound = 0.5_dp
      displacement_height_coefficient%upper_bound = 0.85_dp
    end if
    if (present(momentum_roughness_length_coefficient)) then
      momentum_roughness_length_coefficient%value = ieee_value(momentum_roughness_length_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      momentum_roughness_length_coefficient%optimize = .false.
      momentum_roughness_length_coefficient%lower_bound = ieee_value(momentum_roughness_length_coefficient%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      momentum_roughness_length_coefficient%upper_bound = ieee_value(momentum_roughness_length_coefficient%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      momentum_roughness_length_coefficient%lower_bound = 0.09_dp
      momentum_roughness_length_coefficient%upper_bound = 0.16_dp
    end if
    if (present(heat_roughness_length_coefficient)) then
      heat_roughness_length_coefficient%value = ieee_value(heat_roughness_length_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      heat_roughness_length_coefficient%optimize = .false.
      heat_roughness_length_coefficient%lower_bound = ieee_value(heat_roughness_length_coefficient%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      heat_roughness_length_coefficient%upper_bound = ieee_value(heat_roughness_length_coefficient%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      heat_roughness_length_coefficient%lower_bound = 0.07_dp
      heat_roughness_length_coefficient%upper_bound = 0.13_dp
    end if
    if (present(stomatal_resistance)) then
      stomatal_resistance%value = ieee_value(stomatal_resistance%value, ieee_quiet_nan) ! sentinel for derived component value
      stomatal_resistance%optimize = .false.
      stomatal_resistance%lower_bound = ieee_value(stomatal_resistance%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      stomatal_resistance%upper_bound = ieee_value(stomatal_resistance%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      stomatal_resistance%lower_bound = 10.0_dp
      stomatal_resistance%upper_bound = 200.0_dp
    end if
  end function nml_pet_3_init_type


  !> \brief Read pet_3 namelist from file
  integer function nml_pet_3_from_file(this, file, errmsg) result(status)
    class(nml_pet_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: canopy_height_forest
    type(parameter_t) :: canopy_height_impervious
    type(parameter_t) :: canopy_height_pervious
    type(parameter_t) :: displacement_height_coefficient
    type(parameter_t) :: momentum_roughness_length_coefficient
    type(parameter_t) :: heat_roughness_length_coefficient
    type(parameter_t) :: stomatal_resistance
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet_3/ &
      canopy_height_forest, &
      canopy_height_impervious, &
      canopy_height_pervious, &
      displacement_height_coefficient, &
      momentum_roughness_length_coefficient, &
      heat_roughness_length_coefficient, &
      stomatal_resistance

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    canopy_height_forest = this%canopy_height_forest
    canopy_height_impervious = this%canopy_height_impervious
    canopy_height_pervious = this%canopy_height_pervious
    displacement_height_coefficient = this%displacement_height_coefficient
    momentum_roughness_length_coefficient = this%momentum_roughness_length_coefficient
    heat_roughness_length_coefficient = this%heat_roughness_length_coefficient
    stomatal_resistance = this%stomatal_resistance

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("pet_3", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet_3, iostat=iostat, iomsg=iomsg)
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
    this%canopy_height_forest = canopy_height_forest
    this%canopy_height_impervious = canopy_height_impervious
    this%canopy_height_pervious = canopy_height_pervious
    this%displacement_height_coefficient = displacement_height_coefficient
    this%momentum_roughness_length_coefficient = momentum_roughness_length_coefficient
    this%heat_roughness_length_coefficient = heat_roughness_length_coefficient
    this%stomatal_resistance = stomatal_resistance

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_3_from_file

  !> \brief Set pet_3 values
  integer function nml_pet_3_set(this, &
    canopy_height_forest, &
    canopy_height_impervious, &
    canopy_height_pervious, &
    displacement_height_coefficient, &
    momentum_roughness_length_coefficient, &
    heat_roughness_length_coefficient, &
    stomatal_resistance, &
    errmsg) result(status)

    class(nml_pet_3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: canopy_height_forest !< Canopy height for forest
    type(parameter_t), intent(in) :: canopy_height_impervious !< Canopy height for impervious areas
    type(parameter_t), intent(in) :: canopy_height_pervious !< Canopy height for pervious areas
    type(parameter_t), intent(in) :: displacement_height_coefficient !< Displacement-height coefficient
    type(parameter_t), intent(in) :: momentum_roughness_length_coefficient !< Momentum roughness-length coefficient
    type(parameter_t), intent(in) :: heat_roughness_length_coefficient !< Heat roughness-length coefficient
    type(parameter_t), intent(in) :: stomatal_resistance !< Stomatal resistance

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%canopy_height_forest = canopy_height_forest
    this%canopy_height_impervious = canopy_height_impervious
    this%canopy_height_pervious = canopy_height_pervious
    this%displacement_height_coefficient = displacement_height_coefficient
    this%momentum_roughness_length_coefficient = momentum_roughness_length_coefficient
    this%heat_roughness_length_coefficient = heat_roughness_length_coefficient
    this%stomatal_resistance = stomatal_resistance

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_3_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet_3_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet_3_t), intent(in) :: this !< namelist instance
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
    case ("canopy_height_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_forest'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_forest%value)) status = NML_ERR_NOT_SET
    case ("canopy_height_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_forest'"
        return
      end if
    case ("canopy_height_forest%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_forest'"
        return
      end if
    case ("canopy_height_forest%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_forest'"
        return
      end if
    case ("canopy_height_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_forest'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("canopy_height_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_impervious'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_impervious%value)) status = NML_ERR_NOT_SET
    case ("canopy_height_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_impervious'"
        return
      end if
    case ("canopy_height_impervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_impervious'"
        return
      end if
    case ("canopy_height_impervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_impervious'"
        return
      end if
    case ("canopy_height_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_impervious'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("canopy_height_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_pervious'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_pervious%value)) status = NML_ERR_NOT_SET
    case ("canopy_height_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_pervious'"
        return
      end if
    case ("canopy_height_pervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_pervious'"
        return
      end if
    case ("canopy_height_pervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_pervious'"
        return
      end if
    case ("canopy_height_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'canopy_height_pervious'"
        return
      end if
      if (ieee_is_nan(this%canopy_height_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("displacement_height_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'displacement_height_coefficient'"
        return
      end if
      if (ieee_is_nan(this%displacement_height_coefficient%value)) status = NML_ERR_NOT_SET
    case ("displacement_height_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'displacement_height_coefficient'"
        return
      end if
    case ("displacement_height_coefficient%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'displacement_height_coefficient'"
        return
      end if
    case ("displacement_height_coefficient%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'displacement_height_coefficient'"
        return
      end if
    case ("displacement_height_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'displacement_height_coefficient'"
        return
      end if
      if (ieee_is_nan(this%displacement_height_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("momentum_roughness_length_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'momentum_roughness_length_coefficient'"
        return
      end if
      if (ieee_is_nan(this%momentum_roughness_length_coefficient%value)) status = NML_ERR_NOT_SET
    case ("momentum_roughness_length_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'momentum_roughness_length_coefficient'"
        return
      end if
    case ("momentum_roughness_length_coefficient%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'momentum_roughness_length_coefficient'"
        return
      end if
    case ("momentum_roughness_length_coefficient%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'momentum_roughness_length_coefficient'"
        return
      end if
    case ("momentum_roughness_length_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'momentum_roughness_length_coefficient'"
        return
      end if
      if (ieee_is_nan(this%momentum_roughness_length_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("heat_roughness_length_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'heat_roughness_length_coefficient'"
        return
      end if
      if (ieee_is_nan(this%heat_roughness_length_coefficient%value)) status = NML_ERR_NOT_SET
    case ("heat_roughness_length_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'heat_roughness_length_coefficient'"
        return
      end if
    case ("heat_roughness_length_coefficient%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'heat_roughness_length_coefficient'"
        return
      end if
    case ("heat_roughness_length_coefficient%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'heat_roughness_length_coefficient'"
        return
      end if
    case ("heat_roughness_length_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'heat_roughness_length_coefficient'"
        return
      end if
      if (ieee_is_nan(this%heat_roughness_length_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("stomatal_resistance%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'stomatal_resistance'"
        return
      end if
      if (ieee_is_nan(this%stomatal_resistance%value)) status = NML_ERR_NOT_SET
    case ("stomatal_resistance%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'stomatal_resistance'"
        return
      end if
    case ("stomatal_resistance%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'stomatal_resistance'"
        return
      end if
    case ("stomatal_resistance%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'stomatal_resistance'"
        return
      end if
    case ("stomatal_resistance")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'stomatal_resistance'"
        return
      end if
      if (ieee_is_nan(this%stomatal_resistance%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet_3_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet_3_is_valid(this, errmsg) result(status)
    class(nml_pet_3_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("canopy_height_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: canopy_height_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("canopy_height_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: canopy_height_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("canopy_height_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: canopy_height_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("displacement_height_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: displacement_height_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("momentum_roughness_length_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: momentum_roughness_length_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("heat_roughness_length_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: heat_roughness_length_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("stomatal_resistance", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: stomatal_resistance"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_pet_3_is_valid

end module nml_pet_3
