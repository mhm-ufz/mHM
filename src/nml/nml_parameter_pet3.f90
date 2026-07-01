!> \file nml_parameter_pet3.f90
!> \copydoc nml_pet3

!> \brief PET - Case 3
!> \details Parameters for PET (case 3 - Penman-Monteith).
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet3
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

  !> \class nml_pet3_t
  !> \brief PET - Case 3
  !> \details Parameters for PET (case 3 - Penman-Monteith).
  type, public :: nml_pet3_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: canopyheigth_forest !< Canopy height forest
    real(dp), dimension(5) :: canopyheigth_impervious !< Canopy height impervious
    real(dp), dimension(5) :: canopyheigth_pervious !< Canopy height pervious
    real(dp), dimension(5) :: displacementheight_coeff !< Displacement height coefficient
    real(dp), dimension(5) :: roughnesslength_momentum_coeff !< Roughness length momentum coefficient
    real(dp), dimension(5) :: roughnesslength_heat_coeff !< Roughness length heat coefficient
    real(dp), dimension(5) :: stomatal_resistance !< Stomatal resistance
  contains
    procedure :: init => nml_pet3_init
    procedure :: from_file => nml_pet3_from_file
    procedure :: set => nml_pet3_set
    procedure :: is_set => nml_pet3_is_set
    procedure :: is_valid => nml_pet3_is_valid
  end type nml_pet3_t

contains

  !> \brief Initialize defaults and sentinels for pet3
  integer function nml_pet3_init(this, errmsg) result(status)
    class(nml_pet3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%canopyheigth_forest = ieee_value(this%canopyheigth_forest, ieee_quiet_nan) ! sentinel for required real array
    this%canopyheigth_impervious = ieee_value(this%canopyheigth_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%canopyheigth_pervious = ieee_value(this%canopyheigth_pervious, ieee_quiet_nan) ! sentinel for required real array
    this%displacementheight_coeff = ieee_value(this%displacementheight_coeff, ieee_quiet_nan) ! sentinel for required real array
    this%roughnesslength_momentum_coeff = ieee_value(this%roughnesslength_momentum_coeff, ieee_quiet_nan) ! sentinel for required real array
    this%roughnesslength_heat_coeff = ieee_value(this%roughnesslength_heat_coeff, ieee_quiet_nan) ! sentinel for required real array
    this%stomatal_resistance = ieee_value(this%stomatal_resistance, ieee_quiet_nan) ! sentinel for required real array
  end function nml_pet3_init


  !> \brief Read pet3 namelist from file
  integer function nml_pet3_from_file(this, file, errmsg) result(status)
    class(nml_pet3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: canopyheigth_forest
    real(dp), dimension(5) :: canopyheigth_impervious
    real(dp), dimension(5) :: canopyheigth_pervious
    real(dp), dimension(5) :: displacementheight_coeff
    real(dp), dimension(5) :: roughnesslength_momentum_coeff
    real(dp), dimension(5) :: roughnesslength_heat_coeff
    real(dp), dimension(5) :: stomatal_resistance
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet3/ &
      canopyheigth_forest, &
      canopyheigth_impervious, &
      canopyheigth_pervious, &
      displacementheight_coeff, &
      roughnesslength_momentum_coeff, &
      roughnesslength_heat_coeff, &
      stomatal_resistance

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    canopyheigth_forest = this%canopyheigth_forest
    canopyheigth_impervious = this%canopyheigth_impervious
    canopyheigth_pervious = this%canopyheigth_pervious
    displacementheight_coeff = this%displacementheight_coeff
    roughnesslength_momentum_coeff = this%roughnesslength_momentum_coeff
    roughnesslength_heat_coeff = this%roughnesslength_heat_coeff
    stomatal_resistance = this%stomatal_resistance

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("pet3", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet3, iostat=iostat, iomsg=iomsg)
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
    this%canopyheigth_forest = canopyheigth_forest
    this%canopyheigth_impervious = canopyheigth_impervious
    this%canopyheigth_pervious = canopyheigth_pervious
    this%displacementheight_coeff = displacementheight_coeff
    this%roughnesslength_momentum_coeff = roughnesslength_momentum_coeff
    this%roughnesslength_heat_coeff = roughnesslength_heat_coeff
    this%stomatal_resistance = stomatal_resistance

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet3_from_file

  !> \brief Set pet3 values
  integer function nml_pet3_set(this, &
    canopyheigth_forest, &
    canopyheigth_impervious, &
    canopyheigth_pervious, &
    displacementheight_coeff, &
    roughnesslength_momentum_coeff, &
    roughnesslength_heat_coeff, &
    stomatal_resistance, &
    errmsg) result(status)

    class(nml_pet3_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: canopyheigth_forest !< Canopy height forest
    real(dp), dimension(:), intent(in) :: canopyheigth_impervious !< Canopy height impervious
    real(dp), dimension(:), intent(in) :: canopyheigth_pervious !< Canopy height pervious
    real(dp), dimension(:), intent(in) :: displacementheight_coeff !< Displacement height coefficient
    real(dp), dimension(:), intent(in) :: roughnesslength_momentum_coeff !< Roughness length momentum coefficient
    real(dp), dimension(:), intent(in) :: roughnesslength_heat_coeff !< Roughness length heat coefficient
    real(dp), dimension(:), intent(in) :: stomatal_resistance !< Stomatal resistance
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(canopyheigth_forest, 1) > size(this%canopyheigth_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'canopyheigth_forest'"
      return
    end if
    lb__1 = lbound(this%canopyheigth_forest, 1)
    ub__1 = lb__1 + size(canopyheigth_forest, 1) - 1
    this%canopyheigth_forest(lb__1:ub__1) = canopyheigth_forest
    if (size(canopyheigth_impervious, 1) > size(this%canopyheigth_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'canopyheigth_impervious'"
      return
    end if
    lb__1 = lbound(this%canopyheigth_impervious, 1)
    ub__1 = lb__1 + size(canopyheigth_impervious, 1) - 1
    this%canopyheigth_impervious(lb__1:ub__1) = canopyheigth_impervious
    if (size(canopyheigth_pervious, 1) > size(this%canopyheigth_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'canopyheigth_pervious'"
      return
    end if
    lb__1 = lbound(this%canopyheigth_pervious, 1)
    ub__1 = lb__1 + size(canopyheigth_pervious, 1) - 1
    this%canopyheigth_pervious(lb__1:ub__1) = canopyheigth_pervious
    if (size(displacementheight_coeff, 1) > size(this%displacementheight_coeff, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'displacementheight_coeff'"
      return
    end if
    lb__1 = lbound(this%displacementheight_coeff, 1)
    ub__1 = lb__1 + size(displacementheight_coeff, 1) - 1
    this%displacementheight_coeff(lb__1:ub__1) = displacementheight_coeff
    if (size(roughnesslength_momentum_coeff, 1) > size(this%roughnesslength_momentum_coeff, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'roughnesslength_momentum_coeff'"
      return
    end if
    lb__1 = lbound(this%roughnesslength_momentum_coeff, 1)
    ub__1 = lb__1 + size(roughnesslength_momentum_coeff, 1) - 1
    this%roughnesslength_momentum_coeff(lb__1:ub__1) = roughnesslength_momentum_coeff
    if (size(roughnesslength_heat_coeff, 1) > size(this%roughnesslength_heat_coeff, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'roughnesslength_heat_coeff'"
      return
    end if
    lb__1 = lbound(this%roughnesslength_heat_coeff, 1)
    ub__1 = lb__1 + size(roughnesslength_heat_coeff, 1) - 1
    this%roughnesslength_heat_coeff(lb__1:ub__1) = roughnesslength_heat_coeff
    if (size(stomatal_resistance, 1) > size(this%stomatal_resistance, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'stomatal_resistance'"
      return
    end if
    lb__1 = lbound(this%stomatal_resistance, 1)
    ub__1 = lb__1 + size(stomatal_resistance, 1) - 1
    this%stomatal_resistance(lb__1:ub__1) = stomatal_resistance

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet3_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet3_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet3_t), intent(in) :: this !< namelist instance
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
    case ("canopyheigth_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%canopyheigth_forest), ubound(this%canopyheigth_forest), &
          "canopyheigth_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%canopyheigth_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%canopyheigth_forest))) status = NML_ERR_NOT_SET
      end if
    case ("canopyheigth_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%canopyheigth_impervious), ubound(this%canopyheigth_impervious), &
          "canopyheigth_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%canopyheigth_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%canopyheigth_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("canopyheigth_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%canopyheigth_pervious), ubound(this%canopyheigth_pervious), &
          "canopyheigth_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%canopyheigth_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%canopyheigth_pervious))) status = NML_ERR_NOT_SET
      end if
    case ("displacementheight_coeff")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%displacementheight_coeff), ubound(this%displacementheight_coeff), &
          "displacementheight_coeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%displacementheight_coeff(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%displacementheight_coeff))) status = NML_ERR_NOT_SET
      end if
    case ("roughnesslength_momentum_coeff")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%roughnesslength_momentum_coeff), ubound(this%roughnesslength_momentum_coeff), &
          "roughnesslength_momentum_coeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%roughnesslength_momentum_coeff(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%roughnesslength_momentum_coeff))) status = NML_ERR_NOT_SET
      end if
    case ("roughnesslength_heat_coeff")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%roughnesslength_heat_coeff), ubound(this%roughnesslength_heat_coeff), &
          "roughnesslength_heat_coeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%roughnesslength_heat_coeff(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%roughnesslength_heat_coeff))) status = NML_ERR_NOT_SET
      end if
    case ("stomatal_resistance")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%stomatal_resistance), ubound(this%stomatal_resistance), &
          "stomatal_resistance", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%stomatal_resistance(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%stomatal_resistance))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet3_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet3_is_valid(this, errmsg) result(status)
    class(nml_pet3_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%canopyheigth_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: canopyheigth_forest"
      return
    end if
    if (any(ieee_is_nan(this%canopyheigth_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: canopyheigth_forest"
      return
    end if
    if (all(ieee_is_nan(this%canopyheigth_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: canopyheigth_impervious"
      return
    end if
    if (any(ieee_is_nan(this%canopyheigth_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: canopyheigth_impervious"
      return
    end if
    if (all(ieee_is_nan(this%canopyheigth_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: canopyheigth_pervious"
      return
    end if
    if (any(ieee_is_nan(this%canopyheigth_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: canopyheigth_pervious"
      return
    end if
    if (all(ieee_is_nan(this%displacementheight_coeff(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: displacementheight_coeff"
      return
    end if
    if (any(ieee_is_nan(this%displacementheight_coeff(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: displacementheight_coeff"
      return
    end if
    if (all(ieee_is_nan(this%roughnesslength_momentum_coeff(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: roughnesslength_momentum_coeff"
      return
    end if
    if (any(ieee_is_nan(this%roughnesslength_momentum_coeff(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: roughnesslength_momentum_coeff"
      return
    end if
    if (all(ieee_is_nan(this%roughnesslength_heat_coeff(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: roughnesslength_heat_coeff"
      return
    end if
    if (any(ieee_is_nan(this%roughnesslength_heat_coeff(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: roughnesslength_heat_coeff"
      return
    end if
    if (all(ieee_is_nan(this%stomatal_resistance(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: stomatal_resistance"
      return
    end if
    if (any(ieee_is_nan(this%stomatal_resistance(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: stomatal_resistance"
      return
    end if
  end function nml_pet3_is_valid

end module nml_pet3
