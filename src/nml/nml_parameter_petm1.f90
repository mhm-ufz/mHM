!> \file nml_parameter_petm1.f90
!> \copydoc nml_petm1

!> \brief PET - Case -1
!> \details Parameters for PET (case -1 - LAI correction).
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_petm1
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

  !> \class nml_petm1_t
  !> \brief PET - Case -1
  !> \details Parameters for PET (case -1 - LAI correction).
  type, public :: nml_petm1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: pet_a_forest !< Potential evapotranspiration forest
    real(dp), dimension(5) :: pet_a_impervious !< Potential evapotranspiration impervious
    real(dp), dimension(5) :: pet_a_pervious !< Potential evapotranspiration pervious
    real(dp), dimension(5) :: pet_b !< Potential evapotranspiration b
    real(dp), dimension(5) :: pet_c !< Potential evapotranspiration c
  contains
    procedure :: init => nml_petm1_init
    procedure :: from_file => nml_petm1_from_file
    procedure :: set => nml_petm1_set
    procedure :: is_set => nml_petm1_is_set
    procedure :: is_valid => nml_petm1_is_valid
  end type nml_petm1_t

contains

  !> \brief Initialize defaults and sentinels for petm1
  integer function nml_petm1_init(this, errmsg) result(status)
    class(nml_petm1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%pet_a_forest = ieee_value(this%pet_a_forest, ieee_quiet_nan) ! sentinel for required real array
    this%pet_a_impervious = ieee_value(this%pet_a_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%pet_a_pervious = ieee_value(this%pet_a_pervious, ieee_quiet_nan) ! sentinel for required real array
    this%pet_b = ieee_value(this%pet_b, ieee_quiet_nan) ! sentinel for required real array
    this%pet_c = ieee_value(this%pet_c, ieee_quiet_nan) ! sentinel for required real array
  end function nml_petm1_init


  !> \brief Read petm1 namelist from file
  integer function nml_petm1_from_file(this, file, errmsg) result(status)
    class(nml_petm1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: pet_a_forest
    real(dp), dimension(5) :: pet_a_impervious
    real(dp), dimension(5) :: pet_a_pervious
    real(dp), dimension(5) :: pet_b
    real(dp), dimension(5) :: pet_c
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /petm1/ &
      pet_a_forest, &
      pet_a_impervious, &
      pet_a_pervious, &
      pet_b, &
      pet_c

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    pet_a_forest = this%pet_a_forest
    pet_a_impervious = this%pet_a_impervious
    pet_a_pervious = this%pet_a_pervious
    pet_b = this%pet_b
    pet_c = this%pet_c

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("petm1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=petm1, iostat=iostat, iomsg=iomsg)
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
    this%pet_a_forest = pet_a_forest
    this%pet_a_impervious = pet_a_impervious
    this%pet_a_pervious = pet_a_pervious
    this%pet_b = pet_b
    this%pet_c = pet_c

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_petm1_from_file

  !> \brief Set petm1 values
  integer function nml_petm1_set(this, &
    pet_a_forest, &
    pet_a_impervious, &
    pet_a_pervious, &
    pet_b, &
    pet_c, &
    errmsg) result(status)

    class(nml_petm1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: pet_a_forest !< Potential evapotranspiration forest
    real(dp), dimension(:), intent(in) :: pet_a_impervious !< Potential evapotranspiration impervious
    real(dp), dimension(:), intent(in) :: pet_a_pervious !< Potential evapotranspiration pervious
    real(dp), dimension(:), intent(in) :: pet_b !< Potential evapotranspiration b
    real(dp), dimension(:), intent(in) :: pet_c !< Potential evapotranspiration c
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(pet_a_forest, 1) > size(this%pet_a_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_a_forest'"
      return
    end if
    lb__1 = lbound(this%pet_a_forest, 1)
    ub__1 = lb__1 + size(pet_a_forest, 1) - 1
    this%pet_a_forest(lb__1:ub__1) = pet_a_forest
    if (size(pet_a_impervious, 1) > size(this%pet_a_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_a_impervious'"
      return
    end if
    lb__1 = lbound(this%pet_a_impervious, 1)
    ub__1 = lb__1 + size(pet_a_impervious, 1) - 1
    this%pet_a_impervious(lb__1:ub__1) = pet_a_impervious
    if (size(pet_a_pervious, 1) > size(this%pet_a_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_a_pervious'"
      return
    end if
    lb__1 = lbound(this%pet_a_pervious, 1)
    ub__1 = lb__1 + size(pet_a_pervious, 1) - 1
    this%pet_a_pervious(lb__1:ub__1) = pet_a_pervious
    if (size(pet_b, 1) > size(this%pet_b, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_b'"
      return
    end if
    lb__1 = lbound(this%pet_b, 1)
    ub__1 = lb__1 + size(pet_b, 1) - 1
    this%pet_b(lb__1:ub__1) = pet_b
    if (size(pet_c, 1) > size(this%pet_c, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_c'"
      return
    end if
    lb__1 = lbound(this%pet_c, 1)
    ub__1 = lb__1 + size(pet_c, 1) - 1
    this%pet_c(lb__1:ub__1) = pet_c

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_petm1_set

  !> \brief Check whether a namelist value was set
  integer function nml_petm1_is_set(this, name, idx, errmsg) result(status)
    class(nml_petm1_t), intent(in) :: this !< namelist instance
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
    case ("pet_a_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_a_forest), ubound(this%pet_a_forest), &
          "PET_a_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pet_a_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pet_a_forest))) status = NML_ERR_NOT_SET
      end if
    case ("pet_a_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_a_impervious), ubound(this%pet_a_impervious), &
          "PET_a_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pet_a_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pet_a_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("pet_a_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_a_pervious), ubound(this%pet_a_pervious), &
          "PET_a_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pet_a_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pet_a_pervious))) status = NML_ERR_NOT_SET
      end if
    case ("pet_b")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_b), ubound(this%pet_b), &
          "PET_b", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pet_b(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pet_b))) status = NML_ERR_NOT_SET
      end if
    case ("pet_c")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_c), ubound(this%pet_c), &
          "PET_c", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%pet_c(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%pet_c))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_petm1_is_set

  !> \brief Validate required values and constraints
  integer function nml_petm1_is_valid(this, errmsg) result(status)
    class(nml_petm1_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%pet_a_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PET_a_forest"
      return
    end if
    if (any(ieee_is_nan(this%pet_a_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PET_a_forest"
      return
    end if
    if (all(ieee_is_nan(this%pet_a_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PET_a_impervious"
      return
    end if
    if (any(ieee_is_nan(this%pet_a_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PET_a_impervious"
      return
    end if
    if (all(ieee_is_nan(this%pet_a_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PET_a_pervious"
      return
    end if
    if (any(ieee_is_nan(this%pet_a_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PET_a_pervious"
      return
    end if
    if (all(ieee_is_nan(this%pet_b(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PET_b"
      return
    end if
    if (any(ieee_is_nan(this%pet_b(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PET_b"
      return
    end if
    if (all(ieee_is_nan(this%pet_c(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PET_c"
      return
    end if
    if (any(ieee_is_nan(this%pet_c(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PET_c"
      return
    end if
  end function nml_petm1_is_valid

end module nml_petm1
