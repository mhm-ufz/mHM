!> \file nml_parameter_pet_m1.f90
!> \copydoc nml_pet_m1

!> \brief PET - Case -1
!> \details Parameters for LAI correction of externally supplied PET.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet_m1
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

  !> \class nml_pet_m1_t
  !> \brief PET - Case -1
  !> \details Parameters for LAI correction of externally supplied PET.
  type, public :: nml_pet_m1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: pet_a_forest !< PET correction coefficient a for forest
    type(parameter_t) :: pet_a_impervious !< PET correction coefficient a for impervious areas
    type(parameter_t) :: pet_a_pervious !< PET correction coefficient a for pervious areas
    type(parameter_t) :: pet_b !< PET correction coefficient b
    type(parameter_t) :: pet_c !< PET correction coefficient c
  contains
    procedure :: init => nml_pet_m1_init
    procedure :: init_type => nml_pet_m1_init_type
    procedure :: from_file => nml_pet_m1_from_file
    procedure :: set => nml_pet_m1_set
    procedure :: is_set => nml_pet_m1_is_set
    procedure :: is_valid => nml_pet_m1_is_valid
  end type nml_pet_m1_t

contains

  !> \brief Initialize defaults and sentinels for pet_m1
  integer function nml_pet_m1_init(this, errmsg) result(status)
    class(nml_pet_m1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      pet_a_forest=this%pet_a_forest, &
      pet_a_impervious=this%pet_a_impervious, &
      pet_a_pervious=this%pet_a_pervious, &
      pet_b=this%pet_b, &
      pet_c=this%pet_c, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_pet_m1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_pet_m1_init_type(this, &
    pet_a_forest, &
    pet_a_impervious, &
    pet_a_pervious, &
    pet_b, &
    pet_c, &
    errmsg) result(status)
    class(nml_pet_m1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: pet_a_forest !< PET correction coefficient a for forest
    type(parameter_t), intent(inout), optional :: pet_a_impervious !< PET correction coefficient a for impervious areas
    type(parameter_t), intent(inout), optional :: pet_a_pervious !< PET correction coefficient a for pervious areas
    type(parameter_t), intent(inout), optional :: pet_b !< PET correction coefficient b
    type(parameter_t), intent(inout), optional :: pet_c !< PET correction coefficient c
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(pet_a_forest)) then
      pet_a_forest%value = ieee_value(pet_a_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      pet_a_forest%optimize = .false.
      pet_a_forest%lower_bound = ieee_value(pet_a_forest%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pet_a_forest%upper_bound = ieee_value(pet_a_forest%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pet_a_forest%lower_bound = 0.3_dp
      pet_a_forest%upper_bound = 1.3_dp
    end if
    if (present(pet_a_impervious)) then
      pet_a_impervious%value = ieee_value(pet_a_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      pet_a_impervious%optimize = .false.
      pet_a_impervious%lower_bound = ieee_value(pet_a_impervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pet_a_impervious%upper_bound = ieee_value(pet_a_impervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pet_a_impervious%lower_bound = 0.3_dp
      pet_a_impervious%upper_bound = 1.3_dp
    end if
    if (present(pet_a_pervious)) then
      pet_a_pervious%value = ieee_value(pet_a_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      pet_a_pervious%optimize = .false.
      pet_a_pervious%lower_bound = ieee_value(pet_a_pervious%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pet_a_pervious%upper_bound = ieee_value(pet_a_pervious%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pet_a_pervious%lower_bound = 0.3_dp
      pet_a_pervious%upper_bound = 1.3_dp
    end if
    if (present(pet_b)) then
      pet_b%value = ieee_value(pet_b%value, ieee_quiet_nan) ! sentinel for derived component value
      pet_b%optimize = .false.
      pet_b%lower_bound = ieee_value(pet_b%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pet_b%upper_bound = ieee_value(pet_b%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pet_b%lower_bound = 0.0_dp
      pet_b%upper_bound = 1.5_dp
    end if
    if (present(pet_c)) then
      pet_c%value = ieee_value(pet_c%value, ieee_quiet_nan) ! sentinel for derived component value
      pet_c%optimize = .false.
      pet_c%lower_bound = ieee_value(pet_c%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      pet_c%upper_bound = ieee_value(pet_c%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      pet_c%lower_bound = -2.0_dp
      pet_c%upper_bound = 0.0_dp
    end if
  end function nml_pet_m1_init_type


  !> \brief Read pet_m1 namelist from file
  integer function nml_pet_m1_from_file(this, file, errmsg) result(status)
    class(nml_pet_m1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: pet_a_forest
    type(parameter_t) :: pet_a_impervious
    type(parameter_t) :: pet_a_pervious
    type(parameter_t) :: pet_b
    type(parameter_t) :: pet_c
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet_m1/ &
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

    status = nml%find("pet_m1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet_m1, iostat=iostat, iomsg=iomsg)
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
  end function nml_pet_m1_from_file

  !> \brief Set pet_m1 values
  integer function nml_pet_m1_set(this, &
    pet_a_forest, &
    pet_a_impervious, &
    pet_a_pervious, &
    pet_b, &
    pet_c, &
    errmsg) result(status)

    class(nml_pet_m1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: pet_a_forest !< PET correction coefficient a for forest
    type(parameter_t), intent(in) :: pet_a_impervious !< PET correction coefficient a for impervious areas
    type(parameter_t), intent(in) :: pet_a_pervious !< PET correction coefficient a for pervious areas
    type(parameter_t), intent(in) :: pet_b !< PET correction coefficient b
    type(parameter_t), intent(in) :: pet_c !< PET correction coefficient c

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%pet_a_forest = pet_a_forest
    this%pet_a_impervious = pet_a_impervious
    this%pet_a_pervious = pet_a_pervious
    this%pet_b = pet_b
    this%pet_c = pet_c

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_m1_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet_m1_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet_m1_t), intent(in) :: this !< namelist instance
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
    case ("pet_a_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_forest'"
        return
      end if
      if (ieee_is_nan(this%pet_a_forest%value)) status = NML_ERR_NOT_SET
    case ("pet_a_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_forest'"
        return
      end if
    case ("pet_a_forest%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_forest'"
        return
      end if
    case ("pet_a_forest%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_forest'"
        return
      end if
    case ("pet_a_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_forest'"
        return
      end if
      if (ieee_is_nan(this%pet_a_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("pet_a_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_impervious'"
        return
      end if
      if (ieee_is_nan(this%pet_a_impervious%value)) status = NML_ERR_NOT_SET
    case ("pet_a_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_impervious'"
        return
      end if
    case ("pet_a_impervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_impervious'"
        return
      end if
    case ("pet_a_impervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_impervious'"
        return
      end if
    case ("pet_a_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_impervious'"
        return
      end if
      if (ieee_is_nan(this%pet_a_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("pet_a_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_pervious'"
        return
      end if
      if (ieee_is_nan(this%pet_a_pervious%value)) status = NML_ERR_NOT_SET
    case ("pet_a_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_pervious'"
        return
      end if
    case ("pet_a_pervious%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_pervious'"
        return
      end if
    case ("pet_a_pervious%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_pervious'"
        return
      end if
    case ("pet_a_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_a_pervious'"
        return
      end if
      if (ieee_is_nan(this%pet_a_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("pet_b%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_b'"
        return
      end if
      if (ieee_is_nan(this%pet_b%value)) status = NML_ERR_NOT_SET
    case ("pet_b%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_b'"
        return
      end if
    case ("pet_b%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_b'"
        return
      end if
    case ("pet_b%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_b'"
        return
      end if
    case ("pet_b")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_b'"
        return
      end if
      if (ieee_is_nan(this%pet_b%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("pet_c%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_c'"
        return
      end if
      if (ieee_is_nan(this%pet_c%value)) status = NML_ERR_NOT_SET
    case ("pet_c%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_c'"
        return
      end if
    case ("pet_c%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_c'"
        return
      end if
    case ("pet_c%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_c'"
        return
      end if
    case ("pet_c")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet_c'"
        return
      end if
      if (ieee_is_nan(this%pet_c%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet_m1_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet_m1_is_valid(this, errmsg) result(status)
    class(nml_pet_m1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("pet_a_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pet_a_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("pet_a_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pet_a_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("pet_a_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pet_a_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("pet_b", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pet_b"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("pet_c", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: pet_c"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_pet_m1_is_valid

end module nml_pet_m1
