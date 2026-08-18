!> \file nml_parameter_pet_m2.f90
!> \copydoc nml_pet_m2

!> \brief PET - Case -2
!> \details Parameters for aspect correction of externally supplied PET.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet_m2
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

  !> \class nml_pet_m2_t
  !> \brief PET - Case -2
  !> \details Parameters for aspect correction of externally supplied PET.
  type, public :: nml_pet_m2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: correction_factor_min !< Minimum PET aspect-correction factor
    type(parameter_t) :: correction_factor_delta !< Delta added to the minimum PET aspect-correction factor
    type(parameter_t) :: aspect_threshold !< Aspect threshold
  contains
    procedure :: init => nml_pet_m2_init
    procedure :: init_type => nml_pet_m2_init_type
    procedure :: from_file => nml_pet_m2_from_file
    procedure :: set => nml_pet_m2_set
    procedure :: is_set => nml_pet_m2_is_set
    procedure :: is_valid => nml_pet_m2_is_valid
  end type nml_pet_m2_t

contains

  !> \brief Initialize defaults and sentinels for pet_m2
  integer function nml_pet_m2_init(this, errmsg) result(status)
    class(nml_pet_m2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      correction_factor_min=this%correction_factor_min, &
      correction_factor_delta=this%correction_factor_delta, &
      aspect_threshold=this%aspect_threshold, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_pet_m2_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_pet_m2_init_type(this, &
    correction_factor_min, &
    correction_factor_delta, &
    aspect_threshold, &
    errmsg) result(status)
    class(nml_pet_m2_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: correction_factor_min !< Minimum PET aspect-correction factor
    type(parameter_t), intent(inout), optional :: correction_factor_delta !< Delta added to the minimum PET aspect-correction factor
    type(parameter_t), intent(inout), optional :: aspect_threshold !< Aspect threshold
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(correction_factor_min)) then
      correction_factor_min%value = ieee_value(correction_factor_min%value, ieee_quiet_nan) ! sentinel for derived component value
      correction_factor_min%optimize = .false.
      correction_factor_min%min = ieee_value(correction_factor_min%min, ieee_quiet_nan) ! sentinel for derived component min
      correction_factor_min%max = ieee_value(correction_factor_min%max, ieee_quiet_nan) ! sentinel for derived component max
      correction_factor_min%min = 0.7_dp
      correction_factor_min%max = 1.3_dp
    end if
    if (present(correction_factor_delta)) then
      correction_factor_delta%value = ieee_value(correction_factor_delta%value, ieee_quiet_nan) ! sentinel for derived component value
      correction_factor_delta%optimize = .false.
      correction_factor_delta%min = ieee_value(correction_factor_delta%min, ieee_quiet_nan) ! sentinel for derived component min
      correction_factor_delta%max = ieee_value(correction_factor_delta%max, ieee_quiet_nan) ! sentinel for derived component max
      correction_factor_delta%min = 0.0_dp
      correction_factor_delta%max = 0.2_dp
    end if
    if (present(aspect_threshold)) then
      aspect_threshold%value = ieee_value(aspect_threshold%value, ieee_quiet_nan) ! sentinel for derived component value
      aspect_threshold%optimize = .false.
      aspect_threshold%min = ieee_value(aspect_threshold%min, ieee_quiet_nan) ! sentinel for derived component min
      aspect_threshold%max = ieee_value(aspect_threshold%max, ieee_quiet_nan) ! sentinel for derived component max
      aspect_threshold%min = 160.0_dp
      aspect_threshold%max = 200.0_dp
    end if
  end function nml_pet_m2_init_type


  !> \brief Read pet_m2 namelist from file
  integer function nml_pet_m2_from_file(this, file, errmsg) result(status)
    class(nml_pet_m2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: correction_factor_min
    type(parameter_t) :: correction_factor_delta
    type(parameter_t) :: aspect_threshold
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet_m2/ &
      correction_factor_min, &
      correction_factor_delta, &
      aspect_threshold

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    correction_factor_min = this%correction_factor_min
    correction_factor_delta = this%correction_factor_delta
    aspect_threshold = this%aspect_threshold

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("pet_m2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet_m2, iostat=iostat, iomsg=iomsg)
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
    this%correction_factor_min = correction_factor_min
    this%correction_factor_delta = correction_factor_delta
    this%aspect_threshold = aspect_threshold

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_m2_from_file

  !> \brief Set pet_m2 values
  integer function nml_pet_m2_set(this, &
    correction_factor_min, &
    correction_factor_delta, &
    aspect_threshold, &
    errmsg) result(status)

    class(nml_pet_m2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: correction_factor_min !< Minimum PET aspect-correction factor
    type(parameter_t), intent(in) :: correction_factor_delta !< Delta added to the minimum PET aspect-correction factor
    type(parameter_t), intent(in) :: aspect_threshold !< Aspect threshold

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%correction_factor_min = correction_factor_min
    this%correction_factor_delta = correction_factor_delta
    this%aspect_threshold = aspect_threshold

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_m2_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet_m2_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet_m2_t), intent(in) :: this !< namelist instance
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
    case ("correction_factor_min%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_min'"
        return
      end if
      if (ieee_is_nan(this%correction_factor_min%value)) status = NML_ERR_NOT_SET
    case ("correction_factor_min%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_min'"
        return
      end if
    case ("correction_factor_min%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_min'"
        return
      end if
    case ("correction_factor_min%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_min'"
        return
      end if
    case ("correction_factor_min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_min'"
        return
      end if
      if (ieee_is_nan(this%correction_factor_min%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("correction_factor_delta%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_delta'"
        return
      end if
      if (ieee_is_nan(this%correction_factor_delta%value)) status = NML_ERR_NOT_SET
    case ("correction_factor_delta%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_delta'"
        return
      end if
    case ("correction_factor_delta%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_delta'"
        return
      end if
    case ("correction_factor_delta%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_delta'"
        return
      end if
    case ("correction_factor_delta")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'correction_factor_delta'"
        return
      end if
      if (ieee_is_nan(this%correction_factor_delta%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("aspect_threshold%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'aspect_threshold'"
        return
      end if
      if (ieee_is_nan(this%aspect_threshold%value)) status = NML_ERR_NOT_SET
    case ("aspect_threshold%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'aspect_threshold'"
        return
      end if
    case ("aspect_threshold%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'aspect_threshold'"
        return
      end if
    case ("aspect_threshold%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'aspect_threshold'"
        return
      end if
    case ("aspect_threshold")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'aspect_threshold'"
        return
      end if
      if (ieee_is_nan(this%aspect_threshold%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet_m2_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet_m2_is_valid(this, errmsg) result(status)
    class(nml_pet_m2_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("correction_factor_min", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: correction_factor_min"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("correction_factor_delta", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: correction_factor_delta"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("aspect_threshold", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: aspect_threshold"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_pet_m2_is_valid

end module nml_pet_m2
