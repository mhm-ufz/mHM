!> \file nml_output_mrm.f90
!> \copydoc nml_output_mrm

!> \brief mRM output configuration
!> \details Output configuration for mRM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_output_mrm
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
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  integer(i4), parameter, public :: output_deflate_level__default = 6_i4
  logical, parameter, public :: output_double_precision__default = .false.
  integer(i4), parameter, public :: output_time_reference__default = 2_i4
  integer(i4), parameter, public :: output_frequency__default = -1_i4
  logical, parameter, public :: out_qrouted__default = .false.
  logical, parameter, public :: out_rivtemp__default = .false.

  ! enum values
  integer(i4), parameter, public :: output_time_reference__enum_values(3) = [0_i4, 1_i4, 2_i4]

  ! bounds values
  integer(i4), parameter, public :: output_deflate_level__min = 0_i4
  integer(i4), parameter, public :: output_deflate_level__max = 9_i4
  integer(i4), parameter, public :: output_frequency__min = -3_i4

  !> \class nml_output_mrm_t
  !> \brief mRM output configuration
  !> \details Output configuration for mRM.
  type, public :: nml_output_mrm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer(i4) :: output_deflate_level !< Output deflate level
    logical :: output_double_precision !< Output double precision
    integer(i4) :: output_time_reference !< Output time reference
    integer(i4) :: output_frequency !< Output time step
    logical :: out_qrouted !< Routed Streamflow
    logical :: out_rivtemp !< Routed Temperature
  contains
    procedure :: init => nml_output_mrm_init
    procedure :: from_file => nml_output_mrm_from_file
    procedure :: set => nml_output_mrm_set
    procedure :: is_set => nml_output_mrm_is_set
    procedure :: is_valid => nml_output_mrm_is_valid
  end type nml_output_mrm_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function output_time_reference__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == output_time_reference__enum_values)
  end function output_time_reference__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function output_deflate_level__in_bounds(val, allow_missing) result(in_bounds)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < output_deflate_level__min) in_bounds = .false.
    if (val > output_deflate_level__max) in_bounds = .false.
  end function output_deflate_level__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function output_frequency__in_bounds(val, allow_missing) result(in_bounds)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < output_frequency__min) in_bounds = .false.
  end function output_frequency__in_bounds

  !> \brief Initialize defaults and sentinels for output_mrm
  integer function nml_output_mrm_init(this, errmsg) result(status)
    class(nml_output_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! default values
    this%output_deflate_level = output_deflate_level__default
    this%output_double_precision = output_double_precision__default ! bool values always need a default
    this%output_time_reference = output_time_reference__default
    this%output_frequency = output_frequency__default
    this%out_qrouted = out_qrouted__default ! bool values always need a default
    this%out_rivtemp = out_rivtemp__default ! bool values always need a default
  end function nml_output_mrm_init


  !> \brief Read output_mrm namelist from file
  integer function nml_output_mrm_from_file(this, file, errmsg) result(status)
    class(nml_output_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4) :: output_deflate_level
    logical :: output_double_precision
    integer(i4) :: output_time_reference
    integer(i4) :: output_frequency
    logical :: out_qrouted
    logical :: out_rivtemp
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /output_mrm/ &
      output_deflate_level, &
      output_double_precision, &
      output_time_reference, &
      output_frequency, &
      out_qrouted, &
      out_rivtemp

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    output_deflate_level = this%output_deflate_level
    output_double_precision = this%output_double_precision
    output_time_reference = this%output_time_reference
    output_frequency = this%output_frequency
    out_qrouted = this%out_qrouted
    out_rivtemp = this%out_rivtemp

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("output_mrm", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=output_mrm, iostat=iostat, iomsg=iomsg)
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
    this%output_deflate_level = output_deflate_level
    this%output_double_precision = output_double_precision
    this%output_time_reference = output_time_reference
    this%output_frequency = output_frequency
    this%out_qrouted = out_qrouted
    this%out_rivtemp = out_rivtemp

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_output_mrm_from_file

  !> \brief Set output_mrm values
  integer function nml_output_mrm_set(this, &
    output_deflate_level, &
    output_double_precision, &
    output_time_reference, &
    output_frequency, &
    out_qrouted, &
    out_rivtemp, &
    errmsg) result(status)

    class(nml_output_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), intent(in), optional :: output_deflate_level !< Output deflate level
    logical, intent(in), optional :: output_double_precision !< Output double precision
    integer(i4), intent(in), optional :: output_time_reference !< Output time reference
    integer(i4), intent(in), optional :: output_frequency !< Output time step
    logical, intent(in), optional :: out_qrouted !< Routed Streamflow
    logical, intent(in), optional :: out_rivtemp !< Routed Temperature

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(output_deflate_level)) this%output_deflate_level = output_deflate_level
    if (present(output_double_precision)) this%output_double_precision = output_double_precision
    if (present(output_time_reference)) this%output_time_reference = output_time_reference
    if (present(output_frequency)) this%output_frequency = output_frequency
    if (present(out_qrouted)) this%out_qrouted = out_qrouted
    if (present(out_rivtemp)) this%out_rivtemp = out_rivtemp

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_output_mrm_set

  !> \brief Check whether a namelist value was set
  integer function nml_output_mrm_is_set(this, name, idx, errmsg) result(status)
    class(nml_output_mrm_t), intent(in) :: this !< namelist instance
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
    case ("output_deflate_level")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_deflate_level'"
        return
      end if
    case ("output_double_precision")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_double_precision'"
        return
      end if
    case ("output_time_reference")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_time_reference'"
        return
      end if
    case ("output_frequency")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_frequency'"
        return
      end if
    case ("out_qrouted")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_Qrouted'"
        return
      end if
    case ("out_rivtemp")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_RivTemp'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_output_mrm_is_set

  !> \brief Validate required values and constraints
  integer function nml_output_mrm_is_valid(this, errmsg) result(status)
    class(nml_output_mrm_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! enum constraints
    istat = this%is_set("output_time_reference", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_time_reference__in_enum(this%output_time_reference)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: output_time_reference"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    ! bounds constraints
    istat = this%is_set("output_deflate_level", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_deflate_level__in_bounds(this%output_deflate_level)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: output_deflate_level"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("output_frequency", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_frequency__in_bounds(this%output_frequency)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: output_frequency"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
  end function nml_output_mrm_is_valid

end module nml_output_mrm
