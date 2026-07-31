!> \file nml_parameter_pet_2.f90
!> \copydoc nml_pet_2

!> \brief PET - Case 2
!> \details Parameters for the Priestley-Taylor PET method.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_pet_2
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

  !> \class nml_pet_2_t
  !> \brief PET - Case 2
  !> \details Parameters for the Priestley-Taylor PET method.
  type, public :: nml_pet_2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: priestley_taylor_coefficient !< Priestley-Taylor coefficient
    type(parameter_t) :: priestley_taylor_lai_correction !< Priestley-Taylor LAI correction factor
  contains
    procedure :: init => nml_pet_2_init
    procedure :: init_type => nml_pet_2_init_type
    procedure :: from_file => nml_pet_2_from_file
    procedure :: set => nml_pet_2_set
    procedure :: is_set => nml_pet_2_is_set
    procedure :: is_valid => nml_pet_2_is_valid
  end type nml_pet_2_t

contains

  !> \brief Initialize defaults and sentinels for pet_2
  integer function nml_pet_2_init(this, errmsg) result(status)
    class(nml_pet_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      priestley_taylor_coefficient=this%priestley_taylor_coefficient, &
      priestley_taylor_lai_correction=this%priestley_taylor_lai_correction, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_pet_2_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_pet_2_init_type(this, &
    priestley_taylor_coefficient, &
    priestley_taylor_lai_correction, &
    errmsg) result(status)
    class(nml_pet_2_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: priestley_taylor_coefficient !< Priestley-Taylor coefficient
    type(parameter_t), intent(inout), optional :: priestley_taylor_lai_correction !< Priestley-Taylor LAI correction factor
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(priestley_taylor_coefficient)) then
      priestley_taylor_coefficient%value = ieee_value(priestley_taylor_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      priestley_taylor_coefficient%optimize = .false.
      priestley_taylor_coefficient%min = ieee_value(priestley_taylor_coefficient%min, ieee_quiet_nan) ! sentinel for derived component min
      priestley_taylor_coefficient%max = ieee_value(priestley_taylor_coefficient%max, ieee_quiet_nan) ! sentinel for derived component max
      priestley_taylor_coefficient%min = 0.75_dp
      priestley_taylor_coefficient%max = 1.75_dp
    end if
    if (present(priestley_taylor_lai_correction)) then
      priestley_taylor_lai_correction%value = ieee_value(priestley_taylor_lai_correction%value, ieee_quiet_nan) ! sentinel for derived component value
      priestley_taylor_lai_correction%optimize = .false.
      priestley_taylor_lai_correction%min = ieee_value(priestley_taylor_lai_correction%min, ieee_quiet_nan) ! sentinel for derived component min
      priestley_taylor_lai_correction%max = ieee_value(priestley_taylor_lai_correction%max, ieee_quiet_nan) ! sentinel for derived component max
      priestley_taylor_lai_correction%min = -0.5_dp
      priestley_taylor_lai_correction%max = 0.2_dp
    end if
  end function nml_pet_2_init_type


  !> \brief Read pet_2 namelist from file
  integer function nml_pet_2_from_file(this, file, errmsg) result(status)
    class(nml_pet_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: priestley_taylor_coefficient
    type(parameter_t) :: priestley_taylor_lai_correction
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /pet_2/ &
      priestley_taylor_coefficient, &
      priestley_taylor_lai_correction

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    priestley_taylor_coefficient = this%priestley_taylor_coefficient
    priestley_taylor_lai_correction = this%priestley_taylor_lai_correction

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("pet_2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=pet_2, iostat=iostat, iomsg=iomsg)
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
    this%priestley_taylor_coefficient = priestley_taylor_coefficient
    this%priestley_taylor_lai_correction = priestley_taylor_lai_correction

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_2_from_file

  !> \brief Set pet_2 values
  integer function nml_pet_2_set(this, &
    priestley_taylor_coefficient, &
    priestley_taylor_lai_correction, &
    errmsg) result(status)

    class(nml_pet_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: priestley_taylor_coefficient !< Priestley-Taylor coefficient
    type(parameter_t), intent(in) :: priestley_taylor_lai_correction !< Priestley-Taylor LAI correction factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%priestley_taylor_coefficient = priestley_taylor_coefficient
    this%priestley_taylor_lai_correction = priestley_taylor_lai_correction

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_pet_2_set

  !> \brief Check whether a namelist value was set
  integer function nml_pet_2_is_set(this, name, idx, errmsg) result(status)
    class(nml_pet_2_t), intent(in) :: this !< namelist instance
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
    case ("priestley_taylor_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_coefficient'"
        return
      end if
      if (ieee_is_nan(this%priestley_taylor_coefficient%value)) status = NML_ERR_NOT_SET
    case ("priestley_taylor_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_coefficient'"
        return
      end if
    case ("priestley_taylor_coefficient%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_coefficient'"
        return
      end if
    case ("priestley_taylor_coefficient%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_coefficient'"
        return
      end if
    case ("priestley_taylor_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_coefficient'"
        return
      end if
      if (ieee_is_nan(this%priestley_taylor_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("priestley_taylor_lai_correction%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_lai_correction'"
        return
      end if
      if (ieee_is_nan(this%priestley_taylor_lai_correction%value)) status = NML_ERR_NOT_SET
    case ("priestley_taylor_lai_correction%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_lai_correction'"
        return
      end if
    case ("priestley_taylor_lai_correction%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_lai_correction'"
        return
      end if
    case ("priestley_taylor_lai_correction%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_lai_correction'"
        return
      end if
    case ("priestley_taylor_lai_correction")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'priestley_taylor_lai_correction'"
        return
      end if
      if (ieee_is_nan(this%priestley_taylor_lai_correction%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_pet_2_is_set

  !> \brief Validate required values and constraints
  integer function nml_pet_2_is_valid(this, errmsg) result(status)
    class(nml_pet_2_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("priestley_taylor_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: priestley_taylor_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("priestley_taylor_lai_correction", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: priestley_taylor_lai_correction"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_pet_2_is_valid

end module nml_pet_2
