!> \file nml_parameter_percolation_1.f90
!> \copydoc nml_percolation_1

!> \brief Percolation - Case 1
!> \details Parameters for percolation and karst recharge.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_percolation_1
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

  !> \class nml_percolation_1_t
  !> \brief Percolation - Case 1
  !> \details Parameters for percolation and karst recharge.
  type, public :: nml_percolation_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: recharge_coefficient !< Recharge coefficient
    type(parameter_t) :: karstic_recharge_factor !< Karstic recharge factor
  contains
    procedure :: init => nml_percolation_1_init
    procedure :: init_type => nml_percolation_1_init_type
    procedure :: from_file => nml_percolation_1_from_file
    procedure :: set => nml_percolation_1_set
    procedure :: is_set => nml_percolation_1_is_set
    procedure :: is_valid => nml_percolation_1_is_valid
  end type nml_percolation_1_t

contains

  !> \brief Initialize defaults and sentinels for percolation_1
  integer function nml_percolation_1_init(this, errmsg) result(status)
    class(nml_percolation_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      recharge_coefficient=this%recharge_coefficient, &
      karstic_recharge_factor=this%karstic_recharge_factor, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_percolation_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_percolation_1_init_type(this, &
    recharge_coefficient, &
    karstic_recharge_factor, &
    errmsg) result(status)
    class(nml_percolation_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: recharge_coefficient !< Recharge coefficient
    type(parameter_t), intent(inout), optional :: karstic_recharge_factor !< Karstic recharge factor
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(recharge_coefficient)) then
      recharge_coefficient%value = ieee_value(recharge_coefficient%value, ieee_quiet_nan) ! sentinel for derived component value
      recharge_coefficient%optimize = .false.
      recharge_coefficient%lower_bound = ieee_value(recharge_coefficient%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      recharge_coefficient%upper_bound = ieee_value(recharge_coefficient%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      recharge_coefficient%lower_bound = 0.0_dp
      recharge_coefficient%upper_bound = 50.0_dp
    end if
    if (present(karstic_recharge_factor)) then
      karstic_recharge_factor%value = ieee_value(karstic_recharge_factor%value, ieee_quiet_nan) ! sentinel for derived component value
      karstic_recharge_factor%optimize = .false.
      karstic_recharge_factor%lower_bound = ieee_value(karstic_recharge_factor%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      karstic_recharge_factor%upper_bound = ieee_value(karstic_recharge_factor%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      karstic_recharge_factor%lower_bound = -5.0_dp
      karstic_recharge_factor%upper_bound = 5.0_dp
    end if
  end function nml_percolation_1_init_type


  !> \brief Read percolation_1 namelist from file
  integer function nml_percolation_1_from_file(this, file, errmsg) result(status)
    class(nml_percolation_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: recharge_coefficient
    type(parameter_t) :: karstic_recharge_factor
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /percolation_1/ &
      recharge_coefficient, &
      karstic_recharge_factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    recharge_coefficient = this%recharge_coefficient
    karstic_recharge_factor = this%karstic_recharge_factor

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("percolation_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=percolation_1, iostat=iostat, iomsg=iomsg)
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
    this%recharge_coefficient = recharge_coefficient
    this%karstic_recharge_factor = karstic_recharge_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_percolation_1_from_file

  !> \brief Set percolation_1 values
  integer function nml_percolation_1_set(this, &
    recharge_coefficient, &
    karstic_recharge_factor, &
    errmsg) result(status)

    class(nml_percolation_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: recharge_coefficient !< Recharge coefficient
    type(parameter_t), intent(in) :: karstic_recharge_factor !< Karstic recharge factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%recharge_coefficient = recharge_coefficient
    this%karstic_recharge_factor = karstic_recharge_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_percolation_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_percolation_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_percolation_1_t), intent(in) :: this !< namelist instance
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
    case ("recharge_coefficient%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recharge_coefficient'"
        return
      end if
      if (ieee_is_nan(this%recharge_coefficient%value)) status = NML_ERR_NOT_SET
    case ("recharge_coefficient%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recharge_coefficient'"
        return
      end if
    case ("recharge_coefficient%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recharge_coefficient'"
        return
      end if
    case ("recharge_coefficient%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recharge_coefficient'"
        return
      end if
    case ("recharge_coefficient")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recharge_coefficient'"
        return
      end if
      if (ieee_is_nan(this%recharge_coefficient%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("karstic_recharge_factor%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'karstic_recharge_factor'"
        return
      end if
      if (ieee_is_nan(this%karstic_recharge_factor%value)) status = NML_ERR_NOT_SET
    case ("karstic_recharge_factor%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'karstic_recharge_factor'"
        return
      end if
    case ("karstic_recharge_factor%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'karstic_recharge_factor'"
        return
      end if
    case ("karstic_recharge_factor%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'karstic_recharge_factor'"
        return
      end if
    case ("karstic_recharge_factor")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'karstic_recharge_factor'"
        return
      end if
      if (ieee_is_nan(this%karstic_recharge_factor%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_percolation_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_percolation_1_is_valid(this, errmsg) result(status)
    class(nml_percolation_1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("recharge_coefficient", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: recharge_coefficient"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("karstic_recharge_factor", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: karstic_recharge_factor"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_percolation_1_is_valid

end module nml_percolation_1
