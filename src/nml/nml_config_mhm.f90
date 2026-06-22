!> \file nml_config_mhm.f90
!> \copydoc nml_config_mhm

!> \brief mHM model configuration
!> \details Configuration for the mHM model setup.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_mhm
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
    n_domains__default, &
    buf
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  ! default values
  logical, parameter, public :: read_restart__default = .false.
  logical, parameter, public :: write_restart__default = .false.
  logical, parameter, public :: share_evap_coeff__default = .true.

  !> \class nml_config_mhm_t
  !> \brief mHM model configuration
  !> \details Configuration for the mHM model setup.
  type, public :: nml_config_mhm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    character(len=buf), allocatable, dimension(:) :: output_path !< Output path
    logical, allocatable, dimension(:) :: read_restart !< Read restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path !< Restart input path
    logical, allocatable, dimension(:) :: write_restart !< Write restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path !< Restart output path
    real(dp), allocatable, dimension(:, :) :: evap_coeff !< Evaporation coefficients
    logical :: share_evap_coeff !< Share evaporation coefficients between domains
  contains
    procedure :: init => nml_config_mhm_init
    procedure :: set_dims => nml_config_mhm_set_dims
    procedure :: from_file => nml_config_mhm_from_file
    procedure :: set => nml_config_mhm_set
    procedure :: is_set => nml_config_mhm_is_set
    procedure :: is_valid => nml_config_mhm_is_valid
  end type nml_config_mhm_t

contains

  !> \brief Initialize defaults and sentinels for config_mhm
  integer function nml_config_mhm_init(this, errmsg) result(status)
    class(nml_config_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%output_path)) deallocate(this%output_path)
    allocate(character(len=buf) :: this%output_path(this%n_domains))
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    allocate(this%read_restart(this%n_domains))
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    allocate(character(len=buf) :: this%restart_input_path(this%n_domains))
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    allocate(this%write_restart(this%n_domains))
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    allocate(character(len=buf) :: this%restart_output_path(this%n_domains))
    if (allocated(this%evap_coeff)) deallocate(this%evap_coeff)
    allocate(this%evap_coeff(12, this%n_domains))

    ! sentinel values for required/optional parameters
    this%output_path = achar(0) ! sentinel for optional string array
    this%restart_input_path = achar(0) ! sentinel for optional string array
    this%restart_output_path = achar(0) ! sentinel for optional string array
    this%evap_coeff = ieee_value(this%evap_coeff, ieee_quiet_nan) ! sentinel for optional real array
    ! default values
    this%read_restart = read_restart__default
    this%write_restart = write_restart__default
    this%share_evap_coeff = share_evap_coeff__default ! bool values always need a default
  end function nml_config_mhm_init

  !> \brief Reset runtime dimensions for config_mhm
  integer function nml_config_mhm_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_mhm_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer :: candidate__n_domains
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(n_domains)) then
      candidate__n_domains = n_domains
    else
      candidate__n_domains = n_domains__default
    end if
    if (candidate__n_domains <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_domains' must be positive"
      return
    end if
    this%n_domains = candidate__n_domains

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%output_path)) deallocate(this%output_path)
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    if (allocated(this%evap_coeff)) deallocate(this%evap_coeff)
    this%is_configured = .false.
  end function nml_config_mhm_set_dims


  !> \brief Read config_mhm namelist from file
  integer function nml_config_mhm_from_file(this, file, errmsg) result(status)
    class(nml_config_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf), allocatable, dimension(:) :: output_path
    logical, allocatable, dimension(:) :: read_restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path
    logical, allocatable, dimension(:) :: write_restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path
    real(dp), allocatable, dimension(:, :) :: evap_coeff
    logical :: share_evap_coeff
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_mhm/ &
      output_path, &
      read_restart, &
      restart_input_path, &
      write_restart, &
      restart_output_path, &
      evap_coeff, &
      share_evap_coeff

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(output_path)) deallocate(output_path)
    allocate(character(len=buf) :: output_path(this%n_domains))
    if (allocated(read_restart)) deallocate(read_restart)
    allocate(read_restart(this%n_domains))
    if (allocated(restart_input_path)) deallocate(restart_input_path)
    allocate(character(len=buf) :: restart_input_path(this%n_domains))
    if (allocated(write_restart)) deallocate(write_restart)
    allocate(write_restart(this%n_domains))
    if (allocated(restart_output_path)) deallocate(restart_output_path)
    allocate(character(len=buf) :: restart_output_path(this%n_domains))
    if (allocated(evap_coeff)) deallocate(evap_coeff)
    allocate(evap_coeff(12, this%n_domains))
    output_path = this%output_path
    read_restart = this%read_restart
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_output_path = this%restart_output_path
    evap_coeff = this%evap_coeff
    share_evap_coeff = this%share_evap_coeff

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_mhm", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_mhm, iostat=iostat, iomsg=iomsg)
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
    this%output_path = output_path
    this%read_restart = read_restart
    this%restart_input_path = restart_input_path
    this%write_restart = write_restart
    this%restart_output_path = restart_output_path
    this%evap_coeff = evap_coeff
    this%share_evap_coeff = share_evap_coeff

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mhm_from_file

  !> \brief Set config_mhm values
  integer function nml_config_mhm_set(this, &
    output_path, &
    read_restart, &
    restart_input_path, &
    write_restart, &
    restart_output_path, &
    evap_coeff, &
    share_evap_coeff, &
    errmsg) result(status)

    class(nml_config_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), dimension(:), intent(in), optional :: output_path !< Output path
    logical, dimension(:), intent(in), optional :: read_restart !< Read restart
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< Restart input path
    logical, dimension(:), intent(in), optional :: write_restart !< Write restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< Restart output path
    real(dp), dimension(:, :), intent(in), optional :: evap_coeff !< Evaporation coefficients
    logical, intent(in), optional :: share_evap_coeff !< Share evaporation coefficients between domains
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(output_path)) then
      if (size(output_path, 1) > size(this%output_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'output_path'"
        return
      end if
      lb__1 = lbound(this%output_path, 1)
      ub__1 = lb__1 + size(output_path, 1) - 1
      this%output_path(lb__1:ub__1) = output_path
    end if
    if (present(read_restart)) then
      if (size(read_restart, 1) > size(this%read_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'read_restart'"
        return
      end if
      lb__1 = lbound(this%read_restart, 1)
      ub__1 = lb__1 + size(read_restart, 1) - 1
      this%read_restart(lb__1:ub__1) = read_restart
    end if
    if (present(restart_input_path)) then
      if (size(restart_input_path, 1) > size(this%restart_input_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_input_path'"
        return
      end if
      lb__1 = lbound(this%restart_input_path, 1)
      ub__1 = lb__1 + size(restart_input_path, 1) - 1
      this%restart_input_path(lb__1:ub__1) = restart_input_path
    end if
    if (present(write_restart)) then
      if (size(write_restart, 1) > size(this%write_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'write_restart'"
        return
      end if
      lb__1 = lbound(this%write_restart, 1)
      ub__1 = lb__1 + size(write_restart, 1) - 1
      this%write_restart(lb__1:ub__1) = write_restart
    end if
    if (present(restart_output_path)) then
      if (size(restart_output_path, 1) > size(this%restart_output_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_output_path'"
        return
      end if
      lb__1 = lbound(this%restart_output_path, 1)
      ub__1 = lb__1 + size(restart_output_path, 1) - 1
      this%restart_output_path(lb__1:ub__1) = restart_output_path
    end if
    if (present(evap_coeff)) then
      if (size(evap_coeff, 1) > size(this%evap_coeff, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'evap_coeff'"
        return
      end if
      lb__1 = lbound(this%evap_coeff, 1)
      ub__1 = lb__1 + size(evap_coeff, 1) - 1
      if (size(evap_coeff, 2) > size(this%evap_coeff, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'evap_coeff'"
        return
      end if
      lb__2 = lbound(this%evap_coeff, 2)
      ub__2 = lb__2 + size(evap_coeff, 2) - 1
      this%evap_coeff(lb__1:ub__1, lb__2:ub__2) = evap_coeff
    end if
    if (present(share_evap_coeff)) this%share_evap_coeff = share_evap_coeff

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mhm_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_mhm_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_mhm_t), intent(in) :: this !< namelist instance
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
    case ("output_path")
      if (.not. allocated(this%output_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%output_path), ubound(this%output_path), &
          "output_path", errmsg)
        if (status /= NML_OK) return
        if (this%output_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%output_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("read_restart")
      if (.not. allocated(this%read_restart)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%read_restart), ubound(this%read_restart), &
          "read_restart", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("restart_input_path")
      if (.not. allocated(this%restart_input_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%restart_input_path), ubound(this%restart_input_path), &
          "restart_input_path", errmsg)
        if (status /= NML_OK) return
        if (this%restart_input_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%restart_input_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("write_restart")
      if (.not. allocated(this%write_restart)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%write_restart), ubound(this%write_restart), &
          "write_restart", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("restart_output_path")
      if (.not. allocated(this%restart_output_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%restart_output_path), ubound(this%restart_output_path), &
          "restart_output_path", errmsg)
        if (status /= NML_OK) return
        if (this%restart_output_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%restart_output_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("evap_coeff")
      if (.not. allocated(this%evap_coeff)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%evap_coeff), ubound(this%evap_coeff), &
          "evap_coeff", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%evap_coeff(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%evap_coeff))) status = NML_ERR_NOT_SET
      end if
    case ("share_evap_coeff")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'share_evap_coeff'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_mhm_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_mhm_is_valid(this, errmsg) result(status)
    class(nml_config_mhm_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_mhm_is_valid

end module nml_config_mhm
