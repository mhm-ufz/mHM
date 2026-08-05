!> \file nml_config_time.f90
!> \copydoc nml_config_time

!> \brief Time configuration
!> \details Configuration for simulation and evaluation time periods in mHM.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_time
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
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  logical, parameter, public :: share_time_period__default = .false.
  integer(i4), parameter, public :: time_step__default = 1_i4
  logical, parameter, public :: share_time_step__default = .true.

  !> \class nml_config_time_t
  !> \brief Time configuration
  !> \details Configuration for simulation and evaluation time periods in mHM.
  type, public :: nml_config_time_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    character(len=buf), allocatable, dimension(:) :: sim_start !< Simulation start
    character(len=buf), allocatable, dimension(:) :: eval_start !< Evaluation start
    character(len=buf), allocatable, dimension(:) :: sim_end !< Simulation end
    logical :: share_time_period !< Share time period between domains
    integer(i4), allocatable, dimension(:) :: time_step !< Global model time step
    logical :: share_time_step !< Share time step between domains
  contains
    procedure :: init => nml_config_time_init
    procedure :: set_dims => nml_config_time_set_dims
    procedure :: from_file => nml_config_time_from_file
    procedure :: set => nml_config_time_set
    procedure :: is_set => nml_config_time_is_set
    procedure :: is_valid => nml_config_time_is_valid
  end type nml_config_time_t

contains

  !> \brief Initialize defaults and sentinels for config_time
  integer function nml_config_time_init(this, errmsg) result(status)
    class(nml_config_time_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%sim_start)) deallocate(this%sim_start)
    allocate(character(len=buf) :: this%sim_start(this%n_domains))
    if (allocated(this%eval_start)) deallocate(this%eval_start)
    allocate(character(len=buf) :: this%eval_start(this%n_domains))
    if (allocated(this%sim_end)) deallocate(this%sim_end)
    allocate(character(len=buf) :: this%sim_end(this%n_domains))
    if (allocated(this%time_step)) deallocate(this%time_step)
    allocate(this%time_step(this%n_domains))

    ! sentinel values for required/optional parameters
    this%sim_start = achar(0) ! sentinel for optional string array
    this%eval_start = achar(0) ! sentinel for optional string array
    this%sim_end = achar(0) ! sentinel for optional string array
    ! default values
    this%share_time_period = share_time_period__default ! bool values always need a default
    this%time_step = time_step__default
    this%share_time_step = share_time_step__default ! bool values always need a default
  end function nml_config_time_init

  !> \brief Reset runtime dimensions for config_time
  integer function nml_config_time_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_time_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%sim_start)) deallocate(this%sim_start)
    if (allocated(this%eval_start)) deallocate(this%eval_start)
    if (allocated(this%sim_end)) deallocate(this%sim_end)
    if (allocated(this%time_step)) deallocate(this%time_step)
    this%is_configured = .false.
  end function nml_config_time_set_dims


  !> \brief Read config_time namelist from file
  integer function nml_config_time_from_file(this, file, errmsg) result(status)
    class(nml_config_time_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf), allocatable, dimension(:) :: sim_start
    character(len=buf), allocatable, dimension(:) :: eval_start
    character(len=buf), allocatable, dimension(:) :: sim_end
    logical :: share_time_period
    integer(i4), allocatable, dimension(:) :: time_step
    logical :: share_time_step
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_time/ &
      sim_start, &
      eval_start, &
      sim_end, &
      share_time_period, &
      time_step, &
      share_time_step

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(sim_start)) deallocate(sim_start)
    allocate(character(len=buf) :: sim_start(this%n_domains))
    if (allocated(eval_start)) deallocate(eval_start)
    allocate(character(len=buf) :: eval_start(this%n_domains))
    if (allocated(sim_end)) deallocate(sim_end)
    allocate(character(len=buf) :: sim_end(this%n_domains))
    if (allocated(time_step)) deallocate(time_step)
    allocate(time_step(this%n_domains))
    sim_start = this%sim_start
    eval_start = this%eval_start
    sim_end = this%sim_end
    share_time_period = this%share_time_period
    time_step = this%time_step
    share_time_step = this%share_time_step

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_time", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_time, iostat=iostat, iomsg=iomsg)
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
    this%sim_start = sim_start
    this%eval_start = eval_start
    this%sim_end = sim_end
    this%share_time_period = share_time_period
    this%time_step = time_step
    this%share_time_step = share_time_step

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_time_from_file

  !> \brief Set config_time values
  integer function nml_config_time_set(this, &
    sim_start, &
    eval_start, &
    sim_end, &
    share_time_period, &
    time_step, &
    share_time_step, &
    errmsg) result(status)

    class(nml_config_time_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), dimension(:), intent(in), optional :: sim_start !< Simulation start
    character(len=*), dimension(:), intent(in), optional :: eval_start !< Evaluation start
    character(len=*), dimension(:), intent(in), optional :: sim_end !< Simulation end
    logical, intent(in), optional :: share_time_period !< Share time period between domains
    integer(i4), dimension(:), intent(in), optional :: time_step !< Global model time step
    logical, intent(in), optional :: share_time_step !< Share time step between domains
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(sim_start)) then
      if (size(sim_start, 1) > size(this%sim_start, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'sim_start'"
        return
      end if
      lb__1 = lbound(this%sim_start, 1)
      ub__1 = lb__1 + size(sim_start, 1) - 1
      this%sim_start(lb__1:ub__1) = sim_start
    end if
    if (present(eval_start)) then
      if (size(eval_start, 1) > size(this%eval_start, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'eval_start'"
        return
      end if
      lb__1 = lbound(this%eval_start, 1)
      ub__1 = lb__1 + size(eval_start, 1) - 1
      this%eval_start(lb__1:ub__1) = eval_start
    end if
    if (present(sim_end)) then
      if (size(sim_end, 1) > size(this%sim_end, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'sim_end'"
        return
      end if
      lb__1 = lbound(this%sim_end, 1)
      ub__1 = lb__1 + size(sim_end, 1) - 1
      this%sim_end(lb__1:ub__1) = sim_end
    end if
    if (present(share_time_period)) this%share_time_period = share_time_period
    if (present(time_step)) then
      if (size(time_step, 1) > size(this%time_step, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'time_step'"
        return
      end if
      lb__1 = lbound(this%time_step, 1)
      ub__1 = lb__1 + size(time_step, 1) - 1
      this%time_step(lb__1:ub__1) = time_step
    end if
    if (present(share_time_step)) this%share_time_step = share_time_step

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_time_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_time_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_time_t), intent(in) :: this !< namelist instance
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
    case ("sim_start")
      if (.not. allocated(this%sim_start)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%sim_start), ubound(this%sim_start), &
          "sim_start", errmsg)
        if (status /= NML_OK) return
        if (this%sim_start(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%sim_start == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("eval_start")
      if (.not. allocated(this%eval_start)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%eval_start), ubound(this%eval_start), &
          "eval_start", errmsg)
        if (status /= NML_OK) return
        if (this%eval_start(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%eval_start == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("sim_end")
      if (.not. allocated(this%sim_end)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%sim_end), ubound(this%sim_end), &
          "sim_end", errmsg)
        if (status /= NML_OK) return
        if (this%sim_end(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%sim_end == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("share_time_period")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'share_time_period'"
        return
      end if
    case ("time_step")
      if (.not. allocated(this%time_step)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%time_step), ubound(this%time_step), &
          "time_step", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("share_time_step")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'share_time_step'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_time_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_time_is_valid(this, errmsg) result(status)
    class(nml_config_time_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_time_is_valid

end module nml_config_time
