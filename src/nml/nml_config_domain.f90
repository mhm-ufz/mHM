!> \file nml_config_domain.f90
!> \copydoc nml_config_domain

!> \brief Domain configuration
!> \details Domain-indexed configuration for combined mHM runs with domain files.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_domain
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

  implicit none

  ! default values
  character(len=buf), parameter, public :: domain_nmls__default = "mhm.nml"

  !> \class nml_config_domain_t
  !> \brief Domain configuration
  !> \details Domain-indexed configuration for combined mHM runs with domain files.
  type, public :: nml_config_domain_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    character(len=buf), allocatable, dimension(:) :: domain_dirs !< Domain directories
    character(len=buf), allocatable, dimension(:) :: domain_nmls !< Domain namelists
  contains
    procedure :: init => nml_config_domain_init
    procedure :: set_dims => nml_config_domain_set_dims
    procedure :: from_file => nml_config_domain_from_file
    procedure :: set => nml_config_domain_set
    procedure :: is_set => nml_config_domain_is_set
    procedure :: is_valid => nml_config_domain_is_valid
  end type nml_config_domain_t

contains

  !> \brief Initialize defaults and sentinels for config_domain
  integer function nml_config_domain_init(this, errmsg) result(status)
    class(nml_config_domain_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%domain_dirs)) deallocate(this%domain_dirs)
    allocate(character(len=buf) :: this%domain_dirs(this%n_domains))
    if (allocated(this%domain_nmls)) deallocate(this%domain_nmls)
    allocate(character(len=buf) :: this%domain_nmls(this%n_domains))

    ! sentinel values for required/optional parameters
    this%domain_dirs = achar(0) ! sentinel for optional string array
    ! default values
    this%domain_nmls = domain_nmls__default
  end function nml_config_domain_init

  !> \brief Reset runtime dimensions for config_domain
  integer function nml_config_domain_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_domain_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%domain_dirs)) deallocate(this%domain_dirs)
    if (allocated(this%domain_nmls)) deallocate(this%domain_nmls)
    this%is_configured = .false.
  end function nml_config_domain_set_dims


  !> \brief Read config_domain namelist from file
  integer function nml_config_domain_from_file(this, file, errmsg) result(status)
    class(nml_config_domain_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf), allocatable, dimension(:) :: domain_dirs
    character(len=buf), allocatable, dimension(:) :: domain_nmls
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_domain/ &
      domain_dirs, &
      domain_nmls

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(domain_dirs)) deallocate(domain_dirs)
    allocate(character(len=buf) :: domain_dirs(this%n_domains))
    if (allocated(domain_nmls)) deallocate(domain_nmls)
    allocate(character(len=buf) :: domain_nmls(this%n_domains))
    domain_dirs = this%domain_dirs
    domain_nmls = this%domain_nmls

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_domain", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_domain, iostat=iostat, iomsg=iomsg)
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
    this%domain_dirs = domain_dirs
    this%domain_nmls = domain_nmls

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_domain_from_file

  !> \brief Set config_domain values
  integer function nml_config_domain_set(this, &
    domain_dirs, &
    domain_nmls, &
    errmsg) result(status)

    class(nml_config_domain_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), dimension(:), intent(in), optional :: domain_dirs !< Domain directories
    character(len=*), dimension(:), intent(in), optional :: domain_nmls !< Domain namelists
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(domain_dirs)) then
      if (size(domain_dirs, 1) > size(this%domain_dirs, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'domain_dirs'"
        return
      end if
      lb__1 = lbound(this%domain_dirs, 1)
      ub__1 = lb__1 + size(domain_dirs, 1) - 1
      this%domain_dirs(lb__1:ub__1) = domain_dirs
    end if
    if (present(domain_nmls)) then
      if (size(domain_nmls, 1) > size(this%domain_nmls, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'domain_nmls'"
        return
      end if
      lb__1 = lbound(this%domain_nmls, 1)
      ub__1 = lb__1 + size(domain_nmls, 1) - 1
      this%domain_nmls(lb__1:ub__1) = domain_nmls
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_domain_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_domain_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_domain_t), intent(in) :: this !< namelist instance
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
    case ("domain_dirs")
      if (.not. allocated(this%domain_dirs)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%domain_dirs), ubound(this%domain_dirs), &
          "domain_dirs", errmsg)
        if (status /= NML_OK) return
        if (this%domain_dirs(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%domain_dirs == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("domain_nmls")
      if (.not. allocated(this%domain_nmls)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%domain_nmls), ubound(this%domain_nmls), &
          "domain_nmls", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_domain_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_domain_is_valid(this, errmsg) result(status)
    class(nml_config_domain_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_domain_is_valid

end module nml_config_domain
