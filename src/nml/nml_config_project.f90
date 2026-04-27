!> \file nml_config_project.f90
!> \copydoc nml_config_project

!> \brief Project configuration
!> \details Configuration for the overall project setup in mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_project
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
    buf, &
    max_domains, &
    NML_ERR_PARTLY_SET
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  character(len=buf), parameter, public :: project_details_default = "mHM project"
  character(len=buf), parameter, public :: setup_description_default = "Model run"
  character(len=buf), parameter, public :: simulation_type_default = "Simulation"
  character(len=buf), parameter, public :: conventions_default = "None"
  character(len=buf), parameter, public :: contact_default = "Developer"
  character(len=buf), parameter, public :: mhm_details_default = "Research unit"
  character(len=buf), parameter, public :: history_default = "Model run version 1"
  integer(i4), parameter, public :: n_domains_default = 1_i4
  logical, parameter, public :: read_domains_from_dirs_default = .false.
  character(len=buf), parameter, public :: domain_nmls_default = "mhm.nml"

  !> \class nml_config_project_t
  !> \brief Project configuration
  !> \details Configuration for the overall project setup in mHM.
  type, public :: nml_config_project_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    character(len=buf) :: project_details !< Project name
    character(len=buf) :: setup_description !< Description of the setup
    character(len=buf) :: simulation_type !< Type of simulation
    character(len=buf) :: conventions !< Convention used for dataset
    character(len=buf) :: contact !< Contact details, incl. PI name, modellers
    character(len=buf) :: mhm_details !< Developing institution
    character(len=buf) :: history !< Some details on data/model run version.
    integer(i4) :: n_domains !< Number of domains
    logical :: read_domains_from_dirs !< Flag for separate domains
    character(len=buf), dimension(max_domains) :: domain_dirs !< Domain directories
    character(len=buf), dimension(max_domains) :: domain_nmls !< Domain namelists
  contains
    procedure :: init => nml_config_project_init
    procedure :: from_file => nml_config_project_from_file
    procedure :: set => nml_config_project_set
    procedure :: is_set => nml_config_project_is_set
    procedure :: filled_shape => nml_config_project_filled_shape
    procedure :: is_valid => nml_config_project_is_valid
  end type nml_config_project_t

contains

  !> \brief Initialize defaults and sentinels for config_project
  integer function nml_config_project_init(this, errmsg) result(status)
    class(nml_config_project_t), intent(inout) :: this
    character(len=*), intent(out), optional :: errmsg

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%domain_dirs = repeat(achar(0), len(this%domain_dirs)) ! sentinel for optional string array
    ! default values
    this%project_details = project_details_default
    this%setup_description = setup_description_default
    this%simulation_type = simulation_type_default
    this%conventions = conventions_default
    this%contact = contact_default
    this%mhm_details = mhm_details_default
    this%history = history_default
    this%n_domains = n_domains_default
    this%read_domains_from_dirs = read_domains_from_dirs_default ! bool values always need a default
    this%domain_nmls = domain_nmls_default
  end function nml_config_project_init

  !> \brief Read config_project namelist from file
  integer function nml_config_project_from_file(this, file, errmsg) result(status)
    class(nml_config_project_t), intent(inout) :: this
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg
    ! namelist variables
    character(len=buf) :: project_details
    character(len=buf) :: setup_description
    character(len=buf) :: simulation_type
    character(len=buf) :: conventions
    character(len=buf) :: contact
    character(len=buf) :: mhm_details
    character(len=buf) :: history
    integer(i4) :: n_domains
    logical :: read_domains_from_dirs
    character(len=buf), dimension(max_domains) :: domain_dirs
    character(len=buf), dimension(max_domains) :: domain_nmls
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_project/ &
      project_details, &
      setup_description, &
      simulation_type, &
      conventions, &
      contact, &
      mhm_details, &
      history, &
      n_domains, &
      read_domains_from_dirs, &
      domain_dirs, &
      domain_nmls

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    project_details = this%project_details
    setup_description = this%setup_description
    simulation_type = this%simulation_type
    conventions = this%conventions
    contact = this%contact
    mhm_details = this%mhm_details
    history = this%history
    n_domains = this%n_domains
    read_domains_from_dirs = this%read_domains_from_dirs
    domain_dirs = this%domain_dirs
    domain_nmls = this%domain_nmls

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_project", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_project, iostat=iostat, iomsg=iomsg)
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
    this%project_details = project_details
    this%setup_description = setup_description
    this%simulation_type = simulation_type
    this%conventions = conventions
    this%contact = contact
    this%mhm_details = mhm_details
    this%history = history
    this%n_domains = n_domains
    this%read_domains_from_dirs = read_domains_from_dirs
    this%domain_dirs = domain_dirs
    this%domain_nmls = domain_nmls

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_project_from_file

  !> \brief Set config_project values
  integer function nml_config_project_set(this, &
    project_details, &
    setup_description, &
    simulation_type, &
    conventions, &
    contact, &
    mhm_details, &
    history, &
    n_domains, &
    read_domains_from_dirs, &
    domain_dirs, &
    domain_nmls, &
    errmsg) result(status)

    class(nml_config_project_t), intent(inout) :: this
    character(len=*), intent(out), optional :: errmsg
    character(len=*), intent(in), optional :: project_details
    character(len=*), intent(in), optional :: setup_description
    character(len=*), intent(in), optional :: simulation_type
    character(len=*), intent(in), optional :: conventions
    character(len=*), intent(in), optional :: contact
    character(len=*), intent(in), optional :: mhm_details
    character(len=*), intent(in), optional :: history
    integer(i4), intent(in), optional :: n_domains
    logical, intent(in), optional :: read_domains_from_dirs
    character(len=*), dimension(:), intent(in), optional :: domain_dirs
    character(len=*), dimension(:), intent(in), optional :: domain_nmls
    integer :: &
      lb_1, &
      ub_1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(project_details)) this%project_details = project_details
    if (present(setup_description)) this%setup_description = setup_description
    if (present(simulation_type)) this%simulation_type = simulation_type
    if (present(conventions)) this%conventions = conventions
    if (present(contact)) this%contact = contact
    if (present(mhm_details)) this%mhm_details = mhm_details
    if (present(history)) this%history = history
    if (present(n_domains)) this%n_domains = n_domains
    if (present(read_domains_from_dirs)) this%read_domains_from_dirs = read_domains_from_dirs
    if (present(domain_dirs)) then
      if (size(domain_dirs, 1) > size(this%domain_dirs, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'domain_dirs'"
        return
      end if
      lb_1 = lbound(this%domain_dirs, 1)
      ub_1 = lb_1 + size(domain_dirs, 1) - 1
      this%domain_dirs(lb_1:ub_1) = domain_dirs
    end if
    if (present(domain_nmls)) then
      if (size(domain_nmls, 1) > size(this%domain_nmls, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'domain_nmls'"
        return
      end if
      lb_1 = lbound(this%domain_nmls, 1)
      ub_1 = lb_1 + size(domain_nmls, 1) - 1
      this%domain_nmls(lb_1:ub_1) = domain_nmls
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_project_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_project_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_project_t), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: idx(:)
    character(len=*), intent(out), optional :: errmsg

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    select case (to_lower(trim(name)))
    case ("project_details")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'project_details'"
        return
      end if
    case ("setup_description")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'setup_description'"
        return
      end if
    case ("simulation_type")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'simulation_type'"
        return
      end if
    case ("conventions")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'Conventions'"
        return
      end if
    case ("contact")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'contact'"
        return
      end if
    case ("mhm_details")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'mHM_details'"
        return
      end if
    case ("history")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'history'"
        return
      end if
    case ("n_domains")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'n_domains'"
        return
      end if
    case ("read_domains_from_dirs")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'read_domains_from_dirs'"
        return
      end if
    case ("domain_dirs")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%domain_dirs), ubound(this%domain_dirs), &
          "domain_dirs", errmsg)
        if (status /= NML_OK) return
        if (this%domain_dirs(idx(1)) == repeat(achar(0), len(this%domain_dirs))) status = NML_ERR_NOT_SET
      else
        if (all(this%domain_dirs == repeat(achar(0), len(this%domain_dirs)))) status = NML_ERR_NOT_SET
      end if
    case ("domain_nmls")
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
  end function nml_config_project_is_set

  !> \brief Determine the filled shape along flexible dimensions
  integer function nml_config_project_filled_shape(this, name, filled, errmsg) result(status)
    class(nml_config_project_t), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, intent(out) :: filled(:)
    character(len=*), intent(out), optional :: errmsg
    integer :: idx
    integer :: dim
    integer :: &
      lb_1, &
      ub_1

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    select case (to_lower(trim(name)))
    case ("domain_dirs")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'domain_dirs'"
        return
      end if
      do dim = 1, 1
        filled(dim) = size(this%domain_dirs, dim)
      end do
      filled(1) = 0
      do idx = ubound(this%domain_dirs, 1), &
        lbound(this%domain_dirs, 1), -1
        if (.not. (this%domain_dirs(idx) == repeat(achar(0), len(this%domain_dirs)))) then
          filled(1) = idx - lbound(this%domain_dirs, 1) + 1
          exit
        end if
      end do
      if (minval(filled) > 0) then
        lb_1 = lbound(this%domain_dirs, 1)
        ub_1 = lb_1 + filled(1) - 1
        if (any(this%domain_dirs(lb_1:ub_1) == repeat(achar(0), len(this%domain_dirs)))) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: domain_dirs"
          return
        end if
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "field is not a flexible array: " // trim(name)
    end select
  end function nml_config_project_filled_shape

  !> \brief Validate required values and constraints
  integer function nml_config_project_is_valid(this, errmsg) result(status)
    class(nml_config_project_t), intent(in) :: this
    character(len=*), intent(out), optional :: errmsg
    integer :: istat
    integer, allocatable :: filled(:)

    status = NML_OK
    if (present(errmsg)) errmsg = ""

    ! flexible arrays
    if (allocated(filled)) deallocate(filled)
    allocate(filled(1))
    istat = this%filled_shape("domain_dirs", filled, errmsg=errmsg)
    if (istat == NML_ERR_PARTLY_SET) then
      status = istat
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) errmsg = "array partly set: domain_dirs"
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_config_project_is_valid

end module nml_config_project
