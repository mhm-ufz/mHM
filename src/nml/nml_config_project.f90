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
    buf
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  character(len=buf), parameter, public :: project_details__default = "mHM project"
  character(len=buf), parameter, public :: setup_description__default = "Model run"
  character(len=buf), parameter, public :: simulation_type__default = "Simulation"
  character(len=buf), parameter, public :: conventions__default = "None"
  character(len=buf), parameter, public :: contact__default = "Developer"
  character(len=buf), parameter, public :: mhm_details__default = "Research unit"
  character(len=buf), parameter, public :: history__default = "Model run version 1"
  integer(i4), parameter, public :: n_domains__default = 1_i4
  integer(i4), parameter, public :: n_geo_units__default = 25_i4
  integer(i4), parameter, public :: max_layers__default = 10_i4
  logical, parameter, public :: read_domains_from_dirs__default = .false.

  ! bounds values
  integer(i4), parameter, public :: max_layers__min = 1_i4

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
    integer(i4) :: n_geo_units !< Number of geological units
    integer(i4) :: max_layers !< Maximum number of soil layers
    logical :: read_domains_from_dirs !< Flag for separate domains
  contains
    procedure :: init => nml_config_project_init
    procedure :: from_file => nml_config_project_from_file
    procedure :: set => nml_config_project_set
    procedure :: is_set => nml_config_project_is_set
    procedure :: is_valid => nml_config_project_is_valid
  end type nml_config_project_t

contains

  !> \brief Check whether a value is within bounds
  elemental logical function max_layers__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < max_layers__min) in_bounds = .false.
  end function max_layers__in_bounds

  !> \brief Initialize defaults and sentinels for config_project
  integer function nml_config_project_init(this, errmsg) result(status)
    class(nml_config_project_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! default values
    this%project_details = project_details__default
    this%setup_description = setup_description__default
    this%simulation_type = simulation_type__default
    this%conventions = conventions__default
    this%contact = contact__default
    this%mhm_details = mhm_details__default
    this%history = history__default
    this%n_domains = n_domains__default
    this%n_geo_units = n_geo_units__default
    this%max_layers = max_layers__default
    this%read_domains_from_dirs = read_domains_from_dirs__default ! bool values always need a default
  end function nml_config_project_init


  !> \brief Read config_project namelist from file
  integer function nml_config_project_from_file(this, file, errmsg) result(status)
    class(nml_config_project_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf) :: project_details
    character(len=buf) :: setup_description
    character(len=buf) :: simulation_type
    character(len=buf) :: conventions
    character(len=buf) :: contact
    character(len=buf) :: mhm_details
    character(len=buf) :: history
    integer(i4) :: n_domains
    integer(i4) :: n_geo_units
    integer(i4) :: max_layers
    logical :: read_domains_from_dirs
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
      n_geo_units, &
      max_layers, &
      read_domains_from_dirs

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
    n_geo_units = this%n_geo_units
    max_layers = this%max_layers
    read_domains_from_dirs = this%read_domains_from_dirs

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
    this%n_geo_units = n_geo_units
    this%max_layers = max_layers
    this%read_domains_from_dirs = read_domains_from_dirs

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
    n_geo_units, &
    max_layers, &
    read_domains_from_dirs, &
    errmsg) result(status)

    class(nml_config_project_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), intent(in), optional :: project_details !< Project name
    character(len=*), intent(in), optional :: setup_description !< Description of the setup
    character(len=*), intent(in), optional :: simulation_type !< Type of simulation
    character(len=*), intent(in), optional :: conventions !< Convention used for dataset
    character(len=*), intent(in), optional :: contact !< Contact details, incl. PI name, modellers
    character(len=*), intent(in), optional :: mhm_details !< Developing institution
    character(len=*), intent(in), optional :: history !< Some details on data/model run version.
    integer(i4), intent(in), optional :: n_domains !< Number of domains
    integer(i4), intent(in), optional :: n_geo_units !< Number of geological units
    integer(i4), intent(in), optional :: max_layers !< Maximum number of soil layers
    logical, intent(in), optional :: read_domains_from_dirs !< Flag for separate domains

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
    if (present(n_geo_units)) this%n_geo_units = n_geo_units
    if (present(max_layers)) this%max_layers = max_layers
    if (present(read_domains_from_dirs)) this%read_domains_from_dirs = read_domains_from_dirs

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_project_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_project_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_project_t), intent(in) :: this !< namelist instance
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
    case ("n_geo_units")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'n_geo_units'"
        return
      end if
    case ("max_layers")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'max_layers'"
        return
      end if
    case ("read_domains_from_dirs")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'read_domains_from_dirs'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_project_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_project_is_valid(this, errmsg) result(status)
    class(nml_config_project_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! bounds constraints
    istat = this%is_set("max_layers", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. max_layers__in_bounds(this%max_layers)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: max_layers"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
  end function nml_config_project_is_valid

end module nml_config_project
