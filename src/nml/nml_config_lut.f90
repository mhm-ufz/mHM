!> \file nml_config_lut.f90
!> \copydoc nml_config_lut

!> \brief Look-up table configuration
!> \details Domain-indexed look-up table (LUT) paths provided as metadata independent of the
!! process containers. Attached to the exchange of each domain so that any container
!! (land-surface MPR, Input, mLM) can resolve the shared class-to-property tables
!! without duplicating the path in its own process configuration.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_lut
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

  !> \class nml_config_lut_t
  !> \brief Look-up table configuration
  !> \details Domain-indexed look-up table (LUT) paths provided as metadata independent of the
  !! process containers. Attached to the exchange of each domain so that any container
  !! (land-surface MPR, Input, mLM) can resolve the shared class-to-property tables
  !! without duplicating the path in its own process configuration.
  type, public :: nml_config_lut_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    character(len=buf), allocatable, dimension(:) :: soil_lut_path !< Soil look-up table path
    character(len=buf), allocatable, dimension(:) :: geo_lut_path !< Geology look-up table path
  contains
    procedure :: init => nml_config_lut_init
    procedure :: set_dims => nml_config_lut_set_dims
    procedure :: from_file => nml_config_lut_from_file
    procedure :: set => nml_config_lut_set
    procedure :: is_set => nml_config_lut_is_set
    procedure :: is_valid => nml_config_lut_is_valid
  end type nml_config_lut_t

contains

  !> \brief Initialize defaults and sentinels for config_lut
  integer function nml_config_lut_init(this, errmsg) result(status)
    class(nml_config_lut_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%soil_lut_path)) deallocate(this%soil_lut_path)
    allocate(character(len=buf) :: this%soil_lut_path(this%n_domains))
    if (allocated(this%geo_lut_path)) deallocate(this%geo_lut_path)
    allocate(character(len=buf) :: this%geo_lut_path(this%n_domains))

    ! sentinel values for required/optional parameters
    this%soil_lut_path = achar(0) ! sentinel for optional string array
    this%geo_lut_path = achar(0) ! sentinel for optional string array
  end function nml_config_lut_init

  !> \brief Reset runtime dimensions for config_lut
  integer function nml_config_lut_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_lut_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%soil_lut_path)) deallocate(this%soil_lut_path)
    if (allocated(this%geo_lut_path)) deallocate(this%geo_lut_path)
    this%is_configured = .false.
  end function nml_config_lut_set_dims


  !> \brief Read config_lut namelist from file
  integer function nml_config_lut_from_file(this, file, errmsg) result(status)
    class(nml_config_lut_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf), allocatable, dimension(:) :: soil_lut_path
    character(len=buf), allocatable, dimension(:) :: geo_lut_path
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_lut/ &
      soil_lut_path, &
      geo_lut_path

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(soil_lut_path)) deallocate(soil_lut_path)
    allocate(character(len=buf) :: soil_lut_path(this%n_domains))
    if (allocated(geo_lut_path)) deallocate(geo_lut_path)
    allocate(character(len=buf) :: geo_lut_path(this%n_domains))
    soil_lut_path = this%soil_lut_path
    geo_lut_path = this%geo_lut_path

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_lut", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_lut, iostat=iostat, iomsg=iomsg)
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
    this%soil_lut_path = soil_lut_path
    this%geo_lut_path = geo_lut_path

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_lut_from_file

  !> \brief Set config_lut values
  integer function nml_config_lut_set(this, &
    soil_lut_path, &
    geo_lut_path, &
    errmsg) result(status)

    class(nml_config_lut_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), dimension(:), intent(in), optional :: soil_lut_path !< Soil look-up table path
    character(len=*), dimension(:), intent(in), optional :: geo_lut_path !< Geology look-up table path
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(soil_lut_path)) then
      if (size(soil_lut_path, 1) > size(this%soil_lut_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_lut_path'"
        return
      end if
      lb__1 = lbound(this%soil_lut_path, 1)
      ub__1 = lb__1 + size(soil_lut_path, 1) - 1
      this%soil_lut_path(lb__1:ub__1) = soil_lut_path
    end if
    if (present(geo_lut_path)) then
      if (size(geo_lut_path, 1) > size(this%geo_lut_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'geo_lut_path'"
        return
      end if
      lb__1 = lbound(this%geo_lut_path, 1)
      ub__1 = lb__1 + size(geo_lut_path, 1) - 1
      this%geo_lut_path(lb__1:ub__1) = geo_lut_path
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_lut_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_lut_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_lut_t), intent(in) :: this !< namelist instance
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
    case ("soil_lut_path")
      if (.not. allocated(this%soil_lut_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_lut_path), ubound(this%soil_lut_path), &
          "soil_lut_path", errmsg)
        if (status /= NML_OK) return
        if (this%soil_lut_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%soil_lut_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("geo_lut_path")
      if (.not. allocated(this%geo_lut_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%geo_lut_path), ubound(this%geo_lut_path), &
          "geo_lut_path", errmsg)
        if (status /= NML_OK) return
        if (this%geo_lut_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%geo_lut_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_lut_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_lut_is_valid(this, errmsg) result(status)
    class(nml_config_lut_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_lut_is_valid

end module nml_config_lut
