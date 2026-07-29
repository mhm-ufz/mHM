!> \file nml_config_resolution.f90
!> \copydoc nml_config_resolution

!> \brief Grid resolution configuration
!> \details Domain-indexed target resolutions for internally derived grids.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_resolution
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
    n_domains__default
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  ! default values
  real(dp), parameter, public :: hydro__default = 0.0_dp
  real(dp), parameter, public :: route__default = 0.0_dp

  ! bounds values
  real(dp), parameter, public :: hydro__min = 0.0_dp
  real(dp), parameter, public :: route__min = 0.0_dp

  !> \class nml_config_resolution_t
  !> \brief Grid resolution configuration
  !> \details Domain-indexed target resolutions for internally derived grids.
  type, public :: nml_config_resolution_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    real(dp), allocatable, dimension(:) :: hydro !< Hydrological grid resolution
    real(dp), allocatable, dimension(:) :: route !< Routing grid resolution
  contains
    procedure :: init => nml_config_resolution_init
    procedure :: set_dims => nml_config_resolution_set_dims
    procedure :: from_file => nml_config_resolution_from_file
    procedure :: set => nml_config_resolution_set
    procedure :: is_set => nml_config_resolution_is_set
    procedure :: is_valid => nml_config_resolution_is_valid
  end type nml_config_resolution_t

contains

  !> \brief Check whether a value is within bounds
  elemental logical function hydro__in_bounds(val, allow_missing) result(in_bounds)
    real(dp), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (ieee_is_nan(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < hydro__min) in_bounds = .false.
  end function hydro__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function route__in_bounds(val, allow_missing) result(in_bounds)
    real(dp), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (ieee_is_nan(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < route__min) in_bounds = .false.
  end function route__in_bounds

  !> \brief Initialize defaults and sentinels for config_resolution
  integer function nml_config_resolution_init(this, errmsg) result(status)
    class(nml_config_resolution_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%hydro)) deallocate(this%hydro)
    allocate(this%hydro(this%n_domains))
    if (allocated(this%route)) deallocate(this%route)
    allocate(this%route(this%n_domains))

    ! default values
    this%hydro = hydro__default
    this%route = route__default
  end function nml_config_resolution_init

  !> \brief Reset runtime dimensions for config_resolution
  integer function nml_config_resolution_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_resolution_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%hydro)) deallocate(this%hydro)
    if (allocated(this%route)) deallocate(this%route)
    this%is_configured = .false.
  end function nml_config_resolution_set_dims


  !> \brief Read config_resolution namelist from file
  integer function nml_config_resolution_from_file(this, file, errmsg) result(status)
    class(nml_config_resolution_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), allocatable, dimension(:) :: hydro
    real(dp), allocatable, dimension(:) :: route
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_resolution/ &
      hydro, &
      route

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(hydro)) deallocate(hydro)
    allocate(hydro(this%n_domains))
    if (allocated(route)) deallocate(route)
    allocate(route(this%n_domains))
    hydro = this%hydro
    route = this%route

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_resolution", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_resolution, iostat=iostat, iomsg=iomsg)
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
    this%hydro = hydro
    this%route = route

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_resolution_from_file

  !> \brief Set config_resolution values
  integer function nml_config_resolution_set(this, &
    hydro, &
    route, &
    errmsg) result(status)

    class(nml_config_resolution_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in), optional :: hydro !< Hydrological grid resolution
    real(dp), dimension(:), intent(in), optional :: route !< Routing grid resolution
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(hydro)) then
      if (size(hydro, 1) > size(this%hydro, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro'"
        return
      end if
      lb__1 = lbound(this%hydro, 1)
      ub__1 = lb__1 + size(hydro, 1) - 1
      this%hydro(lb__1:ub__1) = hydro
    end if
    if (present(route)) then
      if (size(route, 1) > size(this%route, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'route'"
        return
      end if
      lb__1 = lbound(this%route, 1)
      ub__1 = lb__1 + size(route, 1) - 1
      this%route(lb__1:ub__1) = route
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_resolution_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_resolution_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_resolution_t), intent(in) :: this !< namelist instance
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
    case ("hydro")
      if (.not. allocated(this%hydro)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro), ubound(this%hydro), &
          "hydro", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("route")
      if (.not. allocated(this%route)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%route), ubound(this%route), &
          "route", errmsg)
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
  end function nml_config_resolution_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_resolution_is_valid(this, errmsg) result(status)
    class(nml_config_resolution_t), intent(in) :: this !< namelist instance
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
    if (allocated(this%hydro)) then
    if (.not. all(hydro__in_bounds(this%hydro, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: hydro"
      return
    end if
    end if
    if (allocated(this%route)) then
    if (.not. all(route__in_bounds(this%route, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: route"
      return
    end if
    end if
  end function nml_config_resolution_is_valid

end module nml_config_resolution
