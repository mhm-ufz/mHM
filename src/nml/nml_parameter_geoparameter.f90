!> \file nml_parameter_geoparameter.f90
!> \copydoc nml_geoparameter

!> \brief Geological parameters
!> \details Parameters for geoparameter.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_geoparameter
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
    n_geo_units__default, &
    NML_ERR_PARTLY_SET
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  !> \class nml_geoparameter_t
  !> \brief Geological parameters
  !> \details Parameters for geoparameter.
  type, public :: nml_geoparameter_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_geo_units = n_geo_units__default !< runtime dimension for n_geo_units
    real(dp), allocatable, dimension(:, :) :: geoparam !< Geological parameters
  contains
    procedure :: init => nml_geoparameter_init
    procedure :: set_dims => nml_geoparameter_set_dims
    procedure :: from_file => nml_geoparameter_from_file
    procedure :: set => nml_geoparameter_set
    procedure :: is_set => nml_geoparameter_is_set
    procedure :: is_valid => nml_geoparameter_is_valid
  end type nml_geoparameter_t

contains

  !> \brief Initialize defaults and sentinels for geoparameter
  integer function nml_geoparameter_init(this, errmsg) result(status)
    class(nml_geoparameter_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%geoparam)) deallocate(this%geoparam)
    allocate(this%geoparam(5, this%n_geo_units))

    ! sentinel values for required/optional parameters
    this%geoparam = ieee_value(this%geoparam, ieee_quiet_nan) ! sentinel for required real array
  end function nml_geoparameter_init

  !> \brief Reset runtime dimensions for geoparameter
  integer function nml_geoparameter_set_dims(this, &
    n_geo_units, &
    errmsg) result(status)
    class(nml_geoparameter_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_geo_units !< runtime dimension override for n_geo_units
    integer :: candidate__n_geo_units
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(n_geo_units)) then
      candidate__n_geo_units = n_geo_units
    else
      candidate__n_geo_units = n_geo_units__default
    end if
    if (candidate__n_geo_units <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_geo_units' must be positive"
      return
    end if
    this%n_geo_units = candidate__n_geo_units

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%geoparam)) deallocate(this%geoparam)
    this%is_configured = .false.
  end function nml_geoparameter_set_dims


  !> \brief Read geoparameter namelist from file
  integer function nml_geoparameter_from_file(this, file, errmsg) result(status)
    class(nml_geoparameter_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), allocatable, dimension(:, :) :: geoparam
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /geoparameter/ &
      geoparam

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(geoparam)) deallocate(geoparam)
    allocate(geoparam(5, this%n_geo_units))
    geoparam = this%geoparam

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("geoparameter", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=geoparameter, iostat=iostat, iomsg=iomsg)
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
    this%geoparam = geoparam

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_geoparameter_from_file

  !> \brief Set geoparameter values
  integer function nml_geoparameter_set(this, &
    geoparam, &
    errmsg) result(status)

    class(nml_geoparameter_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:, :), intent(in) :: geoparam !< Geological parameters
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(geoparam, 1) > size(this%geoparam, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'geoparam'"
      return
    end if
    lb__1 = lbound(this%geoparam, 1)
    ub__1 = lb__1 + size(geoparam, 1) - 1
    if (size(geoparam, 2) > size(this%geoparam, 2)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'geoparam'"
      return
    end if
    lb__2 = lbound(this%geoparam, 2)
    ub__2 = lb__2 + size(geoparam, 2) - 1
    this%geoparam(lb__1:ub__1, lb__2:ub__2) = geoparam

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_geoparameter_set

  !> \brief Check whether a namelist value was set
  integer function nml_geoparameter_is_set(this, name, idx, errmsg) result(status)
    class(nml_geoparameter_t), intent(in) :: this !< namelist instance
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
    case ("geoparam")
      if (.not. allocated(this%geoparam)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%geoparam), ubound(this%geoparam), &
          "GeoParam", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%geoparam(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%geoparam))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_geoparameter_is_set

  !> \brief Validate required values and constraints
  integer function nml_geoparameter_is_valid(this, errmsg) result(status)
    class(nml_geoparameter_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required arrays
    if (.not. allocated(this%geoparam)) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: GeoParam"
      return
    end if
    if (all(ieee_is_nan(this%geoparam(:, :)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: GeoParam"
      return
    end if
    if (any(ieee_is_nan(this%geoparam(:, :)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: GeoParam"
      return
    end if
  end function nml_geoparameter_is_valid

end module nml_geoparameter
