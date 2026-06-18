!> \file nml_config_mpr.f90
!> \copydoc nml_config_mpr

!> \brief MPR configuration
!> \details Configuration for the multiscale parameter regionalization in mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_mpr
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
    max_layers__default, &
    buf
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4, &
    dp

  implicit none

  ! default values
  integer(i4), parameter, public :: soil_db_mode__default = 0_i4
  integer(i4), parameter, public :: soil_depth__default = 0_i4
  character(len=buf), parameter, public :: land_cover_var__default = "land_cover"
  character(len=buf), parameter, public :: lai_var__default = "lai"
  logical, parameter, public :: read_restart__default = .false.
  logical, parameter, public :: write_restart__default = .false.

  ! enum values
  integer(i4), parameter, public :: soil_db_mode__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: lai_time_step__enum_values(5) = [-3_i4, -2_i4, -1_i4, 0_i4, 1_i4]

  ! bounds values
  integer(i4), parameter, public :: n_layers__min = 1_i4

  !> \class nml_config_mpr_t
  !> \brief MPR configuration
  !> \details Configuration for the multiscale parameter regionalization in mHM.
  type, public :: nml_config_mpr_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer :: max_layers = max_layers__default !< runtime dimension for max_layers
    integer(i4), allocatable, dimension(:) :: soil_db_mode !< Soil database mode
    integer(i4), allocatable, dimension(:) :: tillage_depth !< Tillage depth
    integer(i4), allocatable, dimension(:) :: n_layers !< Number of soil layers
    integer(i4), allocatable, dimension(:, :) :: soil_depth !< Soil horizon depth
    real(dp), allocatable, dimension(:) :: fracsealed_cityarea !< Sealed fraction of city area
    character(len=buf), allocatable, dimension(:) :: land_cover_path !< Land cover path
    character(len=buf), allocatable, dimension(:) :: land_cover_var !< Land cover variable
    integer(i4), allocatable, dimension(:) :: lai_time_step !< LAI time step
    character(len=buf), allocatable, dimension(:) :: lai_path !< LAI path
    character(len=buf), allocatable, dimension(:) :: lai_var !< LAI variable
    character(len=buf), allocatable, dimension(:) :: soil_lut_path !< Soil LUT path
    character(len=buf), allocatable, dimension(:) :: geo_lut_path !< Geology LUT path
    character(len=buf), allocatable, dimension(:) :: lai_lut_path !< LAI LUT path
    logical, allocatable, dimension(:) :: read_restart !< Read restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path !< Restart input path
    logical, allocatable, dimension(:) :: write_restart !< Write restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path !< Restart output path
  contains
    procedure :: init => nml_config_mpr_init
    procedure :: set_dims => nml_config_mpr_set_dims
    procedure :: from_file => nml_config_mpr_from_file
    procedure :: set => nml_config_mpr_set
    procedure :: is_set => nml_config_mpr_is_set
    procedure :: is_valid => nml_config_mpr_is_valid
  end type nml_config_mpr_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function soil_db_mode__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == soil_db_mode__enum_values)
  end function soil_db_mode__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function lai_time_step__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == lai_time_step__enum_values)
  end function lai_time_step__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function n_layers__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < n_layers__min) in_bounds = .false.
  end function n_layers__in_bounds

  !> \brief Initialize defaults and sentinels for config_mpr
  integer function nml_config_mpr_init(this, errmsg) result(status)
    class(nml_config_mpr_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%soil_db_mode)) deallocate(this%soil_db_mode)
    allocate(this%soil_db_mode(this%n_domains))
    if (allocated(this%tillage_depth)) deallocate(this%tillage_depth)
    allocate(this%tillage_depth(this%n_domains))
    if (allocated(this%n_layers)) deallocate(this%n_layers)
    allocate(this%n_layers(this%n_domains))
    if (allocated(this%soil_depth)) deallocate(this%soil_depth)
    allocate(this%soil_depth(this%max_layers, this%n_domains))
    if (allocated(this%fracsealed_cityarea)) deallocate(this%fracsealed_cityarea)
    allocate(this%fracsealed_cityarea(this%n_domains))
    if (allocated(this%land_cover_path)) deallocate(this%land_cover_path)
    allocate(character(len=buf) :: this%land_cover_path(this%n_domains))
    if (allocated(this%land_cover_var)) deallocate(this%land_cover_var)
    allocate(character(len=buf) :: this%land_cover_var(this%n_domains))
    if (allocated(this%lai_time_step)) deallocate(this%lai_time_step)
    allocate(this%lai_time_step(this%n_domains))
    if (allocated(this%lai_path)) deallocate(this%lai_path)
    allocate(character(len=buf) :: this%lai_path(this%n_domains))
    if (allocated(this%lai_var)) deallocate(this%lai_var)
    allocate(character(len=buf) :: this%lai_var(this%n_domains))
    if (allocated(this%soil_lut_path)) deallocate(this%soil_lut_path)
    allocate(character(len=buf) :: this%soil_lut_path(this%n_domains))
    if (allocated(this%geo_lut_path)) deallocate(this%geo_lut_path)
    allocate(character(len=buf) :: this%geo_lut_path(this%n_domains))
    if (allocated(this%lai_lut_path)) deallocate(this%lai_lut_path)
    allocate(character(len=buf) :: this%lai_lut_path(this%n_domains))
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    allocate(this%read_restart(this%n_domains))
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    allocate(character(len=buf) :: this%restart_input_path(this%n_domains))
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    allocate(this%write_restart(this%n_domains))
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    allocate(character(len=buf) :: this%restart_output_path(this%n_domains))

    ! sentinel values for required/optional parameters
    this%tillage_depth = -huge(this%tillage_depth) ! sentinel for optional integer array
    this%n_layers = -huge(this%n_layers) ! sentinel for optional integer array
    this%fracsealed_cityarea = ieee_value(this%fracsealed_cityarea, ieee_quiet_nan) ! sentinel for optional real array
    this%land_cover_path = achar(0) ! sentinel for optional string array
    this%lai_time_step = -huge(this%lai_time_step) ! sentinel for optional integer array
    this%lai_path = achar(0) ! sentinel for optional string array
    this%soil_lut_path = achar(0) ! sentinel for optional string array
    this%geo_lut_path = achar(0) ! sentinel for optional string array
    this%lai_lut_path = achar(0) ! sentinel for optional string array
    this%restart_input_path = achar(0) ! sentinel for optional string array
    this%restart_output_path = achar(0) ! sentinel for optional string array
    ! default values
    this%soil_db_mode = soil_db_mode__default
    this%soil_depth = soil_depth__default
    this%land_cover_var = land_cover_var__default
    this%lai_var = lai_var__default
    this%read_restart = read_restart__default
    this%write_restart = write_restart__default
  end function nml_config_mpr_init

  !> \brief Reset runtime dimensions for config_mpr
  integer function nml_config_mpr_set_dims(this, &
    n_domains, &
    max_layers, &
    errmsg) result(status)
    class(nml_config_mpr_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer, intent(in), optional :: max_layers !< runtime dimension override for max_layers
    integer :: candidate__n_domains
    integer :: candidate__max_layers
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
    if (present(max_layers)) then
      candidate__max_layers = max_layers
    else
      candidate__max_layers = max_layers__default
    end if
    if (candidate__max_layers <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'max_layers' must be positive"
      return
    end if
    this%n_domains = candidate__n_domains
    this%max_layers = candidate__max_layers

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%soil_db_mode)) deallocate(this%soil_db_mode)
    if (allocated(this%tillage_depth)) deallocate(this%tillage_depth)
    if (allocated(this%n_layers)) deallocate(this%n_layers)
    if (allocated(this%soil_depth)) deallocate(this%soil_depth)
    if (allocated(this%fracsealed_cityarea)) deallocate(this%fracsealed_cityarea)
    if (allocated(this%land_cover_path)) deallocate(this%land_cover_path)
    if (allocated(this%land_cover_var)) deallocate(this%land_cover_var)
    if (allocated(this%lai_time_step)) deallocate(this%lai_time_step)
    if (allocated(this%lai_path)) deallocate(this%lai_path)
    if (allocated(this%lai_var)) deallocate(this%lai_var)
    if (allocated(this%soil_lut_path)) deallocate(this%soil_lut_path)
    if (allocated(this%geo_lut_path)) deallocate(this%geo_lut_path)
    if (allocated(this%lai_lut_path)) deallocate(this%lai_lut_path)
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    this%is_configured = .false.
  end function nml_config_mpr_set_dims


  !> \brief Read config_mpr namelist from file
  integer function nml_config_mpr_from_file(this, file, errmsg) result(status)
    class(nml_config_mpr_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4), allocatable, dimension(:) :: soil_db_mode
    integer(i4), allocatable, dimension(:) :: tillage_depth
    integer(i4), allocatable, dimension(:) :: n_layers
    integer(i4), allocatable, dimension(:, :) :: soil_depth
    real(dp), allocatable, dimension(:) :: fracsealed_cityarea
    character(len=buf), allocatable, dimension(:) :: land_cover_path
    character(len=buf), allocatable, dimension(:) :: land_cover_var
    integer(i4), allocatable, dimension(:) :: lai_time_step
    character(len=buf), allocatable, dimension(:) :: lai_path
    character(len=buf), allocatable, dimension(:) :: lai_var
    character(len=buf), allocatable, dimension(:) :: soil_lut_path
    character(len=buf), allocatable, dimension(:) :: geo_lut_path
    character(len=buf), allocatable, dimension(:) :: lai_lut_path
    logical, allocatable, dimension(:) :: read_restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path
    logical, allocatable, dimension(:) :: write_restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_mpr/ &
      soil_db_mode, &
      tillage_depth, &
      n_layers, &
      soil_depth, &
      fracsealed_cityarea, &
      land_cover_path, &
      land_cover_var, &
      lai_time_step, &
      lai_path, &
      lai_var, &
      soil_lut_path, &
      geo_lut_path, &
      lai_lut_path, &
      read_restart, &
      restart_input_path, &
      write_restart, &
      restart_output_path

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(soil_db_mode)) deallocate(soil_db_mode)
    allocate(soil_db_mode(this%n_domains))
    if (allocated(tillage_depth)) deallocate(tillage_depth)
    allocate(tillage_depth(this%n_domains))
    if (allocated(n_layers)) deallocate(n_layers)
    allocate(n_layers(this%n_domains))
    if (allocated(soil_depth)) deallocate(soil_depth)
    allocate(soil_depth(this%max_layers, this%n_domains))
    if (allocated(fracsealed_cityarea)) deallocate(fracsealed_cityarea)
    allocate(fracsealed_cityarea(this%n_domains))
    if (allocated(land_cover_path)) deallocate(land_cover_path)
    allocate(character(len=buf) :: land_cover_path(this%n_domains))
    if (allocated(land_cover_var)) deallocate(land_cover_var)
    allocate(character(len=buf) :: land_cover_var(this%n_domains))
    if (allocated(lai_time_step)) deallocate(lai_time_step)
    allocate(lai_time_step(this%n_domains))
    if (allocated(lai_path)) deallocate(lai_path)
    allocate(character(len=buf) :: lai_path(this%n_domains))
    if (allocated(lai_var)) deallocate(lai_var)
    allocate(character(len=buf) :: lai_var(this%n_domains))
    if (allocated(soil_lut_path)) deallocate(soil_lut_path)
    allocate(character(len=buf) :: soil_lut_path(this%n_domains))
    if (allocated(geo_lut_path)) deallocate(geo_lut_path)
    allocate(character(len=buf) :: geo_lut_path(this%n_domains))
    if (allocated(lai_lut_path)) deallocate(lai_lut_path)
    allocate(character(len=buf) :: lai_lut_path(this%n_domains))
    if (allocated(read_restart)) deallocate(read_restart)
    allocate(read_restart(this%n_domains))
    if (allocated(restart_input_path)) deallocate(restart_input_path)
    allocate(character(len=buf) :: restart_input_path(this%n_domains))
    if (allocated(write_restart)) deallocate(write_restart)
    allocate(write_restart(this%n_domains))
    if (allocated(restart_output_path)) deallocate(restart_output_path)
    allocate(character(len=buf) :: restart_output_path(this%n_domains))
    soil_db_mode = this%soil_db_mode
    tillage_depth = this%tillage_depth
    n_layers = this%n_layers
    soil_depth = this%soil_depth
    fracsealed_cityarea = this%fracsealed_cityarea
    land_cover_path = this%land_cover_path
    land_cover_var = this%land_cover_var
    lai_time_step = this%lai_time_step
    lai_path = this%lai_path
    lai_var = this%lai_var
    soil_lut_path = this%soil_lut_path
    geo_lut_path = this%geo_lut_path
    lai_lut_path = this%lai_lut_path
    read_restart = this%read_restart
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_output_path = this%restart_output_path

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_mpr", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_mpr, iostat=iostat, iomsg=iomsg)
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
    this%soil_db_mode = soil_db_mode
    this%tillage_depth = tillage_depth
    this%n_layers = n_layers
    this%soil_depth = soil_depth
    this%fracsealed_cityarea = fracsealed_cityarea
    this%land_cover_path = land_cover_path
    this%land_cover_var = land_cover_var
    this%lai_time_step = lai_time_step
    this%lai_path = lai_path
    this%lai_var = lai_var
    this%soil_lut_path = soil_lut_path
    this%geo_lut_path = geo_lut_path
    this%lai_lut_path = lai_lut_path
    this%read_restart = read_restart
    this%restart_input_path = restart_input_path
    this%write_restart = write_restart
    this%restart_output_path = restart_output_path

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mpr_from_file

  !> \brief Set config_mpr values
  integer function nml_config_mpr_set(this, &
    soil_db_mode, &
    tillage_depth, &
    n_layers, &
    soil_depth, &
    fracsealed_cityarea, &
    land_cover_path, &
    land_cover_var, &
    lai_time_step, &
    lai_path, &
    lai_var, &
    soil_lut_path, &
    geo_lut_path, &
    lai_lut_path, &
    read_restart, &
    restart_input_path, &
    write_restart, &
    restart_output_path, &
    errmsg) result(status)

    class(nml_config_mpr_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), dimension(:), intent(in), optional :: soil_db_mode !< Soil database mode
    integer(i4), dimension(:), intent(in), optional :: tillage_depth !< Tillage depth
    integer(i4), dimension(:), intent(in), optional :: n_layers !< Number of soil layers
    integer(i4), dimension(:, :), intent(in), optional :: soil_depth !< Soil horizon depth
    real(dp), dimension(:), intent(in), optional :: fracsealed_cityarea !< Sealed fraction of city area
    character(len=*), dimension(:), intent(in), optional :: land_cover_path !< Land cover path
    character(len=*), dimension(:), intent(in), optional :: land_cover_var !< Land cover variable
    integer(i4), dimension(:), intent(in), optional :: lai_time_step !< LAI time step
    character(len=*), dimension(:), intent(in), optional :: lai_path !< LAI path
    character(len=*), dimension(:), intent(in), optional :: lai_var !< LAI variable
    character(len=*), dimension(:), intent(in), optional :: soil_lut_path !< Soil LUT path
    character(len=*), dimension(:), intent(in), optional :: geo_lut_path !< Geology LUT path
    character(len=*), dimension(:), intent(in), optional :: lai_lut_path !< LAI LUT path
    logical, dimension(:), intent(in), optional :: read_restart !< Read restart
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< Restart input path
    logical, dimension(:), intent(in), optional :: write_restart !< Write restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< Restart output path
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(soil_db_mode)) then
      if (size(soil_db_mode, 1) > size(this%soil_db_mode, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_db_mode'"
        return
      end if
      lb__1 = lbound(this%soil_db_mode, 1)
      ub__1 = lb__1 + size(soil_db_mode, 1) - 1
      this%soil_db_mode(lb__1:ub__1) = soil_db_mode
    end if
    if (present(tillage_depth)) then
      if (size(tillage_depth, 1) > size(this%tillage_depth, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tillage_depth'"
        return
      end if
      lb__1 = lbound(this%tillage_depth, 1)
      ub__1 = lb__1 + size(tillage_depth, 1) - 1
      this%tillage_depth(lb__1:ub__1) = tillage_depth
    end if
    if (present(n_layers)) then
      if (size(n_layers, 1) > size(this%n_layers, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'n_layers'"
        return
      end if
      lb__1 = lbound(this%n_layers, 1)
      ub__1 = lb__1 + size(n_layers, 1) - 1
      this%n_layers(lb__1:ub__1) = n_layers
    end if
    if (present(soil_depth)) then
      if (size(soil_depth, 1) > size(this%soil_depth, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_depth'"
        return
      end if
      lb__1 = lbound(this%soil_depth, 1)
      ub__1 = lb__1 + size(soil_depth, 1) - 1
      if (size(soil_depth, 2) > size(this%soil_depth, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'soil_depth'"
        return
      end if
      lb__2 = lbound(this%soil_depth, 2)
      ub__2 = lb__2 + size(soil_depth, 2) - 1
      this%soil_depth(lb__1:ub__1, lb__2:ub__2) = soil_depth
    end if
    if (present(fracsealed_cityarea)) then
      if (size(fracsealed_cityarea, 1) > size(this%fracsealed_cityarea, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'fracsealed_cityarea'"
        return
      end if
      lb__1 = lbound(this%fracsealed_cityarea, 1)
      ub__1 = lb__1 + size(fracsealed_cityarea, 1) - 1
      this%fracsealed_cityarea(lb__1:ub__1) = fracsealed_cityarea
    end if
    if (present(land_cover_path)) then
      if (size(land_cover_path, 1) > size(this%land_cover_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'land_cover_path'"
        return
      end if
      lb__1 = lbound(this%land_cover_path, 1)
      ub__1 = lb__1 + size(land_cover_path, 1) - 1
      this%land_cover_path(lb__1:ub__1) = land_cover_path
    end if
    if (present(land_cover_var)) then
      if (size(land_cover_var, 1) > size(this%land_cover_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'land_cover_var'"
        return
      end if
      lb__1 = lbound(this%land_cover_var, 1)
      ub__1 = lb__1 + size(land_cover_var, 1) - 1
      this%land_cover_var(lb__1:ub__1) = land_cover_var
    end if
    if (present(lai_time_step)) then
      if (size(lai_time_step, 1) > size(this%lai_time_step, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_time_step'"
        return
      end if
      lb__1 = lbound(this%lai_time_step, 1)
      ub__1 = lb__1 + size(lai_time_step, 1) - 1
      this%lai_time_step(lb__1:ub__1) = lai_time_step
    end if
    if (present(lai_path)) then
      if (size(lai_path, 1) > size(this%lai_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_path'"
        return
      end if
      lb__1 = lbound(this%lai_path, 1)
      ub__1 = lb__1 + size(lai_path, 1) - 1
      this%lai_path(lb__1:ub__1) = lai_path
    end if
    if (present(lai_var)) then
      if (size(lai_var, 1) > size(this%lai_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_var'"
        return
      end if
      lb__1 = lbound(this%lai_var, 1)
      ub__1 = lb__1 + size(lai_var, 1) - 1
      this%lai_var(lb__1:ub__1) = lai_var
    end if
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
    if (present(lai_lut_path)) then
      if (size(lai_lut_path, 1) > size(this%lai_lut_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_lut_path'"
        return
      end if
      lb__1 = lbound(this%lai_lut_path, 1)
      ub__1 = lb__1 + size(lai_lut_path, 1) - 1
      this%lai_lut_path(lb__1:ub__1) = lai_lut_path
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

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mpr_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_mpr_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_mpr_t), intent(in) :: this !< namelist instance
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
    case ("soil_db_mode")
      if (.not. allocated(this%soil_db_mode)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_db_mode), ubound(this%soil_db_mode), &
          "soil_db_mode", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tillage_depth")
      if (.not. allocated(this%tillage_depth)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tillage_depth), ubound(this%tillage_depth), &
          "tillage_depth", errmsg)
        if (status /= NML_OK) return
        if (this%tillage_depth(idx(1)) == -huge(this%tillage_depth(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%tillage_depth == -huge(this%tillage_depth))) status = NML_ERR_NOT_SET
      end if
    case ("n_layers")
      if (.not. allocated(this%n_layers)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%n_layers), ubound(this%n_layers), &
          "n_layers", errmsg)
        if (status /= NML_OK) return
        if (this%n_layers(idx(1)) == -huge(this%n_layers(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%n_layers == -huge(this%n_layers))) status = NML_ERR_NOT_SET
      end if
    case ("soil_depth")
      if (.not. allocated(this%soil_depth)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_depth), ubound(this%soil_depth), &
          "soil_depth", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("fracsealed_cityarea")
      if (.not. allocated(this%fracsealed_cityarea)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%fracsealed_cityarea), ubound(this%fracsealed_cityarea), &
          "fracSealed_cityArea", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%fracsealed_cityarea(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%fracsealed_cityarea))) status = NML_ERR_NOT_SET
      end if
    case ("land_cover_path")
      if (.not. allocated(this%land_cover_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%land_cover_path), ubound(this%land_cover_path), &
          "land_cover_path", errmsg)
        if (status /= NML_OK) return
        if (this%land_cover_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%land_cover_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("land_cover_var")
      if (.not. allocated(this%land_cover_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%land_cover_var), ubound(this%land_cover_var), &
          "land_cover_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("lai_time_step")
      if (.not. allocated(this%lai_time_step)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_time_step), ubound(this%lai_time_step), &
          "lai_time_step", errmsg)
        if (status /= NML_OK) return
        if (this%lai_time_step(idx(1)) == -huge(this%lai_time_step(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_time_step == -huge(this%lai_time_step))) status = NML_ERR_NOT_SET
      end if
    case ("lai_path")
      if (.not. allocated(this%lai_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_path), ubound(this%lai_path), &
          "lai_path", errmsg)
        if (status /= NML_OK) return
        if (this%lai_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("lai_var")
      if (.not. allocated(this%lai_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_var), ubound(this%lai_var), &
          "lai_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
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
    case ("lai_lut_path")
      if (.not. allocated(this%lai_lut_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_lut_path), ubound(this%lai_lut_path), &
          "lai_lut_path", errmsg)
        if (status /= NML_OK) return
        if (this%lai_lut_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_lut_path == achar(0))) status = NML_ERR_NOT_SET
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
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_mpr_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_mpr_is_valid(this, errmsg) result(status)
    class(nml_config_mpr_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! enum constraints
    if (allocated(this%soil_db_mode)) then
    if (.not. all(soil_db_mode__in_enum(this%soil_db_mode, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: soil_db_mode"
      return
    end if
    end if
    if (allocated(this%lai_time_step)) then
    if (.not. all(lai_time_step__in_enum(this%lai_time_step, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: lai_time_step"
      return
    end if
    end if
    ! bounds constraints
    if (allocated(this%n_layers)) then
    if (.not. all(n_layers__in_bounds(this%n_layers, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: n_layers"
      return
    end if
    end if
  end function nml_config_mpr_is_valid

end module nml_config_mpr
