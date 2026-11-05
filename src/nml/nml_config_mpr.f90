!> \file nml_config_mpr.f90
!> \copydoc nml_config_mpr

!> \brief MPR configuration
!> \details Configuration for the multiscale parameter regionalization in mHM.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
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
    buf, &
    NML_ERR_PARTLY_SET
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
    procedure, private :: set_mask__rank__1 => nml_config_mpr_set_mask__rank__1
    procedure, private :: set_mask__rank__2 => nml_config_mpr_set_mask__rank__2
    generic :: set_mask => &
      set_mask__rank__1, &
      set_mask__rank__2
    procedure :: filled_shape => nml_config_mpr_filled_shape
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
    read_restart = this%read_restart
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_output_path = this%restart_output_path

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_mpr", errmsg=errmsg)
    if (status /= NML_OK) then
      if (status == NML_ERR_NML_NOT_FOUND) then
        close_status = nml%close(errmsg=errmsg)
        if (close_status /= NML_OK) then
          status = close_status
          return
        end if
        this%is_configured = .true.
        status = NML_OK
        return
      end if
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
    logical, dimension(:), intent(in), optional :: read_restart !< Read restart
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< Restart input path
    logical, dimension(:), intent(in), optional :: write_restart !< Write restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< Restart output path

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
      this%soil_db_mode(1:size(soil_db_mode, 1)) = soil_db_mode
    end if
    if (present(tillage_depth)) then
      if (size(tillage_depth, 1) > size(this%tillage_depth, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tillage_depth'"
        return
      end if
      this%tillage_depth(1:size(tillage_depth, 1)) = tillage_depth
    end if
    if (present(n_layers)) then
      if (size(n_layers, 1) > size(this%n_layers, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'n_layers'"
        return
      end if
      this%n_layers(1:size(n_layers, 1)) = n_layers
    end if
    if (present(soil_depth)) then
      if (size(soil_depth, 1) > size(this%soil_depth, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_depth'"
        return
      end if
      if (size(soil_depth, 2) > size(this%soil_depth, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'soil_depth'"
        return
      end if
      this%soil_depth(1:size(soil_depth, 1), 1:size(soil_depth, 2)) = soil_depth
    end if
    if (present(fracsealed_cityarea)) then
      if (size(fracsealed_cityarea, 1) > size(this%fracsealed_cityarea, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'fracsealed_cityarea'"
        return
      end if
      this%fracsealed_cityarea(1:size(fracsealed_cityarea, 1)) = fracsealed_cityarea
    end if
    if (present(land_cover_path)) then
      if (size(land_cover_path, 1) > size(this%land_cover_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'land_cover_path'"
        return
      end if
      this%land_cover_path(1:size(land_cover_path, 1)) = land_cover_path
    end if
    if (present(land_cover_var)) then
      if (size(land_cover_var, 1) > size(this%land_cover_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'land_cover_var'"
        return
      end if
      this%land_cover_var(1:size(land_cover_var, 1)) = land_cover_var
    end if
    if (present(lai_time_step)) then
      if (size(lai_time_step, 1) > size(this%lai_time_step, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_time_step'"
        return
      end if
      this%lai_time_step(1:size(lai_time_step, 1)) = lai_time_step
    end if
    if (present(lai_path)) then
      if (size(lai_path, 1) > size(this%lai_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_path'"
        return
      end if
      this%lai_path(1:size(lai_path, 1)) = lai_path
    end if
    if (present(lai_var)) then
      if (size(lai_var, 1) > size(this%lai_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_var'"
        return
      end if
      this%lai_var(1:size(lai_var, 1)) = lai_var
    end if
    if (present(read_restart)) then
      if (size(read_restart, 1) > size(this%read_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'read_restart'"
        return
      end if
      this%read_restart(1:size(read_restart, 1)) = read_restart
    end if
    if (present(restart_input_path)) then
      if (size(restart_input_path, 1) > size(this%restart_input_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_input_path'"
        return
      end if
      this%restart_input_path(1:size(restart_input_path, 1)) = restart_input_path
    end if
    if (present(write_restart)) then
      if (size(write_restart, 1) > size(this%write_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'write_restart'"
        return
      end if
      this%write_restart(1:size(write_restart, 1)) = write_restart
    end if
    if (present(restart_output_path)) then
      if (size(restart_output_path, 1) > size(this%restart_output_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_output_path'"
        return
      end if
      this%restart_output_path(1:size(restart_output_path, 1)) = restart_output_path
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
    intrinsic :: shape

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
        if (present(errmsg)) errmsg = "array storage unavailable: soil_db_mode"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%soil_db_mode), "soil_db_mode", errmsg)
        if (status /= NML_OK) return
      end if
    case ("tillage_depth")
      if (.not. allocated(this%tillage_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: tillage_depth"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%tillage_depth), "tillage_depth", errmsg)
        if (status /= NML_OK) return
        if (this%tillage_depth(idx(1)) == -huge(this%tillage_depth(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%tillage_depth == -huge(this%tillage_depth))) status = NML_ERR_NOT_SET
      end if
    case ("n_layers")
      if (.not. allocated(this%n_layers)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: n_layers"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%n_layers), "n_layers", errmsg)
        if (status /= NML_OK) return
        if (this%n_layers(idx(1)) == -huge(this%n_layers(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%n_layers == -huge(this%n_layers))) status = NML_ERR_NOT_SET
      end if
    case ("soil_depth")
      if (.not. allocated(this%soil_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: soil_depth"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%soil_depth), "soil_depth", errmsg)
        if (status /= NML_OK) return
      end if
    case ("fracsealed_cityarea")
      if (.not. allocated(this%fracsealed_cityarea)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: fracSealed_cityArea"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%fracsealed_cityarea), "fracSealed_cityArea", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%fracsealed_cityarea(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%fracsealed_cityarea))) status = NML_ERR_NOT_SET
      end if
    case ("land_cover_path")
      if (.not. allocated(this%land_cover_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_path"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%land_cover_path), "land_cover_path", errmsg)
        if (status /= NML_OK) return
        if (this%land_cover_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%land_cover_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("land_cover_var")
      if (.not. allocated(this%land_cover_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_var"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%land_cover_var), "land_cover_var", errmsg)
        if (status /= NML_OK) return
      end if
    case ("lai_time_step")
      if (.not. allocated(this%lai_time_step)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_time_step"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%lai_time_step), "lai_time_step", errmsg)
        if (status /= NML_OK) return
        if (this%lai_time_step(idx(1)) == -huge(this%lai_time_step(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_time_step == -huge(this%lai_time_step))) status = NML_ERR_NOT_SET
      end if
    case ("lai_path")
      if (.not. allocated(this%lai_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_path"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%lai_path), "lai_path", errmsg)
        if (status /= NML_OK) return
        if (this%lai_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("lai_var")
      if (.not. allocated(this%lai_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_var"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%lai_var), "lai_var", errmsg)
        if (status /= NML_OK) return
      end if
    case ("read_restart")
      if (.not. allocated(this%read_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: read_restart"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%read_restart), "read_restart", errmsg)
        if (status /= NML_OK) return
      end if
    case ("restart_input_path")
      if (.not. allocated(this%restart_input_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_input_path"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%restart_input_path), "restart_input_path", errmsg)
        if (status /= NML_OK) return
        if (this%restart_input_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%restart_input_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("write_restart")
      if (.not. allocated(this%write_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: write_restart"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%write_restart), "write_restart", errmsg)
        if (status /= NML_OK) return
      end if
    case ("restart_output_path")
      if (.not. allocated(this%restart_output_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_output_path"
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, shape(this%restart_output_path), "restart_output_path", errmsg)
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

  !> \brief Obtain the item presence mask with one field-name dispatch
  integer function nml_config_mpr_set_mask__rank__1(this, name, mask, errmsg) result(status)
    class(nml_config_mpr_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field or derived component path
    logical, allocatable, intent(out) :: mask(:) !< item presence
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: idx__1

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
        if (present(errmsg)) errmsg = "array storage unavailable: soil_db_mode"
        return
      end if
      allocate(mask( &
        size(this%soil_db_mode, 1)))
      mask = .true.
    case ("tillage_depth")
      if (.not. allocated(this%tillage_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: tillage_depth"
        return
      end if
      allocate(mask( &
        size(this%tillage_depth, 1)))
      do idx__1 = 1, size(this%tillage_depth, 1)
        if (this%tillage_depth(idx__1) == -huge(this%tillage_depth(idx__1))) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("n_layers")
      if (.not. allocated(this%n_layers)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: n_layers"
        return
      end if
      allocate(mask( &
        size(this%n_layers, 1)))
      do idx__1 = 1, size(this%n_layers, 1)
        if (this%n_layers(idx__1) == -huge(this%n_layers(idx__1))) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("soil_depth")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'soil_depth'"
    case ("fracsealed_cityarea")
      if (.not. allocated(this%fracsealed_cityarea)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: fracSealed_cityArea"
        return
      end if
      allocate(mask( &
        size(this%fracsealed_cityarea, 1)))
      do idx__1 = 1, size(this%fracsealed_cityarea, 1)
        if (ieee_is_nan(this%fracsealed_cityarea(idx__1))) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("land_cover_path")
      if (.not. allocated(this%land_cover_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_path"
        return
      end if
      allocate(mask( &
        size(this%land_cover_path, 1)))
      do idx__1 = 1, size(this%land_cover_path, 1)
        if (this%land_cover_path(idx__1) == achar(0)) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("land_cover_var")
      if (.not. allocated(this%land_cover_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_var"
        return
      end if
      allocate(mask( &
        size(this%land_cover_var, 1)))
      mask = .true.
    case ("lai_time_step")
      if (.not. allocated(this%lai_time_step)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_time_step"
        return
      end if
      allocate(mask( &
        size(this%lai_time_step, 1)))
      do idx__1 = 1, size(this%lai_time_step, 1)
        if (this%lai_time_step(idx__1) == -huge(this%lai_time_step(idx__1))) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("lai_path")
      if (.not. allocated(this%lai_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_path"
        return
      end if
      allocate(mask( &
        size(this%lai_path, 1)))
      do idx__1 = 1, size(this%lai_path, 1)
        if (this%lai_path(idx__1) == achar(0)) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("lai_var")
      if (.not. allocated(this%lai_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_var"
        return
      end if
      allocate(mask( &
        size(this%lai_var, 1)))
      mask = .true.
    case ("read_restart")
      if (.not. allocated(this%read_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: read_restart"
        return
      end if
      allocate(mask( &
        size(this%read_restart, 1)))
      mask = .true.
    case ("restart_input_path")
      if (.not. allocated(this%restart_input_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_input_path"
        return
      end if
      allocate(mask( &
        size(this%restart_input_path, 1)))
      do idx__1 = 1, size(this%restart_input_path, 1)
        if (this%restart_input_path(idx__1) == achar(0)) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case ("write_restart")
      if (.not. allocated(this%write_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: write_restart"
        return
      end if
      allocate(mask( &
        size(this%write_restart, 1)))
      mask = .true.
    case ("restart_output_path")
      if (.not. allocated(this%restart_output_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_output_path"
        return
      end if
      allocate(mask( &
        size(this%restart_output_path, 1)))
      do idx__1 = 1, size(this%restart_output_path, 1)
        if (this%restart_output_path(idx__1) == achar(0)) then
          mask(idx__1) = .false.
        else
          mask(idx__1) = .true.
        end if
      end do
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
  end function nml_config_mpr_set_mask__rank__1

  !> \brief Obtain the item presence mask with one field-name dispatch
  integer function nml_config_mpr_set_mask__rank__2(this, name, mask, errmsg) result(status)
    class(nml_config_mpr_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field or derived component path
    logical, allocatable, intent(out) :: mask(:, :) !< item presence
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: idx__1
    integer :: idx__2

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if
    select case (to_lower(trim(name)))
    case ("soil_db_mode")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'soil_db_mode'"
    case ("tillage_depth")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'tillage_depth'"
    case ("n_layers")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'n_layers'"
    case ("soil_depth")
      if (.not. allocated(this%soil_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: soil_depth"
        return
      end if
      allocate(mask( &
        size(this%soil_depth, 1), &
        size(this%soil_depth, 2)))
      mask = .true.
    case ("fracsealed_cityarea")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'fracSealed_cityArea'"
    case ("land_cover_path")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'land_cover_path'"
    case ("land_cover_var")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'land_cover_var'"
    case ("lai_time_step")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'lai_time_step'"
    case ("lai_path")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'lai_path'"
    case ("lai_var")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'lai_var'"
    case ("read_restart")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'read_restart'"
    case ("restart_input_path")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'restart_input_path'"
    case ("write_restart")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'write_restart'"
    case ("restart_output_path")
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "mask rank mismatch for 'restart_output_path'"
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
  end function nml_config_mpr_set_mask__rank__2

  !> \brief Determine the smallest complete populated prefix starting at index 1
  integer function nml_config_mpr_filled_shape(this, name, filled, errmsg) result(status)
    class(nml_config_mpr_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field or derived component path
    integer, intent(out) :: filled(:) !< prefix extents, retained when the prefix has holes
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: idx__1
    integer :: idx__2

    filled = 0
    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if
    select case (to_lower(trim(name)))
    case ("soil_db_mode")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'soil_db_mode'"
        return
      end if
      if (.not. allocated(this%soil_db_mode)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: soil_db_mode"
        return
      end if
      filled(1) = size(this%soil_db_mode, 1)
    case ("tillage_depth")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'tillage_depth'"
        return
      end if
      if (.not. allocated(this%tillage_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: tillage_depth"
        return
      end if
      do idx__1 = 1, size(this%tillage_depth, 1)
        if (this%tillage_depth(idx__1) == -huge(this%tillage_depth(idx__1))) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%tillage_depth(idx__1) == -huge(this%tillage_depth(idx__1))) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: tillage_depth"
          return
        end if
      end do
    case ("n_layers")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'n_layers'"
        return
      end if
      if (.not. allocated(this%n_layers)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: n_layers"
        return
      end if
      do idx__1 = 1, size(this%n_layers, 1)
        if (this%n_layers(idx__1) == -huge(this%n_layers(idx__1))) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%n_layers(idx__1) == -huge(this%n_layers(idx__1))) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: n_layers"
          return
        end if
      end do
    case ("soil_depth")
      if (size(filled) /= 2) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'soil_depth'"
        return
      end if
      if (.not. allocated(this%soil_depth)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: soil_depth"
        return
      end if
      filled(1) = size(this%soil_depth, 1)
      filled(2) = size(this%soil_depth, 2)
    case ("fracsealed_cityarea")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'fracSealed_cityArea'"
        return
      end if
      if (.not. allocated(this%fracsealed_cityarea)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: fracSealed_cityArea"
        return
      end if
      do idx__1 = 1, size(this%fracsealed_cityarea, 1)
        if (ieee_is_nan(this%fracsealed_cityarea(idx__1))) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (ieee_is_nan(this%fracsealed_cityarea(idx__1))) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: fracSealed_cityArea"
          return
        end if
      end do
    case ("land_cover_path")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'land_cover_path'"
        return
      end if
      if (.not. allocated(this%land_cover_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_path"
        return
      end if
      do idx__1 = 1, size(this%land_cover_path, 1)
        if (this%land_cover_path(idx__1) == achar(0)) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%land_cover_path(idx__1) == achar(0)) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: land_cover_path"
          return
        end if
      end do
    case ("land_cover_var")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'land_cover_var'"
        return
      end if
      if (.not. allocated(this%land_cover_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: land_cover_var"
        return
      end if
      filled(1) = size(this%land_cover_var, 1)
    case ("lai_time_step")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'lai_time_step'"
        return
      end if
      if (.not. allocated(this%lai_time_step)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_time_step"
        return
      end if
      do idx__1 = 1, size(this%lai_time_step, 1)
        if (this%lai_time_step(idx__1) == -huge(this%lai_time_step(idx__1))) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%lai_time_step(idx__1) == -huge(this%lai_time_step(idx__1))) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: lai_time_step"
          return
        end if
      end do
    case ("lai_path")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'lai_path'"
        return
      end if
      if (.not. allocated(this%lai_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_path"
        return
      end if
      do idx__1 = 1, size(this%lai_path, 1)
        if (this%lai_path(idx__1) == achar(0)) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%lai_path(idx__1) == achar(0)) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: lai_path"
          return
        end if
      end do
    case ("lai_var")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'lai_var'"
        return
      end if
      if (.not. allocated(this%lai_var)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: lai_var"
        return
      end if
      filled(1) = size(this%lai_var, 1)
    case ("read_restart")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'read_restart'"
        return
      end if
      if (.not. allocated(this%read_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: read_restart"
        return
      end if
      filled(1) = size(this%read_restart, 1)
    case ("restart_input_path")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'restart_input_path'"
        return
      end if
      if (.not. allocated(this%restart_input_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_input_path"
        return
      end if
      do idx__1 = 1, size(this%restart_input_path, 1)
        if (this%restart_input_path(idx__1) == achar(0)) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%restart_input_path(idx__1) == achar(0)) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: restart_input_path"
          return
        end if
      end do
    case ("write_restart")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'write_restart'"
        return
      end if
      if (.not. allocated(this%write_restart)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: write_restart"
        return
      end if
      filled(1) = size(this%write_restart, 1)
    case ("restart_output_path")
      if (size(filled) /= 1) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "shape rank mismatch for 'restart_output_path'"
        return
      end if
      if (.not. allocated(this%restart_output_path)) then
        status = NML_ERR_NOT_SET
        if (present(errmsg)) errmsg = "array storage unavailable: restart_output_path"
        return
      end if
      do idx__1 = 1, size(this%restart_output_path, 1)
        if (this%restart_output_path(idx__1) == achar(0)) cycle
        if (idx__1 > filled(1)) filled(1) = idx__1
      end do
      do idx__1 = 1, filled(1)
        if (this%restart_output_path(idx__1) == achar(0)) then
          status = NML_ERR_PARTLY_SET
          if (present(errmsg)) errmsg = "array partly set: restart_output_path"
          return
        end if
      end do
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
  end function nml_config_mpr_filled_shape

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
