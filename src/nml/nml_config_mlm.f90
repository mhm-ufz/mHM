!> \file nml_config_mlm.f90
!> \copydoc nml_config_mlm

!> \brief mLM configuration
!> \details Configuration for the v6 lake component (mLM).
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_mlm
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
    buf, &
    max_lakes__default
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4, &
    dp

  implicit none

  ! default values
  logical, parameter, public :: read_restart__default = .false.
  logical, parameter, public :: write_restart__default = .false.
  logical, parameter, public :: restart_save_topology__default = .false.
  logical, parameter, public :: force_initialize_lake_at_restart__default = .false.
  logical, parameter, public :: force_initialize_lake_at_a_date__default = .false.
  integer(i4), parameter, public :: evaporation_method__default = 1_i4
  real(dp), parameter, public :: input_windspeed_height__default = 2.0_dp
  real(dp), parameter, public :: climatological_windspeed__default = 2.0_dp
  logical, parameter, public :: varying_albedo__default = .true.
  logical, parameter, public :: varying_lambda__default = .true.
  logical, parameter, public :: varying_psychro__default = .true.
  integer(i4), parameter, public :: albedo_method__default = 2_i4
  integer(i4), parameter, public :: bathymetry_type__default = 0_i4
  logical, parameter, public :: natural__default = .true.

  ! enum values
  integer(i4), parameter, public :: evaporation_method__enum_values(3) = [1_i4, 2_i4, 3_i4]
  integer(i4), parameter, public :: albedo_method__enum_values(2) = [1_i4, 2_i4]
  integer(i4), parameter, public :: bathymetry_type__enum_values(3) = [0_i4, 1_i4, 2_i4]

  ! bounds values
  integer(i4), parameter, public :: n_lakes__min = 1_i4
  real(dp), parameter, public :: input_windspeed_height__min = 0.0_dp
  real(dp), parameter, public :: climatological_windspeed__min = 0.0_dp

  !> \class nml_config_mlm_t
  !> \brief mLM configuration
  !> \details Configuration for the v6 lake component (mLM).
  type, public :: nml_config_mlm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer :: max_lakes = max_lakes__default !< runtime dimension for max_lakes
    integer(i4), allocatable, dimension(:) :: n_lakes !< Number of lakes
    character(len=buf), allocatable, dimension(:) :: output_path !< mLM lake-point output file
    logical, allocatable, dimension(:) :: read_restart !< Read mLM restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path !< mLM restart input file
    logical, allocatable, dimension(:) :: write_restart !< Write mLM restart
    logical, allocatable, dimension(:) :: restart_save_topology !< Include level-0 lake topology in restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path !< mLM restart output file
    logical :: force_initialize_lake_at_restart !< Force lake level initialization at restart
    logical :: force_initialize_lake_at_a_date !< Force lake level initialization at a date
    integer(i4) :: evaporation_method !< Lake evaporation method
    real(dp) :: input_windspeed_height !< Input windspeed height [m]
    real(dp) :: climatological_windspeed !< Climatological 2 m windspeed [m s-1]
    logical :: varying_albedo !< Variable water albedo
    logical :: varying_lambda !< Variable latent heat of vaporisation
    logical :: varying_psychro !< Variable psychrometric constant
    integer(i4) :: albedo_method !< Water albedo calculation method
    integer(i4), allocatable, dimension(:, :) :: bathymetry_type !< Lake hAV curve type
    logical, allocatable, dimension(:, :) :: natural !< Natural lake flag
  contains
    procedure :: init => nml_config_mlm_init
    procedure :: set_dims => nml_config_mlm_set_dims
    procedure :: from_file => nml_config_mlm_from_file
    procedure :: set => nml_config_mlm_set
    procedure :: is_set => nml_config_mlm_is_set
    procedure :: is_valid => nml_config_mlm_is_valid
  end type nml_config_mlm_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function evaporation_method__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == evaporation_method__enum_values)
  end function evaporation_method__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function albedo_method__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == albedo_method__enum_values)
  end function albedo_method__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function bathymetry_type__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == bathymetry_type__enum_values)
  end function bathymetry_type__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function n_lakes__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < n_lakes__min) in_bounds = .false.
  end function n_lakes__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function input_windspeed_height__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < input_windspeed_height__min) in_bounds = .false.
  end function input_windspeed_height__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function climatological_windspeed__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < climatological_windspeed__min) in_bounds = .false.
  end function climatological_windspeed__in_bounds

  !> \brief Initialize defaults and sentinels for config_mlm
  integer function nml_config_mlm_init(this, errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%n_lakes)) deallocate(this%n_lakes)
    allocate(this%n_lakes(this%n_domains))
    if (allocated(this%output_path)) deallocate(this%output_path)
    allocate(character(len=buf) :: this%output_path(this%n_domains))
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    allocate(this%read_restart(this%n_domains))
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    allocate(character(len=buf) :: this%restart_input_path(this%n_domains))
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    allocate(this%write_restart(this%n_domains))
    if (allocated(this%restart_save_topology)) deallocate(this%restart_save_topology)
    allocate(this%restart_save_topology(this%n_domains))
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    allocate(character(len=buf) :: this%restart_output_path(this%n_domains))
    if (allocated(this%bathymetry_type)) deallocate(this%bathymetry_type)
    allocate(this%bathymetry_type(this%n_domains, this%max_lakes))
    if (allocated(this%natural)) deallocate(this%natural)
    allocate(this%natural(this%n_domains, this%max_lakes))

    ! sentinel values for required/optional parameters
    this%n_lakes = -huge(this%n_lakes) ! sentinel for optional integer array
    this%output_path = achar(0) ! sentinel for optional string array
    this%restart_input_path = achar(0) ! sentinel for optional string array
    this%restart_output_path = achar(0) ! sentinel for optional string array
    ! default values
    this%read_restart = read_restart__default
    this%write_restart = write_restart__default
    this%restart_save_topology = restart_save_topology__default
    this%force_initialize_lake_at_restart = force_initialize_lake_at_restart__default ! bool values always need a default
    this%force_initialize_lake_at_a_date = force_initialize_lake_at_a_date__default ! bool values always need a default
    this%evaporation_method = evaporation_method__default
    this%input_windspeed_height = input_windspeed_height__default
    this%climatological_windspeed = climatological_windspeed__default
    this%varying_albedo = varying_albedo__default ! bool values always need a default
    this%varying_lambda = varying_lambda__default ! bool values always need a default
    this%varying_psychro = varying_psychro__default ! bool values always need a default
    this%albedo_method = albedo_method__default
    this%bathymetry_type = bathymetry_type__default
    this%natural = natural__default
  end function nml_config_mlm_init

  !> \brief Reset runtime dimensions for config_mlm
  integer function nml_config_mlm_set_dims(this, &
    n_domains, &
    max_lakes, &
    errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer, intent(in), optional :: max_lakes !< runtime dimension override for max_lakes
    integer :: candidate__n_domains
    integer :: candidate__max_lakes
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
    if (present(max_lakes)) then
      candidate__max_lakes = max_lakes
    else
      candidate__max_lakes = max_lakes__default
    end if
    if (candidate__max_lakes <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'max_lakes' must be positive"
      return
    end if
    this%n_domains = candidate__n_domains
    this%max_lakes = candidate__max_lakes

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%n_lakes)) deallocate(this%n_lakes)
    if (allocated(this%output_path)) deallocate(this%output_path)
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    if (allocated(this%restart_save_topology)) deallocate(this%restart_save_topology)
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    if (allocated(this%bathymetry_type)) deallocate(this%bathymetry_type)
    if (allocated(this%natural)) deallocate(this%natural)
    this%is_configured = .false.
  end function nml_config_mlm_set_dims


  !> \brief Read config_mlm namelist from file
  integer function nml_config_mlm_from_file(this, file, errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4), allocatable, dimension(:) :: n_lakes
    character(len=buf), allocatable, dimension(:) :: output_path
    logical, allocatable, dimension(:) :: read_restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path
    logical, allocatable, dimension(:) :: write_restart
    logical, allocatable, dimension(:) :: restart_save_topology
    character(len=buf), allocatable, dimension(:) :: restart_output_path
    logical :: force_initialize_lake_at_restart
    logical :: force_initialize_lake_at_a_date
    integer(i4) :: evaporation_method
    real(dp) :: input_windspeed_height
    real(dp) :: climatological_windspeed
    logical :: varying_albedo
    logical :: varying_lambda
    logical :: varying_psychro
    integer(i4) :: albedo_method
    integer(i4), allocatable, dimension(:, :) :: bathymetry_type
    logical, allocatable, dimension(:, :) :: natural
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_mlm/ &
      n_lakes, &
      output_path, &
      read_restart, &
      restart_input_path, &
      write_restart, &
      restart_save_topology, &
      restart_output_path, &
      force_initialize_lake_at_restart, &
      force_initialize_lake_at_a_date, &
      evaporation_method, &
      input_windspeed_height, &
      climatological_windspeed, &
      varying_albedo, &
      varying_lambda, &
      varying_psychro, &
      albedo_method, &
      bathymetry_type, &
      natural

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(n_lakes)) deallocate(n_lakes)
    allocate(n_lakes(this%n_domains))
    if (allocated(output_path)) deallocate(output_path)
    allocate(character(len=buf) :: output_path(this%n_domains))
    if (allocated(read_restart)) deallocate(read_restart)
    allocate(read_restart(this%n_domains))
    if (allocated(restart_input_path)) deallocate(restart_input_path)
    allocate(character(len=buf) :: restart_input_path(this%n_domains))
    if (allocated(write_restart)) deallocate(write_restart)
    allocate(write_restart(this%n_domains))
    if (allocated(restart_save_topology)) deallocate(restart_save_topology)
    allocate(restart_save_topology(this%n_domains))
    if (allocated(restart_output_path)) deallocate(restart_output_path)
    allocate(character(len=buf) :: restart_output_path(this%n_domains))
    if (allocated(bathymetry_type)) deallocate(bathymetry_type)
    allocate(bathymetry_type(this%n_domains, this%max_lakes))
    if (allocated(natural)) deallocate(natural)
    allocate(natural(this%n_domains, this%max_lakes))
    n_lakes = this%n_lakes
    output_path = this%output_path
    read_restart = this%read_restart
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_save_topology = this%restart_save_topology
    restart_output_path = this%restart_output_path
    force_initialize_lake_at_restart = this%force_initialize_lake_at_restart
    force_initialize_lake_at_a_date = this%force_initialize_lake_at_a_date
    evaporation_method = this%evaporation_method
    input_windspeed_height = this%input_windspeed_height
    climatological_windspeed = this%climatological_windspeed
    varying_albedo = this%varying_albedo
    varying_lambda = this%varying_lambda
    varying_psychro = this%varying_psychro
    albedo_method = this%albedo_method
    bathymetry_type = this%bathymetry_type
    natural = this%natural

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_mlm", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_mlm, iostat=iostat, iomsg=iomsg)
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
    this%n_lakes = n_lakes
    this%output_path = output_path
    this%read_restart = read_restart
    this%restart_input_path = restart_input_path
    this%write_restart = write_restart
    this%restart_save_topology = restart_save_topology
    this%restart_output_path = restart_output_path
    this%force_initialize_lake_at_restart = force_initialize_lake_at_restart
    this%force_initialize_lake_at_a_date = force_initialize_lake_at_a_date
    this%evaporation_method = evaporation_method
    this%input_windspeed_height = input_windspeed_height
    this%climatological_windspeed = climatological_windspeed
    this%varying_albedo = varying_albedo
    this%varying_lambda = varying_lambda
    this%varying_psychro = varying_psychro
    this%albedo_method = albedo_method
    this%bathymetry_type = bathymetry_type
    this%natural = natural

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mlm_from_file

  !> \brief Set config_mlm values
  integer function nml_config_mlm_set(this, &
    n_lakes, &
    output_path, &
    read_restart, &
    restart_input_path, &
    write_restart, &
    restart_save_topology, &
    restart_output_path, &
    force_initialize_lake_at_restart, &
    force_initialize_lake_at_a_date, &
    evaporation_method, &
    input_windspeed_height, &
    climatological_windspeed, &
    varying_albedo, &
    varying_lambda, &
    varying_psychro, &
    albedo_method, &
    bathymetry_type, &
    natural, &
    errmsg) result(status)

    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), dimension(:), intent(in), optional :: n_lakes !< Number of lakes
    character(len=*), dimension(:), intent(in), optional :: output_path !< mLM lake-point output file
    logical, dimension(:), intent(in), optional :: read_restart !< Read mLM restart
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< mLM restart input file
    logical, dimension(:), intent(in), optional :: write_restart !< Write mLM restart
    logical, dimension(:), intent(in), optional :: restart_save_topology !< Include level-0 lake topology in restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< mLM restart output file
    logical, intent(in), optional :: force_initialize_lake_at_restart !< Force lake level initialization at restart
    logical, intent(in), optional :: force_initialize_lake_at_a_date !< Force lake level initialization at a date
    integer(i4), intent(in), optional :: evaporation_method !< Lake evaporation method
    real(dp), intent(in), optional :: input_windspeed_height !< Input windspeed height [m]
    real(dp), intent(in), optional :: climatological_windspeed !< Climatological 2 m windspeed [m s-1]
    logical, intent(in), optional :: varying_albedo !< Variable water albedo
    logical, intent(in), optional :: varying_lambda !< Variable latent heat of vaporisation
    logical, intent(in), optional :: varying_psychro !< Variable psychrometric constant
    integer(i4), intent(in), optional :: albedo_method !< Water albedo calculation method
    integer(i4), dimension(:, :), intent(in), optional :: bathymetry_type !< Lake hAV curve type
    logical, dimension(:, :), intent(in), optional :: natural !< Natural lake flag
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(n_lakes)) then
      if (size(n_lakes, 1) > size(this%n_lakes, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'n_lakes'"
        return
      end if
      lb__1 = lbound(this%n_lakes, 1)
      ub__1 = lb__1 + size(n_lakes, 1) - 1
      this%n_lakes(lb__1:ub__1) = n_lakes
    end if
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
    if (present(restart_save_topology)) then
      if (size(restart_save_topology, 1) > size(this%restart_save_topology, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_save_topology'"
        return
      end if
      lb__1 = lbound(this%restart_save_topology, 1)
      ub__1 = lb__1 + size(restart_save_topology, 1) - 1
      this%restart_save_topology(lb__1:ub__1) = restart_save_topology
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
    if (present(force_initialize_lake_at_restart)) this%force_initialize_lake_at_restart = force_initialize_lake_at_restart
    if (present(force_initialize_lake_at_a_date)) this%force_initialize_lake_at_a_date = force_initialize_lake_at_a_date
    if (present(evaporation_method)) this%evaporation_method = evaporation_method
    if (present(input_windspeed_height)) this%input_windspeed_height = input_windspeed_height
    if (present(climatological_windspeed)) this%climatological_windspeed = climatological_windspeed
    if (present(varying_albedo)) this%varying_albedo = varying_albedo
    if (present(varying_lambda)) this%varying_lambda = varying_lambda
    if (present(varying_psychro)) this%varying_psychro = varying_psychro
    if (present(albedo_method)) this%albedo_method = albedo_method
    if (present(bathymetry_type)) then
      if (size(bathymetry_type, 1) > size(this%bathymetry_type, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'bathymetry_type'"
        return
      end if
      lb__1 = lbound(this%bathymetry_type, 1)
      ub__1 = lb__1 + size(bathymetry_type, 1) - 1
      if (size(bathymetry_type, 2) > size(this%bathymetry_type, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'bathymetry_type'"
        return
      end if
      lb__2 = lbound(this%bathymetry_type, 2)
      ub__2 = lb__2 + size(bathymetry_type, 2) - 1
      this%bathymetry_type(lb__1:ub__1, lb__2:ub__2) = bathymetry_type
    end if
    if (present(natural)) then
      if (size(natural, 1) > size(this%natural, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'natural'"
        return
      end if
      lb__1 = lbound(this%natural, 1)
      ub__1 = lb__1 + size(natural, 1) - 1
      if (size(natural, 2) > size(this%natural, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'natural'"
        return
      end if
      lb__2 = lbound(this%natural, 2)
      ub__2 = lb__2 + size(natural, 2) - 1
      this%natural(lb__1:ub__1, lb__2:ub__2) = natural
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mlm_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_mlm_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_mlm_t), intent(in) :: this !< namelist instance
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
    case ("n_lakes")
      if (.not. allocated(this%n_lakes)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%n_lakes), ubound(this%n_lakes), &
          "n_lakes", errmsg)
        if (status /= NML_OK) return
        if (this%n_lakes(idx(1)) == -huge(this%n_lakes(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%n_lakes == -huge(this%n_lakes))) status = NML_ERR_NOT_SET
      end if
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
    case ("restart_save_topology")
      if (.not. allocated(this%restart_save_topology)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%restart_save_topology), ubound(this%restart_save_topology), &
          "restart_save_topology", errmsg)
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
    case ("force_initialize_lake_at_restart")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'force_initialize_lake_at_restart'"
        return
      end if
    case ("force_initialize_lake_at_a_date")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'force_initialize_lake_at_a_date'"
        return
      end if
    case ("evaporation_method")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'evaporation_method'"
        return
      end if
    case ("input_windspeed_height")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'input_windspeed_height'"
        return
      end if
    case ("climatological_windspeed")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'climatological_windspeed'"
        return
      end if
    case ("varying_albedo")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'varying_albedo'"
        return
      end if
    case ("varying_lambda")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'varying_lambda'"
        return
      end if
    case ("varying_psychro")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'varying_psychro'"
        return
      end if
    case ("albedo_method")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'albedo_method'"
        return
      end if
    case ("bathymetry_type")
      if (.not. allocated(this%bathymetry_type)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%bathymetry_type), ubound(this%bathymetry_type), &
          "bathymetry_type", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("natural")
      if (.not. allocated(this%natural)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%natural), ubound(this%natural), &
          "natural", errmsg)
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
  end function nml_config_mlm_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_mlm_is_valid(this, errmsg) result(status)
    class(nml_config_mlm_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("evaporation_method", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. evaporation_method__in_enum(this%evaporation_method)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: evaporation_method"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("albedo_method", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. albedo_method__in_enum(this%albedo_method)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: albedo_method"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    if (allocated(this%bathymetry_type)) then
    if (.not. all(bathymetry_type__in_enum(this%bathymetry_type, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: bathymetry_type"
      return
    end if
    end if
    ! bounds constraints
    if (allocated(this%n_lakes)) then
    if (.not. all(n_lakes__in_bounds(this%n_lakes, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: n_lakes"
      return
    end if
    end if
    istat = this%is_set("input_windspeed_height", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. input_windspeed_height__in_bounds(this%input_windspeed_height)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: input_windspeed_height"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("climatological_windspeed", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. climatological_windspeed__in_bounds(this%climatological_windspeed)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: climatological_windspeed"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
  end function nml_config_mlm_is_valid

end module nml_config_mlm
