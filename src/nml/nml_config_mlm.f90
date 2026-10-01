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
    n_lakes__default
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp, &
    i4

  implicit none

  ! default values
  logical, parameter, public :: read_restart__default = .false.
  logical, parameter, public :: write_restart__default = .false.
  logical, parameter, public :: force_initialize_lake_at_restart__default = .false.
  logical, parameter, public :: force_initialize_lake_at_a_date__default = .false.
  real(dp), parameter, public :: input_windspeed_height__default = 2.0_dp
  logical, parameter, public :: varying_albedo__default = .true.
  logical, parameter, public :: varying_lambda__default = .true.
  logical, parameter, public :: varying_psychro__default = .true.
  integer(i4), parameter, public :: albedo_method__default = 2_i4
  logical, parameter, public :: hanasaki_2006__default = .true.
  logical, parameter, public :: point_demand__default = .false.
  integer(i4), parameter, public :: bathymetry_type__default = 0_i4
  logical, parameter, public :: natural__default = .true.
  logical, parameter, public :: use_watersupply__default = .false.
  logical, parameter, public :: use_irrigation__default = .false.
  real(dp), parameter, public :: service_distance__default = 600.0_dp
  logical, parameter, public :: print_lake_maps__default = .true.
  logical, parameter, public :: print_links_table__default = .true.

  ! enum values
  integer(i4), parameter, public :: albedo_method__enum_values(2) = [1_i4, 2_i4]
  integer(i4), parameter, public :: bathymetry_type__enum_values(3) = [0_i4, 1_i4, 2_i4]

  ! bounds values
  real(dp), parameter, public :: input_windspeed_height__min = 0.0_dp

  !> \class nml_config_mlm_t
  !> \brief mLM configuration
  !> \details Configuration for the v6 lake component (mLM).
  type, public :: nml_config_mlm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer :: n_lakes = n_lakes__default !< runtime dimension for n_lakes
    character(len=buf), allocatable, dimension(:) :: input_path !< mLM lake input path
    character(len=buf), allocatable, dimension(:) :: output_path !< mLM lake-point output path
    logical, allocatable, dimension(:) :: read_restart !< Read mLM restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path !< mLM restart input path
    logical, allocatable, dimension(:) :: write_restart !< Write mLM restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path !< mLM restart output path
    logical :: force_initialize_lake_at_restart !< Force lake level initialization at restart
    logical :: force_initialize_lake_at_a_date !< Force lake level initialization at a date
    real(dp) :: input_windspeed_height !< Input windspeed height [m]
    logical :: varying_albedo !< Variable water albedo
    logical :: varying_lambda !< Variable latent heat of vaporisation
    logical :: varying_psychro !< Variable psychrometric constant
    integer(i4) :: albedo_method !< Water albedo calculation method
    logical :: hanasaki_2006 !< Hanasaki et al. (2006) release target
    logical :: point_demand !< Point demand input
    integer(i4), allocatable, dimension(:, :) :: bathymetry_type !< Lake hAV curve type
    logical, allocatable, dimension(:, :) :: natural !< Natural lake flag
    logical, allocatable, dimension(:, :) :: use_watersupply !< Water-supply use flag
    logical, allocatable, dimension(:, :) :: use_irrigation !< Irrigation use flag
    real(dp), allocatable, dimension(:, :) :: service_distance !< Downstream service distance [km]
    character(len=buf), allocatable, dimension(:, :) :: bathymetry_filename !< Lake bathymetry file
    character(len=buf), allocatable, dimension(:, :) :: release_target_filename !< Release target file
    character(len=buf), allocatable, dimension(:, :) :: watersupply_demand_filename !< Water-supply demand file
    character(len=buf), allocatable, dimension(:, :) :: irrigation_demand_filename !< Irrigation demand file
    character(len=buf), allocatable, dimension(:, :) :: ibt_filename !< Inter-basin transfer file
    character(len=buf), allocatable, dimension(:, :) :: level_filename !< Observed water-level file
    character(len=buf), allocatable, dimension(:, :) :: area_filename !< Observed surface-area file
    character(len=buf), allocatable, dimension(:, :) :: spill_filename !< Observed spill file
    logical :: print_lake_maps !< Print lake maps
    logical :: print_links_table !< Print links table
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

  !> \brief Initialize defaults and sentinels for config_mlm
  integer function nml_config_mlm_init(this, errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%input_path)) deallocate(this%input_path)
    allocate(character(len=buf) :: this%input_path(this%n_domains))
    if (allocated(this%output_path)) deallocate(this%output_path)
    allocate(character(len=buf) :: this%output_path(this%n_domains))
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    allocate(this%read_restart(this%n_domains))
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    allocate(character(len=buf) :: this%restart_input_path(this%n_domains))
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    allocate(this%write_restart(this%n_domains))
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    allocate(character(len=buf) :: this%restart_output_path(this%n_domains))
    if (allocated(this%bathymetry_type)) deallocate(this%bathymetry_type)
    allocate(this%bathymetry_type(this%n_domains, this%n_lakes))
    if (allocated(this%natural)) deallocate(this%natural)
    allocate(this%natural(this%n_domains, this%n_lakes))
    if (allocated(this%use_watersupply)) deallocate(this%use_watersupply)
    allocate(this%use_watersupply(this%n_domains, this%n_lakes))
    if (allocated(this%use_irrigation)) deallocate(this%use_irrigation)
    allocate(this%use_irrigation(this%n_domains, this%n_lakes))
    if (allocated(this%service_distance)) deallocate(this%service_distance)
    allocate(this%service_distance(this%n_domains, this%n_lakes))
    if (allocated(this%bathymetry_filename)) deallocate(this%bathymetry_filename)
    allocate(character(len=buf) :: this%bathymetry_filename(this%n_domains, this%n_lakes))
    if (allocated(this%release_target_filename)) deallocate(this%release_target_filename)
    allocate(character(len=buf) :: this%release_target_filename(this%n_domains, this%n_lakes))
    if (allocated(this%watersupply_demand_filename)) deallocate(this%watersupply_demand_filename)
    allocate(character(len=buf) :: this%watersupply_demand_filename(this%n_domains, this%n_lakes))
    if (allocated(this%irrigation_demand_filename)) deallocate(this%irrigation_demand_filename)
    allocate(character(len=buf) :: this%irrigation_demand_filename(this%n_domains, this%n_lakes))
    if (allocated(this%ibt_filename)) deallocate(this%ibt_filename)
    allocate(character(len=buf) :: this%ibt_filename(this%n_domains, this%n_lakes))
    if (allocated(this%level_filename)) deallocate(this%level_filename)
    allocate(character(len=buf) :: this%level_filename(this%n_domains, this%n_lakes))
    if (allocated(this%area_filename)) deallocate(this%area_filename)
    allocate(character(len=buf) :: this%area_filename(this%n_domains, this%n_lakes))
    if (allocated(this%spill_filename)) deallocate(this%spill_filename)
    allocate(character(len=buf) :: this%spill_filename(this%n_domains, this%n_lakes))

    ! sentinel values for required/optional parameters
    this%input_path = achar(0) ! sentinel for optional string array
    this%output_path = achar(0) ! sentinel for optional string array
    this%restart_input_path = achar(0) ! sentinel for optional string array
    this%restart_output_path = achar(0) ! sentinel for optional string array
    this%bathymetry_filename = achar(0) ! sentinel for optional string array
    this%release_target_filename = achar(0) ! sentinel for optional string array
    this%watersupply_demand_filename = achar(0) ! sentinel for optional string array
    this%irrigation_demand_filename = achar(0) ! sentinel for optional string array
    this%ibt_filename = achar(0) ! sentinel for optional string array
    this%level_filename = achar(0) ! sentinel for optional string array
    this%area_filename = achar(0) ! sentinel for optional string array
    this%spill_filename = achar(0) ! sentinel for optional string array
    ! default values
    this%read_restart = read_restart__default
    this%write_restart = write_restart__default
    this%force_initialize_lake_at_restart = force_initialize_lake_at_restart__default ! bool values always need a default
    this%force_initialize_lake_at_a_date = force_initialize_lake_at_a_date__default ! bool values always need a default
    this%input_windspeed_height = input_windspeed_height__default
    this%varying_albedo = varying_albedo__default ! bool values always need a default
    this%varying_lambda = varying_lambda__default ! bool values always need a default
    this%varying_psychro = varying_psychro__default ! bool values always need a default
    this%albedo_method = albedo_method__default
    this%hanasaki_2006 = hanasaki_2006__default ! bool values always need a default
    this%point_demand = point_demand__default ! bool values always need a default
    this%bathymetry_type = bathymetry_type__default
    this%natural = natural__default
    this%use_watersupply = use_watersupply__default
    this%use_irrigation = use_irrigation__default
    this%service_distance = service_distance__default
    this%print_lake_maps = print_lake_maps__default ! bool values always need a default
    this%print_links_table = print_links_table__default ! bool values always need a default
  end function nml_config_mlm_init

  !> \brief Reset runtime dimensions for config_mlm
  integer function nml_config_mlm_set_dims(this, &
    n_domains, &
    n_lakes, &
    errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer, intent(in), optional :: n_lakes !< runtime dimension override for n_lakes
    integer :: candidate__n_domains
    integer :: candidate__n_lakes
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
    if (present(n_lakes)) then
      candidate__n_lakes = n_lakes
    else
      candidate__n_lakes = n_lakes__default
    end if
    if (candidate__n_lakes <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_lakes' must be positive"
      return
    end if
    this%n_domains = candidate__n_domains
    this%n_lakes = candidate__n_lakes

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%input_path)) deallocate(this%input_path)
    if (allocated(this%output_path)) deallocate(this%output_path)
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    if (allocated(this%bathymetry_type)) deallocate(this%bathymetry_type)
    if (allocated(this%natural)) deallocate(this%natural)
    if (allocated(this%use_watersupply)) deallocate(this%use_watersupply)
    if (allocated(this%use_irrigation)) deallocate(this%use_irrigation)
    if (allocated(this%service_distance)) deallocate(this%service_distance)
    if (allocated(this%bathymetry_filename)) deallocate(this%bathymetry_filename)
    if (allocated(this%release_target_filename)) deallocate(this%release_target_filename)
    if (allocated(this%watersupply_demand_filename)) deallocate(this%watersupply_demand_filename)
    if (allocated(this%irrigation_demand_filename)) deallocate(this%irrigation_demand_filename)
    if (allocated(this%ibt_filename)) deallocate(this%ibt_filename)
    if (allocated(this%level_filename)) deallocate(this%level_filename)
    if (allocated(this%area_filename)) deallocate(this%area_filename)
    if (allocated(this%spill_filename)) deallocate(this%spill_filename)
    this%is_configured = .false.
  end function nml_config_mlm_set_dims


  !> \brief Read config_mlm namelist from file
  integer function nml_config_mlm_from_file(this, file, errmsg) result(status)
    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    character(len=buf), allocatable, dimension(:) :: input_path
    character(len=buf), allocatable, dimension(:) :: output_path
    logical, allocatable, dimension(:) :: read_restart
    character(len=buf), allocatable, dimension(:) :: restart_input_path
    logical, allocatable, dimension(:) :: write_restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path
    logical :: force_initialize_lake_at_restart
    logical :: force_initialize_lake_at_a_date
    real(dp) :: input_windspeed_height
    logical :: varying_albedo
    logical :: varying_lambda
    logical :: varying_psychro
    integer(i4) :: albedo_method
    logical :: hanasaki_2006
    logical :: point_demand
    integer(i4), allocatable, dimension(:, :) :: bathymetry_type
    logical, allocatable, dimension(:, :) :: natural
    logical, allocatable, dimension(:, :) :: use_watersupply
    logical, allocatable, dimension(:, :) :: use_irrigation
    real(dp), allocatable, dimension(:, :) :: service_distance
    character(len=buf), allocatable, dimension(:, :) :: bathymetry_filename
    character(len=buf), allocatable, dimension(:, :) :: release_target_filename
    character(len=buf), allocatable, dimension(:, :) :: watersupply_demand_filename
    character(len=buf), allocatable, dimension(:, :) :: irrigation_demand_filename
    character(len=buf), allocatable, dimension(:, :) :: ibt_filename
    character(len=buf), allocatable, dimension(:, :) :: level_filename
    character(len=buf), allocatable, dimension(:, :) :: area_filename
    character(len=buf), allocatable, dimension(:, :) :: spill_filename
    logical :: print_lake_maps
    logical :: print_links_table
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_mlm/ &
      input_path, &
      output_path, &
      read_restart, &
      restart_input_path, &
      write_restart, &
      restart_output_path, &
      force_initialize_lake_at_restart, &
      force_initialize_lake_at_a_date, &
      input_windspeed_height, &
      varying_albedo, &
      varying_lambda, &
      varying_psychro, &
      albedo_method, &
      hanasaki_2006, &
      point_demand, &
      bathymetry_type, &
      natural, &
      use_watersupply, &
      use_irrigation, &
      service_distance, &
      bathymetry_filename, &
      release_target_filename, &
      watersupply_demand_filename, &
      irrigation_demand_filename, &
      ibt_filename, &
      level_filename, &
      area_filename, &
      spill_filename, &
      print_lake_maps, &
      print_links_table

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(input_path)) deallocate(input_path)
    allocate(character(len=buf) :: input_path(this%n_domains))
    if (allocated(output_path)) deallocate(output_path)
    allocate(character(len=buf) :: output_path(this%n_domains))
    if (allocated(read_restart)) deallocate(read_restart)
    allocate(read_restart(this%n_domains))
    if (allocated(restart_input_path)) deallocate(restart_input_path)
    allocate(character(len=buf) :: restart_input_path(this%n_domains))
    if (allocated(write_restart)) deallocate(write_restart)
    allocate(write_restart(this%n_domains))
    if (allocated(restart_output_path)) deallocate(restart_output_path)
    allocate(character(len=buf) :: restart_output_path(this%n_domains))
    if (allocated(bathymetry_type)) deallocate(bathymetry_type)
    allocate(bathymetry_type(this%n_domains, this%n_lakes))
    if (allocated(natural)) deallocate(natural)
    allocate(natural(this%n_domains, this%n_lakes))
    if (allocated(use_watersupply)) deallocate(use_watersupply)
    allocate(use_watersupply(this%n_domains, this%n_lakes))
    if (allocated(use_irrigation)) deallocate(use_irrigation)
    allocate(use_irrigation(this%n_domains, this%n_lakes))
    if (allocated(service_distance)) deallocate(service_distance)
    allocate(service_distance(this%n_domains, this%n_lakes))
    if (allocated(bathymetry_filename)) deallocate(bathymetry_filename)
    allocate(character(len=buf) :: bathymetry_filename(this%n_domains, this%n_lakes))
    if (allocated(release_target_filename)) deallocate(release_target_filename)
    allocate(character(len=buf) :: release_target_filename(this%n_domains, this%n_lakes))
    if (allocated(watersupply_demand_filename)) deallocate(watersupply_demand_filename)
    allocate(character(len=buf) :: watersupply_demand_filename(this%n_domains, this%n_lakes))
    if (allocated(irrigation_demand_filename)) deallocate(irrigation_demand_filename)
    allocate(character(len=buf) :: irrigation_demand_filename(this%n_domains, this%n_lakes))
    if (allocated(ibt_filename)) deallocate(ibt_filename)
    allocate(character(len=buf) :: ibt_filename(this%n_domains, this%n_lakes))
    if (allocated(level_filename)) deallocate(level_filename)
    allocate(character(len=buf) :: level_filename(this%n_domains, this%n_lakes))
    if (allocated(area_filename)) deallocate(area_filename)
    allocate(character(len=buf) :: area_filename(this%n_domains, this%n_lakes))
    if (allocated(spill_filename)) deallocate(spill_filename)
    allocate(character(len=buf) :: spill_filename(this%n_domains, this%n_lakes))
    input_path = this%input_path
    output_path = this%output_path
    read_restart = this%read_restart
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_output_path = this%restart_output_path
    force_initialize_lake_at_restart = this%force_initialize_lake_at_restart
    force_initialize_lake_at_a_date = this%force_initialize_lake_at_a_date
    input_windspeed_height = this%input_windspeed_height
    varying_albedo = this%varying_albedo
    varying_lambda = this%varying_lambda
    varying_psychro = this%varying_psychro
    albedo_method = this%albedo_method
    hanasaki_2006 = this%hanasaki_2006
    point_demand = this%point_demand
    bathymetry_type = this%bathymetry_type
    natural = this%natural
    use_watersupply = this%use_watersupply
    use_irrigation = this%use_irrigation
    service_distance = this%service_distance
    bathymetry_filename = this%bathymetry_filename
    release_target_filename = this%release_target_filename
    watersupply_demand_filename = this%watersupply_demand_filename
    irrigation_demand_filename = this%irrigation_demand_filename
    ibt_filename = this%ibt_filename
    level_filename = this%level_filename
    area_filename = this%area_filename
    spill_filename = this%spill_filename
    print_lake_maps = this%print_lake_maps
    print_links_table = this%print_links_table

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
    this%input_path = input_path
    this%output_path = output_path
    this%read_restart = read_restart
    this%restart_input_path = restart_input_path
    this%write_restart = write_restart
    this%restart_output_path = restart_output_path
    this%force_initialize_lake_at_restart = force_initialize_lake_at_restart
    this%force_initialize_lake_at_a_date = force_initialize_lake_at_a_date
    this%input_windspeed_height = input_windspeed_height
    this%varying_albedo = varying_albedo
    this%varying_lambda = varying_lambda
    this%varying_psychro = varying_psychro
    this%albedo_method = albedo_method
    this%hanasaki_2006 = hanasaki_2006
    this%point_demand = point_demand
    this%bathymetry_type = bathymetry_type
    this%natural = natural
    this%use_watersupply = use_watersupply
    this%use_irrigation = use_irrigation
    this%service_distance = service_distance
    this%bathymetry_filename = bathymetry_filename
    this%release_target_filename = release_target_filename
    this%watersupply_demand_filename = watersupply_demand_filename
    this%irrigation_demand_filename = irrigation_demand_filename
    this%ibt_filename = ibt_filename
    this%level_filename = level_filename
    this%area_filename = area_filename
    this%spill_filename = spill_filename
    this%print_lake_maps = print_lake_maps
    this%print_links_table = print_links_table

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mlm_from_file

  !> \brief Set config_mlm values
  integer function nml_config_mlm_set(this, &
    input_path, &
    output_path, &
    read_restart, &
    restart_input_path, &
    write_restart, &
    restart_output_path, &
    force_initialize_lake_at_restart, &
    force_initialize_lake_at_a_date, &
    input_windspeed_height, &
    varying_albedo, &
    varying_lambda, &
    varying_psychro, &
    albedo_method, &
    hanasaki_2006, &
    point_demand, &
    bathymetry_type, &
    natural, &
    use_watersupply, &
    use_irrigation, &
    service_distance, &
    bathymetry_filename, &
    release_target_filename, &
    watersupply_demand_filename, &
    irrigation_demand_filename, &
    ibt_filename, &
    level_filename, &
    area_filename, &
    spill_filename, &
    print_lake_maps, &
    print_links_table, &
    errmsg) result(status)

    class(nml_config_mlm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    character(len=*), dimension(:), intent(in), optional :: input_path !< mLM lake input path
    character(len=*), dimension(:), intent(in), optional :: output_path !< mLM lake-point output path
    logical, dimension(:), intent(in), optional :: read_restart !< Read mLM restart
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< mLM restart input path
    logical, dimension(:), intent(in), optional :: write_restart !< Write mLM restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< mLM restart output path
    logical, intent(in), optional :: force_initialize_lake_at_restart !< Force lake level initialization at restart
    logical, intent(in), optional :: force_initialize_lake_at_a_date !< Force lake level initialization at a date
    real(dp), intent(in), optional :: input_windspeed_height !< Input windspeed height [m]
    logical, intent(in), optional :: varying_albedo !< Variable water albedo
    logical, intent(in), optional :: varying_lambda !< Variable latent heat of vaporisation
    logical, intent(in), optional :: varying_psychro !< Variable psychrometric constant
    integer(i4), intent(in), optional :: albedo_method !< Water albedo calculation method
    logical, intent(in), optional :: hanasaki_2006 !< Hanasaki et al. (2006) release target
    logical, intent(in), optional :: point_demand !< Point demand input
    integer(i4), dimension(:, :), intent(in), optional :: bathymetry_type !< Lake hAV curve type
    logical, dimension(:, :), intent(in), optional :: natural !< Natural lake flag
    logical, dimension(:, :), intent(in), optional :: use_watersupply !< Water-supply use flag
    logical, dimension(:, :), intent(in), optional :: use_irrigation !< Irrigation use flag
    real(dp), dimension(:, :), intent(in), optional :: service_distance !< Downstream service distance [km]
    character(len=*), dimension(:, :), intent(in), optional :: bathymetry_filename !< Lake bathymetry file
    character(len=*), dimension(:, :), intent(in), optional :: release_target_filename !< Release target file
    character(len=*), dimension(:, :), intent(in), optional :: watersupply_demand_filename !< Water-supply demand file
    character(len=*), dimension(:, :), intent(in), optional :: irrigation_demand_filename !< Irrigation demand file
    character(len=*), dimension(:, :), intent(in), optional :: ibt_filename !< Inter-basin transfer file
    character(len=*), dimension(:, :), intent(in), optional :: level_filename !< Observed water-level file
    character(len=*), dimension(:, :), intent(in), optional :: area_filename !< Observed surface-area file
    character(len=*), dimension(:, :), intent(in), optional :: spill_filename !< Observed spill file
    logical, intent(in), optional :: print_lake_maps !< Print lake maps
    logical, intent(in), optional :: print_links_table !< Print links table
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(input_path)) then
      if (size(input_path, 1) > size(this%input_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'input_path'"
        return
      end if
      lb__1 = lbound(this%input_path, 1)
      ub__1 = lb__1 + size(input_path, 1) - 1
      this%input_path(lb__1:ub__1) = input_path
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
    if (present(input_windspeed_height)) this%input_windspeed_height = input_windspeed_height
    if (present(varying_albedo)) this%varying_albedo = varying_albedo
    if (present(varying_lambda)) this%varying_lambda = varying_lambda
    if (present(varying_psychro)) this%varying_psychro = varying_psychro
    if (present(albedo_method)) this%albedo_method = albedo_method
    if (present(hanasaki_2006)) this%hanasaki_2006 = hanasaki_2006
    if (present(point_demand)) this%point_demand = point_demand
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
    if (present(use_watersupply)) then
      if (size(use_watersupply, 1) > size(this%use_watersupply, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'use_watersupply'"
        return
      end if
      lb__1 = lbound(this%use_watersupply, 1)
      ub__1 = lb__1 + size(use_watersupply, 1) - 1
      if (size(use_watersupply, 2) > size(this%use_watersupply, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'use_watersupply'"
        return
      end if
      lb__2 = lbound(this%use_watersupply, 2)
      ub__2 = lb__2 + size(use_watersupply, 2) - 1
      this%use_watersupply(lb__1:ub__1, lb__2:ub__2) = use_watersupply
    end if
    if (present(use_irrigation)) then
      if (size(use_irrigation, 1) > size(this%use_irrigation, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'use_irrigation'"
        return
      end if
      lb__1 = lbound(this%use_irrigation, 1)
      ub__1 = lb__1 + size(use_irrigation, 1) - 1
      if (size(use_irrigation, 2) > size(this%use_irrigation, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'use_irrigation'"
        return
      end if
      lb__2 = lbound(this%use_irrigation, 2)
      ub__2 = lb__2 + size(use_irrigation, 2) - 1
      this%use_irrigation(lb__1:ub__1, lb__2:ub__2) = use_irrigation
    end if
    if (present(service_distance)) then
      if (size(service_distance, 1) > size(this%service_distance, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'service_distance'"
        return
      end if
      lb__1 = lbound(this%service_distance, 1)
      ub__1 = lb__1 + size(service_distance, 1) - 1
      if (size(service_distance, 2) > size(this%service_distance, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'service_distance'"
        return
      end if
      lb__2 = lbound(this%service_distance, 2)
      ub__2 = lb__2 + size(service_distance, 2) - 1
      this%service_distance(lb__1:ub__1, lb__2:ub__2) = service_distance
    end if
    if (present(bathymetry_filename)) then
      if (size(bathymetry_filename, 1) > size(this%bathymetry_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'bathymetry_filename'"
        return
      end if
      lb__1 = lbound(this%bathymetry_filename, 1)
      ub__1 = lb__1 + size(bathymetry_filename, 1) - 1
      if (size(bathymetry_filename, 2) > size(this%bathymetry_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'bathymetry_filename'"
        return
      end if
      lb__2 = lbound(this%bathymetry_filename, 2)
      ub__2 = lb__2 + size(bathymetry_filename, 2) - 1
      this%bathymetry_filename(lb__1:ub__1, lb__2:ub__2) = bathymetry_filename
    end if
    if (present(release_target_filename)) then
      if (size(release_target_filename, 1) > size(this%release_target_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'release_target_filename'"
        return
      end if
      lb__1 = lbound(this%release_target_filename, 1)
      ub__1 = lb__1 + size(release_target_filename, 1) - 1
      if (size(release_target_filename, 2) > size(this%release_target_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'release_target_filename'"
        return
      end if
      lb__2 = lbound(this%release_target_filename, 2)
      ub__2 = lb__2 + size(release_target_filename, 2) - 1
      this%release_target_filename(lb__1:ub__1, lb__2:ub__2) = release_target_filename
    end if
    if (present(watersupply_demand_filename)) then
      if (size(watersupply_demand_filename, 1) > size(this%watersupply_demand_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'watersupply_demand_filename'"
        return
      end if
      lb__1 = lbound(this%watersupply_demand_filename, 1)
      ub__1 = lb__1 + size(watersupply_demand_filename, 1) - 1
      if (size(watersupply_demand_filename, 2) > size(this%watersupply_demand_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'watersupply_demand_filename'"
        return
      end if
      lb__2 = lbound(this%watersupply_demand_filename, 2)
      ub__2 = lb__2 + size(watersupply_demand_filename, 2) - 1
      this%watersupply_demand_filename(lb__1:ub__1, lb__2:ub__2) = watersupply_demand_filename
    end if
    if (present(irrigation_demand_filename)) then
      if (size(irrigation_demand_filename, 1) > size(this%irrigation_demand_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'irrigation_demand_filename'"
        return
      end if
      lb__1 = lbound(this%irrigation_demand_filename, 1)
      ub__1 = lb__1 + size(irrigation_demand_filename, 1) - 1
      if (size(irrigation_demand_filename, 2) > size(this%irrigation_demand_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'irrigation_demand_filename'"
        return
      end if
      lb__2 = lbound(this%irrigation_demand_filename, 2)
      ub__2 = lb__2 + size(irrigation_demand_filename, 2) - 1
      this%irrigation_demand_filename(lb__1:ub__1, lb__2:ub__2) = irrigation_demand_filename
    end if
    if (present(ibt_filename)) then
      if (size(ibt_filename, 1) > size(this%ibt_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ibt_filename'"
        return
      end if
      lb__1 = lbound(this%ibt_filename, 1)
      ub__1 = lb__1 + size(ibt_filename, 1) - 1
      if (size(ibt_filename, 2) > size(this%ibt_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'ibt_filename'"
        return
      end if
      lb__2 = lbound(this%ibt_filename, 2)
      ub__2 = lb__2 + size(ibt_filename, 2) - 1
      this%ibt_filename(lb__1:ub__1, lb__2:ub__2) = ibt_filename
    end if
    if (present(level_filename)) then
      if (size(level_filename, 1) > size(this%level_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'level_filename'"
        return
      end if
      lb__1 = lbound(this%level_filename, 1)
      ub__1 = lb__1 + size(level_filename, 1) - 1
      if (size(level_filename, 2) > size(this%level_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'level_filename'"
        return
      end if
      lb__2 = lbound(this%level_filename, 2)
      ub__2 = lb__2 + size(level_filename, 2) - 1
      this%level_filename(lb__1:ub__1, lb__2:ub__2) = level_filename
    end if
    if (present(area_filename)) then
      if (size(area_filename, 1) > size(this%area_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'area_filename'"
        return
      end if
      lb__1 = lbound(this%area_filename, 1)
      ub__1 = lb__1 + size(area_filename, 1) - 1
      if (size(area_filename, 2) > size(this%area_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'area_filename'"
        return
      end if
      lb__2 = lbound(this%area_filename, 2)
      ub__2 = lb__2 + size(area_filename, 2) - 1
      this%area_filename(lb__1:ub__1, lb__2:ub__2) = area_filename
    end if
    if (present(spill_filename)) then
      if (size(spill_filename, 1) > size(this%spill_filename, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'spill_filename'"
        return
      end if
      lb__1 = lbound(this%spill_filename, 1)
      ub__1 = lb__1 + size(spill_filename, 1) - 1
      if (size(spill_filename, 2) > size(this%spill_filename, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'spill_filename'"
        return
      end if
      lb__2 = lbound(this%spill_filename, 2)
      ub__2 = lb__2 + size(spill_filename, 2) - 1
      this%spill_filename(lb__1:ub__1, lb__2:ub__2) = spill_filename
    end if
    if (present(print_lake_maps)) this%print_lake_maps = print_lake_maps
    if (present(print_links_table)) this%print_links_table = print_links_table

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
    case ("input_path")
      if (.not. allocated(this%input_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%input_path), ubound(this%input_path), &
          "input_path", errmsg)
        if (status /= NML_OK) return
        if (this%input_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%input_path == achar(0))) status = NML_ERR_NOT_SET
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
    case ("input_windspeed_height")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'input_windspeed_height'"
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
    case ("hanasaki_2006")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'hanasaki_2006'"
        return
      end if
    case ("point_demand")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'point_demand'"
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
    case ("use_watersupply")
      if (.not. allocated(this%use_watersupply)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%use_watersupply), ubound(this%use_watersupply), &
          "use_watersupply", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("use_irrigation")
      if (.not. allocated(this%use_irrigation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%use_irrigation), ubound(this%use_irrigation), &
          "use_irrigation", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("service_distance")
      if (.not. allocated(this%service_distance)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%service_distance), ubound(this%service_distance), &
          "service_distance", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("bathymetry_filename")
      if (.not. allocated(this%bathymetry_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%bathymetry_filename), ubound(this%bathymetry_filename), &
          "bathymetry_filename", errmsg)
        if (status /= NML_OK) return
        if (this%bathymetry_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%bathymetry_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("release_target_filename")
      if (.not. allocated(this%release_target_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%release_target_filename), ubound(this%release_target_filename), &
          "release_target_filename", errmsg)
        if (status /= NML_OK) return
        if (this%release_target_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%release_target_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("watersupply_demand_filename")
      if (.not. allocated(this%watersupply_demand_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%watersupply_demand_filename), ubound(this%watersupply_demand_filename), &
          "watersupply_demand_filename", errmsg)
        if (status /= NML_OK) return
        if (this%watersupply_demand_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%watersupply_demand_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("irrigation_demand_filename")
      if (.not. allocated(this%irrigation_demand_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%irrigation_demand_filename), ubound(this%irrigation_demand_filename), &
          "irrigation_demand_filename", errmsg)
        if (status /= NML_OK) return
        if (this%irrigation_demand_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%irrigation_demand_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("ibt_filename")
      if (.not. allocated(this%ibt_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ibt_filename), ubound(this%ibt_filename), &
          "ibt_filename", errmsg)
        if (status /= NML_OK) return
        if (this%ibt_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%ibt_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("level_filename")
      if (.not. allocated(this%level_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%level_filename), ubound(this%level_filename), &
          "level_filename", errmsg)
        if (status /= NML_OK) return
        if (this%level_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%level_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("area_filename")
      if (.not. allocated(this%area_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%area_filename), ubound(this%area_filename), &
          "area_filename", errmsg)
        if (status /= NML_OK) return
        if (this%area_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%area_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("spill_filename")
      if (.not. allocated(this%spill_filename)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_filename), ubound(this%spill_filename), &
          "spill_filename", errmsg)
        if (status /= NML_OK) return
        if (this%spill_filename(idx(1), idx(2)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%spill_filename == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("print_lake_maps")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'print_lake_maps'"
        return
      end if
    case ("print_links_table")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'print_links_table'"
        return
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
  end function nml_config_mlm_is_valid

end module nml_config_mlm
