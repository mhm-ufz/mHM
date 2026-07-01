!> \file nml_config_input.f90
!> \copydoc nml_config_input

!> \brief Input configuration
!> \details Paths and variable names for input data used by mHM.
!! Arrays are indexed by domain (dimension 1). Most paths are optional.
!! Variable name entries define the NetCDF variable names to read.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_input
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
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  integer(i4), parameter, public :: chunking__default = 0_i4
  integer(i4), parameter, public :: time_stamp_location__default = 0_i4
  logical, parameter, public :: morph_latlon__default = .false.
  character(len=buf), parameter, public :: pre_var__default = "pre"
  character(len=buf), parameter, public :: pet_var__default = "pet"
  character(len=buf), parameter, public :: temp_var__default = "tavg"
  character(len=buf), parameter, public :: tann_var__default = "tann"
  character(len=buf), parameter, public :: tmin_var__default = "tmin"
  character(len=buf), parameter, public :: tmax_var__default = "tmax"
  character(len=buf), parameter, public :: ssrd_var__default = "ssrd"
  character(len=buf), parameter, public :: strd_var__default = "strd"
  character(len=buf), parameter, public :: netrad_var__default = "net_rad"
  character(len=buf), parameter, public :: eabs_var__default = "eabs"
  character(len=buf), parameter, public :: wind_var__default = "windspeed"
  character(len=buf), parameter, public :: meteo_mask_var__default = "mask"
  character(len=buf), parameter, public :: runoff_var__default = "runoff"
  character(len=buf), parameter, public :: runoff_sealed_var__default = "runoff_sealed"
  character(len=buf), parameter, public :: interflow_fast_var__default = "interflow_fast"
  character(len=buf), parameter, public :: interflow_slow_var__default = "interflow_slow"
  character(len=buf), parameter, public :: baseflow_var__default = "baseflow"
  character(len=buf), parameter, public :: hydro_mask_var__default = "mask"
  character(len=buf), parameter, public :: dem_var__default = "dem"
  character(len=buf), parameter, public :: slope_var__default = "slope"
  character(len=buf), parameter, public :: aspect_var__default = "aspect"
  character(len=buf), parameter, public :: fdir_var__default = "fdir"
  character(len=buf), parameter, public :: facc_var__default = "facc"
  character(len=buf), parameter, public :: geo_class_var__default = "geology_class"
  character(len=buf), parameter, public :: soil_class_var__default = "soil_class"
  character(len=buf), parameter, public :: lai_class_var__default = "LAI_class"
  character(len=buf), parameter, public :: river_width_var__default = "P_bkfl"
  character(len=buf), parameter, public :: morph_mask_var__default = "mask"
  character(len=buf), parameter, public :: hydro_lat_var__default = "lat"
  character(len=buf), parameter, public :: hydro_lon_var__default = "lon"
  character(len=buf), parameter, public :: morph_lat_var__default = "lat_l0"
  character(len=buf), parameter, public :: morph_lon_var__default = "lon_l0"
  character(len=buf), parameter, public :: route_lat_var__default = "lat_l11"
  character(len=buf), parameter, public :: route_lon_var__default = "lon_l11"

  ! enum values
  integer(i4), parameter, public :: time_stamp_location__enum_values(3) = [0_i4, 1_i4, 2_i4]

  ! bounds values
  integer(i4), parameter, public :: chunking__min = -3_i4

  !> \class nml_config_input_t
  !> \brief Input configuration
  !> \details Paths and variable names for input data used by mHM.
  !! Arrays are indexed by domain (dimension 1). Most paths are optional.
  !! Variable name entries define the NetCDF variable names to read.
  type, public :: nml_config_input_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer(i4), allocatable, dimension(:) :: chunking !< Chunking for input data
    integer(i4), allocatable, dimension(:) :: time_stamp_location !< NetCDF time-stamp location
    character(len=buf), allocatable, dimension(:) :: latlon_path !< Latlon specification file path
    logical, allocatable, dimension(:) :: morph_latlon !< DEM in latlon coordinates
    character(len=buf), allocatable, dimension(:) :: pre_path !< Precipitation input
    character(len=buf), allocatable, dimension(:) :: pet_path !< Potential evapotranspiration input
    character(len=buf), allocatable, dimension(:) :: temp_path !< Air temperature input
    character(len=buf), allocatable, dimension(:) :: tann_path !< Air temperature annual mean input
    character(len=buf), allocatable, dimension(:) :: tmin_path !< Air temperature daily minimum input
    character(len=buf), allocatable, dimension(:) :: tmax_path !< Air temperature daily maximum input
    character(len=buf), allocatable, dimension(:) :: ssrd_path !< Surface shortwave radiation downwards input
    character(len=buf), allocatable, dimension(:) :: strd_path !< Surface thermal radiation downwards input
    character(len=buf), allocatable, dimension(:) :: netrad_path !< Net radiation input
    character(len=buf), allocatable, dimension(:) :: eabs_path !< Vapor pressure input
    character(len=buf), allocatable, dimension(:) :: wind_path !< Wind speed input
    character(len=buf), allocatable, dimension(:) :: meteo_mask_path !< Meteorological mask file path
    character(len=buf), allocatable, dimension(:) :: runoff_path !< Runoff input
    character(len=buf), allocatable, dimension(:) :: runoff_sealed_path !< Sealed runoff input
    character(len=buf), allocatable, dimension(:) :: interflow_fast_path !< Fast interflow input
    character(len=buf), allocatable, dimension(:) :: interflow_slow_path !< Slow interflow input
    character(len=buf), allocatable, dimension(:) :: baseflow_path !< Baseflow input
    character(len=buf), allocatable, dimension(:) :: hydro_mask_path !< Hydrological mask file path
    character(len=buf), allocatable, dimension(:) :: dem_path !< DEM input
    character(len=buf), allocatable, dimension(:) :: slope_path !< Slope input
    character(len=buf), allocatable, dimension(:) :: aspect_path !< Aspect input
    character(len=buf), allocatable, dimension(:) :: fdir_path !< Flow direction input
    character(len=buf), allocatable, dimension(:) :: facc_path !< Flow accumulation input
    character(len=buf), allocatable, dimension(:) :: geo_class_path !< Geology class input
    character(len=buf), allocatable, dimension(:) :: soil_class_path !< Soil class input
    character(len=buf), allocatable, dimension(:) :: soil_horizon_class_path !< Soil horizon class input
    character(len=buf), allocatable, dimension(:) :: lai_class_path !< LAI class input
    character(len=buf), allocatable, dimension(:) :: river_width_path !< River width input
    character(len=buf), allocatable, dimension(:) :: morph_mask_path !< Morphology mask file path
    character(len=buf), allocatable, dimension(:) :: pre_var !< Precipitation variable name
    character(len=buf), allocatable, dimension(:) :: pet_var !< Potential evapotranspiration variable name
    character(len=buf), allocatable, dimension(:) :: temp_var !< Air temperature variable name
    character(len=buf), allocatable, dimension(:) :: tann_var !< Air temperature annual mean variable name
    character(len=buf), allocatable, dimension(:) :: tmin_var !< Air temperature daily minimum variable name
    character(len=buf), allocatable, dimension(:) :: tmax_var !< Air temperature daily maximum variable name
    character(len=buf), allocatable, dimension(:) :: ssrd_var !< Surface shortwave radiation variable name
    character(len=buf), allocatable, dimension(:) :: strd_var !< Surface thermal radiation variable name
    character(len=buf), allocatable, dimension(:) :: netrad_var !< Net radiation variable name
    character(len=buf), allocatable, dimension(:) :: eabs_var !< Vapor pressure variable name
    character(len=buf), allocatable, dimension(:) :: wind_var !< Wind speed variable name
    character(len=buf), allocatable, dimension(:) :: meteo_mask_var !< Meteorological mask variable name
    character(len=buf), allocatable, dimension(:) :: runoff_var !< Runoff variable name
    character(len=buf), allocatable, dimension(:) :: runoff_sealed_var !< Sealed runoff variable name
    character(len=buf), allocatable, dimension(:) :: interflow_fast_var !< Fast interflow variable name
    character(len=buf), allocatable, dimension(:) :: interflow_slow_var !< Slow interflow variable name
    character(len=buf), allocatable, dimension(:) :: baseflow_var !< Baseflow variable name
    character(len=buf), allocatable, dimension(:) :: hydro_mask_var !< Hydrological mask variable name
    character(len=buf), allocatable, dimension(:) :: dem_var !< DEM variable name
    character(len=buf), allocatable, dimension(:) :: slope_var !< Slope variable name
    character(len=buf), allocatable, dimension(:) :: aspect_var !< Aspect variable name
    character(len=buf), allocatable, dimension(:) :: fdir_var !< Flow direction variable name
    character(len=buf), allocatable, dimension(:) :: facc_var !< Flow accumulation variable name
    character(len=buf), allocatable, dimension(:) :: geo_class_var !< Geology class variable name
    character(len=buf), allocatable, dimension(:) :: soil_class_var !< Soil class variable name
    character(len=buf), allocatable, dimension(:) :: lai_class_var !< LAI class variable name
    character(len=buf), allocatable, dimension(:) :: river_width_var !< River width variable name
    character(len=buf), allocatable, dimension(:) :: morph_mask_var !< Morphology mask variable name
    character(len=buf), allocatable, dimension(:) :: hydro_lat_var !< Hydrological latitude variable name
    character(len=buf), allocatable, dimension(:) :: hydro_lon_var !< Hydrological longitude variable name
    character(len=buf), allocatable, dimension(:) :: morph_lat_var !< Morphology latitude variable name
    character(len=buf), allocatable, dimension(:) :: morph_lon_var !< Morphology longitude variable name
    character(len=buf), allocatable, dimension(:) :: route_lat_var !< Routing latitude variable name
    character(len=buf), allocatable, dimension(:) :: route_lon_var !< Routing longitude variable name
  contains
    procedure :: init => nml_config_input_init
    procedure :: set_dims => nml_config_input_set_dims
    procedure :: from_file => nml_config_input_from_file
    procedure :: set => nml_config_input_set
    procedure :: is_set => nml_config_input_is_set
    procedure :: is_valid => nml_config_input_is_valid
  end type nml_config_input_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function time_stamp_location__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == time_stamp_location__enum_values)
  end function time_stamp_location__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function chunking__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < chunking__min) in_bounds = .false.
  end function chunking__in_bounds

  !> \brief Initialize defaults and sentinels for config_input
  integer function nml_config_input_init(this, errmsg) result(status)
    class(nml_config_input_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%chunking)) deallocate(this%chunking)
    allocate(this%chunking(this%n_domains))
    if (allocated(this%time_stamp_location)) deallocate(this%time_stamp_location)
    allocate(this%time_stamp_location(this%n_domains))
    if (allocated(this%latlon_path)) deallocate(this%latlon_path)
    allocate(character(len=buf) :: this%latlon_path(this%n_domains))
    if (allocated(this%morph_latlon)) deallocate(this%morph_latlon)
    allocate(this%morph_latlon(this%n_domains))
    if (allocated(this%pre_path)) deallocate(this%pre_path)
    allocate(character(len=buf) :: this%pre_path(this%n_domains))
    if (allocated(this%pet_path)) deallocate(this%pet_path)
    allocate(character(len=buf) :: this%pet_path(this%n_domains))
    if (allocated(this%temp_path)) deallocate(this%temp_path)
    allocate(character(len=buf) :: this%temp_path(this%n_domains))
    if (allocated(this%tann_path)) deallocate(this%tann_path)
    allocate(character(len=buf) :: this%tann_path(this%n_domains))
    if (allocated(this%tmin_path)) deallocate(this%tmin_path)
    allocate(character(len=buf) :: this%tmin_path(this%n_domains))
    if (allocated(this%tmax_path)) deallocate(this%tmax_path)
    allocate(character(len=buf) :: this%tmax_path(this%n_domains))
    if (allocated(this%ssrd_path)) deallocate(this%ssrd_path)
    allocate(character(len=buf) :: this%ssrd_path(this%n_domains))
    if (allocated(this%strd_path)) deallocate(this%strd_path)
    allocate(character(len=buf) :: this%strd_path(this%n_domains))
    if (allocated(this%netrad_path)) deallocate(this%netrad_path)
    allocate(character(len=buf) :: this%netrad_path(this%n_domains))
    if (allocated(this%eabs_path)) deallocate(this%eabs_path)
    allocate(character(len=buf) :: this%eabs_path(this%n_domains))
    if (allocated(this%wind_path)) deallocate(this%wind_path)
    allocate(character(len=buf) :: this%wind_path(this%n_domains))
    if (allocated(this%meteo_mask_path)) deallocate(this%meteo_mask_path)
    allocate(character(len=buf) :: this%meteo_mask_path(this%n_domains))
    if (allocated(this%runoff_path)) deallocate(this%runoff_path)
    allocate(character(len=buf) :: this%runoff_path(this%n_domains))
    if (allocated(this%runoff_sealed_path)) deallocate(this%runoff_sealed_path)
    allocate(character(len=buf) :: this%runoff_sealed_path(this%n_domains))
    if (allocated(this%interflow_fast_path)) deallocate(this%interflow_fast_path)
    allocate(character(len=buf) :: this%interflow_fast_path(this%n_domains))
    if (allocated(this%interflow_slow_path)) deallocate(this%interflow_slow_path)
    allocate(character(len=buf) :: this%interflow_slow_path(this%n_domains))
    if (allocated(this%baseflow_path)) deallocate(this%baseflow_path)
    allocate(character(len=buf) :: this%baseflow_path(this%n_domains))
    if (allocated(this%hydro_mask_path)) deallocate(this%hydro_mask_path)
    allocate(character(len=buf) :: this%hydro_mask_path(this%n_domains))
    if (allocated(this%dem_path)) deallocate(this%dem_path)
    allocate(character(len=buf) :: this%dem_path(this%n_domains))
    if (allocated(this%slope_path)) deallocate(this%slope_path)
    allocate(character(len=buf) :: this%slope_path(this%n_domains))
    if (allocated(this%aspect_path)) deallocate(this%aspect_path)
    allocate(character(len=buf) :: this%aspect_path(this%n_domains))
    if (allocated(this%fdir_path)) deallocate(this%fdir_path)
    allocate(character(len=buf) :: this%fdir_path(this%n_domains))
    if (allocated(this%facc_path)) deallocate(this%facc_path)
    allocate(character(len=buf) :: this%facc_path(this%n_domains))
    if (allocated(this%geo_class_path)) deallocate(this%geo_class_path)
    allocate(character(len=buf) :: this%geo_class_path(this%n_domains))
    if (allocated(this%soil_class_path)) deallocate(this%soil_class_path)
    allocate(character(len=buf) :: this%soil_class_path(this%n_domains))
    if (allocated(this%soil_horizon_class_path)) deallocate(this%soil_horizon_class_path)
    allocate(character(len=buf) :: this%soil_horizon_class_path(this%n_domains))
    if (allocated(this%lai_class_path)) deallocate(this%lai_class_path)
    allocate(character(len=buf) :: this%lai_class_path(this%n_domains))
    if (allocated(this%river_width_path)) deallocate(this%river_width_path)
    allocate(character(len=buf) :: this%river_width_path(this%n_domains))
    if (allocated(this%morph_mask_path)) deallocate(this%morph_mask_path)
    allocate(character(len=buf) :: this%morph_mask_path(this%n_domains))
    if (allocated(this%pre_var)) deallocate(this%pre_var)
    allocate(character(len=buf) :: this%pre_var(this%n_domains))
    if (allocated(this%pet_var)) deallocate(this%pet_var)
    allocate(character(len=buf) :: this%pet_var(this%n_domains))
    if (allocated(this%temp_var)) deallocate(this%temp_var)
    allocate(character(len=buf) :: this%temp_var(this%n_domains))
    if (allocated(this%tann_var)) deallocate(this%tann_var)
    allocate(character(len=buf) :: this%tann_var(this%n_domains))
    if (allocated(this%tmin_var)) deallocate(this%tmin_var)
    allocate(character(len=buf) :: this%tmin_var(this%n_domains))
    if (allocated(this%tmax_var)) deallocate(this%tmax_var)
    allocate(character(len=buf) :: this%tmax_var(this%n_domains))
    if (allocated(this%ssrd_var)) deallocate(this%ssrd_var)
    allocate(character(len=buf) :: this%ssrd_var(this%n_domains))
    if (allocated(this%strd_var)) deallocate(this%strd_var)
    allocate(character(len=buf) :: this%strd_var(this%n_domains))
    if (allocated(this%netrad_var)) deallocate(this%netrad_var)
    allocate(character(len=buf) :: this%netrad_var(this%n_domains))
    if (allocated(this%eabs_var)) deallocate(this%eabs_var)
    allocate(character(len=buf) :: this%eabs_var(this%n_domains))
    if (allocated(this%wind_var)) deallocate(this%wind_var)
    allocate(character(len=buf) :: this%wind_var(this%n_domains))
    if (allocated(this%meteo_mask_var)) deallocate(this%meteo_mask_var)
    allocate(character(len=buf) :: this%meteo_mask_var(this%n_domains))
    if (allocated(this%runoff_var)) deallocate(this%runoff_var)
    allocate(character(len=buf) :: this%runoff_var(this%n_domains))
    if (allocated(this%runoff_sealed_var)) deallocate(this%runoff_sealed_var)
    allocate(character(len=buf) :: this%runoff_sealed_var(this%n_domains))
    if (allocated(this%interflow_fast_var)) deallocate(this%interflow_fast_var)
    allocate(character(len=buf) :: this%interflow_fast_var(this%n_domains))
    if (allocated(this%interflow_slow_var)) deallocate(this%interflow_slow_var)
    allocate(character(len=buf) :: this%interflow_slow_var(this%n_domains))
    if (allocated(this%baseflow_var)) deallocate(this%baseflow_var)
    allocate(character(len=buf) :: this%baseflow_var(this%n_domains))
    if (allocated(this%hydro_mask_var)) deallocate(this%hydro_mask_var)
    allocate(character(len=buf) :: this%hydro_mask_var(this%n_domains))
    if (allocated(this%dem_var)) deallocate(this%dem_var)
    allocate(character(len=buf) :: this%dem_var(this%n_domains))
    if (allocated(this%slope_var)) deallocate(this%slope_var)
    allocate(character(len=buf) :: this%slope_var(this%n_domains))
    if (allocated(this%aspect_var)) deallocate(this%aspect_var)
    allocate(character(len=buf) :: this%aspect_var(this%n_domains))
    if (allocated(this%fdir_var)) deallocate(this%fdir_var)
    allocate(character(len=buf) :: this%fdir_var(this%n_domains))
    if (allocated(this%facc_var)) deallocate(this%facc_var)
    allocate(character(len=buf) :: this%facc_var(this%n_domains))
    if (allocated(this%geo_class_var)) deallocate(this%geo_class_var)
    allocate(character(len=buf) :: this%geo_class_var(this%n_domains))
    if (allocated(this%soil_class_var)) deallocate(this%soil_class_var)
    allocate(character(len=buf) :: this%soil_class_var(this%n_domains))
    if (allocated(this%lai_class_var)) deallocate(this%lai_class_var)
    allocate(character(len=buf) :: this%lai_class_var(this%n_domains))
    if (allocated(this%river_width_var)) deallocate(this%river_width_var)
    allocate(character(len=buf) :: this%river_width_var(this%n_domains))
    if (allocated(this%morph_mask_var)) deallocate(this%morph_mask_var)
    allocate(character(len=buf) :: this%morph_mask_var(this%n_domains))
    if (allocated(this%hydro_lat_var)) deallocate(this%hydro_lat_var)
    allocate(character(len=buf) :: this%hydro_lat_var(this%n_domains))
    if (allocated(this%hydro_lon_var)) deallocate(this%hydro_lon_var)
    allocate(character(len=buf) :: this%hydro_lon_var(this%n_domains))
    if (allocated(this%morph_lat_var)) deallocate(this%morph_lat_var)
    allocate(character(len=buf) :: this%morph_lat_var(this%n_domains))
    if (allocated(this%morph_lon_var)) deallocate(this%morph_lon_var)
    allocate(character(len=buf) :: this%morph_lon_var(this%n_domains))
    if (allocated(this%route_lat_var)) deallocate(this%route_lat_var)
    allocate(character(len=buf) :: this%route_lat_var(this%n_domains))
    if (allocated(this%route_lon_var)) deallocate(this%route_lon_var)
    allocate(character(len=buf) :: this%route_lon_var(this%n_domains))

    ! sentinel values for required/optional parameters
    this%latlon_path = achar(0) ! sentinel for optional string array
    this%pre_path = achar(0) ! sentinel for optional string array
    this%pet_path = achar(0) ! sentinel for optional string array
    this%temp_path = achar(0) ! sentinel for optional string array
    this%tann_path = achar(0) ! sentinel for optional string array
    this%tmin_path = achar(0) ! sentinel for optional string array
    this%tmax_path = achar(0) ! sentinel for optional string array
    this%ssrd_path = achar(0) ! sentinel for optional string array
    this%strd_path = achar(0) ! sentinel for optional string array
    this%netrad_path = achar(0) ! sentinel for optional string array
    this%eabs_path = achar(0) ! sentinel for optional string array
    this%wind_path = achar(0) ! sentinel for optional string array
    this%meteo_mask_path = achar(0) ! sentinel for optional string array
    this%runoff_path = achar(0) ! sentinel for optional string array
    this%runoff_sealed_path = achar(0) ! sentinel for optional string array
    this%interflow_fast_path = achar(0) ! sentinel for optional string array
    this%interflow_slow_path = achar(0) ! sentinel for optional string array
    this%baseflow_path = achar(0) ! sentinel for optional string array
    this%hydro_mask_path = achar(0) ! sentinel for optional string array
    this%dem_path = achar(0) ! sentinel for optional string array
    this%slope_path = achar(0) ! sentinel for optional string array
    this%aspect_path = achar(0) ! sentinel for optional string array
    this%fdir_path = achar(0) ! sentinel for optional string array
    this%facc_path = achar(0) ! sentinel for optional string array
    this%geo_class_path = achar(0) ! sentinel for optional string array
    this%soil_class_path = achar(0) ! sentinel for optional string array
    this%soil_horizon_class_path = achar(0) ! sentinel for optional string array
    this%lai_class_path = achar(0) ! sentinel for optional string array
    this%river_width_path = achar(0) ! sentinel for optional string array
    this%morph_mask_path = achar(0) ! sentinel for optional string array
    ! default values
    this%chunking = chunking__default
    this%time_stamp_location = time_stamp_location__default
    this%morph_latlon = morph_latlon__default
    this%pre_var = pre_var__default
    this%pet_var = pet_var__default
    this%temp_var = temp_var__default
    this%tann_var = tann_var__default
    this%tmin_var = tmin_var__default
    this%tmax_var = tmax_var__default
    this%ssrd_var = ssrd_var__default
    this%strd_var = strd_var__default
    this%netrad_var = netrad_var__default
    this%eabs_var = eabs_var__default
    this%wind_var = wind_var__default
    this%meteo_mask_var = meteo_mask_var__default
    this%runoff_var = runoff_var__default
    this%runoff_sealed_var = runoff_sealed_var__default
    this%interflow_fast_var = interflow_fast_var__default
    this%interflow_slow_var = interflow_slow_var__default
    this%baseflow_var = baseflow_var__default
    this%hydro_mask_var = hydro_mask_var__default
    this%dem_var = dem_var__default
    this%slope_var = slope_var__default
    this%aspect_var = aspect_var__default
    this%fdir_var = fdir_var__default
    this%facc_var = facc_var__default
    this%geo_class_var = geo_class_var__default
    this%soil_class_var = soil_class_var__default
    this%lai_class_var = lai_class_var__default
    this%river_width_var = river_width_var__default
    this%morph_mask_var = morph_mask_var__default
    this%hydro_lat_var = hydro_lat_var__default
    this%hydro_lon_var = hydro_lon_var__default
    this%morph_lat_var = morph_lat_var__default
    this%morph_lon_var = morph_lon_var__default
    this%route_lat_var = route_lat_var__default
    this%route_lon_var = route_lon_var__default
  end function nml_config_input_init

  !> \brief Reset runtime dimensions for config_input
  integer function nml_config_input_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_input_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%chunking)) deallocate(this%chunking)
    if (allocated(this%time_stamp_location)) deallocate(this%time_stamp_location)
    if (allocated(this%latlon_path)) deallocate(this%latlon_path)
    if (allocated(this%morph_latlon)) deallocate(this%morph_latlon)
    if (allocated(this%pre_path)) deallocate(this%pre_path)
    if (allocated(this%pet_path)) deallocate(this%pet_path)
    if (allocated(this%temp_path)) deallocate(this%temp_path)
    if (allocated(this%tann_path)) deallocate(this%tann_path)
    if (allocated(this%tmin_path)) deallocate(this%tmin_path)
    if (allocated(this%tmax_path)) deallocate(this%tmax_path)
    if (allocated(this%ssrd_path)) deallocate(this%ssrd_path)
    if (allocated(this%strd_path)) deallocate(this%strd_path)
    if (allocated(this%netrad_path)) deallocate(this%netrad_path)
    if (allocated(this%eabs_path)) deallocate(this%eabs_path)
    if (allocated(this%wind_path)) deallocate(this%wind_path)
    if (allocated(this%meteo_mask_path)) deallocate(this%meteo_mask_path)
    if (allocated(this%runoff_path)) deallocate(this%runoff_path)
    if (allocated(this%runoff_sealed_path)) deallocate(this%runoff_sealed_path)
    if (allocated(this%interflow_fast_path)) deallocate(this%interflow_fast_path)
    if (allocated(this%interflow_slow_path)) deallocate(this%interflow_slow_path)
    if (allocated(this%baseflow_path)) deallocate(this%baseflow_path)
    if (allocated(this%hydro_mask_path)) deallocate(this%hydro_mask_path)
    if (allocated(this%dem_path)) deallocate(this%dem_path)
    if (allocated(this%slope_path)) deallocate(this%slope_path)
    if (allocated(this%aspect_path)) deallocate(this%aspect_path)
    if (allocated(this%fdir_path)) deallocate(this%fdir_path)
    if (allocated(this%facc_path)) deallocate(this%facc_path)
    if (allocated(this%geo_class_path)) deallocate(this%geo_class_path)
    if (allocated(this%soil_class_path)) deallocate(this%soil_class_path)
    if (allocated(this%soil_horizon_class_path)) deallocate(this%soil_horizon_class_path)
    if (allocated(this%lai_class_path)) deallocate(this%lai_class_path)
    if (allocated(this%river_width_path)) deallocate(this%river_width_path)
    if (allocated(this%morph_mask_path)) deallocate(this%morph_mask_path)
    if (allocated(this%pre_var)) deallocate(this%pre_var)
    if (allocated(this%pet_var)) deallocate(this%pet_var)
    if (allocated(this%temp_var)) deallocate(this%temp_var)
    if (allocated(this%tann_var)) deallocate(this%tann_var)
    if (allocated(this%tmin_var)) deallocate(this%tmin_var)
    if (allocated(this%tmax_var)) deallocate(this%tmax_var)
    if (allocated(this%ssrd_var)) deallocate(this%ssrd_var)
    if (allocated(this%strd_var)) deallocate(this%strd_var)
    if (allocated(this%netrad_var)) deallocate(this%netrad_var)
    if (allocated(this%eabs_var)) deallocate(this%eabs_var)
    if (allocated(this%wind_var)) deallocate(this%wind_var)
    if (allocated(this%meteo_mask_var)) deallocate(this%meteo_mask_var)
    if (allocated(this%runoff_var)) deallocate(this%runoff_var)
    if (allocated(this%runoff_sealed_var)) deallocate(this%runoff_sealed_var)
    if (allocated(this%interflow_fast_var)) deallocate(this%interflow_fast_var)
    if (allocated(this%interflow_slow_var)) deallocate(this%interflow_slow_var)
    if (allocated(this%baseflow_var)) deallocate(this%baseflow_var)
    if (allocated(this%hydro_mask_var)) deallocate(this%hydro_mask_var)
    if (allocated(this%dem_var)) deallocate(this%dem_var)
    if (allocated(this%slope_var)) deallocate(this%slope_var)
    if (allocated(this%aspect_var)) deallocate(this%aspect_var)
    if (allocated(this%fdir_var)) deallocate(this%fdir_var)
    if (allocated(this%facc_var)) deallocate(this%facc_var)
    if (allocated(this%geo_class_var)) deallocate(this%geo_class_var)
    if (allocated(this%soil_class_var)) deallocate(this%soil_class_var)
    if (allocated(this%lai_class_var)) deallocate(this%lai_class_var)
    if (allocated(this%river_width_var)) deallocate(this%river_width_var)
    if (allocated(this%morph_mask_var)) deallocate(this%morph_mask_var)
    if (allocated(this%hydro_lat_var)) deallocate(this%hydro_lat_var)
    if (allocated(this%hydro_lon_var)) deallocate(this%hydro_lon_var)
    if (allocated(this%morph_lat_var)) deallocate(this%morph_lat_var)
    if (allocated(this%morph_lon_var)) deallocate(this%morph_lon_var)
    if (allocated(this%route_lat_var)) deallocate(this%route_lat_var)
    if (allocated(this%route_lon_var)) deallocate(this%route_lon_var)
    this%is_configured = .false.
  end function nml_config_input_set_dims


  !> \brief Read config_input namelist from file
  integer function nml_config_input_from_file(this, file, errmsg) result(status)
    class(nml_config_input_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4), allocatable, dimension(:) :: chunking
    integer(i4), allocatable, dimension(:) :: time_stamp_location
    character(len=buf), allocatable, dimension(:) :: latlon_path
    logical, allocatable, dimension(:) :: morph_latlon
    character(len=buf), allocatable, dimension(:) :: pre_path
    character(len=buf), allocatable, dimension(:) :: pet_path
    character(len=buf), allocatable, dimension(:) :: temp_path
    character(len=buf), allocatable, dimension(:) :: tann_path
    character(len=buf), allocatable, dimension(:) :: tmin_path
    character(len=buf), allocatable, dimension(:) :: tmax_path
    character(len=buf), allocatable, dimension(:) :: ssrd_path
    character(len=buf), allocatable, dimension(:) :: strd_path
    character(len=buf), allocatable, dimension(:) :: netrad_path
    character(len=buf), allocatable, dimension(:) :: eabs_path
    character(len=buf), allocatable, dimension(:) :: wind_path
    character(len=buf), allocatable, dimension(:) :: meteo_mask_path
    character(len=buf), allocatable, dimension(:) :: runoff_path
    character(len=buf), allocatable, dimension(:) :: runoff_sealed_path
    character(len=buf), allocatable, dimension(:) :: interflow_fast_path
    character(len=buf), allocatable, dimension(:) :: interflow_slow_path
    character(len=buf), allocatable, dimension(:) :: baseflow_path
    character(len=buf), allocatable, dimension(:) :: hydro_mask_path
    character(len=buf), allocatable, dimension(:) :: dem_path
    character(len=buf), allocatable, dimension(:) :: slope_path
    character(len=buf), allocatable, dimension(:) :: aspect_path
    character(len=buf), allocatable, dimension(:) :: fdir_path
    character(len=buf), allocatable, dimension(:) :: facc_path
    character(len=buf), allocatable, dimension(:) :: geo_class_path
    character(len=buf), allocatable, dimension(:) :: soil_class_path
    character(len=buf), allocatable, dimension(:) :: soil_horizon_class_path
    character(len=buf), allocatable, dimension(:) :: lai_class_path
    character(len=buf), allocatable, dimension(:) :: river_width_path
    character(len=buf), allocatable, dimension(:) :: morph_mask_path
    character(len=buf), allocatable, dimension(:) :: pre_var
    character(len=buf), allocatable, dimension(:) :: pet_var
    character(len=buf), allocatable, dimension(:) :: temp_var
    character(len=buf), allocatable, dimension(:) :: tann_var
    character(len=buf), allocatable, dimension(:) :: tmin_var
    character(len=buf), allocatable, dimension(:) :: tmax_var
    character(len=buf), allocatable, dimension(:) :: ssrd_var
    character(len=buf), allocatable, dimension(:) :: strd_var
    character(len=buf), allocatable, dimension(:) :: netrad_var
    character(len=buf), allocatable, dimension(:) :: eabs_var
    character(len=buf), allocatable, dimension(:) :: wind_var
    character(len=buf), allocatable, dimension(:) :: meteo_mask_var
    character(len=buf), allocatable, dimension(:) :: runoff_var
    character(len=buf), allocatable, dimension(:) :: runoff_sealed_var
    character(len=buf), allocatable, dimension(:) :: interflow_fast_var
    character(len=buf), allocatable, dimension(:) :: interflow_slow_var
    character(len=buf), allocatable, dimension(:) :: baseflow_var
    character(len=buf), allocatable, dimension(:) :: hydro_mask_var
    character(len=buf), allocatable, dimension(:) :: dem_var
    character(len=buf), allocatable, dimension(:) :: slope_var
    character(len=buf), allocatable, dimension(:) :: aspect_var
    character(len=buf), allocatable, dimension(:) :: fdir_var
    character(len=buf), allocatable, dimension(:) :: facc_var
    character(len=buf), allocatable, dimension(:) :: geo_class_var
    character(len=buf), allocatable, dimension(:) :: soil_class_var
    character(len=buf), allocatable, dimension(:) :: lai_class_var
    character(len=buf), allocatable, dimension(:) :: river_width_var
    character(len=buf), allocatable, dimension(:) :: morph_mask_var
    character(len=buf), allocatable, dimension(:) :: hydro_lat_var
    character(len=buf), allocatable, dimension(:) :: hydro_lon_var
    character(len=buf), allocatable, dimension(:) :: morph_lat_var
    character(len=buf), allocatable, dimension(:) :: morph_lon_var
    character(len=buf), allocatable, dimension(:) :: route_lat_var
    character(len=buf), allocatable, dimension(:) :: route_lon_var
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_input/ &
      chunking, &
      time_stamp_location, &
      latlon_path, &
      morph_latlon, &
      pre_path, &
      pet_path, &
      temp_path, &
      tann_path, &
      tmin_path, &
      tmax_path, &
      ssrd_path, &
      strd_path, &
      netrad_path, &
      eabs_path, &
      wind_path, &
      meteo_mask_path, &
      runoff_path, &
      runoff_sealed_path, &
      interflow_fast_path, &
      interflow_slow_path, &
      baseflow_path, &
      hydro_mask_path, &
      dem_path, &
      slope_path, &
      aspect_path, &
      fdir_path, &
      facc_path, &
      geo_class_path, &
      soil_class_path, &
      soil_horizon_class_path, &
      lai_class_path, &
      river_width_path, &
      morph_mask_path, &
      pre_var, &
      pet_var, &
      temp_var, &
      tann_var, &
      tmin_var, &
      tmax_var, &
      ssrd_var, &
      strd_var, &
      netrad_var, &
      eabs_var, &
      wind_var, &
      meteo_mask_var, &
      runoff_var, &
      runoff_sealed_var, &
      interflow_fast_var, &
      interflow_slow_var, &
      baseflow_var, &
      hydro_mask_var, &
      dem_var, &
      slope_var, &
      aspect_var, &
      fdir_var, &
      facc_var, &
      geo_class_var, &
      soil_class_var, &
      lai_class_var, &
      river_width_var, &
      morph_mask_var, &
      hydro_lat_var, &
      hydro_lon_var, &
      morph_lat_var, &
      morph_lon_var, &
      route_lat_var, &
      route_lon_var

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(chunking)) deallocate(chunking)
    allocate(chunking(this%n_domains))
    if (allocated(time_stamp_location)) deallocate(time_stamp_location)
    allocate(time_stamp_location(this%n_domains))
    if (allocated(latlon_path)) deallocate(latlon_path)
    allocate(character(len=buf) :: latlon_path(this%n_domains))
    if (allocated(morph_latlon)) deallocate(morph_latlon)
    allocate(morph_latlon(this%n_domains))
    if (allocated(pre_path)) deallocate(pre_path)
    allocate(character(len=buf) :: pre_path(this%n_domains))
    if (allocated(pet_path)) deallocate(pet_path)
    allocate(character(len=buf) :: pet_path(this%n_domains))
    if (allocated(temp_path)) deallocate(temp_path)
    allocate(character(len=buf) :: temp_path(this%n_domains))
    if (allocated(tann_path)) deallocate(tann_path)
    allocate(character(len=buf) :: tann_path(this%n_domains))
    if (allocated(tmin_path)) deallocate(tmin_path)
    allocate(character(len=buf) :: tmin_path(this%n_domains))
    if (allocated(tmax_path)) deallocate(tmax_path)
    allocate(character(len=buf) :: tmax_path(this%n_domains))
    if (allocated(ssrd_path)) deallocate(ssrd_path)
    allocate(character(len=buf) :: ssrd_path(this%n_domains))
    if (allocated(strd_path)) deallocate(strd_path)
    allocate(character(len=buf) :: strd_path(this%n_domains))
    if (allocated(netrad_path)) deallocate(netrad_path)
    allocate(character(len=buf) :: netrad_path(this%n_domains))
    if (allocated(eabs_path)) deallocate(eabs_path)
    allocate(character(len=buf) :: eabs_path(this%n_domains))
    if (allocated(wind_path)) deallocate(wind_path)
    allocate(character(len=buf) :: wind_path(this%n_domains))
    if (allocated(meteo_mask_path)) deallocate(meteo_mask_path)
    allocate(character(len=buf) :: meteo_mask_path(this%n_domains))
    if (allocated(runoff_path)) deallocate(runoff_path)
    allocate(character(len=buf) :: runoff_path(this%n_domains))
    if (allocated(runoff_sealed_path)) deallocate(runoff_sealed_path)
    allocate(character(len=buf) :: runoff_sealed_path(this%n_domains))
    if (allocated(interflow_fast_path)) deallocate(interflow_fast_path)
    allocate(character(len=buf) :: interflow_fast_path(this%n_domains))
    if (allocated(interflow_slow_path)) deallocate(interflow_slow_path)
    allocate(character(len=buf) :: interflow_slow_path(this%n_domains))
    if (allocated(baseflow_path)) deallocate(baseflow_path)
    allocate(character(len=buf) :: baseflow_path(this%n_domains))
    if (allocated(hydro_mask_path)) deallocate(hydro_mask_path)
    allocate(character(len=buf) :: hydro_mask_path(this%n_domains))
    if (allocated(dem_path)) deallocate(dem_path)
    allocate(character(len=buf) :: dem_path(this%n_domains))
    if (allocated(slope_path)) deallocate(slope_path)
    allocate(character(len=buf) :: slope_path(this%n_domains))
    if (allocated(aspect_path)) deallocate(aspect_path)
    allocate(character(len=buf) :: aspect_path(this%n_domains))
    if (allocated(fdir_path)) deallocate(fdir_path)
    allocate(character(len=buf) :: fdir_path(this%n_domains))
    if (allocated(facc_path)) deallocate(facc_path)
    allocate(character(len=buf) :: facc_path(this%n_domains))
    if (allocated(geo_class_path)) deallocate(geo_class_path)
    allocate(character(len=buf) :: geo_class_path(this%n_domains))
    if (allocated(soil_class_path)) deallocate(soil_class_path)
    allocate(character(len=buf) :: soil_class_path(this%n_domains))
    if (allocated(soil_horizon_class_path)) deallocate(soil_horizon_class_path)
    allocate(character(len=buf) :: soil_horizon_class_path(this%n_domains))
    if (allocated(lai_class_path)) deallocate(lai_class_path)
    allocate(character(len=buf) :: lai_class_path(this%n_domains))
    if (allocated(river_width_path)) deallocate(river_width_path)
    allocate(character(len=buf) :: river_width_path(this%n_domains))
    if (allocated(morph_mask_path)) deallocate(morph_mask_path)
    allocate(character(len=buf) :: morph_mask_path(this%n_domains))
    if (allocated(pre_var)) deallocate(pre_var)
    allocate(character(len=buf) :: pre_var(this%n_domains))
    if (allocated(pet_var)) deallocate(pet_var)
    allocate(character(len=buf) :: pet_var(this%n_domains))
    if (allocated(temp_var)) deallocate(temp_var)
    allocate(character(len=buf) :: temp_var(this%n_domains))
    if (allocated(tann_var)) deallocate(tann_var)
    allocate(character(len=buf) :: tann_var(this%n_domains))
    if (allocated(tmin_var)) deallocate(tmin_var)
    allocate(character(len=buf) :: tmin_var(this%n_domains))
    if (allocated(tmax_var)) deallocate(tmax_var)
    allocate(character(len=buf) :: tmax_var(this%n_domains))
    if (allocated(ssrd_var)) deallocate(ssrd_var)
    allocate(character(len=buf) :: ssrd_var(this%n_domains))
    if (allocated(strd_var)) deallocate(strd_var)
    allocate(character(len=buf) :: strd_var(this%n_domains))
    if (allocated(netrad_var)) deallocate(netrad_var)
    allocate(character(len=buf) :: netrad_var(this%n_domains))
    if (allocated(eabs_var)) deallocate(eabs_var)
    allocate(character(len=buf) :: eabs_var(this%n_domains))
    if (allocated(wind_var)) deallocate(wind_var)
    allocate(character(len=buf) :: wind_var(this%n_domains))
    if (allocated(meteo_mask_var)) deallocate(meteo_mask_var)
    allocate(character(len=buf) :: meteo_mask_var(this%n_domains))
    if (allocated(runoff_var)) deallocate(runoff_var)
    allocate(character(len=buf) :: runoff_var(this%n_domains))
    if (allocated(runoff_sealed_var)) deallocate(runoff_sealed_var)
    allocate(character(len=buf) :: runoff_sealed_var(this%n_domains))
    if (allocated(interflow_fast_var)) deallocate(interflow_fast_var)
    allocate(character(len=buf) :: interflow_fast_var(this%n_domains))
    if (allocated(interflow_slow_var)) deallocate(interflow_slow_var)
    allocate(character(len=buf) :: interflow_slow_var(this%n_domains))
    if (allocated(baseflow_var)) deallocate(baseflow_var)
    allocate(character(len=buf) :: baseflow_var(this%n_domains))
    if (allocated(hydro_mask_var)) deallocate(hydro_mask_var)
    allocate(character(len=buf) :: hydro_mask_var(this%n_domains))
    if (allocated(dem_var)) deallocate(dem_var)
    allocate(character(len=buf) :: dem_var(this%n_domains))
    if (allocated(slope_var)) deallocate(slope_var)
    allocate(character(len=buf) :: slope_var(this%n_domains))
    if (allocated(aspect_var)) deallocate(aspect_var)
    allocate(character(len=buf) :: aspect_var(this%n_domains))
    if (allocated(fdir_var)) deallocate(fdir_var)
    allocate(character(len=buf) :: fdir_var(this%n_domains))
    if (allocated(facc_var)) deallocate(facc_var)
    allocate(character(len=buf) :: facc_var(this%n_domains))
    if (allocated(geo_class_var)) deallocate(geo_class_var)
    allocate(character(len=buf) :: geo_class_var(this%n_domains))
    if (allocated(soil_class_var)) deallocate(soil_class_var)
    allocate(character(len=buf) :: soil_class_var(this%n_domains))
    if (allocated(lai_class_var)) deallocate(lai_class_var)
    allocate(character(len=buf) :: lai_class_var(this%n_domains))
    if (allocated(river_width_var)) deallocate(river_width_var)
    allocate(character(len=buf) :: river_width_var(this%n_domains))
    if (allocated(morph_mask_var)) deallocate(morph_mask_var)
    allocate(character(len=buf) :: morph_mask_var(this%n_domains))
    if (allocated(hydro_lat_var)) deallocate(hydro_lat_var)
    allocate(character(len=buf) :: hydro_lat_var(this%n_domains))
    if (allocated(hydro_lon_var)) deallocate(hydro_lon_var)
    allocate(character(len=buf) :: hydro_lon_var(this%n_domains))
    if (allocated(morph_lat_var)) deallocate(morph_lat_var)
    allocate(character(len=buf) :: morph_lat_var(this%n_domains))
    if (allocated(morph_lon_var)) deallocate(morph_lon_var)
    allocate(character(len=buf) :: morph_lon_var(this%n_domains))
    if (allocated(route_lat_var)) deallocate(route_lat_var)
    allocate(character(len=buf) :: route_lat_var(this%n_domains))
    if (allocated(route_lon_var)) deallocate(route_lon_var)
    allocate(character(len=buf) :: route_lon_var(this%n_domains))
    chunking = this%chunking
    time_stamp_location = this%time_stamp_location
    latlon_path = this%latlon_path
    morph_latlon = this%morph_latlon
    pre_path = this%pre_path
    pet_path = this%pet_path
    temp_path = this%temp_path
    tann_path = this%tann_path
    tmin_path = this%tmin_path
    tmax_path = this%tmax_path
    ssrd_path = this%ssrd_path
    strd_path = this%strd_path
    netrad_path = this%netrad_path
    eabs_path = this%eabs_path
    wind_path = this%wind_path
    meteo_mask_path = this%meteo_mask_path
    runoff_path = this%runoff_path
    runoff_sealed_path = this%runoff_sealed_path
    interflow_fast_path = this%interflow_fast_path
    interflow_slow_path = this%interflow_slow_path
    baseflow_path = this%baseflow_path
    hydro_mask_path = this%hydro_mask_path
    dem_path = this%dem_path
    slope_path = this%slope_path
    aspect_path = this%aspect_path
    fdir_path = this%fdir_path
    facc_path = this%facc_path
    geo_class_path = this%geo_class_path
    soil_class_path = this%soil_class_path
    soil_horizon_class_path = this%soil_horizon_class_path
    lai_class_path = this%lai_class_path
    river_width_path = this%river_width_path
    morph_mask_path = this%morph_mask_path
    pre_var = this%pre_var
    pet_var = this%pet_var
    temp_var = this%temp_var
    tann_var = this%tann_var
    tmin_var = this%tmin_var
    tmax_var = this%tmax_var
    ssrd_var = this%ssrd_var
    strd_var = this%strd_var
    netrad_var = this%netrad_var
    eabs_var = this%eabs_var
    wind_var = this%wind_var
    meteo_mask_var = this%meteo_mask_var
    runoff_var = this%runoff_var
    runoff_sealed_var = this%runoff_sealed_var
    interflow_fast_var = this%interflow_fast_var
    interflow_slow_var = this%interflow_slow_var
    baseflow_var = this%baseflow_var
    hydro_mask_var = this%hydro_mask_var
    dem_var = this%dem_var
    slope_var = this%slope_var
    aspect_var = this%aspect_var
    fdir_var = this%fdir_var
    facc_var = this%facc_var
    geo_class_var = this%geo_class_var
    soil_class_var = this%soil_class_var
    lai_class_var = this%lai_class_var
    river_width_var = this%river_width_var
    morph_mask_var = this%morph_mask_var
    hydro_lat_var = this%hydro_lat_var
    hydro_lon_var = this%hydro_lon_var
    morph_lat_var = this%morph_lat_var
    morph_lon_var = this%morph_lon_var
    route_lat_var = this%route_lat_var
    route_lon_var = this%route_lon_var

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_input", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_input, iostat=iostat, iomsg=iomsg)
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
    this%chunking = chunking
    this%time_stamp_location = time_stamp_location
    this%latlon_path = latlon_path
    this%morph_latlon = morph_latlon
    this%pre_path = pre_path
    this%pet_path = pet_path
    this%temp_path = temp_path
    this%tann_path = tann_path
    this%tmin_path = tmin_path
    this%tmax_path = tmax_path
    this%ssrd_path = ssrd_path
    this%strd_path = strd_path
    this%netrad_path = netrad_path
    this%eabs_path = eabs_path
    this%wind_path = wind_path
    this%meteo_mask_path = meteo_mask_path
    this%runoff_path = runoff_path
    this%runoff_sealed_path = runoff_sealed_path
    this%interflow_fast_path = interflow_fast_path
    this%interflow_slow_path = interflow_slow_path
    this%baseflow_path = baseflow_path
    this%hydro_mask_path = hydro_mask_path
    this%dem_path = dem_path
    this%slope_path = slope_path
    this%aspect_path = aspect_path
    this%fdir_path = fdir_path
    this%facc_path = facc_path
    this%geo_class_path = geo_class_path
    this%soil_class_path = soil_class_path
    this%soil_horizon_class_path = soil_horizon_class_path
    this%lai_class_path = lai_class_path
    this%river_width_path = river_width_path
    this%morph_mask_path = morph_mask_path
    this%pre_var = pre_var
    this%pet_var = pet_var
    this%temp_var = temp_var
    this%tann_var = tann_var
    this%tmin_var = tmin_var
    this%tmax_var = tmax_var
    this%ssrd_var = ssrd_var
    this%strd_var = strd_var
    this%netrad_var = netrad_var
    this%eabs_var = eabs_var
    this%wind_var = wind_var
    this%meteo_mask_var = meteo_mask_var
    this%runoff_var = runoff_var
    this%runoff_sealed_var = runoff_sealed_var
    this%interflow_fast_var = interflow_fast_var
    this%interflow_slow_var = interflow_slow_var
    this%baseflow_var = baseflow_var
    this%hydro_mask_var = hydro_mask_var
    this%dem_var = dem_var
    this%slope_var = slope_var
    this%aspect_var = aspect_var
    this%fdir_var = fdir_var
    this%facc_var = facc_var
    this%geo_class_var = geo_class_var
    this%soil_class_var = soil_class_var
    this%lai_class_var = lai_class_var
    this%river_width_var = river_width_var
    this%morph_mask_var = morph_mask_var
    this%hydro_lat_var = hydro_lat_var
    this%hydro_lon_var = hydro_lon_var
    this%morph_lat_var = morph_lat_var
    this%morph_lon_var = morph_lon_var
    this%route_lat_var = route_lat_var
    this%route_lon_var = route_lon_var

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_input_from_file

  !> \brief Set config_input values
  integer function nml_config_input_set(this, &
    chunking, &
    time_stamp_location, &
    latlon_path, &
    morph_latlon, &
    pre_path, &
    pet_path, &
    temp_path, &
    tann_path, &
    tmin_path, &
    tmax_path, &
    ssrd_path, &
    strd_path, &
    netrad_path, &
    eabs_path, &
    wind_path, &
    meteo_mask_path, &
    runoff_path, &
    runoff_sealed_path, &
    interflow_fast_path, &
    interflow_slow_path, &
    baseflow_path, &
    hydro_mask_path, &
    dem_path, &
    slope_path, &
    aspect_path, &
    fdir_path, &
    facc_path, &
    geo_class_path, &
    soil_class_path, &
    soil_horizon_class_path, &
    lai_class_path, &
    river_width_path, &
    morph_mask_path, &
    pre_var, &
    pet_var, &
    temp_var, &
    tann_var, &
    tmin_var, &
    tmax_var, &
    ssrd_var, &
    strd_var, &
    netrad_var, &
    eabs_var, &
    wind_var, &
    meteo_mask_var, &
    runoff_var, &
    runoff_sealed_var, &
    interflow_fast_var, &
    interflow_slow_var, &
    baseflow_var, &
    hydro_mask_var, &
    dem_var, &
    slope_var, &
    aspect_var, &
    fdir_var, &
    facc_var, &
    geo_class_var, &
    soil_class_var, &
    lai_class_var, &
    river_width_var, &
    morph_mask_var, &
    hydro_lat_var, &
    hydro_lon_var, &
    morph_lat_var, &
    morph_lon_var, &
    route_lat_var, &
    route_lon_var, &
    errmsg) result(status)

    class(nml_config_input_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), dimension(:), intent(in), optional :: chunking !< Chunking for input data
    integer(i4), dimension(:), intent(in), optional :: time_stamp_location !< NetCDF time-stamp location
    character(len=*), dimension(:), intent(in), optional :: latlon_path !< Latlon specification file path
    logical, dimension(:), intent(in), optional :: morph_latlon !< DEM in latlon coordinates
    character(len=*), dimension(:), intent(in), optional :: pre_path !< Precipitation input
    character(len=*), dimension(:), intent(in), optional :: pet_path !< Potential evapotranspiration input
    character(len=*), dimension(:), intent(in), optional :: temp_path !< Air temperature input
    character(len=*), dimension(:), intent(in), optional :: tann_path !< Air temperature annual mean input
    character(len=*), dimension(:), intent(in), optional :: tmin_path !< Air temperature daily minimum input
    character(len=*), dimension(:), intent(in), optional :: tmax_path !< Air temperature daily maximum input
    character(len=*), dimension(:), intent(in), optional :: ssrd_path !< Surface shortwave radiation downwards input
    character(len=*), dimension(:), intent(in), optional :: strd_path !< Surface thermal radiation downwards input
    character(len=*), dimension(:), intent(in), optional :: netrad_path !< Net radiation input
    character(len=*), dimension(:), intent(in), optional :: eabs_path !< Vapor pressure input
    character(len=*), dimension(:), intent(in), optional :: wind_path !< Wind speed input
    character(len=*), dimension(:), intent(in), optional :: meteo_mask_path !< Meteorological mask file path
    character(len=*), dimension(:), intent(in), optional :: runoff_path !< Runoff input
    character(len=*), dimension(:), intent(in), optional :: runoff_sealed_path !< Sealed runoff input
    character(len=*), dimension(:), intent(in), optional :: interflow_fast_path !< Fast interflow input
    character(len=*), dimension(:), intent(in), optional :: interflow_slow_path !< Slow interflow input
    character(len=*), dimension(:), intent(in), optional :: baseflow_path !< Baseflow input
    character(len=*), dimension(:), intent(in), optional :: hydro_mask_path !< Hydrological mask file path
    character(len=*), dimension(:), intent(in), optional :: dem_path !< DEM input
    character(len=*), dimension(:), intent(in), optional :: slope_path !< Slope input
    character(len=*), dimension(:), intent(in), optional :: aspect_path !< Aspect input
    character(len=*), dimension(:), intent(in), optional :: fdir_path !< Flow direction input
    character(len=*), dimension(:), intent(in), optional :: facc_path !< Flow accumulation input
    character(len=*), dimension(:), intent(in), optional :: geo_class_path !< Geology class input
    character(len=*), dimension(:), intent(in), optional :: soil_class_path !< Soil class input
    character(len=*), dimension(:), intent(in), optional :: soil_horizon_class_path !< Soil horizon class input
    character(len=*), dimension(:), intent(in), optional :: lai_class_path !< LAI class input
    character(len=*), dimension(:), intent(in), optional :: river_width_path !< River width input
    character(len=*), dimension(:), intent(in), optional :: morph_mask_path !< Morphology mask file path
    character(len=*), dimension(:), intent(in), optional :: pre_var !< Precipitation variable name
    character(len=*), dimension(:), intent(in), optional :: pet_var !< Potential evapotranspiration variable name
    character(len=*), dimension(:), intent(in), optional :: temp_var !< Air temperature variable name
    character(len=*), dimension(:), intent(in), optional :: tann_var !< Air temperature annual mean variable name
    character(len=*), dimension(:), intent(in), optional :: tmin_var !< Air temperature daily minimum variable name
    character(len=*), dimension(:), intent(in), optional :: tmax_var !< Air temperature daily maximum variable name
    character(len=*), dimension(:), intent(in), optional :: ssrd_var !< Surface shortwave radiation variable name
    character(len=*), dimension(:), intent(in), optional :: strd_var !< Surface thermal radiation variable name
    character(len=*), dimension(:), intent(in), optional :: netrad_var !< Net radiation variable name
    character(len=*), dimension(:), intent(in), optional :: eabs_var !< Vapor pressure variable name
    character(len=*), dimension(:), intent(in), optional :: wind_var !< Wind speed variable name
    character(len=*), dimension(:), intent(in), optional :: meteo_mask_var !< Meteorological mask variable name
    character(len=*), dimension(:), intent(in), optional :: runoff_var !< Runoff variable name
    character(len=*), dimension(:), intent(in), optional :: runoff_sealed_var !< Sealed runoff variable name
    character(len=*), dimension(:), intent(in), optional :: interflow_fast_var !< Fast interflow variable name
    character(len=*), dimension(:), intent(in), optional :: interflow_slow_var !< Slow interflow variable name
    character(len=*), dimension(:), intent(in), optional :: baseflow_var !< Baseflow variable name
    character(len=*), dimension(:), intent(in), optional :: hydro_mask_var !< Hydrological mask variable name
    character(len=*), dimension(:), intent(in), optional :: dem_var !< DEM variable name
    character(len=*), dimension(:), intent(in), optional :: slope_var !< Slope variable name
    character(len=*), dimension(:), intent(in), optional :: aspect_var !< Aspect variable name
    character(len=*), dimension(:), intent(in), optional :: fdir_var !< Flow direction variable name
    character(len=*), dimension(:), intent(in), optional :: facc_var !< Flow accumulation variable name
    character(len=*), dimension(:), intent(in), optional :: geo_class_var !< Geology class variable name
    character(len=*), dimension(:), intent(in), optional :: soil_class_var !< Soil class variable name
    character(len=*), dimension(:), intent(in), optional :: lai_class_var !< LAI class variable name
    character(len=*), dimension(:), intent(in), optional :: river_width_var !< River width variable name
    character(len=*), dimension(:), intent(in), optional :: morph_mask_var !< Morphology mask variable name
    character(len=*), dimension(:), intent(in), optional :: hydro_lat_var !< Hydrological latitude variable name
    character(len=*), dimension(:), intent(in), optional :: hydro_lon_var !< Hydrological longitude variable name
    character(len=*), dimension(:), intent(in), optional :: morph_lat_var !< Morphology latitude variable name
    character(len=*), dimension(:), intent(in), optional :: morph_lon_var !< Morphology longitude variable name
    character(len=*), dimension(:), intent(in), optional :: route_lat_var !< Routing latitude variable name
    character(len=*), dimension(:), intent(in), optional :: route_lon_var !< Routing longitude variable name
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(chunking)) then
      if (size(chunking, 1) > size(this%chunking, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'chunking'"
        return
      end if
      lb__1 = lbound(this%chunking, 1)
      ub__1 = lb__1 + size(chunking, 1) - 1
      this%chunking(lb__1:ub__1) = chunking
    end if
    if (present(time_stamp_location)) then
      if (size(time_stamp_location, 1) > size(this%time_stamp_location, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'time_stamp_location'"
        return
      end if
      lb__1 = lbound(this%time_stamp_location, 1)
      ub__1 = lb__1 + size(time_stamp_location, 1) - 1
      this%time_stamp_location(lb__1:ub__1) = time_stamp_location
    end if
    if (present(latlon_path)) then
      if (size(latlon_path, 1) > size(this%latlon_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'latlon_path'"
        return
      end if
      lb__1 = lbound(this%latlon_path, 1)
      ub__1 = lb__1 + size(latlon_path, 1) - 1
      this%latlon_path(lb__1:ub__1) = latlon_path
    end if
    if (present(morph_latlon)) then
      if (size(morph_latlon, 1) > size(this%morph_latlon, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_latlon'"
        return
      end if
      lb__1 = lbound(this%morph_latlon, 1)
      ub__1 = lb__1 + size(morph_latlon, 1) - 1
      this%morph_latlon(lb__1:ub__1) = morph_latlon
    end if
    if (present(pre_path)) then
      if (size(pre_path, 1) > size(this%pre_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pre_path'"
        return
      end if
      lb__1 = lbound(this%pre_path, 1)
      ub__1 = lb__1 + size(pre_path, 1) - 1
      this%pre_path(lb__1:ub__1) = pre_path
    end if
    if (present(pet_path)) then
      if (size(pet_path, 1) > size(this%pet_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_path'"
        return
      end if
      lb__1 = lbound(this%pet_path, 1)
      ub__1 = lb__1 + size(pet_path, 1) - 1
      this%pet_path(lb__1:ub__1) = pet_path
    end if
    if (present(temp_path)) then
      if (size(temp_path, 1) > size(this%temp_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'temp_path'"
        return
      end if
      lb__1 = lbound(this%temp_path, 1)
      ub__1 = lb__1 + size(temp_path, 1) - 1
      this%temp_path(lb__1:ub__1) = temp_path
    end if
    if (present(tann_path)) then
      if (size(tann_path, 1) > size(this%tann_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tann_path'"
        return
      end if
      lb__1 = lbound(this%tann_path, 1)
      ub__1 = lb__1 + size(tann_path, 1) - 1
      this%tann_path(lb__1:ub__1) = tann_path
    end if
    if (present(tmin_path)) then
      if (size(tmin_path, 1) > size(this%tmin_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmin_path'"
        return
      end if
      lb__1 = lbound(this%tmin_path, 1)
      ub__1 = lb__1 + size(tmin_path, 1) - 1
      this%tmin_path(lb__1:ub__1) = tmin_path
    end if
    if (present(tmax_path)) then
      if (size(tmax_path, 1) > size(this%tmax_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmax_path'"
        return
      end if
      lb__1 = lbound(this%tmax_path, 1)
      ub__1 = lb__1 + size(tmax_path, 1) - 1
      this%tmax_path(lb__1:ub__1) = tmax_path
    end if
    if (present(ssrd_path)) then
      if (size(ssrd_path, 1) > size(this%ssrd_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ssrd_path'"
        return
      end if
      lb__1 = lbound(this%ssrd_path, 1)
      ub__1 = lb__1 + size(ssrd_path, 1) - 1
      this%ssrd_path(lb__1:ub__1) = ssrd_path
    end if
    if (present(strd_path)) then
      if (size(strd_path, 1) > size(this%strd_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'strd_path'"
        return
      end if
      lb__1 = lbound(this%strd_path, 1)
      ub__1 = lb__1 + size(strd_path, 1) - 1
      this%strd_path(lb__1:ub__1) = strd_path
    end if
    if (present(netrad_path)) then
      if (size(netrad_path, 1) > size(this%netrad_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'netrad_path'"
        return
      end if
      lb__1 = lbound(this%netrad_path, 1)
      ub__1 = lb__1 + size(netrad_path, 1) - 1
      this%netrad_path(lb__1:ub__1) = netrad_path
    end if
    if (present(eabs_path)) then
      if (size(eabs_path, 1) > size(this%eabs_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'eabs_path'"
        return
      end if
      lb__1 = lbound(this%eabs_path, 1)
      ub__1 = lb__1 + size(eabs_path, 1) - 1
      this%eabs_path(lb__1:ub__1) = eabs_path
    end if
    if (present(wind_path)) then
      if (size(wind_path, 1) > size(this%wind_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'wind_path'"
        return
      end if
      lb__1 = lbound(this%wind_path, 1)
      ub__1 = lb__1 + size(wind_path, 1) - 1
      this%wind_path(lb__1:ub__1) = wind_path
    end if
    if (present(meteo_mask_path)) then
      if (size(meteo_mask_path, 1) > size(this%meteo_mask_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_mask_path'"
        return
      end if
      lb__1 = lbound(this%meteo_mask_path, 1)
      ub__1 = lb__1 + size(meteo_mask_path, 1) - 1
      this%meteo_mask_path(lb__1:ub__1) = meteo_mask_path
    end if
    if (present(runoff_path)) then
      if (size(runoff_path, 1) > size(this%runoff_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_path'"
        return
      end if
      lb__1 = lbound(this%runoff_path, 1)
      ub__1 = lb__1 + size(runoff_path, 1) - 1
      this%runoff_path(lb__1:ub__1) = runoff_path
    end if
    if (present(runoff_sealed_path)) then
      if (size(runoff_sealed_path, 1) > size(this%runoff_sealed_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_sealed_path'"
        return
      end if
      lb__1 = lbound(this%runoff_sealed_path, 1)
      ub__1 = lb__1 + size(runoff_sealed_path, 1) - 1
      this%runoff_sealed_path(lb__1:ub__1) = runoff_sealed_path
    end if
    if (present(interflow_fast_path)) then
      if (size(interflow_fast_path, 1) > size(this%interflow_fast_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_fast_path'"
        return
      end if
      lb__1 = lbound(this%interflow_fast_path, 1)
      ub__1 = lb__1 + size(interflow_fast_path, 1) - 1
      this%interflow_fast_path(lb__1:ub__1) = interflow_fast_path
    end if
    if (present(interflow_slow_path)) then
      if (size(interflow_slow_path, 1) > size(this%interflow_slow_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_slow_path'"
        return
      end if
      lb__1 = lbound(this%interflow_slow_path, 1)
      ub__1 = lb__1 + size(interflow_slow_path, 1) - 1
      this%interflow_slow_path(lb__1:ub__1) = interflow_slow_path
    end if
    if (present(baseflow_path)) then
      if (size(baseflow_path, 1) > size(this%baseflow_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'baseflow_path'"
        return
      end if
      lb__1 = lbound(this%baseflow_path, 1)
      ub__1 = lb__1 + size(baseflow_path, 1) - 1
      this%baseflow_path(lb__1:ub__1) = baseflow_path
    end if
    if (present(hydro_mask_path)) then
      if (size(hydro_mask_path, 1) > size(this%hydro_mask_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_mask_path'"
        return
      end if
      lb__1 = lbound(this%hydro_mask_path, 1)
      ub__1 = lb__1 + size(hydro_mask_path, 1) - 1
      this%hydro_mask_path(lb__1:ub__1) = hydro_mask_path
    end if
    if (present(dem_path)) then
      if (size(dem_path, 1) > size(this%dem_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'dem_path'"
        return
      end if
      lb__1 = lbound(this%dem_path, 1)
      ub__1 = lb__1 + size(dem_path, 1) - 1
      this%dem_path(lb__1:ub__1) = dem_path
    end if
    if (present(slope_path)) then
      if (size(slope_path, 1) > size(this%slope_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'slope_path'"
        return
      end if
      lb__1 = lbound(this%slope_path, 1)
      ub__1 = lb__1 + size(slope_path, 1) - 1
      this%slope_path(lb__1:ub__1) = slope_path
    end if
    if (present(aspect_path)) then
      if (size(aspect_path, 1) > size(this%aspect_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'aspect_path'"
        return
      end if
      lb__1 = lbound(this%aspect_path, 1)
      ub__1 = lb__1 + size(aspect_path, 1) - 1
      this%aspect_path(lb__1:ub__1) = aspect_path
    end if
    if (present(fdir_path)) then
      if (size(fdir_path, 1) > size(this%fdir_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'fdir_path'"
        return
      end if
      lb__1 = lbound(this%fdir_path, 1)
      ub__1 = lb__1 + size(fdir_path, 1) - 1
      this%fdir_path(lb__1:ub__1) = fdir_path
    end if
    if (present(facc_path)) then
      if (size(facc_path, 1) > size(this%facc_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'facc_path'"
        return
      end if
      lb__1 = lbound(this%facc_path, 1)
      ub__1 = lb__1 + size(facc_path, 1) - 1
      this%facc_path(lb__1:ub__1) = facc_path
    end if
    if (present(geo_class_path)) then
      if (size(geo_class_path, 1) > size(this%geo_class_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'geo_class_path'"
        return
      end if
      lb__1 = lbound(this%geo_class_path, 1)
      ub__1 = lb__1 + size(geo_class_path, 1) - 1
      this%geo_class_path(lb__1:ub__1) = geo_class_path
    end if
    if (present(soil_class_path)) then
      if (size(soil_class_path, 1) > size(this%soil_class_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_class_path'"
        return
      end if
      lb__1 = lbound(this%soil_class_path, 1)
      ub__1 = lb__1 + size(soil_class_path, 1) - 1
      this%soil_class_path(lb__1:ub__1) = soil_class_path
    end if
    if (present(soil_horizon_class_path)) then
      if (size(soil_horizon_class_path, 1) > size(this%soil_horizon_class_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_horizon_class_path'"
        return
      end if
      lb__1 = lbound(this%soil_horizon_class_path, 1)
      ub__1 = lb__1 + size(soil_horizon_class_path, 1) - 1
      this%soil_horizon_class_path(lb__1:ub__1) = soil_horizon_class_path
    end if
    if (present(lai_class_path)) then
      if (size(lai_class_path, 1) > size(this%lai_class_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_class_path'"
        return
      end if
      lb__1 = lbound(this%lai_class_path, 1)
      ub__1 = lb__1 + size(lai_class_path, 1) - 1
      this%lai_class_path(lb__1:ub__1) = lai_class_path
    end if
    if (present(river_width_path)) then
      if (size(river_width_path, 1) > size(this%river_width_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'river_width_path'"
        return
      end if
      lb__1 = lbound(this%river_width_path, 1)
      ub__1 = lb__1 + size(river_width_path, 1) - 1
      this%river_width_path(lb__1:ub__1) = river_width_path
    end if
    if (present(morph_mask_path)) then
      if (size(morph_mask_path, 1) > size(this%morph_mask_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_mask_path'"
        return
      end if
      lb__1 = lbound(this%morph_mask_path, 1)
      ub__1 = lb__1 + size(morph_mask_path, 1) - 1
      this%morph_mask_path(lb__1:ub__1) = morph_mask_path
    end if
    if (present(pre_var)) then
      if (size(pre_var, 1) > size(this%pre_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pre_var'"
        return
      end if
      lb__1 = lbound(this%pre_var, 1)
      ub__1 = lb__1 + size(pre_var, 1) - 1
      this%pre_var(lb__1:ub__1) = pre_var
    end if
    if (present(pet_var)) then
      if (size(pet_var, 1) > size(this%pet_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_var'"
        return
      end if
      lb__1 = lbound(this%pet_var, 1)
      ub__1 = lb__1 + size(pet_var, 1) - 1
      this%pet_var(lb__1:ub__1) = pet_var
    end if
    if (present(temp_var)) then
      if (size(temp_var, 1) > size(this%temp_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'temp_var'"
        return
      end if
      lb__1 = lbound(this%temp_var, 1)
      ub__1 = lb__1 + size(temp_var, 1) - 1
      this%temp_var(lb__1:ub__1) = temp_var
    end if
    if (present(tann_var)) then
      if (size(tann_var, 1) > size(this%tann_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tann_var'"
        return
      end if
      lb__1 = lbound(this%tann_var, 1)
      ub__1 = lb__1 + size(tann_var, 1) - 1
      this%tann_var(lb__1:ub__1) = tann_var
    end if
    if (present(tmin_var)) then
      if (size(tmin_var, 1) > size(this%tmin_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmin_var'"
        return
      end if
      lb__1 = lbound(this%tmin_var, 1)
      ub__1 = lb__1 + size(tmin_var, 1) - 1
      this%tmin_var(lb__1:ub__1) = tmin_var
    end if
    if (present(tmax_var)) then
      if (size(tmax_var, 1) > size(this%tmax_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmax_var'"
        return
      end if
      lb__1 = lbound(this%tmax_var, 1)
      ub__1 = lb__1 + size(tmax_var, 1) - 1
      this%tmax_var(lb__1:ub__1) = tmax_var
    end if
    if (present(ssrd_var)) then
      if (size(ssrd_var, 1) > size(this%ssrd_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ssrd_var'"
        return
      end if
      lb__1 = lbound(this%ssrd_var, 1)
      ub__1 = lb__1 + size(ssrd_var, 1) - 1
      this%ssrd_var(lb__1:ub__1) = ssrd_var
    end if
    if (present(strd_var)) then
      if (size(strd_var, 1) > size(this%strd_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'strd_var'"
        return
      end if
      lb__1 = lbound(this%strd_var, 1)
      ub__1 = lb__1 + size(strd_var, 1) - 1
      this%strd_var(lb__1:ub__1) = strd_var
    end if
    if (present(netrad_var)) then
      if (size(netrad_var, 1) > size(this%netrad_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'netrad_var'"
        return
      end if
      lb__1 = lbound(this%netrad_var, 1)
      ub__1 = lb__1 + size(netrad_var, 1) - 1
      this%netrad_var(lb__1:ub__1) = netrad_var
    end if
    if (present(eabs_var)) then
      if (size(eabs_var, 1) > size(this%eabs_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'eabs_var'"
        return
      end if
      lb__1 = lbound(this%eabs_var, 1)
      ub__1 = lb__1 + size(eabs_var, 1) - 1
      this%eabs_var(lb__1:ub__1) = eabs_var
    end if
    if (present(wind_var)) then
      if (size(wind_var, 1) > size(this%wind_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'wind_var'"
        return
      end if
      lb__1 = lbound(this%wind_var, 1)
      ub__1 = lb__1 + size(wind_var, 1) - 1
      this%wind_var(lb__1:ub__1) = wind_var
    end if
    if (present(meteo_mask_var)) then
      if (size(meteo_mask_var, 1) > size(this%meteo_mask_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_mask_var'"
        return
      end if
      lb__1 = lbound(this%meteo_mask_var, 1)
      ub__1 = lb__1 + size(meteo_mask_var, 1) - 1
      this%meteo_mask_var(lb__1:ub__1) = meteo_mask_var
    end if
    if (present(runoff_var)) then
      if (size(runoff_var, 1) > size(this%runoff_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_var'"
        return
      end if
      lb__1 = lbound(this%runoff_var, 1)
      ub__1 = lb__1 + size(runoff_var, 1) - 1
      this%runoff_var(lb__1:ub__1) = runoff_var
    end if
    if (present(runoff_sealed_var)) then
      if (size(runoff_sealed_var, 1) > size(this%runoff_sealed_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_sealed_var'"
        return
      end if
      lb__1 = lbound(this%runoff_sealed_var, 1)
      ub__1 = lb__1 + size(runoff_sealed_var, 1) - 1
      this%runoff_sealed_var(lb__1:ub__1) = runoff_sealed_var
    end if
    if (present(interflow_fast_var)) then
      if (size(interflow_fast_var, 1) > size(this%interflow_fast_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_fast_var'"
        return
      end if
      lb__1 = lbound(this%interflow_fast_var, 1)
      ub__1 = lb__1 + size(interflow_fast_var, 1) - 1
      this%interflow_fast_var(lb__1:ub__1) = interflow_fast_var
    end if
    if (present(interflow_slow_var)) then
      if (size(interflow_slow_var, 1) > size(this%interflow_slow_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_slow_var'"
        return
      end if
      lb__1 = lbound(this%interflow_slow_var, 1)
      ub__1 = lb__1 + size(interflow_slow_var, 1) - 1
      this%interflow_slow_var(lb__1:ub__1) = interflow_slow_var
    end if
    if (present(baseflow_var)) then
      if (size(baseflow_var, 1) > size(this%baseflow_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'baseflow_var'"
        return
      end if
      lb__1 = lbound(this%baseflow_var, 1)
      ub__1 = lb__1 + size(baseflow_var, 1) - 1
      this%baseflow_var(lb__1:ub__1) = baseflow_var
    end if
    if (present(hydro_mask_var)) then
      if (size(hydro_mask_var, 1) > size(this%hydro_mask_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_mask_var'"
        return
      end if
      lb__1 = lbound(this%hydro_mask_var, 1)
      ub__1 = lb__1 + size(hydro_mask_var, 1) - 1
      this%hydro_mask_var(lb__1:ub__1) = hydro_mask_var
    end if
    if (present(dem_var)) then
      if (size(dem_var, 1) > size(this%dem_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'dem_var'"
        return
      end if
      lb__1 = lbound(this%dem_var, 1)
      ub__1 = lb__1 + size(dem_var, 1) - 1
      this%dem_var(lb__1:ub__1) = dem_var
    end if
    if (present(slope_var)) then
      if (size(slope_var, 1) > size(this%slope_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'slope_var'"
        return
      end if
      lb__1 = lbound(this%slope_var, 1)
      ub__1 = lb__1 + size(slope_var, 1) - 1
      this%slope_var(lb__1:ub__1) = slope_var
    end if
    if (present(aspect_var)) then
      if (size(aspect_var, 1) > size(this%aspect_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'aspect_var'"
        return
      end if
      lb__1 = lbound(this%aspect_var, 1)
      ub__1 = lb__1 + size(aspect_var, 1) - 1
      this%aspect_var(lb__1:ub__1) = aspect_var
    end if
    if (present(fdir_var)) then
      if (size(fdir_var, 1) > size(this%fdir_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'fdir_var'"
        return
      end if
      lb__1 = lbound(this%fdir_var, 1)
      ub__1 = lb__1 + size(fdir_var, 1) - 1
      this%fdir_var(lb__1:ub__1) = fdir_var
    end if
    if (present(facc_var)) then
      if (size(facc_var, 1) > size(this%facc_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'facc_var'"
        return
      end if
      lb__1 = lbound(this%facc_var, 1)
      ub__1 = lb__1 + size(facc_var, 1) - 1
      this%facc_var(lb__1:ub__1) = facc_var
    end if
    if (present(geo_class_var)) then
      if (size(geo_class_var, 1) > size(this%geo_class_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'geo_class_var'"
        return
      end if
      lb__1 = lbound(this%geo_class_var, 1)
      ub__1 = lb__1 + size(geo_class_var, 1) - 1
      this%geo_class_var(lb__1:ub__1) = geo_class_var
    end if
    if (present(soil_class_var)) then
      if (size(soil_class_var, 1) > size(this%soil_class_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_class_var'"
        return
      end if
      lb__1 = lbound(this%soil_class_var, 1)
      ub__1 = lb__1 + size(soil_class_var, 1) - 1
      this%soil_class_var(lb__1:ub__1) = soil_class_var
    end if
    if (present(lai_class_var)) then
      if (size(lai_class_var, 1) > size(this%lai_class_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_class_var'"
        return
      end if
      lb__1 = lbound(this%lai_class_var, 1)
      ub__1 = lb__1 + size(lai_class_var, 1) - 1
      this%lai_class_var(lb__1:ub__1) = lai_class_var
    end if
    if (present(river_width_var)) then
      if (size(river_width_var, 1) > size(this%river_width_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'river_width_var'"
        return
      end if
      lb__1 = lbound(this%river_width_var, 1)
      ub__1 = lb__1 + size(river_width_var, 1) - 1
      this%river_width_var(lb__1:ub__1) = river_width_var
    end if
    if (present(morph_mask_var)) then
      if (size(morph_mask_var, 1) > size(this%morph_mask_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_mask_var'"
        return
      end if
      lb__1 = lbound(this%morph_mask_var, 1)
      ub__1 = lb__1 + size(morph_mask_var, 1) - 1
      this%morph_mask_var(lb__1:ub__1) = morph_mask_var
    end if
    if (present(hydro_lat_var)) then
      if (size(hydro_lat_var, 1) > size(this%hydro_lat_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_lat_var'"
        return
      end if
      lb__1 = lbound(this%hydro_lat_var, 1)
      ub__1 = lb__1 + size(hydro_lat_var, 1) - 1
      this%hydro_lat_var(lb__1:ub__1) = hydro_lat_var
    end if
    if (present(hydro_lon_var)) then
      if (size(hydro_lon_var, 1) > size(this%hydro_lon_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_lon_var'"
        return
      end if
      lb__1 = lbound(this%hydro_lon_var, 1)
      ub__1 = lb__1 + size(hydro_lon_var, 1) - 1
      this%hydro_lon_var(lb__1:ub__1) = hydro_lon_var
    end if
    if (present(morph_lat_var)) then
      if (size(morph_lat_var, 1) > size(this%morph_lat_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_lat_var'"
        return
      end if
      lb__1 = lbound(this%morph_lat_var, 1)
      ub__1 = lb__1 + size(morph_lat_var, 1) - 1
      this%morph_lat_var(lb__1:ub__1) = morph_lat_var
    end if
    if (present(morph_lon_var)) then
      if (size(morph_lon_var, 1) > size(this%morph_lon_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_lon_var'"
        return
      end if
      lb__1 = lbound(this%morph_lon_var, 1)
      ub__1 = lb__1 + size(morph_lon_var, 1) - 1
      this%morph_lon_var(lb__1:ub__1) = morph_lon_var
    end if
    if (present(route_lat_var)) then
      if (size(route_lat_var, 1) > size(this%route_lat_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'route_lat_var'"
        return
      end if
      lb__1 = lbound(this%route_lat_var, 1)
      ub__1 = lb__1 + size(route_lat_var, 1) - 1
      this%route_lat_var(lb__1:ub__1) = route_lat_var
    end if
    if (present(route_lon_var)) then
      if (size(route_lon_var, 1) > size(this%route_lon_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'route_lon_var'"
        return
      end if
      lb__1 = lbound(this%route_lon_var, 1)
      ub__1 = lb__1 + size(route_lon_var, 1) - 1
      this%route_lon_var(lb__1:ub__1) = route_lon_var
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_input_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_input_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_input_t), intent(in) :: this !< namelist instance
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
    case ("chunking")
      if (.not. allocated(this%chunking)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%chunking), ubound(this%chunking), &
          "chunking", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("time_stamp_location")
      if (.not. allocated(this%time_stamp_location)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%time_stamp_location), ubound(this%time_stamp_location), &
          "time_stamp_location", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("latlon_path")
      if (.not. allocated(this%latlon_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%latlon_path), ubound(this%latlon_path), &
          "latlon_path", errmsg)
        if (status /= NML_OK) return
        if (this%latlon_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%latlon_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("morph_latlon")
      if (.not. allocated(this%morph_latlon)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_latlon), ubound(this%morph_latlon), &
          "morph_latlon", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pre_path")
      if (.not. allocated(this%pre_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pre_path), ubound(this%pre_path), &
          "pre_path", errmsg)
        if (status /= NML_OK) return
        if (this%pre_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%pre_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("pet_path")
      if (.not. allocated(this%pet_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_path), ubound(this%pet_path), &
          "pet_path", errmsg)
        if (status /= NML_OK) return
        if (this%pet_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%pet_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("temp_path")
      if (.not. allocated(this%temp_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%temp_path), ubound(this%temp_path), &
          "temp_path", errmsg)
        if (status /= NML_OK) return
        if (this%temp_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%temp_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("tann_path")
      if (.not. allocated(this%tann_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tann_path), ubound(this%tann_path), &
          "tann_path", errmsg)
        if (status /= NML_OK) return
        if (this%tann_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%tann_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("tmin_path")
      if (.not. allocated(this%tmin_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmin_path), ubound(this%tmin_path), &
          "tmin_path", errmsg)
        if (status /= NML_OK) return
        if (this%tmin_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%tmin_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("tmax_path")
      if (.not. allocated(this%tmax_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmax_path), ubound(this%tmax_path), &
          "tmax_path", errmsg)
        if (status /= NML_OK) return
        if (this%tmax_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%tmax_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("ssrd_path")
      if (.not. allocated(this%ssrd_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ssrd_path), ubound(this%ssrd_path), &
          "ssrd_path", errmsg)
        if (status /= NML_OK) return
        if (this%ssrd_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%ssrd_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("strd_path")
      if (.not. allocated(this%strd_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%strd_path), ubound(this%strd_path), &
          "strd_path", errmsg)
        if (status /= NML_OK) return
        if (this%strd_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%strd_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("netrad_path")
      if (.not. allocated(this%netrad_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%netrad_path), ubound(this%netrad_path), &
          "netrad_path", errmsg)
        if (status /= NML_OK) return
        if (this%netrad_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%netrad_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("eabs_path")
      if (.not. allocated(this%eabs_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%eabs_path), ubound(this%eabs_path), &
          "eabs_path", errmsg)
        if (status /= NML_OK) return
        if (this%eabs_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%eabs_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("wind_path")
      if (.not. allocated(this%wind_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%wind_path), ubound(this%wind_path), &
          "wind_path", errmsg)
        if (status /= NML_OK) return
        if (this%wind_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%wind_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_mask_path")
      if (.not. allocated(this%meteo_mask_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_mask_path), ubound(this%meteo_mask_path), &
          "meteo_mask_path", errmsg)
        if (status /= NML_OK) return
        if (this%meteo_mask_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%meteo_mask_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("runoff_path")
      if (.not. allocated(this%runoff_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_path), ubound(this%runoff_path), &
          "runoff_path", errmsg)
        if (status /= NML_OK) return
        if (this%runoff_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%runoff_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("runoff_sealed_path")
      if (.not. allocated(this%runoff_sealed_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_sealed_path), ubound(this%runoff_sealed_path), &
          "runoff_sealed_path", errmsg)
        if (status /= NML_OK) return
        if (this%runoff_sealed_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%runoff_sealed_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("interflow_fast_path")
      if (.not. allocated(this%interflow_fast_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_fast_path), ubound(this%interflow_fast_path), &
          "interflow_fast_path", errmsg)
        if (status /= NML_OK) return
        if (this%interflow_fast_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%interflow_fast_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("interflow_slow_path")
      if (.not. allocated(this%interflow_slow_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_slow_path), ubound(this%interflow_slow_path), &
          "interflow_slow_path", errmsg)
        if (status /= NML_OK) return
        if (this%interflow_slow_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%interflow_slow_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("baseflow_path")
      if (.not. allocated(this%baseflow_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_path), ubound(this%baseflow_path), &
          "baseflow_path", errmsg)
        if (status /= NML_OK) return
        if (this%baseflow_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%baseflow_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_mask_path")
      if (.not. allocated(this%hydro_mask_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_mask_path), ubound(this%hydro_mask_path), &
          "hydro_mask_path", errmsg)
        if (status /= NML_OK) return
        if (this%hydro_mask_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%hydro_mask_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("dem_path")
      if (.not. allocated(this%dem_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%dem_path), ubound(this%dem_path), &
          "dem_path", errmsg)
        if (status /= NML_OK) return
        if (this%dem_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%dem_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("slope_path")
      if (.not. allocated(this%slope_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%slope_path), ubound(this%slope_path), &
          "slope_path", errmsg)
        if (status /= NML_OK) return
        if (this%slope_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%slope_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("aspect_path")
      if (.not. allocated(this%aspect_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%aspect_path), ubound(this%aspect_path), &
          "aspect_path", errmsg)
        if (status /= NML_OK) return
        if (this%aspect_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%aspect_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("fdir_path")
      if (.not. allocated(this%fdir_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%fdir_path), ubound(this%fdir_path), &
          "fdir_path", errmsg)
        if (status /= NML_OK) return
        if (this%fdir_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%fdir_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("facc_path")
      if (.not. allocated(this%facc_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%facc_path), ubound(this%facc_path), &
          "facc_path", errmsg)
        if (status /= NML_OK) return
        if (this%facc_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%facc_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("geo_class_path")
      if (.not. allocated(this%geo_class_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%geo_class_path), ubound(this%geo_class_path), &
          "geo_class_path", errmsg)
        if (status /= NML_OK) return
        if (this%geo_class_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%geo_class_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("soil_class_path")
      if (.not. allocated(this%soil_class_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_class_path), ubound(this%soil_class_path), &
          "soil_class_path", errmsg)
        if (status /= NML_OK) return
        if (this%soil_class_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%soil_class_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("soil_horizon_class_path")
      if (.not. allocated(this%soil_horizon_class_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_horizon_class_path), ubound(this%soil_horizon_class_path), &
          "soil_horizon_class_path", errmsg)
        if (status /= NML_OK) return
        if (this%soil_horizon_class_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%soil_horizon_class_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("lai_class_path")
      if (.not. allocated(this%lai_class_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_class_path), ubound(this%lai_class_path), &
          "lai_class_path", errmsg)
        if (status /= NML_OK) return
        if (this%lai_class_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%lai_class_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("river_width_path")
      if (.not. allocated(this%river_width_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%river_width_path), ubound(this%river_width_path), &
          "river_width_path", errmsg)
        if (status /= NML_OK) return
        if (this%river_width_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%river_width_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("morph_mask_path")
      if (.not. allocated(this%morph_mask_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_mask_path), ubound(this%morph_mask_path), &
          "morph_mask_path", errmsg)
        if (status /= NML_OK) return
        if (this%morph_mask_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%morph_mask_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("pre_var")
      if (.not. allocated(this%pre_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pre_var), ubound(this%pre_var), &
          "pre_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pet_var")
      if (.not. allocated(this%pet_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_var), ubound(this%pet_var), &
          "pet_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("temp_var")
      if (.not. allocated(this%temp_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%temp_var), ubound(this%temp_var), &
          "temp_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tann_var")
      if (.not. allocated(this%tann_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tann_var), ubound(this%tann_var), &
          "tann_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tmin_var")
      if (.not. allocated(this%tmin_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmin_var), ubound(this%tmin_var), &
          "tmin_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tmax_var")
      if (.not. allocated(this%tmax_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmax_var), ubound(this%tmax_var), &
          "tmax_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("ssrd_var")
      if (.not. allocated(this%ssrd_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ssrd_var), ubound(this%ssrd_var), &
          "ssrd_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("strd_var")
      if (.not. allocated(this%strd_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%strd_var), ubound(this%strd_var), &
          "strd_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("netrad_var")
      if (.not. allocated(this%netrad_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%netrad_var), ubound(this%netrad_var), &
          "netrad_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("eabs_var")
      if (.not. allocated(this%eabs_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%eabs_var), ubound(this%eabs_var), &
          "eabs_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("wind_var")
      if (.not. allocated(this%wind_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%wind_var), ubound(this%wind_var), &
          "wind_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("meteo_mask_var")
      if (.not. allocated(this%meteo_mask_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_mask_var), ubound(this%meteo_mask_var), &
          "meteo_mask_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("runoff_var")
      if (.not. allocated(this%runoff_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_var), ubound(this%runoff_var), &
          "runoff_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("runoff_sealed_var")
      if (.not. allocated(this%runoff_sealed_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_sealed_var), ubound(this%runoff_sealed_var), &
          "runoff_sealed_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("interflow_fast_var")
      if (.not. allocated(this%interflow_fast_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_fast_var), ubound(this%interflow_fast_var), &
          "interflow_fast_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("interflow_slow_var")
      if (.not. allocated(this%interflow_slow_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_slow_var), ubound(this%interflow_slow_var), &
          "interflow_slow_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("baseflow_var")
      if (.not. allocated(this%baseflow_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_var), ubound(this%baseflow_var), &
          "baseflow_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_mask_var")
      if (.not. allocated(this%hydro_mask_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_mask_var), ubound(this%hydro_mask_var), &
          "hydro_mask_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("dem_var")
      if (.not. allocated(this%dem_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%dem_var), ubound(this%dem_var), &
          "dem_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("slope_var")
      if (.not. allocated(this%slope_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%slope_var), ubound(this%slope_var), &
          "slope_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("aspect_var")
      if (.not. allocated(this%aspect_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%aspect_var), ubound(this%aspect_var), &
          "aspect_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("fdir_var")
      if (.not. allocated(this%fdir_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%fdir_var), ubound(this%fdir_var), &
          "fdir_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("facc_var")
      if (.not. allocated(this%facc_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%facc_var), ubound(this%facc_var), &
          "facc_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("geo_class_var")
      if (.not. allocated(this%geo_class_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%geo_class_var), ubound(this%geo_class_var), &
          "geo_class_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("soil_class_var")
      if (.not. allocated(this%soil_class_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_class_var), ubound(this%soil_class_var), &
          "soil_class_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("lai_class_var")
      if (.not. allocated(this%lai_class_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_class_var), ubound(this%lai_class_var), &
          "lai_class_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("river_width_var")
      if (.not. allocated(this%river_width_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%river_width_var), ubound(this%river_width_var), &
          "river_width_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_mask_var")
      if (.not. allocated(this%morph_mask_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_mask_var), ubound(this%morph_mask_var), &
          "morph_mask_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_lat_var")
      if (.not. allocated(this%hydro_lat_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_lat_var), ubound(this%hydro_lat_var), &
          "hydro_lat_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_lon_var")
      if (.not. allocated(this%hydro_lon_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_lon_var), ubound(this%hydro_lon_var), &
          "hydro_lon_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_lat_var")
      if (.not. allocated(this%morph_lat_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_lat_var), ubound(this%morph_lat_var), &
          "morph_lat_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_lon_var")
      if (.not. allocated(this%morph_lon_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_lon_var), ubound(this%morph_lon_var), &
          "morph_lon_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("route_lat_var")
      if (.not. allocated(this%route_lat_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%route_lat_var), ubound(this%route_lat_var), &
          "route_lat_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("route_lon_var")
      if (.not. allocated(this%route_lon_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%route_lon_var), ubound(this%route_lon_var), &
          "route_lon_var", errmsg)
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
  end function nml_config_input_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_input_is_valid(this, errmsg) result(status)
    class(nml_config_input_t), intent(in) :: this !< namelist instance
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
    if (allocated(this%time_stamp_location)) then
    if (.not. all(time_stamp_location__in_enum(this%time_stamp_location, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: time_stamp_location"
      return
    end if
    end if
    ! bounds constraints
    if (allocated(this%chunking)) then
    if (.not. all(chunking__in_bounds(this%chunking, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: chunking"
      return
    end if
    end if
  end function nml_config_input_is_valid

end module nml_config_input
