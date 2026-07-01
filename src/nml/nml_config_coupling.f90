!> \file nml_config_coupling.f90
!> \copydoc nml_config_coupling

!> \brief Coupling configuration
!> \details Coupling grid and variable settings for meteorological, hydrological, and morphology components.
!! Arrays are indexed by domain (dimension 1).
!! Grid parameters define resolutions and coordinate systems for coupled runs.
!! Coupled flags indicate whether inputs are provided by an external model.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_coupling
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
    i4, &
    dp

  implicit none

  ! default values
  integer(i4), parameter, public :: meteo_grid_ydir__default = 0_i4
  integer(i4), parameter, public :: meteo_grid_coordsys__default = 0_i4
  integer(i4), parameter, public :: hydro_grid_ydir__default = 0_i4
  integer(i4), parameter, public :: hydro_grid_coordsys__default = 0_i4
  integer(i4), parameter, public :: morph_grid_ydir__default = 0_i4
  integer(i4), parameter, public :: morph_grid_coordsys__default = 0_i4
  logical, parameter, public :: pre_coupled__default = .false.
  logical, parameter, public :: pet_coupled__default = .false.
  logical, parameter, public :: temp_coupled__default = .false.
  logical, parameter, public :: tann_coupled__default = .false.
  logical, parameter, public :: tmin_coupled__default = .false.
  logical, parameter, public :: tmax_coupled__default = .false.
  logical, parameter, public :: ssrd_coupled__default = .false.
  logical, parameter, public :: strd_coupled__default = .false.
  logical, parameter, public :: netrad_coupled__default = .false.
  logical, parameter, public :: eabs_coupled__default = .false.
  logical, parameter, public :: wind_coupled__default = .false.
  logical, parameter, public :: runoff_coupled__default = .false.
  logical, parameter, public :: runoff_sealed_coupled__default = .false.
  logical, parameter, public :: interflow_fast_coupled__default = .false.
  logical, parameter, public :: interflow_slow_coupled__default = .false.
  logical, parameter, public :: baseflow_coupled__default = .false.
  logical, parameter, public :: dem_coupled__default = .false.
  logical, parameter, public :: slope_coupled__default = .false.
  logical, parameter, public :: aspect_coupled__default = .false.
  logical, parameter, public :: geo_class_coupled__default = .false.
  logical, parameter, public :: soil_class_coupled__default = .false.
  logical, parameter, public :: lai_class_coupled__default = .false.
  logical, parameter, public :: river_width_coupled__default = .false.
  logical, parameter, public :: meteo_mask_coupled__default = .false.
  logical, parameter, public :: hydro_mask_coupled__default = .false.
  logical, parameter, public :: morph_mask_coupled__default = .false.
  logical, parameter, public :: hydro_latlon_coupled__default = .false.
  logical, parameter, public :: morph_latlon_coupled__default = .false.
  logical, parameter, public :: route_latlon_coupled__default = .false.

  ! enum values
  integer(i4), parameter, public :: meteo_grid_ydir__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: meteo_grid_coordsys__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: hydro_grid_ydir__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: hydro_grid_coordsys__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: morph_grid_ydir__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: morph_grid_coordsys__enum_values(2) = [0_i4, 1_i4]

  ! bounds values
  integer(i4), parameter, public :: meteo_grid_nx__min = 1_i4
  integer(i4), parameter, public :: meteo_grid_ny__min = 1_i4
  real(dp), parameter, public :: meteo_grid_cellsize__min_excl = 0.0_dp
  integer(i4), parameter, public :: hydro_grid_nx__min = 1_i4
  integer(i4), parameter, public :: hydro_grid_ny__min = 1_i4
  real(dp), parameter, public :: hydro_grid_cellsize__min_excl = 0.0_dp
  integer(i4), parameter, public :: morph_grid_nx__min = 1_i4
  integer(i4), parameter, public :: morph_grid_ny__min = 1_i4
  real(dp), parameter, public :: morph_grid_cellsize__min_excl = 0.0_dp

  !> \class nml_config_coupling_t
  !> \brief Coupling configuration
  !> \details Coupling grid and variable settings for meteorological, hydrological, and morphology components.
  !! Arrays are indexed by domain (dimension 1).
  !! Grid parameters define resolutions and coordinate systems for coupled runs.
  !! Coupled flags indicate whether inputs are provided by an external model.
  type, public :: nml_config_coupling_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer(i4), allocatable, dimension(:) :: meteo_grid_nx !< Meteo grid size in x-direction
    integer(i4), allocatable, dimension(:) :: meteo_grid_ny !< Meteo grid size in y-direction
    real(dp), allocatable, dimension(:) :: meteo_grid_xll !< Meteo grid x origin
    real(dp), allocatable, dimension(:) :: meteo_grid_yll !< Meteo grid y origin
    real(dp), allocatable, dimension(:) :: meteo_grid_cellsize !< Meteo grid cell size
    integer(i4), allocatable, dimension(:) :: meteo_grid_ydir !< Meteo grid y direction
    integer(i4), allocatable, dimension(:) :: meteo_grid_coordsys !< Meteo grid coordinate system
    integer(i4), allocatable, dimension(:) :: hydro_grid_nx !< Hydro grid size in x-direction
    integer(i4), allocatable, dimension(:) :: hydro_grid_ny !< Hydro grid size in y-direction
    real(dp), allocatable, dimension(:) :: hydro_grid_xll !< Hydro grid x origin
    real(dp), allocatable, dimension(:) :: hydro_grid_yll !< Hydro grid y origin
    real(dp), allocatable, dimension(:) :: hydro_grid_cellsize !< Hydro grid cell size
    integer(i4), allocatable, dimension(:) :: hydro_grid_ydir !< Hydro grid y direction
    integer(i4), allocatable, dimension(:) :: hydro_grid_coordsys !< Hydro grid coordinate system
    integer(i4), allocatable, dimension(:) :: morph_grid_nx !< Morph grid size in x-direction
    integer(i4), allocatable, dimension(:) :: morph_grid_ny !< Morph grid size in y-direction
    real(dp), allocatable, dimension(:) :: morph_grid_xll !< Morph grid x origin
    real(dp), allocatable, dimension(:) :: morph_grid_yll !< Morph grid y origin
    real(dp), allocatable, dimension(:) :: morph_grid_cellsize !< Morph grid cell size
    integer(i4), allocatable, dimension(:) :: morph_grid_ydir !< Morph grid y direction
    integer(i4), allocatable, dimension(:) :: morph_grid_coordsys !< Morph grid coordinate system
    logical, allocatable, dimension(:) :: pre_coupled !< Precipitation coupled
    logical, allocatable, dimension(:) :: pet_coupled !< Potential evapotranspiration coupled
    logical, allocatable, dimension(:) :: temp_coupled !< Air temperature coupled
    logical, allocatable, dimension(:) :: tann_coupled !< Air temperature annual mean coupled
    logical, allocatable, dimension(:) :: tmin_coupled !< Air temperature daily minimum coupled
    logical, allocatable, dimension(:) :: tmax_coupled !< Air temperature daily maximum coupled
    logical, allocatable, dimension(:) :: ssrd_coupled !< Surface shortwave radiation coupled
    logical, allocatable, dimension(:) :: strd_coupled !< Surface thermal radiation coupled
    logical, allocatable, dimension(:) :: netrad_coupled !< Net radiation coupled
    logical, allocatable, dimension(:) :: eabs_coupled !< Vapor pressure coupled
    logical, allocatable, dimension(:) :: wind_coupled !< Wind speed coupled
    logical, allocatable, dimension(:) :: runoff_coupled !< Runoff coupled
    logical, allocatable, dimension(:) :: runoff_sealed_coupled !< Sealed runoff coupled
    logical, allocatable, dimension(:) :: interflow_fast_coupled !< Fast interflow coupled
    logical, allocatable, dimension(:) :: interflow_slow_coupled !< Slow interflow coupled
    logical, allocatable, dimension(:) :: baseflow_coupled !< Baseflow coupled
    logical, allocatable, dimension(:) :: dem_coupled !< DEM coupled
    logical, allocatable, dimension(:) :: slope_coupled !< Slope coupled
    logical, allocatable, dimension(:) :: aspect_coupled !< Aspect coupled
    logical, allocatable, dimension(:) :: geo_class_coupled !< Geology class coupled
    logical, allocatable, dimension(:) :: soil_class_coupled !< Soil class coupled
    logical, allocatable, dimension(:) :: lai_class_coupled !< LAI class coupled
    logical, allocatable, dimension(:) :: river_width_coupled !< River width coupled
    logical, allocatable, dimension(:) :: meteo_mask_coupled !< Meteorological mask coupled
    logical, allocatable, dimension(:) :: hydro_mask_coupled !< Hydrological mask coupled
    logical, allocatable, dimension(:) :: morph_mask_coupled !< Morphology mask coupled
    logical, allocatable, dimension(:) :: hydro_latlon_coupled !< Hydrological latlon coupled
    logical, allocatable, dimension(:) :: morph_latlon_coupled !< Morphological latlon coupled
    logical, allocatable, dimension(:) :: route_latlon_coupled !< Routing latlon coupled
  contains
    procedure :: init => nml_config_coupling_init
    procedure :: set_dims => nml_config_coupling_set_dims
    procedure :: from_file => nml_config_coupling_from_file
    procedure :: set => nml_config_coupling_set
    procedure :: is_set => nml_config_coupling_is_set
    procedure :: is_valid => nml_config_coupling_is_valid
  end type nml_config_coupling_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function meteo_grid_ydir__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == meteo_grid_ydir__enum_values)
  end function meteo_grid_ydir__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function meteo_grid_coordsys__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == meteo_grid_coordsys__enum_values)
  end function meteo_grid_coordsys__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function hydro_grid_ydir__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == hydro_grid_ydir__enum_values)
  end function hydro_grid_ydir__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function hydro_grid_coordsys__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == hydro_grid_coordsys__enum_values)
  end function hydro_grid_coordsys__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function morph_grid_ydir__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == morph_grid_ydir__enum_values)
  end function morph_grid_ydir__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function morph_grid_coordsys__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == morph_grid_coordsys__enum_values)
  end function morph_grid_coordsys__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function meteo_grid_nx__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < meteo_grid_nx__min) in_bounds = .false.
  end function meteo_grid_nx__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function meteo_grid_ny__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < meteo_grid_ny__min) in_bounds = .false.
  end function meteo_grid_ny__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function meteo_grid_cellsize__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val <= meteo_grid_cellsize__min_excl) in_bounds = .false.
  end function meteo_grid_cellsize__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function hydro_grid_nx__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < hydro_grid_nx__min) in_bounds = .false.
  end function hydro_grid_nx__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function hydro_grid_ny__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < hydro_grid_ny__min) in_bounds = .false.
  end function hydro_grid_ny__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function hydro_grid_cellsize__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val <= hydro_grid_cellsize__min_excl) in_bounds = .false.
  end function hydro_grid_cellsize__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function morph_grid_nx__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < morph_grid_nx__min) in_bounds = .false.
  end function morph_grid_nx__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function morph_grid_ny__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < morph_grid_ny__min) in_bounds = .false.
  end function morph_grid_ny__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function morph_grid_cellsize__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val <= morph_grid_cellsize__min_excl) in_bounds = .false.
  end function morph_grid_cellsize__in_bounds

  !> \brief Initialize defaults and sentinels for config_coupling
  integer function nml_config_coupling_init(this, errmsg) result(status)
    class(nml_config_coupling_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%meteo_grid_nx)) deallocate(this%meteo_grid_nx)
    allocate(this%meteo_grid_nx(this%n_domains))
    if (allocated(this%meteo_grid_ny)) deallocate(this%meteo_grid_ny)
    allocate(this%meteo_grid_ny(this%n_domains))
    if (allocated(this%meteo_grid_xll)) deallocate(this%meteo_grid_xll)
    allocate(this%meteo_grid_xll(this%n_domains))
    if (allocated(this%meteo_grid_yll)) deallocate(this%meteo_grid_yll)
    allocate(this%meteo_grid_yll(this%n_domains))
    if (allocated(this%meteo_grid_cellsize)) deallocate(this%meteo_grid_cellsize)
    allocate(this%meteo_grid_cellsize(this%n_domains))
    if (allocated(this%meteo_grid_ydir)) deallocate(this%meteo_grid_ydir)
    allocate(this%meteo_grid_ydir(this%n_domains))
    if (allocated(this%meteo_grid_coordsys)) deallocate(this%meteo_grid_coordsys)
    allocate(this%meteo_grid_coordsys(this%n_domains))
    if (allocated(this%hydro_grid_nx)) deallocate(this%hydro_grid_nx)
    allocate(this%hydro_grid_nx(this%n_domains))
    if (allocated(this%hydro_grid_ny)) deallocate(this%hydro_grid_ny)
    allocate(this%hydro_grid_ny(this%n_domains))
    if (allocated(this%hydro_grid_xll)) deallocate(this%hydro_grid_xll)
    allocate(this%hydro_grid_xll(this%n_domains))
    if (allocated(this%hydro_grid_yll)) deallocate(this%hydro_grid_yll)
    allocate(this%hydro_grid_yll(this%n_domains))
    if (allocated(this%hydro_grid_cellsize)) deallocate(this%hydro_grid_cellsize)
    allocate(this%hydro_grid_cellsize(this%n_domains))
    if (allocated(this%hydro_grid_ydir)) deallocate(this%hydro_grid_ydir)
    allocate(this%hydro_grid_ydir(this%n_domains))
    if (allocated(this%hydro_grid_coordsys)) deallocate(this%hydro_grid_coordsys)
    allocate(this%hydro_grid_coordsys(this%n_domains))
    if (allocated(this%morph_grid_nx)) deallocate(this%morph_grid_nx)
    allocate(this%morph_grid_nx(this%n_domains))
    if (allocated(this%morph_grid_ny)) deallocate(this%morph_grid_ny)
    allocate(this%morph_grid_ny(this%n_domains))
    if (allocated(this%morph_grid_xll)) deallocate(this%morph_grid_xll)
    allocate(this%morph_grid_xll(this%n_domains))
    if (allocated(this%morph_grid_yll)) deallocate(this%morph_grid_yll)
    allocate(this%morph_grid_yll(this%n_domains))
    if (allocated(this%morph_grid_cellsize)) deallocate(this%morph_grid_cellsize)
    allocate(this%morph_grid_cellsize(this%n_domains))
    if (allocated(this%morph_grid_ydir)) deallocate(this%morph_grid_ydir)
    allocate(this%morph_grid_ydir(this%n_domains))
    if (allocated(this%morph_grid_coordsys)) deallocate(this%morph_grid_coordsys)
    allocate(this%morph_grid_coordsys(this%n_domains))
    if (allocated(this%pre_coupled)) deallocate(this%pre_coupled)
    allocate(this%pre_coupled(this%n_domains))
    if (allocated(this%pet_coupled)) deallocate(this%pet_coupled)
    allocate(this%pet_coupled(this%n_domains))
    if (allocated(this%temp_coupled)) deallocate(this%temp_coupled)
    allocate(this%temp_coupled(this%n_domains))
    if (allocated(this%tann_coupled)) deallocate(this%tann_coupled)
    allocate(this%tann_coupled(this%n_domains))
    if (allocated(this%tmin_coupled)) deallocate(this%tmin_coupled)
    allocate(this%tmin_coupled(this%n_domains))
    if (allocated(this%tmax_coupled)) deallocate(this%tmax_coupled)
    allocate(this%tmax_coupled(this%n_domains))
    if (allocated(this%ssrd_coupled)) deallocate(this%ssrd_coupled)
    allocate(this%ssrd_coupled(this%n_domains))
    if (allocated(this%strd_coupled)) deallocate(this%strd_coupled)
    allocate(this%strd_coupled(this%n_domains))
    if (allocated(this%netrad_coupled)) deallocate(this%netrad_coupled)
    allocate(this%netrad_coupled(this%n_domains))
    if (allocated(this%eabs_coupled)) deallocate(this%eabs_coupled)
    allocate(this%eabs_coupled(this%n_domains))
    if (allocated(this%wind_coupled)) deallocate(this%wind_coupled)
    allocate(this%wind_coupled(this%n_domains))
    if (allocated(this%runoff_coupled)) deallocate(this%runoff_coupled)
    allocate(this%runoff_coupled(this%n_domains))
    if (allocated(this%runoff_sealed_coupled)) deallocate(this%runoff_sealed_coupled)
    allocate(this%runoff_sealed_coupled(this%n_domains))
    if (allocated(this%interflow_fast_coupled)) deallocate(this%interflow_fast_coupled)
    allocate(this%interflow_fast_coupled(this%n_domains))
    if (allocated(this%interflow_slow_coupled)) deallocate(this%interflow_slow_coupled)
    allocate(this%interflow_slow_coupled(this%n_domains))
    if (allocated(this%baseflow_coupled)) deallocate(this%baseflow_coupled)
    allocate(this%baseflow_coupled(this%n_domains))
    if (allocated(this%dem_coupled)) deallocate(this%dem_coupled)
    allocate(this%dem_coupled(this%n_domains))
    if (allocated(this%slope_coupled)) deallocate(this%slope_coupled)
    allocate(this%slope_coupled(this%n_domains))
    if (allocated(this%aspect_coupled)) deallocate(this%aspect_coupled)
    allocate(this%aspect_coupled(this%n_domains))
    if (allocated(this%geo_class_coupled)) deallocate(this%geo_class_coupled)
    allocate(this%geo_class_coupled(this%n_domains))
    if (allocated(this%soil_class_coupled)) deallocate(this%soil_class_coupled)
    allocate(this%soil_class_coupled(this%n_domains))
    if (allocated(this%lai_class_coupled)) deallocate(this%lai_class_coupled)
    allocate(this%lai_class_coupled(this%n_domains))
    if (allocated(this%river_width_coupled)) deallocate(this%river_width_coupled)
    allocate(this%river_width_coupled(this%n_domains))
    if (allocated(this%meteo_mask_coupled)) deallocate(this%meteo_mask_coupled)
    allocate(this%meteo_mask_coupled(this%n_domains))
    if (allocated(this%hydro_mask_coupled)) deallocate(this%hydro_mask_coupled)
    allocate(this%hydro_mask_coupled(this%n_domains))
    if (allocated(this%morph_mask_coupled)) deallocate(this%morph_mask_coupled)
    allocate(this%morph_mask_coupled(this%n_domains))
    if (allocated(this%hydro_latlon_coupled)) deallocate(this%hydro_latlon_coupled)
    allocate(this%hydro_latlon_coupled(this%n_domains))
    if (allocated(this%morph_latlon_coupled)) deallocate(this%morph_latlon_coupled)
    allocate(this%morph_latlon_coupled(this%n_domains))
    if (allocated(this%route_latlon_coupled)) deallocate(this%route_latlon_coupled)
    allocate(this%route_latlon_coupled(this%n_domains))

    ! sentinel values for required/optional parameters
    this%meteo_grid_nx = -huge(this%meteo_grid_nx) ! sentinel for optional integer array
    this%meteo_grid_ny = -huge(this%meteo_grid_ny) ! sentinel for optional integer array
    this%meteo_grid_xll = ieee_value(this%meteo_grid_xll, ieee_quiet_nan) ! sentinel for optional real array
    this%meteo_grid_yll = ieee_value(this%meteo_grid_yll, ieee_quiet_nan) ! sentinel for optional real array
    this%meteo_grid_cellsize = ieee_value(this%meteo_grid_cellsize, ieee_quiet_nan) ! sentinel for optional real array
    this%hydro_grid_nx = -huge(this%hydro_grid_nx) ! sentinel for optional integer array
    this%hydro_grid_ny = -huge(this%hydro_grid_ny) ! sentinel for optional integer array
    this%hydro_grid_xll = ieee_value(this%hydro_grid_xll, ieee_quiet_nan) ! sentinel for optional real array
    this%hydro_grid_yll = ieee_value(this%hydro_grid_yll, ieee_quiet_nan) ! sentinel for optional real array
    this%hydro_grid_cellsize = ieee_value(this%hydro_grid_cellsize, ieee_quiet_nan) ! sentinel for optional real array
    this%morph_grid_nx = -huge(this%morph_grid_nx) ! sentinel for optional integer array
    this%morph_grid_ny = -huge(this%morph_grid_ny) ! sentinel for optional integer array
    this%morph_grid_xll = ieee_value(this%morph_grid_xll, ieee_quiet_nan) ! sentinel for optional real array
    this%morph_grid_yll = ieee_value(this%morph_grid_yll, ieee_quiet_nan) ! sentinel for optional real array
    this%morph_grid_cellsize = ieee_value(this%morph_grid_cellsize, ieee_quiet_nan) ! sentinel for optional real array
    ! default values
    this%meteo_grid_ydir = meteo_grid_ydir__default
    this%meteo_grid_coordsys = meteo_grid_coordsys__default
    this%hydro_grid_ydir = hydro_grid_ydir__default
    this%hydro_grid_coordsys = hydro_grid_coordsys__default
    this%morph_grid_ydir = morph_grid_ydir__default
    this%morph_grid_coordsys = morph_grid_coordsys__default
    this%pre_coupled = pre_coupled__default
    this%pet_coupled = pet_coupled__default
    this%temp_coupled = temp_coupled__default
    this%tann_coupled = tann_coupled__default
    this%tmin_coupled = tmin_coupled__default
    this%tmax_coupled = tmax_coupled__default
    this%ssrd_coupled = ssrd_coupled__default
    this%strd_coupled = strd_coupled__default
    this%netrad_coupled = netrad_coupled__default
    this%eabs_coupled = eabs_coupled__default
    this%wind_coupled = wind_coupled__default
    this%runoff_coupled = runoff_coupled__default
    this%runoff_sealed_coupled = runoff_sealed_coupled__default
    this%interflow_fast_coupled = interflow_fast_coupled__default
    this%interflow_slow_coupled = interflow_slow_coupled__default
    this%baseflow_coupled = baseflow_coupled__default
    this%dem_coupled = dem_coupled__default
    this%slope_coupled = slope_coupled__default
    this%aspect_coupled = aspect_coupled__default
    this%geo_class_coupled = geo_class_coupled__default
    this%soil_class_coupled = soil_class_coupled__default
    this%lai_class_coupled = lai_class_coupled__default
    this%river_width_coupled = river_width_coupled__default
    this%meteo_mask_coupled = meteo_mask_coupled__default
    this%hydro_mask_coupled = hydro_mask_coupled__default
    this%morph_mask_coupled = morph_mask_coupled__default
    this%hydro_latlon_coupled = hydro_latlon_coupled__default
    this%morph_latlon_coupled = morph_latlon_coupled__default
    this%route_latlon_coupled = route_latlon_coupled__default
  end function nml_config_coupling_init

  !> \brief Reset runtime dimensions for config_coupling
  integer function nml_config_coupling_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_coupling_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%meteo_grid_nx)) deallocate(this%meteo_grid_nx)
    if (allocated(this%meteo_grid_ny)) deallocate(this%meteo_grid_ny)
    if (allocated(this%meteo_grid_xll)) deallocate(this%meteo_grid_xll)
    if (allocated(this%meteo_grid_yll)) deallocate(this%meteo_grid_yll)
    if (allocated(this%meteo_grid_cellsize)) deallocate(this%meteo_grid_cellsize)
    if (allocated(this%meteo_grid_ydir)) deallocate(this%meteo_grid_ydir)
    if (allocated(this%meteo_grid_coordsys)) deallocate(this%meteo_grid_coordsys)
    if (allocated(this%hydro_grid_nx)) deallocate(this%hydro_grid_nx)
    if (allocated(this%hydro_grid_ny)) deallocate(this%hydro_grid_ny)
    if (allocated(this%hydro_grid_xll)) deallocate(this%hydro_grid_xll)
    if (allocated(this%hydro_grid_yll)) deallocate(this%hydro_grid_yll)
    if (allocated(this%hydro_grid_cellsize)) deallocate(this%hydro_grid_cellsize)
    if (allocated(this%hydro_grid_ydir)) deallocate(this%hydro_grid_ydir)
    if (allocated(this%hydro_grid_coordsys)) deallocate(this%hydro_grid_coordsys)
    if (allocated(this%morph_grid_nx)) deallocate(this%morph_grid_nx)
    if (allocated(this%morph_grid_ny)) deallocate(this%morph_grid_ny)
    if (allocated(this%morph_grid_xll)) deallocate(this%morph_grid_xll)
    if (allocated(this%morph_grid_yll)) deallocate(this%morph_grid_yll)
    if (allocated(this%morph_grid_cellsize)) deallocate(this%morph_grid_cellsize)
    if (allocated(this%morph_grid_ydir)) deallocate(this%morph_grid_ydir)
    if (allocated(this%morph_grid_coordsys)) deallocate(this%morph_grid_coordsys)
    if (allocated(this%pre_coupled)) deallocate(this%pre_coupled)
    if (allocated(this%pet_coupled)) deallocate(this%pet_coupled)
    if (allocated(this%temp_coupled)) deallocate(this%temp_coupled)
    if (allocated(this%tann_coupled)) deallocate(this%tann_coupled)
    if (allocated(this%tmin_coupled)) deallocate(this%tmin_coupled)
    if (allocated(this%tmax_coupled)) deallocate(this%tmax_coupled)
    if (allocated(this%ssrd_coupled)) deallocate(this%ssrd_coupled)
    if (allocated(this%strd_coupled)) deallocate(this%strd_coupled)
    if (allocated(this%netrad_coupled)) deallocate(this%netrad_coupled)
    if (allocated(this%eabs_coupled)) deallocate(this%eabs_coupled)
    if (allocated(this%wind_coupled)) deallocate(this%wind_coupled)
    if (allocated(this%runoff_coupled)) deallocate(this%runoff_coupled)
    if (allocated(this%runoff_sealed_coupled)) deallocate(this%runoff_sealed_coupled)
    if (allocated(this%interflow_fast_coupled)) deallocate(this%interflow_fast_coupled)
    if (allocated(this%interflow_slow_coupled)) deallocate(this%interflow_slow_coupled)
    if (allocated(this%baseflow_coupled)) deallocate(this%baseflow_coupled)
    if (allocated(this%dem_coupled)) deallocate(this%dem_coupled)
    if (allocated(this%slope_coupled)) deallocate(this%slope_coupled)
    if (allocated(this%aspect_coupled)) deallocate(this%aspect_coupled)
    if (allocated(this%geo_class_coupled)) deallocate(this%geo_class_coupled)
    if (allocated(this%soil_class_coupled)) deallocate(this%soil_class_coupled)
    if (allocated(this%lai_class_coupled)) deallocate(this%lai_class_coupled)
    if (allocated(this%river_width_coupled)) deallocate(this%river_width_coupled)
    if (allocated(this%meteo_mask_coupled)) deallocate(this%meteo_mask_coupled)
    if (allocated(this%hydro_mask_coupled)) deallocate(this%hydro_mask_coupled)
    if (allocated(this%morph_mask_coupled)) deallocate(this%morph_mask_coupled)
    if (allocated(this%hydro_latlon_coupled)) deallocate(this%hydro_latlon_coupled)
    if (allocated(this%morph_latlon_coupled)) deallocate(this%morph_latlon_coupled)
    if (allocated(this%route_latlon_coupled)) deallocate(this%route_latlon_coupled)
    this%is_configured = .false.
  end function nml_config_coupling_set_dims


  !> \brief Read config_coupling namelist from file
  integer function nml_config_coupling_from_file(this, file, errmsg) result(status)
    class(nml_config_coupling_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4), allocatable, dimension(:) :: meteo_grid_nx
    integer(i4), allocatable, dimension(:) :: meteo_grid_ny
    real(dp), allocatable, dimension(:) :: meteo_grid_xll
    real(dp), allocatable, dimension(:) :: meteo_grid_yll
    real(dp), allocatable, dimension(:) :: meteo_grid_cellsize
    integer(i4), allocatable, dimension(:) :: meteo_grid_ydir
    integer(i4), allocatable, dimension(:) :: meteo_grid_coordsys
    integer(i4), allocatable, dimension(:) :: hydro_grid_nx
    integer(i4), allocatable, dimension(:) :: hydro_grid_ny
    real(dp), allocatable, dimension(:) :: hydro_grid_xll
    real(dp), allocatable, dimension(:) :: hydro_grid_yll
    real(dp), allocatable, dimension(:) :: hydro_grid_cellsize
    integer(i4), allocatable, dimension(:) :: hydro_grid_ydir
    integer(i4), allocatable, dimension(:) :: hydro_grid_coordsys
    integer(i4), allocatable, dimension(:) :: morph_grid_nx
    integer(i4), allocatable, dimension(:) :: morph_grid_ny
    real(dp), allocatable, dimension(:) :: morph_grid_xll
    real(dp), allocatable, dimension(:) :: morph_grid_yll
    real(dp), allocatable, dimension(:) :: morph_grid_cellsize
    integer(i4), allocatable, dimension(:) :: morph_grid_ydir
    integer(i4), allocatable, dimension(:) :: morph_grid_coordsys
    logical, allocatable, dimension(:) :: pre_coupled
    logical, allocatable, dimension(:) :: pet_coupled
    logical, allocatable, dimension(:) :: temp_coupled
    logical, allocatable, dimension(:) :: tann_coupled
    logical, allocatable, dimension(:) :: tmin_coupled
    logical, allocatable, dimension(:) :: tmax_coupled
    logical, allocatable, dimension(:) :: ssrd_coupled
    logical, allocatable, dimension(:) :: strd_coupled
    logical, allocatable, dimension(:) :: netrad_coupled
    logical, allocatable, dimension(:) :: eabs_coupled
    logical, allocatable, dimension(:) :: wind_coupled
    logical, allocatable, dimension(:) :: runoff_coupled
    logical, allocatable, dimension(:) :: runoff_sealed_coupled
    logical, allocatable, dimension(:) :: interflow_fast_coupled
    logical, allocatable, dimension(:) :: interflow_slow_coupled
    logical, allocatable, dimension(:) :: baseflow_coupled
    logical, allocatable, dimension(:) :: dem_coupled
    logical, allocatable, dimension(:) :: slope_coupled
    logical, allocatable, dimension(:) :: aspect_coupled
    logical, allocatable, dimension(:) :: geo_class_coupled
    logical, allocatable, dimension(:) :: soil_class_coupled
    logical, allocatable, dimension(:) :: lai_class_coupled
    logical, allocatable, dimension(:) :: river_width_coupled
    logical, allocatable, dimension(:) :: meteo_mask_coupled
    logical, allocatable, dimension(:) :: hydro_mask_coupled
    logical, allocatable, dimension(:) :: morph_mask_coupled
    logical, allocatable, dimension(:) :: hydro_latlon_coupled
    logical, allocatable, dimension(:) :: morph_latlon_coupled
    logical, allocatable, dimension(:) :: route_latlon_coupled
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_coupling/ &
      meteo_grid_nx, &
      meteo_grid_ny, &
      meteo_grid_xll, &
      meteo_grid_yll, &
      meteo_grid_cellsize, &
      meteo_grid_ydir, &
      meteo_grid_coordsys, &
      hydro_grid_nx, &
      hydro_grid_ny, &
      hydro_grid_xll, &
      hydro_grid_yll, &
      hydro_grid_cellsize, &
      hydro_grid_ydir, &
      hydro_grid_coordsys, &
      morph_grid_nx, &
      morph_grid_ny, &
      morph_grid_xll, &
      morph_grid_yll, &
      morph_grid_cellsize, &
      morph_grid_ydir, &
      morph_grid_coordsys, &
      pre_coupled, &
      pet_coupled, &
      temp_coupled, &
      tann_coupled, &
      tmin_coupled, &
      tmax_coupled, &
      ssrd_coupled, &
      strd_coupled, &
      netrad_coupled, &
      eabs_coupled, &
      wind_coupled, &
      runoff_coupled, &
      runoff_sealed_coupled, &
      interflow_fast_coupled, &
      interflow_slow_coupled, &
      baseflow_coupled, &
      dem_coupled, &
      slope_coupled, &
      aspect_coupled, &
      geo_class_coupled, &
      soil_class_coupled, &
      lai_class_coupled, &
      river_width_coupled, &
      meteo_mask_coupled, &
      hydro_mask_coupled, &
      morph_mask_coupled, &
      hydro_latlon_coupled, &
      morph_latlon_coupled, &
      route_latlon_coupled

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(meteo_grid_nx)) deallocate(meteo_grid_nx)
    allocate(meteo_grid_nx(this%n_domains))
    if (allocated(meteo_grid_ny)) deallocate(meteo_grid_ny)
    allocate(meteo_grid_ny(this%n_domains))
    if (allocated(meteo_grid_xll)) deallocate(meteo_grid_xll)
    allocate(meteo_grid_xll(this%n_domains))
    if (allocated(meteo_grid_yll)) deallocate(meteo_grid_yll)
    allocate(meteo_grid_yll(this%n_domains))
    if (allocated(meteo_grid_cellsize)) deallocate(meteo_grid_cellsize)
    allocate(meteo_grid_cellsize(this%n_domains))
    if (allocated(meteo_grid_ydir)) deallocate(meteo_grid_ydir)
    allocate(meteo_grid_ydir(this%n_domains))
    if (allocated(meteo_grid_coordsys)) deallocate(meteo_grid_coordsys)
    allocate(meteo_grid_coordsys(this%n_domains))
    if (allocated(hydro_grid_nx)) deallocate(hydro_grid_nx)
    allocate(hydro_grid_nx(this%n_domains))
    if (allocated(hydro_grid_ny)) deallocate(hydro_grid_ny)
    allocate(hydro_grid_ny(this%n_domains))
    if (allocated(hydro_grid_xll)) deallocate(hydro_grid_xll)
    allocate(hydro_grid_xll(this%n_domains))
    if (allocated(hydro_grid_yll)) deallocate(hydro_grid_yll)
    allocate(hydro_grid_yll(this%n_domains))
    if (allocated(hydro_grid_cellsize)) deallocate(hydro_grid_cellsize)
    allocate(hydro_grid_cellsize(this%n_domains))
    if (allocated(hydro_grid_ydir)) deallocate(hydro_grid_ydir)
    allocate(hydro_grid_ydir(this%n_domains))
    if (allocated(hydro_grid_coordsys)) deallocate(hydro_grid_coordsys)
    allocate(hydro_grid_coordsys(this%n_domains))
    if (allocated(morph_grid_nx)) deallocate(morph_grid_nx)
    allocate(morph_grid_nx(this%n_domains))
    if (allocated(morph_grid_ny)) deallocate(morph_grid_ny)
    allocate(morph_grid_ny(this%n_domains))
    if (allocated(morph_grid_xll)) deallocate(morph_grid_xll)
    allocate(morph_grid_xll(this%n_domains))
    if (allocated(morph_grid_yll)) deallocate(morph_grid_yll)
    allocate(morph_grid_yll(this%n_domains))
    if (allocated(morph_grid_cellsize)) deallocate(morph_grid_cellsize)
    allocate(morph_grid_cellsize(this%n_domains))
    if (allocated(morph_grid_ydir)) deallocate(morph_grid_ydir)
    allocate(morph_grid_ydir(this%n_domains))
    if (allocated(morph_grid_coordsys)) deallocate(morph_grid_coordsys)
    allocate(morph_grid_coordsys(this%n_domains))
    if (allocated(pre_coupled)) deallocate(pre_coupled)
    allocate(pre_coupled(this%n_domains))
    if (allocated(pet_coupled)) deallocate(pet_coupled)
    allocate(pet_coupled(this%n_domains))
    if (allocated(temp_coupled)) deallocate(temp_coupled)
    allocate(temp_coupled(this%n_domains))
    if (allocated(tann_coupled)) deallocate(tann_coupled)
    allocate(tann_coupled(this%n_domains))
    if (allocated(tmin_coupled)) deallocate(tmin_coupled)
    allocate(tmin_coupled(this%n_domains))
    if (allocated(tmax_coupled)) deallocate(tmax_coupled)
    allocate(tmax_coupled(this%n_domains))
    if (allocated(ssrd_coupled)) deallocate(ssrd_coupled)
    allocate(ssrd_coupled(this%n_domains))
    if (allocated(strd_coupled)) deallocate(strd_coupled)
    allocate(strd_coupled(this%n_domains))
    if (allocated(netrad_coupled)) deallocate(netrad_coupled)
    allocate(netrad_coupled(this%n_domains))
    if (allocated(eabs_coupled)) deallocate(eabs_coupled)
    allocate(eabs_coupled(this%n_domains))
    if (allocated(wind_coupled)) deallocate(wind_coupled)
    allocate(wind_coupled(this%n_domains))
    if (allocated(runoff_coupled)) deallocate(runoff_coupled)
    allocate(runoff_coupled(this%n_domains))
    if (allocated(runoff_sealed_coupled)) deallocate(runoff_sealed_coupled)
    allocate(runoff_sealed_coupled(this%n_domains))
    if (allocated(interflow_fast_coupled)) deallocate(interflow_fast_coupled)
    allocate(interflow_fast_coupled(this%n_domains))
    if (allocated(interflow_slow_coupled)) deallocate(interflow_slow_coupled)
    allocate(interflow_slow_coupled(this%n_domains))
    if (allocated(baseflow_coupled)) deallocate(baseflow_coupled)
    allocate(baseflow_coupled(this%n_domains))
    if (allocated(dem_coupled)) deallocate(dem_coupled)
    allocate(dem_coupled(this%n_domains))
    if (allocated(slope_coupled)) deallocate(slope_coupled)
    allocate(slope_coupled(this%n_domains))
    if (allocated(aspect_coupled)) deallocate(aspect_coupled)
    allocate(aspect_coupled(this%n_domains))
    if (allocated(geo_class_coupled)) deallocate(geo_class_coupled)
    allocate(geo_class_coupled(this%n_domains))
    if (allocated(soil_class_coupled)) deallocate(soil_class_coupled)
    allocate(soil_class_coupled(this%n_domains))
    if (allocated(lai_class_coupled)) deallocate(lai_class_coupled)
    allocate(lai_class_coupled(this%n_domains))
    if (allocated(river_width_coupled)) deallocate(river_width_coupled)
    allocate(river_width_coupled(this%n_domains))
    if (allocated(meteo_mask_coupled)) deallocate(meteo_mask_coupled)
    allocate(meteo_mask_coupled(this%n_domains))
    if (allocated(hydro_mask_coupled)) deallocate(hydro_mask_coupled)
    allocate(hydro_mask_coupled(this%n_domains))
    if (allocated(morph_mask_coupled)) deallocate(morph_mask_coupled)
    allocate(morph_mask_coupled(this%n_domains))
    if (allocated(hydro_latlon_coupled)) deallocate(hydro_latlon_coupled)
    allocate(hydro_latlon_coupled(this%n_domains))
    if (allocated(morph_latlon_coupled)) deallocate(morph_latlon_coupled)
    allocate(morph_latlon_coupled(this%n_domains))
    if (allocated(route_latlon_coupled)) deallocate(route_latlon_coupled)
    allocate(route_latlon_coupled(this%n_domains))
    meteo_grid_nx = this%meteo_grid_nx
    meteo_grid_ny = this%meteo_grid_ny
    meteo_grid_xll = this%meteo_grid_xll
    meteo_grid_yll = this%meteo_grid_yll
    meteo_grid_cellsize = this%meteo_grid_cellsize
    meteo_grid_ydir = this%meteo_grid_ydir
    meteo_grid_coordsys = this%meteo_grid_coordsys
    hydro_grid_nx = this%hydro_grid_nx
    hydro_grid_ny = this%hydro_grid_ny
    hydro_grid_xll = this%hydro_grid_xll
    hydro_grid_yll = this%hydro_grid_yll
    hydro_grid_cellsize = this%hydro_grid_cellsize
    hydro_grid_ydir = this%hydro_grid_ydir
    hydro_grid_coordsys = this%hydro_grid_coordsys
    morph_grid_nx = this%morph_grid_nx
    morph_grid_ny = this%morph_grid_ny
    morph_grid_xll = this%morph_grid_xll
    morph_grid_yll = this%morph_grid_yll
    morph_grid_cellsize = this%morph_grid_cellsize
    morph_grid_ydir = this%morph_grid_ydir
    morph_grid_coordsys = this%morph_grid_coordsys
    pre_coupled = this%pre_coupled
    pet_coupled = this%pet_coupled
    temp_coupled = this%temp_coupled
    tann_coupled = this%tann_coupled
    tmin_coupled = this%tmin_coupled
    tmax_coupled = this%tmax_coupled
    ssrd_coupled = this%ssrd_coupled
    strd_coupled = this%strd_coupled
    netrad_coupled = this%netrad_coupled
    eabs_coupled = this%eabs_coupled
    wind_coupled = this%wind_coupled
    runoff_coupled = this%runoff_coupled
    runoff_sealed_coupled = this%runoff_sealed_coupled
    interflow_fast_coupled = this%interflow_fast_coupled
    interflow_slow_coupled = this%interflow_slow_coupled
    baseflow_coupled = this%baseflow_coupled
    dem_coupled = this%dem_coupled
    slope_coupled = this%slope_coupled
    aspect_coupled = this%aspect_coupled
    geo_class_coupled = this%geo_class_coupled
    soil_class_coupled = this%soil_class_coupled
    lai_class_coupled = this%lai_class_coupled
    river_width_coupled = this%river_width_coupled
    meteo_mask_coupled = this%meteo_mask_coupled
    hydro_mask_coupled = this%hydro_mask_coupled
    morph_mask_coupled = this%morph_mask_coupled
    hydro_latlon_coupled = this%hydro_latlon_coupled
    morph_latlon_coupled = this%morph_latlon_coupled
    route_latlon_coupled = this%route_latlon_coupled

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_coupling", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_coupling, iostat=iostat, iomsg=iomsg)
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
    this%meteo_grid_nx = meteo_grid_nx
    this%meteo_grid_ny = meteo_grid_ny
    this%meteo_grid_xll = meteo_grid_xll
    this%meteo_grid_yll = meteo_grid_yll
    this%meteo_grid_cellsize = meteo_grid_cellsize
    this%meteo_grid_ydir = meteo_grid_ydir
    this%meteo_grid_coordsys = meteo_grid_coordsys
    this%hydro_grid_nx = hydro_grid_nx
    this%hydro_grid_ny = hydro_grid_ny
    this%hydro_grid_xll = hydro_grid_xll
    this%hydro_grid_yll = hydro_grid_yll
    this%hydro_grid_cellsize = hydro_grid_cellsize
    this%hydro_grid_ydir = hydro_grid_ydir
    this%hydro_grid_coordsys = hydro_grid_coordsys
    this%morph_grid_nx = morph_grid_nx
    this%morph_grid_ny = morph_grid_ny
    this%morph_grid_xll = morph_grid_xll
    this%morph_grid_yll = morph_grid_yll
    this%morph_grid_cellsize = morph_grid_cellsize
    this%morph_grid_ydir = morph_grid_ydir
    this%morph_grid_coordsys = morph_grid_coordsys
    this%pre_coupled = pre_coupled
    this%pet_coupled = pet_coupled
    this%temp_coupled = temp_coupled
    this%tann_coupled = tann_coupled
    this%tmin_coupled = tmin_coupled
    this%tmax_coupled = tmax_coupled
    this%ssrd_coupled = ssrd_coupled
    this%strd_coupled = strd_coupled
    this%netrad_coupled = netrad_coupled
    this%eabs_coupled = eabs_coupled
    this%wind_coupled = wind_coupled
    this%runoff_coupled = runoff_coupled
    this%runoff_sealed_coupled = runoff_sealed_coupled
    this%interflow_fast_coupled = interflow_fast_coupled
    this%interflow_slow_coupled = interflow_slow_coupled
    this%baseflow_coupled = baseflow_coupled
    this%dem_coupled = dem_coupled
    this%slope_coupled = slope_coupled
    this%aspect_coupled = aspect_coupled
    this%geo_class_coupled = geo_class_coupled
    this%soil_class_coupled = soil_class_coupled
    this%lai_class_coupled = lai_class_coupled
    this%river_width_coupled = river_width_coupled
    this%meteo_mask_coupled = meteo_mask_coupled
    this%hydro_mask_coupled = hydro_mask_coupled
    this%morph_mask_coupled = morph_mask_coupled
    this%hydro_latlon_coupled = hydro_latlon_coupled
    this%morph_latlon_coupled = morph_latlon_coupled
    this%route_latlon_coupled = route_latlon_coupled

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_coupling_from_file

  !> \brief Set config_coupling values
  integer function nml_config_coupling_set(this, &
    meteo_grid_nx, &
    meteo_grid_ny, &
    meteo_grid_xll, &
    meteo_grid_yll, &
    meteo_grid_cellsize, &
    meteo_grid_ydir, &
    meteo_grid_coordsys, &
    hydro_grid_nx, &
    hydro_grid_ny, &
    hydro_grid_xll, &
    hydro_grid_yll, &
    hydro_grid_cellsize, &
    hydro_grid_ydir, &
    hydro_grid_coordsys, &
    morph_grid_nx, &
    morph_grid_ny, &
    morph_grid_xll, &
    morph_grid_yll, &
    morph_grid_cellsize, &
    morph_grid_ydir, &
    morph_grid_coordsys, &
    pre_coupled, &
    pet_coupled, &
    temp_coupled, &
    tann_coupled, &
    tmin_coupled, &
    tmax_coupled, &
    ssrd_coupled, &
    strd_coupled, &
    netrad_coupled, &
    eabs_coupled, &
    wind_coupled, &
    runoff_coupled, &
    runoff_sealed_coupled, &
    interflow_fast_coupled, &
    interflow_slow_coupled, &
    baseflow_coupled, &
    dem_coupled, &
    slope_coupled, &
    aspect_coupled, &
    geo_class_coupled, &
    soil_class_coupled, &
    lai_class_coupled, &
    river_width_coupled, &
    meteo_mask_coupled, &
    hydro_mask_coupled, &
    morph_mask_coupled, &
    hydro_latlon_coupled, &
    morph_latlon_coupled, &
    route_latlon_coupled, &
    errmsg) result(status)

    class(nml_config_coupling_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), dimension(:), intent(in), optional :: meteo_grid_nx !< Meteo grid size in x-direction
    integer(i4), dimension(:), intent(in), optional :: meteo_grid_ny !< Meteo grid size in y-direction
    real(dp), dimension(:), intent(in), optional :: meteo_grid_xll !< Meteo grid x origin
    real(dp), dimension(:), intent(in), optional :: meteo_grid_yll !< Meteo grid y origin
    real(dp), dimension(:), intent(in), optional :: meteo_grid_cellsize !< Meteo grid cell size
    integer(i4), dimension(:), intent(in), optional :: meteo_grid_ydir !< Meteo grid y direction
    integer(i4), dimension(:), intent(in), optional :: meteo_grid_coordsys !< Meteo grid coordinate system
    integer(i4), dimension(:), intent(in), optional :: hydro_grid_nx !< Hydro grid size in x-direction
    integer(i4), dimension(:), intent(in), optional :: hydro_grid_ny !< Hydro grid size in y-direction
    real(dp), dimension(:), intent(in), optional :: hydro_grid_xll !< Hydro grid x origin
    real(dp), dimension(:), intent(in), optional :: hydro_grid_yll !< Hydro grid y origin
    real(dp), dimension(:), intent(in), optional :: hydro_grid_cellsize !< Hydro grid cell size
    integer(i4), dimension(:), intent(in), optional :: hydro_grid_ydir !< Hydro grid y direction
    integer(i4), dimension(:), intent(in), optional :: hydro_grid_coordsys !< Hydro grid coordinate system
    integer(i4), dimension(:), intent(in), optional :: morph_grid_nx !< Morph grid size in x-direction
    integer(i4), dimension(:), intent(in), optional :: morph_grid_ny !< Morph grid size in y-direction
    real(dp), dimension(:), intent(in), optional :: morph_grid_xll !< Morph grid x origin
    real(dp), dimension(:), intent(in), optional :: morph_grid_yll !< Morph grid y origin
    real(dp), dimension(:), intent(in), optional :: morph_grid_cellsize !< Morph grid cell size
    integer(i4), dimension(:), intent(in), optional :: morph_grid_ydir !< Morph grid y direction
    integer(i4), dimension(:), intent(in), optional :: morph_grid_coordsys !< Morph grid coordinate system
    logical, dimension(:), intent(in), optional :: pre_coupled !< Precipitation coupled
    logical, dimension(:), intent(in), optional :: pet_coupled !< Potential evapotranspiration coupled
    logical, dimension(:), intent(in), optional :: temp_coupled !< Air temperature coupled
    logical, dimension(:), intent(in), optional :: tann_coupled !< Air temperature annual mean coupled
    logical, dimension(:), intent(in), optional :: tmin_coupled !< Air temperature daily minimum coupled
    logical, dimension(:), intent(in), optional :: tmax_coupled !< Air temperature daily maximum coupled
    logical, dimension(:), intent(in), optional :: ssrd_coupled !< Surface shortwave radiation coupled
    logical, dimension(:), intent(in), optional :: strd_coupled !< Surface thermal radiation coupled
    logical, dimension(:), intent(in), optional :: netrad_coupled !< Net radiation coupled
    logical, dimension(:), intent(in), optional :: eabs_coupled !< Vapor pressure coupled
    logical, dimension(:), intent(in), optional :: wind_coupled !< Wind speed coupled
    logical, dimension(:), intent(in), optional :: runoff_coupled !< Runoff coupled
    logical, dimension(:), intent(in), optional :: runoff_sealed_coupled !< Sealed runoff coupled
    logical, dimension(:), intent(in), optional :: interflow_fast_coupled !< Fast interflow coupled
    logical, dimension(:), intent(in), optional :: interflow_slow_coupled !< Slow interflow coupled
    logical, dimension(:), intent(in), optional :: baseflow_coupled !< Baseflow coupled
    logical, dimension(:), intent(in), optional :: dem_coupled !< DEM coupled
    logical, dimension(:), intent(in), optional :: slope_coupled !< Slope coupled
    logical, dimension(:), intent(in), optional :: aspect_coupled !< Aspect coupled
    logical, dimension(:), intent(in), optional :: geo_class_coupled !< Geology class coupled
    logical, dimension(:), intent(in), optional :: soil_class_coupled !< Soil class coupled
    logical, dimension(:), intent(in), optional :: lai_class_coupled !< LAI class coupled
    logical, dimension(:), intent(in), optional :: river_width_coupled !< River width coupled
    logical, dimension(:), intent(in), optional :: meteo_mask_coupled !< Meteorological mask coupled
    logical, dimension(:), intent(in), optional :: hydro_mask_coupled !< Hydrological mask coupled
    logical, dimension(:), intent(in), optional :: morph_mask_coupled !< Morphology mask coupled
    logical, dimension(:), intent(in), optional :: hydro_latlon_coupled !< Hydrological latlon coupled
    logical, dimension(:), intent(in), optional :: morph_latlon_coupled !< Morphological latlon coupled
    logical, dimension(:), intent(in), optional :: route_latlon_coupled !< Routing latlon coupled
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(meteo_grid_nx)) then
      if (size(meteo_grid_nx, 1) > size(this%meteo_grid_nx, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_nx'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_nx, 1)
      ub__1 = lb__1 + size(meteo_grid_nx, 1) - 1
      this%meteo_grid_nx(lb__1:ub__1) = meteo_grid_nx
    end if
    if (present(meteo_grid_ny)) then
      if (size(meteo_grid_ny, 1) > size(this%meteo_grid_ny, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_ny'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_ny, 1)
      ub__1 = lb__1 + size(meteo_grid_ny, 1) - 1
      this%meteo_grid_ny(lb__1:ub__1) = meteo_grid_ny
    end if
    if (present(meteo_grid_xll)) then
      if (size(meteo_grid_xll, 1) > size(this%meteo_grid_xll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_xll'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_xll, 1)
      ub__1 = lb__1 + size(meteo_grid_xll, 1) - 1
      this%meteo_grid_xll(lb__1:ub__1) = meteo_grid_xll
    end if
    if (present(meteo_grid_yll)) then
      if (size(meteo_grid_yll, 1) > size(this%meteo_grid_yll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_yll'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_yll, 1)
      ub__1 = lb__1 + size(meteo_grid_yll, 1) - 1
      this%meteo_grid_yll(lb__1:ub__1) = meteo_grid_yll
    end if
    if (present(meteo_grid_cellsize)) then
      if (size(meteo_grid_cellsize, 1) > size(this%meteo_grid_cellsize, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_cellsize'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_cellsize, 1)
      ub__1 = lb__1 + size(meteo_grid_cellsize, 1) - 1
      this%meteo_grid_cellsize(lb__1:ub__1) = meteo_grid_cellsize
    end if
    if (present(meteo_grid_ydir)) then
      if (size(meteo_grid_ydir, 1) > size(this%meteo_grid_ydir, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_ydir'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_ydir, 1)
      ub__1 = lb__1 + size(meteo_grid_ydir, 1) - 1
      this%meteo_grid_ydir(lb__1:ub__1) = meteo_grid_ydir
    end if
    if (present(meteo_grid_coordsys)) then
      if (size(meteo_grid_coordsys, 1) > size(this%meteo_grid_coordsys, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_grid_coordsys'"
        return
      end if
      lb__1 = lbound(this%meteo_grid_coordsys, 1)
      ub__1 = lb__1 + size(meteo_grid_coordsys, 1) - 1
      this%meteo_grid_coordsys(lb__1:ub__1) = meteo_grid_coordsys
    end if
    if (present(hydro_grid_nx)) then
      if (size(hydro_grid_nx, 1) > size(this%hydro_grid_nx, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_nx'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_nx, 1)
      ub__1 = lb__1 + size(hydro_grid_nx, 1) - 1
      this%hydro_grid_nx(lb__1:ub__1) = hydro_grid_nx
    end if
    if (present(hydro_grid_ny)) then
      if (size(hydro_grid_ny, 1) > size(this%hydro_grid_ny, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_ny'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_ny, 1)
      ub__1 = lb__1 + size(hydro_grid_ny, 1) - 1
      this%hydro_grid_ny(lb__1:ub__1) = hydro_grid_ny
    end if
    if (present(hydro_grid_xll)) then
      if (size(hydro_grid_xll, 1) > size(this%hydro_grid_xll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_xll'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_xll, 1)
      ub__1 = lb__1 + size(hydro_grid_xll, 1) - 1
      this%hydro_grid_xll(lb__1:ub__1) = hydro_grid_xll
    end if
    if (present(hydro_grid_yll)) then
      if (size(hydro_grid_yll, 1) > size(this%hydro_grid_yll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_yll'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_yll, 1)
      ub__1 = lb__1 + size(hydro_grid_yll, 1) - 1
      this%hydro_grid_yll(lb__1:ub__1) = hydro_grid_yll
    end if
    if (present(hydro_grid_cellsize)) then
      if (size(hydro_grid_cellsize, 1) > size(this%hydro_grid_cellsize, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_cellsize'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_cellsize, 1)
      ub__1 = lb__1 + size(hydro_grid_cellsize, 1) - 1
      this%hydro_grid_cellsize(lb__1:ub__1) = hydro_grid_cellsize
    end if
    if (present(hydro_grid_ydir)) then
      if (size(hydro_grid_ydir, 1) > size(this%hydro_grid_ydir, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_ydir'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_ydir, 1)
      ub__1 = lb__1 + size(hydro_grid_ydir, 1) - 1
      this%hydro_grid_ydir(lb__1:ub__1) = hydro_grid_ydir
    end if
    if (present(hydro_grid_coordsys)) then
      if (size(hydro_grid_coordsys, 1) > size(this%hydro_grid_coordsys, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_grid_coordsys'"
        return
      end if
      lb__1 = lbound(this%hydro_grid_coordsys, 1)
      ub__1 = lb__1 + size(hydro_grid_coordsys, 1) - 1
      this%hydro_grid_coordsys(lb__1:ub__1) = hydro_grid_coordsys
    end if
    if (present(morph_grid_nx)) then
      if (size(morph_grid_nx, 1) > size(this%morph_grid_nx, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_nx'"
        return
      end if
      lb__1 = lbound(this%morph_grid_nx, 1)
      ub__1 = lb__1 + size(morph_grid_nx, 1) - 1
      this%morph_grid_nx(lb__1:ub__1) = morph_grid_nx
    end if
    if (present(morph_grid_ny)) then
      if (size(morph_grid_ny, 1) > size(this%morph_grid_ny, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_ny'"
        return
      end if
      lb__1 = lbound(this%morph_grid_ny, 1)
      ub__1 = lb__1 + size(morph_grid_ny, 1) - 1
      this%morph_grid_ny(lb__1:ub__1) = morph_grid_ny
    end if
    if (present(morph_grid_xll)) then
      if (size(morph_grid_xll, 1) > size(this%morph_grid_xll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_xll'"
        return
      end if
      lb__1 = lbound(this%morph_grid_xll, 1)
      ub__1 = lb__1 + size(morph_grid_xll, 1) - 1
      this%morph_grid_xll(lb__1:ub__1) = morph_grid_xll
    end if
    if (present(morph_grid_yll)) then
      if (size(morph_grid_yll, 1) > size(this%morph_grid_yll, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_yll'"
        return
      end if
      lb__1 = lbound(this%morph_grid_yll, 1)
      ub__1 = lb__1 + size(morph_grid_yll, 1) - 1
      this%morph_grid_yll(lb__1:ub__1) = morph_grid_yll
    end if
    if (present(morph_grid_cellsize)) then
      if (size(morph_grid_cellsize, 1) > size(this%morph_grid_cellsize, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_cellsize'"
        return
      end if
      lb__1 = lbound(this%morph_grid_cellsize, 1)
      ub__1 = lb__1 + size(morph_grid_cellsize, 1) - 1
      this%morph_grid_cellsize(lb__1:ub__1) = morph_grid_cellsize
    end if
    if (present(morph_grid_ydir)) then
      if (size(morph_grid_ydir, 1) > size(this%morph_grid_ydir, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_ydir'"
        return
      end if
      lb__1 = lbound(this%morph_grid_ydir, 1)
      ub__1 = lb__1 + size(morph_grid_ydir, 1) - 1
      this%morph_grid_ydir(lb__1:ub__1) = morph_grid_ydir
    end if
    if (present(morph_grid_coordsys)) then
      if (size(morph_grid_coordsys, 1) > size(this%morph_grid_coordsys, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_grid_coordsys'"
        return
      end if
      lb__1 = lbound(this%morph_grid_coordsys, 1)
      ub__1 = lb__1 + size(morph_grid_coordsys, 1) - 1
      this%morph_grid_coordsys(lb__1:ub__1) = morph_grid_coordsys
    end if
    if (present(pre_coupled)) then
      if (size(pre_coupled, 1) > size(this%pre_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pre_coupled'"
        return
      end if
      lb__1 = lbound(this%pre_coupled, 1)
      ub__1 = lb__1 + size(pre_coupled, 1) - 1
      this%pre_coupled(lb__1:ub__1) = pre_coupled
    end if
    if (present(pet_coupled)) then
      if (size(pet_coupled, 1) > size(this%pet_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_coupled'"
        return
      end if
      lb__1 = lbound(this%pet_coupled, 1)
      ub__1 = lb__1 + size(pet_coupled, 1) - 1
      this%pet_coupled(lb__1:ub__1) = pet_coupled
    end if
    if (present(temp_coupled)) then
      if (size(temp_coupled, 1) > size(this%temp_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'temp_coupled'"
        return
      end if
      lb__1 = lbound(this%temp_coupled, 1)
      ub__1 = lb__1 + size(temp_coupled, 1) - 1
      this%temp_coupled(lb__1:ub__1) = temp_coupled
    end if
    if (present(tann_coupled)) then
      if (size(tann_coupled, 1) > size(this%tann_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tann_coupled'"
        return
      end if
      lb__1 = lbound(this%tann_coupled, 1)
      ub__1 = lb__1 + size(tann_coupled, 1) - 1
      this%tann_coupled(lb__1:ub__1) = tann_coupled
    end if
    if (present(tmin_coupled)) then
      if (size(tmin_coupled, 1) > size(this%tmin_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmin_coupled'"
        return
      end if
      lb__1 = lbound(this%tmin_coupled, 1)
      ub__1 = lb__1 + size(tmin_coupled, 1) - 1
      this%tmin_coupled(lb__1:ub__1) = tmin_coupled
    end if
    if (present(tmax_coupled)) then
      if (size(tmax_coupled, 1) > size(this%tmax_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'tmax_coupled'"
        return
      end if
      lb__1 = lbound(this%tmax_coupled, 1)
      ub__1 = lb__1 + size(tmax_coupled, 1) - 1
      this%tmax_coupled(lb__1:ub__1) = tmax_coupled
    end if
    if (present(ssrd_coupled)) then
      if (size(ssrd_coupled, 1) > size(this%ssrd_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ssrd_coupled'"
        return
      end if
      lb__1 = lbound(this%ssrd_coupled, 1)
      ub__1 = lb__1 + size(ssrd_coupled, 1) - 1
      this%ssrd_coupled(lb__1:ub__1) = ssrd_coupled
    end if
    if (present(strd_coupled)) then
      if (size(strd_coupled, 1) > size(this%strd_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'strd_coupled'"
        return
      end if
      lb__1 = lbound(this%strd_coupled, 1)
      ub__1 = lb__1 + size(strd_coupled, 1) - 1
      this%strd_coupled(lb__1:ub__1) = strd_coupled
    end if
    if (present(netrad_coupled)) then
      if (size(netrad_coupled, 1) > size(this%netrad_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'netrad_coupled'"
        return
      end if
      lb__1 = lbound(this%netrad_coupled, 1)
      ub__1 = lb__1 + size(netrad_coupled, 1) - 1
      this%netrad_coupled(lb__1:ub__1) = netrad_coupled
    end if
    if (present(eabs_coupled)) then
      if (size(eabs_coupled, 1) > size(this%eabs_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'eabs_coupled'"
        return
      end if
      lb__1 = lbound(this%eabs_coupled, 1)
      ub__1 = lb__1 + size(eabs_coupled, 1) - 1
      this%eabs_coupled(lb__1:ub__1) = eabs_coupled
    end if
    if (present(wind_coupled)) then
      if (size(wind_coupled, 1) > size(this%wind_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'wind_coupled'"
        return
      end if
      lb__1 = lbound(this%wind_coupled, 1)
      ub__1 = lb__1 + size(wind_coupled, 1) - 1
      this%wind_coupled(lb__1:ub__1) = wind_coupled
    end if
    if (present(runoff_coupled)) then
      if (size(runoff_coupled, 1) > size(this%runoff_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_coupled'"
        return
      end if
      lb__1 = lbound(this%runoff_coupled, 1)
      ub__1 = lb__1 + size(runoff_coupled, 1) - 1
      this%runoff_coupled(lb__1:ub__1) = runoff_coupled
    end if
    if (present(runoff_sealed_coupled)) then
      if (size(runoff_sealed_coupled, 1) > size(this%runoff_sealed_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'runoff_sealed_coupled'"
        return
      end if
      lb__1 = lbound(this%runoff_sealed_coupled, 1)
      ub__1 = lb__1 + size(runoff_sealed_coupled, 1) - 1
      this%runoff_sealed_coupled(lb__1:ub__1) = runoff_sealed_coupled
    end if
    if (present(interflow_fast_coupled)) then
      if (size(interflow_fast_coupled, 1) > size(this%interflow_fast_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_fast_coupled'"
        return
      end if
      lb__1 = lbound(this%interflow_fast_coupled, 1)
      ub__1 = lb__1 + size(interflow_fast_coupled, 1) - 1
      this%interflow_fast_coupled(lb__1:ub__1) = interflow_fast_coupled
    end if
    if (present(interflow_slow_coupled)) then
      if (size(interflow_slow_coupled, 1) > size(this%interflow_slow_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'interflow_slow_coupled'"
        return
      end if
      lb__1 = lbound(this%interflow_slow_coupled, 1)
      ub__1 = lb__1 + size(interflow_slow_coupled, 1) - 1
      this%interflow_slow_coupled(lb__1:ub__1) = interflow_slow_coupled
    end if
    if (present(baseflow_coupled)) then
      if (size(baseflow_coupled, 1) > size(this%baseflow_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'baseflow_coupled'"
        return
      end if
      lb__1 = lbound(this%baseflow_coupled, 1)
      ub__1 = lb__1 + size(baseflow_coupled, 1) - 1
      this%baseflow_coupled(lb__1:ub__1) = baseflow_coupled
    end if
    if (present(dem_coupled)) then
      if (size(dem_coupled, 1) > size(this%dem_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'dem_coupled'"
        return
      end if
      lb__1 = lbound(this%dem_coupled, 1)
      ub__1 = lb__1 + size(dem_coupled, 1) - 1
      this%dem_coupled(lb__1:ub__1) = dem_coupled
    end if
    if (present(slope_coupled)) then
      if (size(slope_coupled, 1) > size(this%slope_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'slope_coupled'"
        return
      end if
      lb__1 = lbound(this%slope_coupled, 1)
      ub__1 = lb__1 + size(slope_coupled, 1) - 1
      this%slope_coupled(lb__1:ub__1) = slope_coupled
    end if
    if (present(aspect_coupled)) then
      if (size(aspect_coupled, 1) > size(this%aspect_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'aspect_coupled'"
        return
      end if
      lb__1 = lbound(this%aspect_coupled, 1)
      ub__1 = lb__1 + size(aspect_coupled, 1) - 1
      this%aspect_coupled(lb__1:ub__1) = aspect_coupled
    end if
    if (present(geo_class_coupled)) then
      if (size(geo_class_coupled, 1) > size(this%geo_class_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'geo_class_coupled'"
        return
      end if
      lb__1 = lbound(this%geo_class_coupled, 1)
      ub__1 = lb__1 + size(geo_class_coupled, 1) - 1
      this%geo_class_coupled(lb__1:ub__1) = geo_class_coupled
    end if
    if (present(soil_class_coupled)) then
      if (size(soil_class_coupled, 1) > size(this%soil_class_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'soil_class_coupled'"
        return
      end if
      lb__1 = lbound(this%soil_class_coupled, 1)
      ub__1 = lb__1 + size(soil_class_coupled, 1) - 1
      this%soil_class_coupled(lb__1:ub__1) = soil_class_coupled
    end if
    if (present(lai_class_coupled)) then
      if (size(lai_class_coupled, 1) > size(this%lai_class_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'lai_class_coupled'"
        return
      end if
      lb__1 = lbound(this%lai_class_coupled, 1)
      ub__1 = lb__1 + size(lai_class_coupled, 1) - 1
      this%lai_class_coupled(lb__1:ub__1) = lai_class_coupled
    end if
    if (present(river_width_coupled)) then
      if (size(river_width_coupled, 1) > size(this%river_width_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'river_width_coupled'"
        return
      end if
      lb__1 = lbound(this%river_width_coupled, 1)
      ub__1 = lb__1 + size(river_width_coupled, 1) - 1
      this%river_width_coupled(lb__1:ub__1) = river_width_coupled
    end if
    if (present(meteo_mask_coupled)) then
      if (size(meteo_mask_coupled, 1) > size(this%meteo_mask_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'meteo_mask_coupled'"
        return
      end if
      lb__1 = lbound(this%meteo_mask_coupled, 1)
      ub__1 = lb__1 + size(meteo_mask_coupled, 1) - 1
      this%meteo_mask_coupled(lb__1:ub__1) = meteo_mask_coupled
    end if
    if (present(hydro_mask_coupled)) then
      if (size(hydro_mask_coupled, 1) > size(this%hydro_mask_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_mask_coupled'"
        return
      end if
      lb__1 = lbound(this%hydro_mask_coupled, 1)
      ub__1 = lb__1 + size(hydro_mask_coupled, 1) - 1
      this%hydro_mask_coupled(lb__1:ub__1) = hydro_mask_coupled
    end if
    if (present(morph_mask_coupled)) then
      if (size(morph_mask_coupled, 1) > size(this%morph_mask_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_mask_coupled'"
        return
      end if
      lb__1 = lbound(this%morph_mask_coupled, 1)
      ub__1 = lb__1 + size(morph_mask_coupled, 1) - 1
      this%morph_mask_coupled(lb__1:ub__1) = morph_mask_coupled
    end if
    if (present(hydro_latlon_coupled)) then
      if (size(hydro_latlon_coupled, 1) > size(this%hydro_latlon_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'hydro_latlon_coupled'"
        return
      end if
      lb__1 = lbound(this%hydro_latlon_coupled, 1)
      ub__1 = lb__1 + size(hydro_latlon_coupled, 1) - 1
      this%hydro_latlon_coupled(lb__1:ub__1) = hydro_latlon_coupled
    end if
    if (present(morph_latlon_coupled)) then
      if (size(morph_latlon_coupled, 1) > size(this%morph_latlon_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'morph_latlon_coupled'"
        return
      end if
      lb__1 = lbound(this%morph_latlon_coupled, 1)
      ub__1 = lb__1 + size(morph_latlon_coupled, 1) - 1
      this%morph_latlon_coupled(lb__1:ub__1) = morph_latlon_coupled
    end if
    if (present(route_latlon_coupled)) then
      if (size(route_latlon_coupled, 1) > size(this%route_latlon_coupled, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'route_latlon_coupled'"
        return
      end if
      lb__1 = lbound(this%route_latlon_coupled, 1)
      ub__1 = lb__1 + size(route_latlon_coupled, 1) - 1
      this%route_latlon_coupled(lb__1:ub__1) = route_latlon_coupled
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_coupling_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_coupling_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_coupling_t), intent(in) :: this !< namelist instance
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
    case ("meteo_grid_nx")
      if (.not. allocated(this%meteo_grid_nx)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_nx), ubound(this%meteo_grid_nx), &
          "meteo_grid_nx", errmsg)
        if (status /= NML_OK) return
        if (this%meteo_grid_nx(idx(1)) == -huge(this%meteo_grid_nx(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%meteo_grid_nx == -huge(this%meteo_grid_nx))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_grid_ny")
      if (.not. allocated(this%meteo_grid_ny)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_ny), ubound(this%meteo_grid_ny), &
          "meteo_grid_ny", errmsg)
        if (status /= NML_OK) return
        if (this%meteo_grid_ny(idx(1)) == -huge(this%meteo_grid_ny(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%meteo_grid_ny == -huge(this%meteo_grid_ny))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_grid_xll")
      if (.not. allocated(this%meteo_grid_xll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_xll), ubound(this%meteo_grid_xll), &
          "meteo_grid_xll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%meteo_grid_xll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%meteo_grid_xll))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_grid_yll")
      if (.not. allocated(this%meteo_grid_yll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_yll), ubound(this%meteo_grid_yll), &
          "meteo_grid_yll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%meteo_grid_yll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%meteo_grid_yll))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_grid_cellsize")
      if (.not. allocated(this%meteo_grid_cellsize)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_cellsize), ubound(this%meteo_grid_cellsize), &
          "meteo_grid_cellsize", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%meteo_grid_cellsize(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%meteo_grid_cellsize))) status = NML_ERR_NOT_SET
      end if
    case ("meteo_grid_ydir")
      if (.not. allocated(this%meteo_grid_ydir)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_ydir), ubound(this%meteo_grid_ydir), &
          "meteo_grid_ydir", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("meteo_grid_coordsys")
      if (.not. allocated(this%meteo_grid_coordsys)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_grid_coordsys), ubound(this%meteo_grid_coordsys), &
          "meteo_grid_coordsys", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_grid_nx")
      if (.not. allocated(this%hydro_grid_nx)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_nx), ubound(this%hydro_grid_nx), &
          "hydro_grid_nx", errmsg)
        if (status /= NML_OK) return
        if (this%hydro_grid_nx(idx(1)) == -huge(this%hydro_grid_nx(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%hydro_grid_nx == -huge(this%hydro_grid_nx))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_grid_ny")
      if (.not. allocated(this%hydro_grid_ny)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_ny), ubound(this%hydro_grid_ny), &
          "hydro_grid_ny", errmsg)
        if (status /= NML_OK) return
        if (this%hydro_grid_ny(idx(1)) == -huge(this%hydro_grid_ny(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%hydro_grid_ny == -huge(this%hydro_grid_ny))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_grid_xll")
      if (.not. allocated(this%hydro_grid_xll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_xll), ubound(this%hydro_grid_xll), &
          "hydro_grid_xll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%hydro_grid_xll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%hydro_grid_xll))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_grid_yll")
      if (.not. allocated(this%hydro_grid_yll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_yll), ubound(this%hydro_grid_yll), &
          "hydro_grid_yll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%hydro_grid_yll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%hydro_grid_yll))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_grid_cellsize")
      if (.not. allocated(this%hydro_grid_cellsize)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_cellsize), ubound(this%hydro_grid_cellsize), &
          "hydro_grid_cellsize", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%hydro_grid_cellsize(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%hydro_grid_cellsize))) status = NML_ERR_NOT_SET
      end if
    case ("hydro_grid_ydir")
      if (.not. allocated(this%hydro_grid_ydir)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_ydir), ubound(this%hydro_grid_ydir), &
          "hydro_grid_ydir", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_grid_coordsys")
      if (.not. allocated(this%hydro_grid_coordsys)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_grid_coordsys), ubound(this%hydro_grid_coordsys), &
          "hydro_grid_coordsys", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_grid_nx")
      if (.not. allocated(this%morph_grid_nx)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_nx), ubound(this%morph_grid_nx), &
          "morph_grid_nx", errmsg)
        if (status /= NML_OK) return
        if (this%morph_grid_nx(idx(1)) == -huge(this%morph_grid_nx(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%morph_grid_nx == -huge(this%morph_grid_nx))) status = NML_ERR_NOT_SET
      end if
    case ("morph_grid_ny")
      if (.not. allocated(this%morph_grid_ny)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_ny), ubound(this%morph_grid_ny), &
          "morph_grid_ny", errmsg)
        if (status /= NML_OK) return
        if (this%morph_grid_ny(idx(1)) == -huge(this%morph_grid_ny(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(this%morph_grid_ny == -huge(this%morph_grid_ny))) status = NML_ERR_NOT_SET
      end if
    case ("morph_grid_xll")
      if (.not. allocated(this%morph_grid_xll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_xll), ubound(this%morph_grid_xll), &
          "morph_grid_xll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%morph_grid_xll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%morph_grid_xll))) status = NML_ERR_NOT_SET
      end if
    case ("morph_grid_yll")
      if (.not. allocated(this%morph_grid_yll)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_yll), ubound(this%morph_grid_yll), &
          "morph_grid_yll", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%morph_grid_yll(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%morph_grid_yll))) status = NML_ERR_NOT_SET
      end if
    case ("morph_grid_cellsize")
      if (.not. allocated(this%morph_grid_cellsize)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_cellsize), ubound(this%morph_grid_cellsize), &
          "morph_grid_cellsize", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%morph_grid_cellsize(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%morph_grid_cellsize))) status = NML_ERR_NOT_SET
      end if
    case ("morph_grid_ydir")
      if (.not. allocated(this%morph_grid_ydir)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_ydir), ubound(this%morph_grid_ydir), &
          "morph_grid_ydir", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_grid_coordsys")
      if (.not. allocated(this%morph_grid_coordsys)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_grid_coordsys), ubound(this%morph_grid_coordsys), &
          "morph_grid_coordsys", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pre_coupled")
      if (.not. allocated(this%pre_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pre_coupled), ubound(this%pre_coupled), &
          "pre_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pet_coupled")
      if (.not. allocated(this%pet_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_coupled), ubound(this%pet_coupled), &
          "pet_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("temp_coupled")
      if (.not. allocated(this%temp_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%temp_coupled), ubound(this%temp_coupled), &
          "temp_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tann_coupled")
      if (.not. allocated(this%tann_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tann_coupled), ubound(this%tann_coupled), &
          "tann_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tmin_coupled")
      if (.not. allocated(this%tmin_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmin_coupled), ubound(this%tmin_coupled), &
          "tmin_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("tmax_coupled")
      if (.not. allocated(this%tmax_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%tmax_coupled), ubound(this%tmax_coupled), &
          "tmax_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("ssrd_coupled")
      if (.not. allocated(this%ssrd_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ssrd_coupled), ubound(this%ssrd_coupled), &
          "ssrd_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("strd_coupled")
      if (.not. allocated(this%strd_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%strd_coupled), ubound(this%strd_coupled), &
          "strd_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("netrad_coupled")
      if (.not. allocated(this%netrad_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%netrad_coupled), ubound(this%netrad_coupled), &
          "netrad_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("eabs_coupled")
      if (.not. allocated(this%eabs_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%eabs_coupled), ubound(this%eabs_coupled), &
          "eabs_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("wind_coupled")
      if (.not. allocated(this%wind_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%wind_coupled), ubound(this%wind_coupled), &
          "wind_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("runoff_coupled")
      if (.not. allocated(this%runoff_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_coupled), ubound(this%runoff_coupled), &
          "runoff_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("runoff_sealed_coupled")
      if (.not. allocated(this%runoff_sealed_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%runoff_sealed_coupled), ubound(this%runoff_sealed_coupled), &
          "runoff_sealed_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("interflow_fast_coupled")
      if (.not. allocated(this%interflow_fast_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_fast_coupled), ubound(this%interflow_fast_coupled), &
          "interflow_fast_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("interflow_slow_coupled")
      if (.not. allocated(this%interflow_slow_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%interflow_slow_coupled), ubound(this%interflow_slow_coupled), &
          "interflow_slow_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("baseflow_coupled")
      if (.not. allocated(this%baseflow_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_coupled), ubound(this%baseflow_coupled), &
          "baseflow_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("dem_coupled")
      if (.not. allocated(this%dem_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%dem_coupled), ubound(this%dem_coupled), &
          "dem_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("slope_coupled")
      if (.not. allocated(this%slope_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%slope_coupled), ubound(this%slope_coupled), &
          "slope_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("aspect_coupled")
      if (.not. allocated(this%aspect_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%aspect_coupled), ubound(this%aspect_coupled), &
          "aspect_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("geo_class_coupled")
      if (.not. allocated(this%geo_class_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%geo_class_coupled), ubound(this%geo_class_coupled), &
          "geo_class_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("soil_class_coupled")
      if (.not. allocated(this%soil_class_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%soil_class_coupled), ubound(this%soil_class_coupled), &
          "soil_class_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("lai_class_coupled")
      if (.not. allocated(this%lai_class_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%lai_class_coupled), ubound(this%lai_class_coupled), &
          "lai_class_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("river_width_coupled")
      if (.not. allocated(this%river_width_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%river_width_coupled), ubound(this%river_width_coupled), &
          "river_width_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("meteo_mask_coupled")
      if (.not. allocated(this%meteo_mask_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%meteo_mask_coupled), ubound(this%meteo_mask_coupled), &
          "meteo_mask_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_mask_coupled")
      if (.not. allocated(this%hydro_mask_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_mask_coupled), ubound(this%hydro_mask_coupled), &
          "hydro_mask_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_mask_coupled")
      if (.not. allocated(this%morph_mask_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_mask_coupled), ubound(this%morph_mask_coupled), &
          "morph_mask_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("hydro_latlon_coupled")
      if (.not. allocated(this%hydro_latlon_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%hydro_latlon_coupled), ubound(this%hydro_latlon_coupled), &
          "hydro_latlon_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("morph_latlon_coupled")
      if (.not. allocated(this%morph_latlon_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%morph_latlon_coupled), ubound(this%morph_latlon_coupled), &
          "morph_latlon_coupled", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("route_latlon_coupled")
      if (.not. allocated(this%route_latlon_coupled)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%route_latlon_coupled), ubound(this%route_latlon_coupled), &
          "route_latlon_coupled", errmsg)
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
  end function nml_config_coupling_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_coupling_is_valid(this, errmsg) result(status)
    class(nml_config_coupling_t), intent(in) :: this !< namelist instance
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
    if (allocated(this%meteo_grid_ydir)) then
    if (.not. all(meteo_grid_ydir__in_enum(this%meteo_grid_ydir, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: meteo_grid_ydir"
      return
    end if
    end if
    if (allocated(this%meteo_grid_coordsys)) then
    if (.not. all(meteo_grid_coordsys__in_enum(this%meteo_grid_coordsys, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: meteo_grid_coordsys"
      return
    end if
    end if
    if (allocated(this%hydro_grid_ydir)) then
    if (.not. all(hydro_grid_ydir__in_enum(this%hydro_grid_ydir, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: hydro_grid_ydir"
      return
    end if
    end if
    if (allocated(this%hydro_grid_coordsys)) then
    if (.not. all(hydro_grid_coordsys__in_enum(this%hydro_grid_coordsys, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: hydro_grid_coordsys"
      return
    end if
    end if
    if (allocated(this%morph_grid_ydir)) then
    if (.not. all(morph_grid_ydir__in_enum(this%morph_grid_ydir, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: morph_grid_ydir"
      return
    end if
    end if
    if (allocated(this%morph_grid_coordsys)) then
    if (.not. all(morph_grid_coordsys__in_enum(this%morph_grid_coordsys, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: morph_grid_coordsys"
      return
    end if
    end if
    ! bounds constraints
    if (allocated(this%meteo_grid_nx)) then
    if (.not. all(meteo_grid_nx__in_bounds(this%meteo_grid_nx, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: meteo_grid_nx"
      return
    end if
    end if
    if (allocated(this%meteo_grid_ny)) then
    if (.not. all(meteo_grid_ny__in_bounds(this%meteo_grid_ny, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: meteo_grid_ny"
      return
    end if
    end if
    if (allocated(this%meteo_grid_cellsize)) then
    if (.not. all(meteo_grid_cellsize__in_bounds(this%meteo_grid_cellsize, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: meteo_grid_cellsize"
      return
    end if
    end if
    if (allocated(this%hydro_grid_nx)) then
    if (.not. all(hydro_grid_nx__in_bounds(this%hydro_grid_nx, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: hydro_grid_nx"
      return
    end if
    end if
    if (allocated(this%hydro_grid_ny)) then
    if (.not. all(hydro_grid_ny__in_bounds(this%hydro_grid_ny, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: hydro_grid_ny"
      return
    end if
    end if
    if (allocated(this%hydro_grid_cellsize)) then
    if (.not. all(hydro_grid_cellsize__in_bounds(this%hydro_grid_cellsize, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: hydro_grid_cellsize"
      return
    end if
    end if
    if (allocated(this%morph_grid_nx)) then
    if (.not. all(morph_grid_nx__in_bounds(this%morph_grid_nx, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: morph_grid_nx"
      return
    end if
    end if
    if (allocated(this%morph_grid_ny)) then
    if (.not. all(morph_grid_ny__in_bounds(this%morph_grid_ny, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: morph_grid_ny"
      return
    end if
    end if
    if (allocated(this%morph_grid_cellsize)) then
    if (.not. all(morph_grid_cellsize__in_bounds(this%morph_grid_cellsize, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: morph_grid_cellsize"
      return
    end if
    end if
  end function nml_config_coupling_is_valid

end module nml_config_coupling
