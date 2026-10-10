!> \file mo_mlm_container.F90
!> \brief \copybrief mo_mlm_container
!> \details \copydetails mo_mlm_container

!> \brief Minimal v6 lake component (mHM Lake Module, mLM).
!> \details
!!   Provides the process--1 delayed pass-through contract between mRM updates.
!!
!!   This module is the hub of the mLM v5.13 -> v6 migration: the v5 `mo_mlm_*` process routines
!!   were hand-ported into standalone v6 leaf modules (`src/mLM/mo_lake_*`), and this container
!!   drives them, grown by analogy to `mo_mhm_container`.
!!
!!   Porting principles (kept throughout the migration):
!!   - The leaf routines keep the v5 LOOK (syntax, comment style, NAME/PURPOSE/INTENT/HISTORY
!!     blocks) and their original authorship/history; each ported file carries a new MODIFICATIONS
!!     entry for the migration. The physics/equations are preserved (no logic restructuring).
!!   - Three deltas only: consistent variable/argument nomenclature, SI units internally, and no
!!     global-variable connections (everything passed by argument).
!!   - Units: SI everywhere internally. Input NetCDF keeps conventional reported units (lake area
!!     km2, volume mcm = 1e6 m3); the Input container converts to SI (m2, m3) at the read boundary.
!!   - Scope so far: the minimal natural lake only (no reservoir regulation, demand or SCC).
!!
!> \changelog
!! - Pallav Shrestha (2018-2023): v5.13 mLM and SCC lake/reservoir concepts (origin of the ported physics).
!! - Sebastian Mueller (2026): v6 architecture and the domain-owned mLM pass-through container with its plugin points.
!! - Pallav Shrestha (2026-09/10): v5.13 -> v6 migration of the minimal natural lake.
!!     * Phase 1: topology, footprint delineation and the hAV bathymetry generator.
!!     * Phase 2: natural spill and the iterative lake water balance; case-1 loop and cold-start.
!!     * Phase 3: lake seepage, lake groundwater storage and baseflow; outflow = spill + baseflow.
!!     * Phase 4: Penman (1948) open-water evaporation (daily, disaggregated to the model step).
!!     * Phase 5: lake-point state/flux output and restart persistence (states + static metadata).
!!     * Priestley-Taylor (1972) evaporation as a second method (evaporation_method selector).
!!     * Varying water albedo on projected grids via geographic lake coordinates (level-0 latlon aux).
!> \authors Pallav Shrestha
!> \authors Sebastian Mueller
!> \date 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
#include "logging.h"
module mo_mlm_container
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_kind, only: i4, i8, dp
  use mo_constants, only: nodata_dp
  use mo_datetime, only: datetime, one_hour
  use mo_exchange_type, only: exchange_t
  use mo_grid, only: grid_t
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_lake_balance, only: lake_balance_outflow
  use mo_lake_topology, only: lake_derive_max_surface_area_dem, lake_read_topology, lake_write_topology
  use mo_lake_bathymetry, only: generate_lake_bathymetry, lake_hav_t
  use mo_lake_spill, only: lake_spill
  use mo_lake_water_balance, only: lake_balance, interpolate_lake_states
  use mo_lake_groundwaterstorage, only: lake_seepage, lake_gws_baseflow
  use mo_mlm_mpr, only: lake_derive_ks
  use mo_lake_evaporation, only: lake_evaporation
  use mo_meteo_temporal_tools, only: temporal_disagg_flux_daynight
  use mo_grid_io, only: var, daily, monthly, yearly, no_time, time_units_delta
  use mo_netcdf_utils, only: add_var
  use mo_points, only: points_t, cartesian, spherical
  use mo_points_io, only: points_output_dataset
  use mo_utils, only: is_close
  use mo_string_utils, only: n2s => num2str
  use nml_config_mlm, only: nml_config_mlm_t
  use nml_output_mlm, only: nml_output_mlm_t
  use nml_helper, only: NML_OK
  implicit none

  character(*), parameter :: s = "mlm" !< logging scope

  !> \class mlm_t
  !> \brief Domain-owned lake pass-through component.
  !> \authors Pallav Shrestha
  !> \authors Sebastian Mueller
  type, public :: mlm_t
    type(nml_config_mlm_t) :: config !< mLM configuration
    type(nml_output_mlm_t) :: output_config !< lake-point output configuration
    type(exchange_t), pointer :: exchange => null() !< owning domain exchange
    integer(i8), allocatable :: lake_ids(:) !< stable IDs in point-set order
    integer(i8), allocatable :: lake_map(:) !< stable IDs on packed level-0 lake cells
    real(dp), allocatable :: lake_sand(:) !< top-horizon sand content [%] on packed level-0 lake cells (lake MPR)
    real(dp), allocatable :: lake_clay(:) !< top-horizon clay content [%] on packed level-0 lake cells (lake MPR)
    real(dp), allocatable :: lake_max_surface_area_dem_derived(:) !< maximum lake surface area from DEM footprint [m2]
    type(grid_t) :: static_lake_grid !< lake footprint geometry retained for restart
    type(points_t) :: static_lake_points !< static point metadata retained for restart output
    real(dp), allocatable :: static_lake_max_level_reported(:) !< maximum levels retained for restart output
    real(dp), allocatable :: static_lake_max_surface_area_reported(:) !< reported maximum surface areas [m2] retained for the lake points
    real(dp), allocatable :: static_lake_max_volume_reported(:) !< reported maximum volumes [m3] retained for the lake points
    real(dp), allocatable :: static_lake_spillway_length(:) !< effective spillway crest length [m] retained for the lake points
    real(dp), allocatable :: static_lake_initial_level(:) !< initial lake level [m] retained for the lake points
    real(dp), allocatable :: static_lake_ks(:) !< per-lake regionalized lakebed Ks [cm d-1] (lake MPR; restart-persisted)
    real(dp), allocatable :: static_lake_latitude(:) !< geographic latitude [degree] retained for the lake points
    real(dp), allocatable :: static_lake_longitude(:) !< geographic longitude [degree] retained for the lake points
    type(lake_hav_t), allocatable :: hav(:) !< per-lake height-area-volume curves built at initialize
    real(dp), allocatable :: lake_level(:) !< current lake surface level [masl]
    real(dp), allocatable :: lake_surface_area(:) !< current lake surface area [m2]
    real(dp), allocatable :: lake_volume(:) !< current lake storage volume [m3]
    real(dp), allocatable :: lake_gw_storage(:) !< current lake groundwater storage [m3]
    real(dp), allocatable :: lake_level_previous(:) !< previous-step lake surface level [masl]
    real(dp), allocatable :: lake_surface_area_previous(:) !< previous-step lake surface area [m2]
    real(dp), allocatable :: lake_volume_previous(:) !< previous-step lake storage volume [m3]
    real(dp), allocatable :: lake_gw_storage_previous(:) !< previous-step lake groundwater storage [m3]
    type(points_t) :: restart_lake_points !< point set restored from restart metadata
    integer(i8), allocatable :: restart_lake_ids(:) !< stable IDs in restart-file order
    real(dp), allocatable :: restart_lake_max_level_reported(:) !< maximum levels restored from restart metadata
    real(dp), allocatable :: restart_lake_max_surface_area_reported(:) !< reported maximum surface areas [m2] restored from restart metadata
    real(dp), allocatable :: restart_lake_max_volume_reported(:) !< reported maximum volumes [m3] restored from restart metadata
    real(dp), allocatable :: restart_lake_spillway_length(:) !< effective spillway crest length [m] restored from restart metadata
    real(dp), allocatable :: restart_lake_initial_level(:) !< initial lake level [m] restored from restart metadata
    real(dp), allocatable :: restart_lake_ks(:) !< per-lake lakebed Ks [cm d-1] restored from restart metadata
    real(dp), allocatable :: restart_lake_latitude(:) !< geographic latitude [degree] restored from restart metadata
    real(dp), allocatable :: restart_lake_longitude(:) !< geographic longitude [degree] restored from restart metadata
    real(dp), allocatable :: outflow(:) !< current lake outflow [m3 s-1]
    real(dp), allocatable :: lake_evaporation(:) !< current lake evaporation [mm per step] (diagnostic, for output)
    real(dp), allocatable :: lake_seepage(:) !< current lake seepage/percolation [mm per step] (diagnostic, for output)
    real(dp), allocatable :: lake_baseflow(:) !< current lake baseflow [m3 s-1] (diagnostic, for output)
    real(dp), allocatable :: lake_spill(:) !< current lake crest spill [m3 s-1] (diagnostic, for output; outflow = spill + baseflow)
    logical :: active = .false. !< whether process -1 is selected
    logical :: read_restart = .false. !< read mLM restart
    logical :: write_restart = .false. !< write mLM restart
    logical :: output_active = .false. !< whether lake-point output is enabled
    logical :: evaporation_active = .false. !< whether the Penman lake-evaporation forcings are provided
    logical :: owns_restart_lake_points = .false. !< whether mLM published the restart point set
    logical :: owns_restart_lake_metadata = .false. !< whether mLM published restart IDs and maximum levels
    logical :: owns_restart_lake_grid = .false. !< whether mLM restored and published lake topology
    character(:), allocatable :: restart_input_path !< resolved input restart path
    character(:), allocatable :: restart_output_path !< resolved output restart path
    character(:), allocatable :: output_path !< resolved lake-point output path
    type(points_output_dataset) :: ds_out !< lake-point output dataset
  contains
    procedure :: set_dims => mlm_set_dims
    procedure :: configure => mlm_configure
    procedure :: prepare_restart => mlm_prepare_restart
    procedure :: connect => mlm_connect
    procedure :: initialize => mlm_initialize
    procedure :: update => mlm_update
    procedure :: finalize => mlm_finalize
    procedure :: destroy => mlm_destroy
    procedure, private :: create_restart => mlm_create_restart
    procedure, private :: read_restart_metadata => mlm_read_restart_metadata
    procedure, private :: derive_lake_ks => mlm_derive_lake_ks
    procedure, private :: read_restart_topology => mlm_read_restart_topology
    procedure, private :: read_restart_state => mlm_read_restart_state
    procedure, private :: validate_restart_metadata => mlm_validate_restart_metadata
    procedure, private :: create_output => mlm_create_output
    procedure, private :: validate_output_timing => mlm_validate_output_timing
    procedure, private :: at_output_boundary => mlm_at_output_boundary
    procedure, private :: update_output => mlm_update_output
  end type mlm_t

contains

  !> \brief Load restart-owned lake definition before field connection.
  subroutine mlm_prepare_restart(self)
    class(mlm_t), target, intent(inout) :: self

    if (self%read_restart) call self%read_restart_metadata()
  end subroutine mlm_prepare_restart

  !> \brief Set generated lake-configuration dimensions.
  subroutine mlm_set_dims(self)
    class(mlm_t), target, intent(inout) :: self
    character(1024) :: errmsg
    integer :: status
    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM: error setting config_mlm dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mlm_set_dims

  !> \brief Configure process -1, restart paths, and optional lake-point output.
  subroutine mlm_configure(self, file, out_file)
    class(mlm_t), target, intent(inout) :: self
    character(*), optional, intent(in) :: file
    character(*), optional, intent(in) :: out_file
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer(i4) :: id(1), lake_case
    integer :: status

    lake_case = self%exchange%config%processes%lake
    self%active = lake_case /= 0_i4
    if (lake_case /= 0_i4 .and. lake_case /= -1_i4 .and. lake_case /= -2_i4 .and. lake_case /= 1_i4) then
      log_fatal(*) "mLM: unsupported lake process case: ", n2s(lake_case)
      error stop 1
    end if
    if (.not.self%active) return
    if (self%exchange%config%processes%routing == 0_i4) then
      log_fatal(*) "mLM lake processes require active mRM routing."
      error stop 1
    end if
    id(1) = self%exchange%nml_domain_id
    if (present(file)) then
      path = self%exchange%get_path(file)
      status = self%config%from_file(path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mLM: error reading config_mlm: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%is_configured) then
      log_fatal(*) "mLM configuration not set."
      error stop 1
    end if
    status = self%config%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM config not valid: ", trim(errmsg)
      error stop 1
    end if

    self%output_active = .true.
    if (present(out_file)) then
      path = self%exchange%get_path(out_file, root=.true.)
      log_info(*) "Read mLM output config: ", path
      status = self%output_config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        self%output_active = .false.
        log_warn(*) "mLM output disabled, config not found: ", trim(errmsg)
      end if
    end if
    if (self%output_config%is_configured) then
      status = self%output_config%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mLM output config invalid: ", trim(errmsg)
        error stop 1
      end if
    else
      self%output_active = .false.
      log_warn(*) "mLM output disabled, config not set."
    end if
    if (self%output_active) then
      status = self%config%is_set("output_path", idx=id, errmsg=errmsg)
      self%output_active = status == NML_OK
      if (status /= NML_OK) then
        log_warn(*) "mLM output disabled, path not set for domain ", n2s(id(1)), ": ", trim(errmsg)
      else
        self%output_path = self%exchange%get_path(self%config%output_path(id(1)))
      end if
    end if

    self%read_restart = self%config%read_restart(id(1))
    self%write_restart = self%config%write_restart(id(1))
    if (self%read_restart) then
      if (self%config%is_set("restart_input_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart input path required when read_restart is true."
        error stop 1
      end if
      self%restart_input_path = self%exchange%get_path(self%config%restart_input_path(id(1)))
    end if
    if (self%write_restart) then
      if (self%config%is_set("restart_output_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart output path required when write_restart is true."
        error stop 1
      end if
      self%restart_output_path = self%exchange%get_path(self%config%restart_output_path(id(1)))
    end if
    if (self%read_restart .and. .not.self%exchange%lake_ids%provided) then
      call self%exchange%lake_ids%provide("mLM")
      call self%exchange%lake_max_level_reported%provide("mLM")
      if (lake_case == -2_i4) call self%exchange%lake_map%provide("mLM")
    end if
    if (lake_case == -2_i4) call self%exchange%lake_max_surface_area_dem_derived%provide("mLM")
    if (lake_case == 1_i4) then
      call self%exchange%lake_level%provide("mLM")
      call self%exchange%lake_surface_area%provide("mLM")
      call self%exchange%lake_volume%provide("mLM")
      call self%exchange%lake_gw_storage%provide("mLM")
      ! Read the lake process parameters (lakes_1) and register them in the process parameter
      ! registry (like routing and the MPR land processes). Registration happens here in the
      ! configure phase, before mo_domain seals the registry, so the sealed-time bounds check
      ! applies to the lake parameters and they become visible to the optimizer. All lakes_1
      ! parameters are global (pan-domain) scalars, so no runtime dimensioning is needed.
      associate(parameter_config => self%exchange%config%parameters)
        if (allocated(self%exchange%parameter_file)) then
          status = parameter_config%lakes_1%from_file(file=self%exchange%parameter_file, errmsg=errmsg)
          if (status /= NML_OK) then
            log_fatal(*) "mLM: error reading lakes_1 parameters: ", trim(errmsg)
            error stop 1
          end if
        end if
        status = parameter_config%lakes_1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "mLM: lakes_1 parameters not valid: ", trim(errmsg)
          error stop 1
        end if
        call self%exchange%parameters%add_process("lake", [ &
          parameter_config%lakes_1%spill_coefficient_of_discharge, &
          parameter_config%lakes_1%spill_exponent, &
          parameter_config%lakes_1%b, &
          parameter_config%lakes_1%ptf_ks_constant, &
          parameter_config%lakes_1%ptf_ks_sand, &
          parameter_config%lakes_1%ptf_ks_clay, &
          parameter_config%lakes_1%baseflow_recession_coefficient], [character(64) :: &
          "spill_coefficient_of_discharge", "spill_exponent", "b", &
          "ptf_ks_constant", "ptf_ks_sand", "ptf_ks_clay", &
          "baseflow_recession_coefficient"], group="lakes_1")
      end associate
    end if
    call self%exchange%lake_outflow%provide("mLM")
  end subroutine mlm_configure

  !> \brief Resolve static lake metadata and publish hourly outflow.
  subroutine mlm_connect(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: n_lakes
    integer(i4) :: lake_case
    logical :: has_points, has_ids, has_levels, has_area, has_volume, has_metadata

    lake_case = self%exchange%config%processes%lake
    has_points = associated(self%exchange%lake_points)
    has_ids = self%exchange%lake_ids%provided .and. self%exchange%lake_ids%has_data()
    has_levels = self%exchange%lake_max_level_reported%provided .and. self%exchange%lake_max_level_reported%has_data()
    has_area = self%exchange%lake_max_surface_area_reported%provided .and. self%exchange%lake_max_surface_area_reported%has_data()
    has_volume = self%exchange%lake_max_volume_reported%provided .and. self%exchange%lake_max_volume_reported%has_data()
    has_metadata = has_points .and. has_ids .and. has_levels
    if ((has_points .neqv. has_ids) .or. (has_points .neqv. has_levels)) then
      log_fatal(*) "mLM: lake points, stable IDs, and maximum levels must be published together."
      error stop 1
    end if

    if (has_metadata) then
      n_lakes = self%exchange%lake_points%n_points
      call self%exchange%check_data(self%exchange%lake_ids, "mLM")
      call self%exchange%check_data(self%exchange%lake_max_level_reported, "mLM")
      if (self%read_restart) call self%validate_restart_metadata()
      self%lake_ids = self%exchange%lake_ids%data
      call self%static_lake_points%init( &
        self%exchange%lake_points%x, self%exchange%lake_points%y, coordsys=self%exchange%lake_points%coordsys)
      self%static_lake_max_level_reported = self%exchange%lake_max_level_reported%data
      if (has_area) self%static_lake_max_surface_area_reported = self%exchange%lake_max_surface_area_reported%data
      if (has_volume) self%static_lake_max_volume_reported = self%exchange%lake_max_volume_reported%data
    else
      if (.not.self%read_restart) then
        log_fatal(*) "mLM: lake process requires lake points."
        error stop 1
      end if
      n_lakes = self%restart_lake_points%n_points
      self%lake_ids = self%restart_lake_ids
      call self%static_lake_points%init( &
        self%restart_lake_points%x, self%restart_lake_points%y, coordsys=self%restart_lake_points%coordsys)
      self%static_lake_max_level_reported = self%restart_lake_max_level_reported
      if (allocated(self%restart_lake_max_surface_area_reported)) &
        self%static_lake_max_surface_area_reported = self%restart_lake_max_surface_area_reported
      if (allocated(self%restart_lake_max_volume_reported)) &
        self%static_lake_max_volume_reported = self%restart_lake_max_volume_reported
      if (allocated(self%restart_lake_spillway_length)) &
        self%static_lake_spillway_length = self%restart_lake_spillway_length
      if (allocated(self%restart_lake_initial_level)) &
        self%static_lake_initial_level = self%restart_lake_initial_level
      if (allocated(self%restart_lake_ks)) &
        self%static_lake_ks = self%restart_lake_ks
      if (allocated(self%restart_lake_latitude)) &
        self%static_lake_latitude = self%restart_lake_latitude
      if (allocated(self%restart_lake_longitude)) &
        self%static_lake_longitude = self%restart_lake_longitude
      self%exchange%lake_points => self%static_lake_points
      call self%exchange%lake_ids%prepare_data("mLM", no_time)
      self%exchange%lake_ids%data => self%lake_ids
      call self%exchange%lake_max_level_reported%prepare_data("mLM", no_time)
      self%exchange%lake_max_level_reported%data => self%static_lake_max_level_reported
      self%owns_restart_lake_points = .true.
      self%owns_restart_lake_metadata = .true.
    end if
    if (n_lakes < 1_i8) then
      log_fatal(*) "mLM: lake process requires at least one lake."
      error stop 1
    end if
    if (lake_case == -2_i4) then
      if (associated(self%exchange%level0_lake) .and. self%exchange%lake_map%provided) then
        call self%exchange%check_data(self%exchange%lake_map, "mLM")
        self%static_lake_grid = self%exchange%level0_lake
        self%lake_map = self%exchange%lake_map%data
      else if (self%read_restart) then
        call self%read_restart_topology()
        self%exchange%level0_lake => self%static_lake_grid
        call self%exchange%lake_map%prepare_data("mLM", no_time)
        self%exchange%lake_map%data => self%lake_map
        self%owns_restart_lake_grid = .true.
      else
        log_fatal(*) "mLM process -2 requires a level-0 lake grid and lake map."
        error stop 1
      end if
      call lake_derive_max_surface_area_dem(self%static_lake_grid, self%lake_map, self%lake_ids, self%lake_max_surface_area_dem_derived)
      call self%exchange%lake_max_surface_area_dem_derived%prepare_data("mLM", no_time)
      self%exchange%lake_max_surface_area_dem_derived%data => self%lake_max_surface_area_dem_derived
      call self%exchange%lake_pre%check_provided("mLM")
      call self%exchange%lake_pet%check_provided("mLM")
    end if
    if (lake_case == 1_i4) then
      ! The minimal natural lake needs the reported crest surface area and volume to build
      ! its height-area-volume (hAV) curves at initialize.
      if (.not.allocated(self%static_lake_max_surface_area_reported) .or. &
          .not.allocated(self%static_lake_max_volume_reported)) then
        log_fatal(*) "mLM process 1 requires reported maximum lake surface area and volume."
        error stop 1
      end if
      ! Capture the level-0 lake footprint, its stable lake map and its top-horizon soil class.
      ! The footprint grid + map drive the detached lake-MPR regionalization of the lakebed Ks at
      ! initialize (soil texture over the footprint -> Cosby Ks -> area-weighted mean per lake),
      ! and are retained for restart. On a restart-only run Input provides none of these; the
      ! lakebed Ks is then restored from the restart metadata instead (see the restart branch above).
      if (associated(self%exchange%level0_lake) .and. self%exchange%lake_map%provided) then
        call self%exchange%check_data(self%exchange%lake_map, "mLM")
        self%static_lake_grid = self%exchange%level0_lake
        self%lake_map = self%exchange%lake_map%data
        if (self%exchange%lake_sand%provided .and. self%exchange%lake_sand%has_data() .and. &
            self%exchange%lake_clay%provided .and. self%exchange%lake_clay%has_data()) then
          call self%exchange%check_data(self%exchange%lake_sand, "mLM")
          call self%exchange%check_data(self%exchange%lake_clay, "mLM")
          self%lake_sand = self%exchange%lake_sand%data
          self%lake_clay = self%exchange%lake_clay%data
        end if
        ! (diagnostic) log the DEM-footprint area against the reported area per lake. The DEM
        ! footprint area is not used by the case-1 physics (the hAV curve is reported-anchored),
        ! but this records how far the level-0 lake delineation lies from the reported value.
        block
          real(dp), allocatable :: dem_area(:)
          integer(i8) :: kk
          call lake_derive_max_surface_area_dem(self%static_lake_grid, self%lake_map, &
            self%lake_ids, dem_area)
          do kk = 1_i8, size(self%lake_ids, kind=i8)
            log_info(*) "mLM lake ", n2s(self%lake_ids(kk)), ": reported area [m2] = ", &
              n2s(self%static_lake_max_surface_area_reported(kk)), " , DEM-footprint area [m2] = ", &
              n2s(dem_area(kk)), " ( ratio DEM/reported = ", &
              n2s(dem_area(kk) / self%static_lake_max_surface_area_reported(kk)), " )"
          end do
        end block
      end if
      ! Allocate the lake states carried by the component (built at initialize, restart-persisted in c4).
      allocate(self%lake_level(n_lakes), source=0.0_dp)
      allocate(self%lake_surface_area(n_lakes), source=0.0_dp)
      allocate(self%lake_volume(n_lakes), source=0.0_dp)
      allocate(self%lake_gw_storage(n_lakes), source=0.0_dp)
      call self%exchange%lake_level%prepare_data("mLM", 1_i4)
      self%exchange%lake_level%data => self%lake_level
      call self%exchange%lake_surface_area%prepare_data("mLM", 1_i4)
      self%exchange%lake_surface_area%data => self%lake_surface_area
      call self%exchange%lake_volume%prepare_data("mLM", 1_i4)
      self%exchange%lake_volume%data => self%lake_volume
      call self%exchange%lake_gw_storage%prepare_data("mLM", 1_i4)
      self%exchange%lake_gw_storage%data => self%lake_gw_storage
      ! The minimal natural lake also needs the spillway length (crest spill) and an
      ! initial level (cold-start). Normally the Input container reads them from the lake
      ! definition file; on a restart-only run (no Input) they are restored from the restart
      ! metadata instead (see read_restart_metadata / the restart branch above).
      ! Spillway length is required (from Input or restart metadata).
      if (self%exchange%lake_spillway_length%provided .and. self%exchange%lake_spillway_length%has_data()) then
        call self%exchange%check_data(self%exchange%lake_spillway_length, "mLM")
        self%static_lake_spillway_length = self%exchange%lake_spillway_length%data
      else if (.not.allocated(self%static_lake_spillway_length)) then
        log_fatal(*) "mLM process 1 requires lake spillway length from inputs or restart metadata."
        error stop 1
      end if
      ! The initial level is optional: if it is absent (or carries the missing-value sentinel
      ! nodata_dp = -9999), the lake is cold-started at 50% of its reported maximum volume (see initialize).
      if (self%exchange%lake_initial_level%provided .and. self%exchange%lake_initial_level%has_data()) then
        self%static_lake_initial_level = self%exchange%lake_initial_level%data
      else if (.not.allocated(self%static_lake_initial_level)) then
        allocate(self%static_lake_initial_level(n_lakes), source=nodata_dp)
      end if
      ! Geographic latitude/longitude of the lake points, used by the varying-albedo solar geometry.
      ! The Input container derives them (from the point coordinates when spherical, or the level-0
      ! latlon auxiliary coordinates when projected); on a restart-only run they come from the restart
      ! metadata. Left as nodata when unavailable; varying_albedo then errors in the guard below.
      if (self%exchange%lake_latitude%provided .and. self%exchange%lake_latitude%has_data()) then
        self%static_lake_latitude = self%exchange%lake_latitude%data
      else if (.not.allocated(self%static_lake_latitude)) then
        allocate(self%static_lake_latitude(n_lakes), source=nodata_dp)
      end if
      if (self%exchange%lake_longitude%provided .and. self%exchange%lake_longitude%has_data()) then
        self%static_lake_longitude = self%exchange%lake_longitude%data
      else if (.not.allocated(self%static_lake_longitude)) then
        allocate(self%static_lake_longitude(n_lakes), source=nodata_dp)
      end if
      ! Minimal scope supports natural lakes only; reservoir regulation is not yet ported.
      if (any(.not. self%config%natural(self%exchange%nml_domain_id, 1:int(n_lakes, i4)))) then
        log_fatal(*) "mLM process 1 supports natural lakes only; reservoir regulation is not yet implemented."
        error stop 1
      end if
      ! Lake evaporation is active when its forcings are provided (kept daily by Meteo). All
      ! methods need temperature and downward short- and long-wave radiation; Penman (method 1)
      ! and the calibration-free Priestley-Taylor (method 3) additionally need windspeed and dew-point
      ! temperature, the classic Priestley-Taylor (method 2) does not.
      self%evaporation_active = self%exchange%lake_temp%provided
      if (self%evaporation_active) then
        if (.not.(self%exchange%lake_ssrd%provided .and. self%exchange%lake_strd%provided)) then
          log_fatal(*) "mLM process 1 lake evaporation requires lake temperature, ssrd and strd together."
          error stop 1
        end if
        if ((self%config%evaporation_method == 1_i4 .or. self%config%evaporation_method == 3_i4) .and. &
            .not.(self%exchange%lake_wind%provided .and. self%exchange%lake_tdew%provided)) then
          log_fatal(*) "mLM process 1 Penman (method 1) and calibration-free Priestley-Taylor (method 3) evaporation also require lake wind and dew-point temperature."
          error stop 1
        end if
        ! A varying water albedo needs the geographic latitude/longitude of each lake (solar geometry).
        ! Spherical lake points carry them directly; projected (cartesian) points get them from the
        ! level-0 latlon auxiliary coordinates (latlon_path). Either way the Input/restart-restored
        ! static coordinates must be present and finite, otherwise the albedo would be computed from a
        ! missing-value coordinate. Require them when varying_albedo is on.
        if (self%config%varying_albedo) then
          if (.not.allocated(self%static_lake_latitude) .or. .not.allocated(self%static_lake_longitude)) then
            log_fatal(*) "mLM process 1 varying_albedo requires lake geographic coordinates (set latlon_path for a projected grid)."
            error stop 1
          end if
          if (any(self%static_lake_latitude < -90.0_dp) .or. any(self%static_lake_latitude > 90.0_dp) .or. &
              any(self%static_lake_longitude < -360.0_dp) .or. any(self%static_lake_longitude > 360.0_dp)) then
            log_fatal(*) "mLM process 1 varying_albedo: lake geographic coordinates are missing or out of range (set latlon_path for a projected grid)."
            error stop 1
          end if
        end if
      end if
      ! Previous-step state storage for the water-balance integration.
      allocate(self%lake_level_previous(n_lakes), source=0.0_dp)
      allocate(self%lake_surface_area_previous(n_lakes), source=0.0_dp)
      allocate(self%lake_volume_previous(n_lakes), source=0.0_dp)
      allocate(self%lake_gw_storage_previous(n_lakes), source=0.0_dp)
      ! Diagnostic per-lake fluxes retained for the lake-point output.
      allocate(self%lake_evaporation(n_lakes), source=0.0_dp)
      allocate(self%lake_seepage(n_lakes), source=0.0_dp)
      allocate(self%lake_baseflow(n_lakes), source=0.0_dp)
      allocate(self%lake_spill(n_lakes), source=0.0_dp)
    end if
    allocate(self%outflow(n_lakes), source=0.0_dp)
    call self%exchange%lake_outflow%prepare_data("mLM", 1_i4)
    self%exchange%lake_outflow%data => self%outflow
  end subroutine mlm_connect

  !> \brief Restore or initialize pass-through outflow after mRM has published inflow.
  subroutine mlm_initialize(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: k, n_lakes
    integer(i4) :: bathymetry_type
    integer(i4) :: domain
    real(dp) :: level_query
    real(dp) :: tsHydroSecs
    real(dp) :: leakance, perc0, kbf0
    if (self%exchange%step_hours /= 1_i4) then
      log_fatal(*) "mLM lake processes require a one-hour model step."
      error stop 1
    end if
    call self%exchange%check_data(self%exchange%lake_inflow, "mLM")
    if (self%exchange%lake_inflow%stepping /= 1_i4) then
      log_fatal(*) "mLM: lake inflow must have one-hour support."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == -2_i4) then
      call self%exchange%check_data(self%exchange%lake_max_surface_area_dem_derived, "mLM")
    end if
    if (self%exchange%config%processes%lake == 1_i4) then
      ! Build the per-lake height-area-volume (hAV) curves from a built-in analytical shape.
      n_lakes = size(self%lake_ids, kind=i8)
      domain = self%exchange%nml_domain_id
      tsHydroSecs = real(self%exchange%step_hours, dp) * 3600.0_dp
      allocate(self%hav(n_lakes))
      do k = 1_i8, n_lakes
        bathymetry_type = self%config%bathymetry_type(self%exchange%nml_domain_id, k)
        if (bathymetry_type /= 1_i4 .and. bathymetry_type /= 2_i4) then
          log_fatal(*) "mLM process 1: bathymetry_type ", n2s(bathymetry_type), &
            " for lake ", n2s(k), " is not supported (0 = external bathymetry file is not yet implemented)."
          error stop 1
        end if
        call generate_lake_bathymetry(bathymetry_type, &
          self%static_lake_max_level_reported(k), self%static_lake_max_surface_area_reported(k), &
          self%static_lake_max_volume_reported(k), &
          self%hav(k)%level, self%hav(k)%surface_area, self%hav(k)%volume)
      end do
      ! Regionalize the per-lake lakebed saturated hydraulic conductivity (detached lake MPR):
      ! the Cosby pedotransfer function over the lake footprint soil texture, aggregated to one
      ! value per lake, from the lake's own pedotransfer parameters (lakes_1) and the lake-footprint
      ! soil texture surfaced by Input. On a restart-only run static_lake_ks was restored from the
      ! restart metadata at connect, so only derive it on a fresh run.
      if (.not. self%read_restart) call self%derive_lake_ks()
      ! Cold-start the lake states from the reported initial level; the matching surface
      ! area and volume are read off the per-lake hAV curve. A restart overrides these below.
      if (.not. self%read_restart) then
        do k = 1_i8, n_lakes
          ! When the initial level was not provided (missing-value sentinel), cold-start the
          ! lake at the level holding 50% of its reported maximum volume (v5.13 behaviour):
          ! invert the hAV curve (volume -> level) at half the reported maximum volume.
          if (self%static_lake_initial_level(k) < -9000.0_dp) then
            block
              real(dp) :: target_volume
              target_volume = 0.5_dp * self%static_lake_max_volume_reported(k)
              call interpolate_lake_states(.false., self%hav(k)%volume, self%hav(k)%level, &
                self%static_lake_max_volume_reported(k), self%static_lake_max_level_reported(k), &
                target_volume, self%static_lake_initial_level(k))
            end block
            log_info(*) "mLM lake ", n2s(self%lake_ids(k)), &
              ": initial level not provided; cold-starting at 50% of reported max volume -> level ", &
              n2s(self%static_lake_initial_level(k)), " masl"
          end if
          self%lake_level(k) = self%static_lake_initial_level(k)
          level_query = self%static_lake_initial_level(k)
          call interpolate_lake_states(.false., self%hav(k)%level, self%hav(k)%surface_area, &
            self%static_lake_max_level_reported(k), self%hav(k)%surface_area(size(self%hav(k)%surface_area)), &
            level_query, self%lake_surface_area(k))
          level_query = self%static_lake_initial_level(k)
          call interpolate_lake_states(.false., self%hav(k)%level, self%hav(k)%volume, &
            self%static_lake_max_level_reported(k), self%hav(k)%volume(size(self%hav(k)%volume)), &
            level_query, self%lake_volume(k))
          ! Cold-start the lake groundwater storage self-consistently: the initial baseflow is
          ! set equal to the initial percolation, so the groundwater store starts at its
          ! own equilibrium and there is no spin-up drift. Evaluate the Darcy seepage at the initial
          ! level -> perc0 [mm/TS]; the per-step recession kbf0 = (tsHydroSecs/86400) /
          ! baseflow_recession_coefficient; then gw0 = (perc0/1000 * A0) / kbf0 gives
          ! baseflow0 = kbf0 * gw0 = perc0 volume. leakance = Ks/b, with Ks the per-lake
          ! regionalized lakebed conductivity (lake MPR) and b the global lakebed thickness.
          leakance = (self%static_lake_ks(k) * 0.01_dp &
            * tsHydroSecs / 86400.0_dp) / self%exchange%config%parameters%lakes_1%b%value
          call lake_seepage(.true., self%hav(k)%level, self%hav(k)%surface_area, &
            self%static_lake_max_level_reported(k), &
            self%hav(k)%surface_area(size(self%hav(k)%surface_area)), &
            self%static_lake_initial_level(k), leakance, perc0)
          kbf0 = (tsHydroSecs / 86400.0_dp) &
            / self%exchange%config%parameters%lakes_1%baseflow_recession_coefficient%value
          self%lake_gw_storage(k) = (perc0 / 1000.0_dp * self%lake_surface_area(k)) / kbf0
        end do
      end if
    end if
    self%outflow(:) = 0.0_dp
    if (self%read_restart) call self%read_restart_state()
    if (self%exchange%config%processes%lake == 1_i4) then
      ! Seed the previous-step states for the water-balance integration.
      self%lake_level_previous = self%lake_level
      self%lake_surface_area_previous = self%lake_surface_area
      self%lake_volume_previous = self%lake_volume
      self%lake_gw_storage_previous = self%lake_gw_storage
    end if
    call self%validate_output_timing()
    call self%create_output()
  end subroutine mlm_initialize

  !> \brief Regionalize the per-lake lakebed saturated hydraulic conductivity (detached lake MPR).
  !> \details Runs the Cosby (1984) pedotransfer function over the lake footprint soil texture and
  !!          aggregates it to one lakebed Ks per lake by an area-weighted arithmetic mean, mirroring
  !!          the detached mRM regionalization (own parameters applied to a level-0 geo-field at
  !!          initialize). Self-contained: the Cosby pedotransfer parameters are lake parameters
  !!          (lakes_1) owned by the lake module, and the lake footprint sand/clay content is surfaced
  !!          by the Input container from the soil LUT. It does not depend on the land-surface
  !!          soil-moisture process or the MPR soil database. The derived field is used in the Darcy
  !!          leakance Ks/b of the lake seepage.
  subroutine mlm_derive_lake_ks(self)
    class(mlm_t), target, intent(inout) :: self
    real(dp) :: ptf_ks(4)  ! Cosby Ks pedotransfer parameters (constant, sand, clay, unit conversion)
    real(dp), parameter :: inch_per_hour_to_centimetre_per_day = 60.96_dp ! Cosby unit conversion

    ! The lakebed Ks regionalization is self-contained: the Cosby pedotransfer parameters are lake
    ! parameters (lakes_1) owned by the lake module, and the lake footprint soil texture (sand/clay)
    ! is surfaced by the Input container from the soil LUT. It does NOT depend on the land-surface
    ! soil-moisture process or the MPR soil database. The lake runs independently of mHM, exactly as
    ! mRM does.
    if (.not.allocated(self%lake_sand) .or. .not.allocated(self%lake_clay)) then
      log_fatal(*) "mLM process 1 lakebed conductivity requires the lake footprint soil texture &
        &(set soil_class_path / soil_horizon_class_path and soil_lut_path so Input can read it)."
      error stop 1
    end if

    ! Cosby et al. (1984) pedotransfer parameters from the lake parameter set, plus the fixed
    ! inch/h -> cm/d unit conversion as the fourth element.
    ptf_ks(1) = self%exchange%config%parameters%lakes_1%ptf_ks_constant%value
    ptf_ks(2) = self%exchange%config%parameters%lakes_1%ptf_ks_sand%value
    ptf_ks(3) = self%exchange%config%parameters%lakes_1%ptf_ks_clay%value
    ptf_ks(4) = inch_per_hour_to_centimetre_per_day

    call lake_derive_ks(self%static_lake_grid, self%lake_map, self%lake_ids, &
      self%lake_sand, self%lake_clay, ptf_ks, self%static_lake_ks)

    block
      integer(i8) :: kk
      do kk = 1_i8, size(self%lake_ids, kind=i8)
        log_info(*) "mLM lake ", n2s(self%lake_ids(kk)), ": regionalized lakebed Ks [cm d-1] = ", &
          n2s(self%static_lake_ks(kk))
      end do
    end block

  end subroutine mlm_derive_lake_ks

  !> \brief Publish the preceding completed hourly inflow as current outflow.
  subroutine mlm_update(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: k, n_lakes
    integer(i4) :: domain
    integer(i4) :: hour, doy
    real(dp) :: spill, volume_adjustment, tsHydroSecs, perc, baseflow, leakance, kbf
    real(dp) :: ntimesteps_day, evap_daily, evap_step, pr_step, qin_d, qout_d, pre_d, lat, lon
    real(dp) :: wind_k, tdew_k
    real(dp) :: evap_components(3)
    logical :: isday
    real(dp) :: input_windspeed_height, climatological_windspeed
    integer(i4) :: albedo_method
    logical :: varying_albedo, varying_lambda, varying_psychro
    select case (self%exchange%config%processes%lake)
      case(-1_i4)
        self%outflow = self%exchange%lake_inflow%data
      case(-2_i4)
      ! if (self%exchange%config%processes%lake == -2_i4) then
        self%outflow = lake_balance_outflow(self%exchange%lake_inflow%data, self%lake_max_surface_area_dem_derived, self%exchange%lake_pre%data, &
          self%exchange%lake_pet%data, 3600.0_dp)
      case(1_i4)
        ! Minimal natural lake: seepage and crest spill, the iterative water balance, then
        ! the groundwater-storage baseflow. Phase 3 scope -> evaporation, water use and
        ! inter-basin transfer are still zero; outflow is the crest spill plus the baseflow.
        n_lakes = size(self%lake_ids, kind=i8)
        domain = self%exchange%nml_domain_id
        tsHydroSecs = real(self%exchange%step_hours, dp) * 3600.0_dp
        ! Scalars for the daily -> step evaporation disaggregation (0.5/0.5 day-night split; the
        ! diurnal signal is carried by the hourly-varying albedo inside the Penman call).
        hour = self%exchange%time_step_start%hour
        doy = self%exchange%time_step_start%doy()
        ntimesteps_day = 24.0_dp / real(self%exchange%step_hours, dp)
        isday = (hour >= 7_i4) .and. (hour < 19_i4)
        input_windspeed_height = self%config%input_windspeed_height
        climatological_windspeed = self%config%climatological_windspeed
        albedo_method = self%config%albedo_method
        varying_albedo = self%config%varying_albedo
        varying_lambda = self%config%varying_lambda
        varying_psychro = self%config%varying_psychro
        !$omp parallel do default(shared) &
        !$omp   private(k, spill, volume_adjustment, perc, baseflow, leakance, kbf, &
        !$omp           evap_daily, evap_step, pr_step, qin_d, qout_d, pre_d, lat, lon, wind_k, tdew_k, evap_components) &
        !$omp   schedule(dynamic)
        do k = 1_i8, n_lakes
          ! Darcy lakebed leakance Ks/b [TS-1 per unit head] and per-step baseflow recession.
          ! leakance replaces the former kp (recession rate) + exponent_for_percolation
          ! (depth-normaliser). Ks is the per-lake regionalized lakebed conductivity (lake MPR) and b
          ! is the global pan-domain lakebed thickness; Ks is cm/d, converted to m/TS then divided by b [m].
          leakance = (self%static_lake_ks(k) * 0.01_dp &
            * tsHydroSecs / 86400.0_dp) / self%exchange%config%parameters%lakes_1%b%value
          kbf = (tsHydroSecs / 86400.0_dp) / &
            self%exchange%config%parameters%lakes_1%baseflow_recession_coefficient%value
          ! Lake seepage from the previous-step level -> percolation to the groundwater store.
          call lake_seepage(.true., self%hav(k)%level, self%hav(k)%surface_area, &
            self%static_lake_max_level_reported(k), &
            self%hav(k)%surface_area(size(self%hav(k)%surface_area)), &
            self%lake_level_previous(k), leakance, &
            perc)
          ! Crest-controlled spill from the previous-step level.
          call lake_spill(self%lake_level_previous(k), self%static_lake_max_level_reported(k), &
            self%exchange%config%parameters%lakes_1%spill_coefficient_of_discharge%value, &
            self%static_lake_spillway_length(k), &
            self%exchange%config%parameters%lakes_1%spill_exponent%value, &
            spill)
          ! Lake precipitation onto the surface (already disaggregated to the step by Meteo).
          pr_step = 0.0_dp
          if (self%exchange%lake_pre%provided) pr_step = self%exchange%lake_pre%data(k)
          ! Lake evaporation (Penman, daily) -> disaggregated to the step. Zero when its forcings
          ! are not provided (keeps earlier-phase runs unchanged).
          !
          ! Design note - why the Penman CALL lives here in mLM and not in mo_meteo_container
          ! (which computes and disaggregates the LAND PET, incl. Penman-Monteith, and publishes it):
          ! land PET is a pure function of meteorological inputs and static parameters, so it belongs
          ! with the meteo handler. Open-water lake evaporation is not: it depends on the lake's own
          ! evolving state - the surface level feeds the atmospheric-pressure / psychrometric-constant
          ! term (varying_psychro), and the advected-heat term uses the inflow/outflow rates over the
          ! lake surface area and the previous-step outflow. Those states (lake_level,
          ! lake_surface_area_previous, previous outflow) are owned by mLM. Computing evaporation here
          ! uses them directly and in the correct water-balance sequence; doing it in meteo would need
          ! meteo to read lagged lake states back from the exchange (a meteo<->mLM mutual dependency)
          ! or to drop the advected-heat term. So the forcing AGGREGATION stays in meteo (the lake
          ! temp/ssrd/strd/wind/tdew fields), but the Penman evaluation + daily->step disaggregation
          ! stay with the lake component that owns the state. (Revisit if/when the advection term is
          ! dropped or river-temperature coupling moves that state elsewhere.)
          evap_step = 0.0_dp
          if (self%evaporation_active) then
            ! Advection rates [mm d-1] over the lake surface (previous-step area and outflow);
            ! Tin/Tout/Tpre are set to the lake air temperature until river-temperature coupling.
            if (self%lake_surface_area_previous(k) > 0.0_dp) then
              qin_d  = self%exchange%lake_inflow%data(k) / self%lake_surface_area_previous(k) * 1000.0_dp * 86400.0_dp
              qout_d = self%outflow(k) / self%lake_surface_area_previous(k) * 1000.0_dp * 86400.0_dp
            else
              qin_d  = 0.0_dp
              qout_d = 0.0_dp
            end if
            pre_d = pr_step * ntimesteps_day
            ! Geographic coordinates for the solar geometry (used only when varying_albedo); these are
            ! true lat/lon even on a projected grid (from the level-0 latlon auxiliary coordinates),
            ! unlike static_lake_points%x/y which hold projected eastings/northings on a cartesian grid.
            lat = self%static_lake_latitude(k)
            lon = self%static_lake_longitude(k)
            ! Windspeed and dew-point temperature are used by Penman only; they may be absent for
            ! Priestley-Taylor, so read them only when provided (zero otherwise, where they are unused).
            wind_k = 0.0_dp
            tdew_k = 0.0_dp
            if (self%exchange%lake_wind%provided) wind_k = self%exchange%lake_wind%data(k)
            if (self%exchange%lake_tdew%provided) tdew_k = self%exchange%lake_tdew%data(k)
            call lake_evaporation(self%config%evaporation_method, hour, doy, &
              self%exchange%lake_temp%data(k), tdew_k, wind_k, self%exchange%lake_ssrd%data(k), &
              self%exchange%lake_strd%data(k), &
              self%exchange%lake_temp%data(k), self%exchange%lake_temp%data(k), self%exchange%lake_temp%data(k), &
              self%lake_level_previous(k), qin_d, qout_d, pre_d, lat, lon, &
              input_windspeed_height, climatological_windspeed, varying_albedo, albedo_method, varying_lambda, varying_psychro, &
              evap_daily, evap_components)
            call temporal_disagg_flux_daynight(isday, ntimesteps_day, evap_daily, 0.5_dp, 0.5_dp, evap_step)
          end if
          ! Iterative water balance -> updates level, surface area and volume (consumes perc, evap, spill).
          call lake_balance(.true., pr_step, evap_step, self%exchange%lake_inflow%data(k), perc, 0.0_dp, spill, 0.0_dp, &
            self%hav(k)%level, self%hav(k)%surface_area, self%hav(k)%volume, &
            self%static_lake_max_level_reported(k), &
            self%hav(k)%surface_area(size(self%hav(k)%surface_area)), &
            self%hav(k)%volume(size(self%hav(k)%volume)), &
            self%lake_level_previous(k), self%lake_surface_area_previous(k), self%lake_volume_previous(k), &
            self%lake_level(k), self%lake_surface_area(k), self%lake_volume(k), &
            volume_adjustment, tsHydroSecs)
          ! Lake groundwater storage update -> baseflow (uses the previous-step surface area and GW storage).
          call lake_gws_baseflow(perc, self%lake_surface_area_previous(k), &
            self%lake_gw_storage_previous(k), self%lake_gw_storage(k), &
            kbf, tsHydroSecs, baseflow)
          ! Natural-lake outflow is the crest spill plus the groundwater baseflow.
          self%outflow(k) = spill + baseflow
          ! Retain the per-lake fluxes for the lake-point output.
          self%lake_evaporation(k) = evap_step
          self%lake_seepage(k) = perc
          self%lake_baseflow(k) = baseflow
          self%lake_spill(k) = spill
        end do
        !$omp end parallel do
        ! Roll the current states into the previous-step slots for the next update.
        self%lake_level_previous = self%lake_level
        self%lake_surface_area_previous = self%lake_surface_area
        self%lake_volume_previous = self%lake_volume
        self%lake_gw_storage_previous = self%lake_gw_storage
    end select
    call self%update_output()
  end subroutine mlm_update

  !> \brief Write optional restart and close the lake output before destroy.
  subroutine mlm_finalize(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%write_restart) call self%create_restart()
    if (self%output_active) then
      call self%ds_out%close()
      log_info(*) "Close mLM output file: ", self%output_path
    else
      log_info(*) "No mLM output file will be written"
    end if
  end subroutine mlm_finalize

  !> \brief Release lake model definition and state after exchange publications are detached.
  subroutine mlm_destroy(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%owns_restart_lake_metadata) then
      self%owns_restart_lake_metadata = .false.
    end if
    if (self%owns_restart_lake_points) then
      self%owns_restart_lake_points = .false.
    end if
    if (self%owns_restart_lake_grid) then
      self%owns_restart_lake_grid = .false.
    end if
    if (allocated(self%lake_ids)) deallocate(self%lake_ids)
    if (allocated(self%lake_map)) deallocate(self%lake_map)
    if (allocated(self%lake_sand)) deallocate(self%lake_sand)
    if (allocated(self%lake_clay)) deallocate(self%lake_clay)
    if (allocated(self%lake_max_surface_area_dem_derived)) deallocate(self%lake_max_surface_area_dem_derived)
    call self%static_lake_points%destroy()
    if (allocated(self%static_lake_max_level_reported)) deallocate(self%static_lake_max_level_reported)
    if (allocated(self%static_lake_max_surface_area_reported)) deallocate(self%static_lake_max_surface_area_reported)
    if (allocated(self%static_lake_max_volume_reported)) deallocate(self%static_lake_max_volume_reported)
    if (allocated(self%static_lake_spillway_length)) deallocate(self%static_lake_spillway_length)
    if (allocated(self%static_lake_initial_level)) deallocate(self%static_lake_initial_level)
    if (allocated(self%static_lake_ks)) deallocate(self%static_lake_ks)
    if (allocated(self%static_lake_latitude)) deallocate(self%static_lake_latitude)
    if (allocated(self%static_lake_longitude)) deallocate(self%static_lake_longitude)
    if (allocated(self%lake_level)) deallocate(self%lake_level)
    if (allocated(self%lake_surface_area)) deallocate(self%lake_surface_area)
    if (allocated(self%lake_volume)) deallocate(self%lake_volume)
    if (allocated(self%lake_gw_storage)) deallocate(self%lake_gw_storage)
    if (allocated(self%lake_level_previous)) deallocate(self%lake_level_previous)
    if (allocated(self%lake_surface_area_previous)) deallocate(self%lake_surface_area_previous)
    if (allocated(self%lake_volume_previous)) deallocate(self%lake_volume_previous)
    if (allocated(self%lake_gw_storage_previous)) deallocate(self%lake_gw_storage_previous)
    if (allocated(self%lake_evaporation)) deallocate(self%lake_evaporation)
    if (allocated(self%lake_seepage)) deallocate(self%lake_seepage)
    if (allocated(self%lake_baseflow)) deallocate(self%lake_baseflow)
    if (allocated(self%lake_spill)) deallocate(self%lake_spill)
    if (allocated(self%hav)) deallocate(self%hav)
    call self%restart_lake_points%destroy()
    if (allocated(self%restart_lake_ids)) deallocate(self%restart_lake_ids)
    if (allocated(self%restart_lake_max_level_reported)) deallocate(self%restart_lake_max_level_reported)
    if (allocated(self%restart_lake_max_surface_area_reported)) deallocate(self%restart_lake_max_surface_area_reported)
    if (allocated(self%restart_lake_max_volume_reported)) deallocate(self%restart_lake_max_volume_reported)
    if (allocated(self%restart_lake_spillway_length)) deallocate(self%restart_lake_spillway_length)
    if (allocated(self%restart_lake_initial_level)) deallocate(self%restart_lake_initial_level)
    if (allocated(self%restart_lake_ks)) deallocate(self%restart_lake_ks)
    if (allocated(self%restart_lake_latitude)) deallocate(self%restart_lake_latitude)
    if (allocated(self%restart_lake_longitude)) deallocate(self%restart_lake_longitude)
    if (allocated(self%outflow)) deallocate(self%outflow)
    call self%static_lake_grid%destroy()
  end subroutine mlm_destroy

  !> \brief Create the lake-point output dataset in stable exchange order.
  subroutine mlm_create_output(self)
    class(mlm_t), target, intent(inout) :: self
    type(var), allocatable :: vars(:)
    character(:), allocatable :: delta, dtype
    integer(i4) :: timestamp

    if (.not.self%output_active) return

    timestamp = self%output_config%output_time_reference
    delta = time_units_delta(self%output_config%output_frequency, timestamp)
    dtype = "f64"
    if (.not.self%output_config%output_double_precision) dtype = "f32"

    call add_var(vars, var(name="lake_id", long_name="stable lake ID", dtype="i64", kind="i8", static=.true.))
    ! Lake states and the component flux diagnostics exist only for the minimal natural lake
    ! (process 1); the states are published exchange fields (metadata via as_output_var) while
    ! evaporation/seepage/baseflow are component members output as their own vars.
    if (self%exchange%config%processes%lake == 1_i4) then
      if (self%output_config%out_lake_level) &
        call add_var(vars, self%exchange%lake_level%as_output_var(dtype=dtype, avg=.true.))
      ! Area and the two volumes are converted from SI to conventional units at the output
      ! boundary (km2, mcm = 1e6 m3) for readability and to mirror the lake-definition input
      ! units; the internal states stay SI (m2, m3). Level stays in m a.s.l.
      if (self%output_config%out_lake_surface_area) &
        call add_var(vars, var(name="lake_surface_area", long_name="current lake surface area", &
          units="km2", dtype=dtype, avg=.true.))
      if (self%output_config%out_lake_volume) &
        call add_var(vars, var(name="lake_volume", long_name="current lake volume", &
          units="1e6 m3", dtype=dtype, avg=.true.))
      if (self%output_config%out_lake_gw_storage) &
        call add_var(vars, var(name="lake_gw_storage", long_name="current lake groundwater storage", &
          units="1e6 m3", dtype=dtype, avg=.true.))
      if (self%output_config%out_lake_evaporation) &
        call add_var(vars, var(name="lake_evaporation", long_name="lake evaporation", units="mm", dtype=dtype))
      if (self%output_config%out_lake_seepage) &
        call add_var(vars, var(name="lake_seepage", long_name="lake seepage to groundwater storage", units="mm", dtype=dtype))
      if (self%output_config%out_lake_baseflow) &
        call add_var(vars, var(name="lake_baseflow", long_name="lake baseflow from groundwater storage", units="m3 s-1", dtype=dtype, avg=.true.))
      if (self%output_config%out_lake_spill) &
        call add_var(vars, var(name="lake_spill", long_name="lake crest spill", units="m3 s-1", dtype=dtype, avg=.true.))
    end if
    ! Lake inflow (exchange field) and precipitation (the lake_pre field, when provided) are
    ! available beyond the minimal natural lake.
    if (self%output_config%out_lake_inflow) &
      call add_var(vars, self%exchange%lake_inflow%as_output_var(name="lake_inflow", long_name="lake inflow", dtype=dtype, avg=.true.))
    if (self%output_config%out_lake_precipitation .and. self%exchange%lake_pre%provided) &
      call add_var(vars, self%exchange%lake_pre%as_output_var(name="lake_precipitation", &
        long_name="lake precipitation", dtype=dtype))
    if (self%output_config%out_lake_outflow) then
      call add_var(vars, self%exchange%lake_outflow%as_output_var(dtype=dtype, avg=.true.))
    end if

    log_info(*) "Create mLM lake-point output file: ", self%output_path
    call self%ds_out%init( &
      path          = self%output_path, &
      points        = self%exchange%lake_points, &
      vars          = vars, &
      start_time    = self%exchange%start_time, &
      delta         = delta, &
      timestamp     = timestamp, &
      deflate_level = self%output_config%output_deflate_level, &
      point_dim_name = "lake", &
      time_series   = self%output_config%output_timeseries)
    call self%ds_out%update("lake_id", self%lake_ids)
    call self%ds_out%write_static()
  end subroutine mlm_create_output

  !> \brief Validate mLM output cadence and restart/output boundary alignment.
  subroutine mlm_validate_output_timing(self)
    class(mlm_t), target, intent(in) :: self
    integer(i4) :: frequency

    if (.not.self%output_active) return
    frequency = self%output_config%output_frequency
    select case (frequency)
      case (daily, monthly, yearly, no_time)
        continue
      case default
        if (frequency < 1_i4 .or. mod(frequency, self%exchange%step_hours) /= 0_i4) then
          log_fatal(*) "mLM output_frequency=", frequency, &
            "h must be a positive whole multiple of model_step=", self%exchange%step_hours, "h."
          error stop 1
        end if
    end select

    if (self%write_restart .and. .not.self%at_output_boundary(self%exchange%end_time)) then
      log_fatal(*) "mLM restart output time must coincide with an mLM output boundary."
      error stop 1
    end if
  end subroutine mlm_validate_output_timing

  !> \brief Report whether a timestamp is a configured mLM output boundary.
  logical function mlm_at_output_boundary(self, time) result(boundary)
    class(mlm_t), target, intent(in) :: self
    type(datetime), intent(in) :: time
    integer(i4) :: elapsed_hours

    select case (self%output_config%output_frequency)
      case (daily)
        boundary = time%is_new_day()
      case (monthly)
        boundary = time%is_new_month()
      case (yearly)
        boundary = time%is_new_year()
      case (no_time)
        boundary = time == self%exchange%end_time
      case default
        elapsed_hours = nint((time - self%exchange%start_time) / one_hour(), i4)
        boundary = elapsed_hours >= 0_i4 .and. &
          mod(elapsed_hours, self%output_config%output_frequency) == 0_i4
    end select
  end function mlm_at_output_boundary

  !> \brief Buffer the current lake release and write completed output intervals.
  subroutine mlm_update_output(self)
    class(mlm_t), target, intent(inout) :: self

    if (.not.self%output_active) return
    if (self%exchange%config%processes%lake == 1_i4) then
      if (self%output_config%out_lake_level) call self%ds_out%update("lake_level", self%lake_level)
      if (self%output_config%out_lake_surface_area) &
        call self%ds_out%update("lake_surface_area", self%lake_surface_area * 1.0e-6_dp)
      if (self%output_config%out_lake_volume) &
        call self%ds_out%update("lake_volume", self%lake_volume * 1.0e-6_dp)
      if (self%output_config%out_lake_gw_storage) &
        call self%ds_out%update("lake_gw_storage", self%lake_gw_storage * 1.0e-6_dp)
      if (self%output_config%out_lake_evaporation) call self%ds_out%update("lake_evaporation", self%lake_evaporation)
      if (self%output_config%out_lake_seepage) call self%ds_out%update("lake_seepage", self%lake_seepage)
      if (self%output_config%out_lake_baseflow) call self%ds_out%update("lake_baseflow", self%lake_baseflow)
      if (self%output_config%out_lake_spill) call self%ds_out%update("lake_spill", self%lake_spill)
    end if
    if (self%output_config%out_lake_inflow) call self%ds_out%update("lake_inflow", self%exchange%lake_inflow%data)
    if (self%output_config%out_lake_precipitation .and. self%exchange%lake_pre%provided) &
      call self%ds_out%update("lake_precipitation", self%exchange%lake_pre%data)
    if (self%output_config%out_lake_outflow) call self%ds_out%update("lake_outflow", self%outflow)
    if (self%at_output_boundary(self%exchange%time)) call self%ds_out%write(self%exchange%time)
  end subroutine mlm_update_output

  !> \brief Persist static lake metadata and current outflow with stable IDs.
  subroutine mlm_create_restart(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcDimension) :: dims(0), lake_dim
    type(NcVariable) :: var

    if (self%static_lake_points%n_points /= size(self%lake_ids, kind=i8) .or. &
        .not.allocated(self%static_lake_max_level_reported)) then
      log_fatal(*) "mLM: cannot write restart without static lake metadata."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == 1_i4 .and. &
        (.not.allocated(self%static_lake_max_surface_area_reported) .or. &
         .not.allocated(self%static_lake_max_volume_reported))) then
      log_fatal(*) "mLM: cannot write process 1 restart without reported maximum surface area and volume."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == 1_i4 .and. &
        (.not.allocated(self%static_lake_spillway_length) .or. &
         .not.allocated(self%static_lake_initial_level))) then
      log_fatal(*) "mLM: cannot write process 1 restart without spillway length and initial level."
      error stop 1
    end if
    nc = NcDataset(self%restart_output_path, "w")
    lake_dim = nc%setDimension("lake", int(size(self%lake_ids), i4))
    var = nc%setVariable("lake_id", "i64", [lake_dim])
    call var%setData(self%lake_ids)
    var = nc%setVariable("lake_outlet_x", "f64", [lake_dim])
    call var%setData(self%static_lake_points%x)
    var = nc%setVariable("lake_outlet_y", "f64", [lake_dim])
    call var%setData(self%static_lake_points%y)
    var = nc%setVariable("lake_max_level_reported", "f64", [lake_dim])
    call var%setAttribute("units", "m")
    call var%setData(self%static_lake_max_level_reported)
    if (allocated(self%static_lake_max_surface_area_reported)) then
      var = nc%setVariable("lake_max_surface_area_reported", "f64", [lake_dim])
      call var%setAttribute("units", "m2")
      call var%setData(self%static_lake_max_surface_area_reported)
    end if
    if (allocated(self%static_lake_max_volume_reported)) then
      var = nc%setVariable("lake_max_volume_reported", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%static_lake_max_volume_reported)
    end if
    if (allocated(self%static_lake_spillway_length)) then
      var = nc%setVariable("lake_spillway_length", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%static_lake_spillway_length)
    end if
    if (allocated(self%static_lake_initial_level)) then
      var = nc%setVariable("lake_initial_level", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%static_lake_initial_level)
    end if
    ! The regionalized lakebed Ks is written so a restart-only run (no Input, no soil database)
    ! can restore it instead of re-deriving it from soil.
    if (allocated(self%static_lake_ks)) then
      var = nc%setVariable("lake_ks", "f64", [lake_dim])
      call var%setAttribute("units", "cm d-1")
      call var%setData(self%static_lake_ks)
    end if
    ! Geographic coordinates are written only when they were derived (i.e. a latlon file was present
    ! for a projected grid, or the points were spherical); they are required at read time only when
    ! varying_albedo is configured (checked in the connect guard).
    if (allocated(self%static_lake_latitude)) then
      var = nc%setVariable("lake_latitude", "f64", [lake_dim])
      call var%setAttribute("units", "degree")
      call var%setData(self%static_lake_latitude)
    end if
    if (allocated(self%static_lake_longitude)) then
      var = nc%setVariable("lake_longitude", "f64", [lake_dim])
      call var%setAttribute("units", "degree")
      call var%setData(self%static_lake_longitude)
    end if
    var = nc%setVariable("lake_coordsys", "i32", dims(:0))
    call var%setData(self%static_lake_points%coordsys)
    var = nc%setVariable("lake_outflow", "f64", [lake_dim])
    call var%setAttribute("units", "m3 s-1")
    call var%setData(self%outflow)
    if (self%exchange%config%processes%lake == 1_i4) then
      var = nc%setVariable("lake_level", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%lake_level)
      var = nc%setVariable("lake_surface_area", "f64", [lake_dim])
      call var%setAttribute("units", "m2")
      call var%setData(self%lake_surface_area)
      var = nc%setVariable("lake_volume", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%lake_volume)
      var = nc%setVariable("lake_gw_storage", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%lake_gw_storage)
    end if
    if (self%exchange%config%processes%lake == -2_i4) then
      call lake_write_topology(nc, self%static_lake_grid, self%lake_map)
    end if
    var = nc%setVariable("mlm_meta", "i32", dims(:0))
    call var%setData(0_i4)
    call var%setAttribute("time_stamp", self%exchange%time%str())
    call nc%close()
  end subroutine mlm_create_restart

  !> \brief Read and validate static lake metadata before resolving the exchange provider.
  subroutine mlm_read_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: var, meta_var
    real(dp), allocatable :: outlet_x(:), outlet_y(:)
    integer(i4), allocatable :: coordsys_data
    integer(i8) :: i
    character(64) :: restart_time
    character(:), allocatable :: expected_time

    nc = NcDataset(self%restart_input_path, "r")
    if (.not.nc%hasVariable("lake_id") .or. .not.nc%hasVariable("lake_outflow")) then
      log_fatal(*) "mLM restart is missing lake IDs or lake outflow: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_id")
    call var%getData(self%restart_lake_ids)
    if (size(self%restart_lake_ids, kind=i8) < 1_i8) then
      log_fatal(*) "mLM restart contains no lakes: ", self%restart_input_path
      error stop 1
    end if
    if (any(self%restart_lake_ids <= 0_i8)) then
      log_fatal(*) "mLM restart contains a non-positive stable lake ID."
      error stop 1
    end if
    do i = 1_i8, size(self%restart_lake_ids, kind=i8) - 1_i8
      if (any(self%restart_lake_ids(i + 1_i8:) == self%restart_lake_ids(i))) then
        log_fatal(*) "mLM restart contains duplicate stable lake IDs."
        error stop 1
      end if
    end do

    if (.not.nc%hasVariable("lake_outlet_x") .or. .not.nc%hasVariable("lake_outlet_y") .or. &
        .not.nc%hasVariable("lake_max_level_reported") .or. .not.nc%hasVariable("lake_coordsys")) then
      log_fatal(*) "mLM restart is missing required static lake metadata: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_outlet_x")
    call var%getData(outlet_x)
    var = nc%getVariable("lake_outlet_y")
    call var%getData(outlet_y)
    var = nc%getVariable("lake_max_level_reported")
    call var%getData(self%restart_lake_max_level_reported)
    if (self%exchange%config%processes%lake == 1_i4) then
      if (.not.nc%hasVariable("lake_max_surface_area_reported") .or. &
          .not.nc%hasVariable("lake_max_volume_reported")) then
        log_fatal(*) "mLM process 1 restart is missing reported maximum surface area or volume: ", self%restart_input_path
        error stop 1
      end if
      if (.not.nc%hasVariable("lake_spillway_length") .or. &
          .not.nc%hasVariable("lake_initial_level")) then
        log_fatal(*) "mLM process 1 restart is missing spillway length or initial level: ", self%restart_input_path
        error stop 1
      end if
    end if
    if (nc%hasVariable("lake_max_surface_area_reported")) then
      var = nc%getVariable("lake_max_surface_area_reported")
      call var%getData(self%restart_lake_max_surface_area_reported)
    end if
    if (nc%hasVariable("lake_max_volume_reported")) then
      var = nc%getVariable("lake_max_volume_reported")
      call var%getData(self%restart_lake_max_volume_reported)
    end if
    if (nc%hasVariable("lake_spillway_length")) then
      var = nc%getVariable("lake_spillway_length")
      call var%getData(self%restart_lake_spillway_length)
    end if
    if (nc%hasVariable("lake_initial_level")) then
      var = nc%getVariable("lake_initial_level")
      call var%getData(self%restart_lake_initial_level)
    end if
    if (nc%hasVariable("lake_ks")) then
      var = nc%getVariable("lake_ks")
      call var%getData(self%restart_lake_ks)
    end if
    ! Geographic coordinates are optional in the restart (back-compat and non-albedo runs); when a
    ! varying_albedo run restarts without them, the connect guard reports the missing coordinates.
    if (nc%hasVariable("lake_latitude")) then
      var = nc%getVariable("lake_latitude")
      call var%getData(self%restart_lake_latitude)
    end if
    if (nc%hasVariable("lake_longitude")) then
      var = nc%getVariable("lake_longitude")
      call var%getData(self%restart_lake_longitude)
    end if
    var = nc%getVariable("lake_coordsys")
    call var%getData(coordsys_data)
    if (size(outlet_x, kind=i8) /= size(self%restart_lake_ids, kind=i8) .or. &
        size(outlet_y, kind=i8) /= size(self%restart_lake_ids, kind=i8) .or. &
        size(self%restart_lake_max_level_reported, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
      log_fatal(*) "mLM restart static lake metadata does not align with stable lake IDs."
      error stop 1
    end if
    if (allocated(self%restart_lake_max_surface_area_reported)) then
      if (size(self%restart_lake_max_surface_area_reported, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart reported maximum surface area does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_spillway_length)) then
      if (size(self%restart_lake_spillway_length, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart spillway length does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_initial_level)) then
      if (size(self%restart_lake_initial_level, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart initial level does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_ks)) then
      if (size(self%restart_lake_ks, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart lakebed Ks does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_latitude)) then
      if (size(self%restart_lake_latitude, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart latitude does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_longitude)) then
      if (size(self%restart_lake_longitude, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart longitude does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (allocated(self%restart_lake_max_volume_reported)) then
      if (size(self%restart_lake_max_volume_reported, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart reported maximum volume does not align with stable lake IDs."
        error stop 1
      end if
    end if
    if (coordsys_data /= cartesian .and. coordsys_data /= spherical) then
      log_fatal(*) "mLM restart contains an unsupported lake coordinate system."
      error stop 1
    end if
    if (.not.all(ieee_is_finite(outlet_x)) .or. .not.all(ieee_is_finite(outlet_y)) .or. &
        .not.all(ieee_is_finite(self%restart_lake_max_level_reported))) then
      log_fatal(*) "mLM restart static lake metadata must be finite."
      error stop 1
    end if
    call self%restart_lake_points%init(outlet_x, outlet_y, coordsys=coordsys_data)

    if (.not.nc%hasVariable("mlm_meta")) then
      log_fatal(*) "mLM restart metadata variable mlm_meta is missing: ", self%restart_input_path
      error stop 1
    end if
    meta_var = nc%getVariable("mlm_meta")
    if (.not.meta_var%hasAttribute("time_stamp")) then
      log_fatal(*) "mLM restart metadata has no time_stamp: ", self%restart_input_path
      error stop 1
    end if
    call meta_var%getAttribute("time_stamp", restart_time)
    expected_time = self%exchange%start_time%str()
    if (trim(restart_time) /= trim(expected_time)) then
      log_fatal(*) "mLM restart timestamp ", trim(restart_time), " does not match domain restart time ", trim(expected_time), "."
      error stop 1
    end if
    call nc%close()
  end subroutine mlm_read_restart_metadata

  !> \brief Restore the static level-0 lake footprint required by lake case -2.
  subroutine mlm_read_restart_topology(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    nc = NcDataset(self%restart_input_path, "r")
    call lake_read_topology(nc, self%static_lake_grid, self%lake_map)
    call nc%close()
  end subroutine mlm_read_restart_topology

  !> \brief Verify restart metadata against Input metadata without imposing restart-file ordering.
  subroutine mlm_validate_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: i, j, n_lakes

    n_lakes = self%exchange%lake_points%n_points
    if (size(self%restart_lake_ids, kind=i8) /= n_lakes) then
      log_fatal(*) "mLM restart lake IDs do not match Input lake metadata."
      error stop 1
    end if
    if (self%restart_lake_points%coordsys /= self%exchange%lake_points%coordsys) then
      log_fatal(*) "mLM restart and Input lake coordinate systems differ."
      error stop 1
    end if
    do i = 1_i8, n_lakes - 1_i8
      if (any(self%exchange%lake_ids%data(i + 1_i8:) == self%exchange%lake_ids%data(i))) then
        log_fatal(*) "Input lake metadata contains duplicate stable lake IDs."
        error stop 1
      end if
    end do
    do i = 1_i8, n_lakes
      do j = 1_i8, size(self%restart_lake_ids, kind=i8)
        if (self%restart_lake_ids(j) == self%exchange%lake_ids%data(i)) exit
      end do
      if (j > size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart and Input lake metadata have different stable lake IDs."
        error stop 1
      end if
      if (.not.is_close(self%restart_lake_points%x(j), self%exchange%lake_points%x(i)) .or. &
          .not.is_close(self%restart_lake_points%y(j), self%exchange%lake_points%y(i)) .or. &
          .not.is_close(self%restart_lake_max_level_reported(j), self%exchange%lake_max_level_reported%data(i))) then
        log_fatal(*) "mLM restart and Input lake metadata differ for stable lake ID."
        error stop 1
      end if
      if (allocated(self%restart_lake_max_surface_area_reported) .and. &
          self%exchange%lake_max_surface_area_reported%has_data()) then
        if (.not.is_close(self%restart_lake_max_surface_area_reported(j), &
                          self%exchange%lake_max_surface_area_reported%data(i))) then
          log_fatal(*) "mLM restart and Input reported maximum surface area differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart_lake_max_volume_reported) .and. &
          self%exchange%lake_max_volume_reported%has_data()) then
        if (.not.is_close(self%restart_lake_max_volume_reported(j), &
                          self%exchange%lake_max_volume_reported%data(i))) then
          log_fatal(*) "mLM restart and Input reported maximum volume differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart_lake_spillway_length) .and. &
          self%exchange%lake_spillway_length%has_data()) then
        if (.not.is_close(self%restart_lake_spillway_length(j), &
                          self%exchange%lake_spillway_length%data(i))) then
          log_fatal(*) "mLM restart and Input spillway length differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart_lake_initial_level) .and. &
          self%exchange%lake_initial_level%has_data()) then
        if (.not.is_close(self%restart_lake_initial_level(j), &
                          self%exchange%lake_initial_level%data(i))) then
          log_fatal(*) "mLM restart and Input initial level differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart_lake_latitude) .and. &
          self%exchange%lake_latitude%has_data()) then
        if (.not.is_close(self%restart_lake_latitude(j), &
                          self%exchange%lake_latitude%data(i))) then
          log_fatal(*) "mLM restart and Input latitude differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart_lake_longitude) .and. &
          self%exchange%lake_longitude%has_data()) then
        if (.not.is_close(self%restart_lake_longitude(j), &
                          self%exchange%lake_longitude%data(i))) then
          log_fatal(*) "mLM restart and Input longitude differ for stable lake ID."
          error stop 1
        end if
      end if
    end do
  end subroutine mlm_validate_restart_metadata

  !> \brief Restore outflow by stable ID rather than restart-file order.
  subroutine mlm_read_restart_state(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: var
    integer(i8), allocatable :: ids(:)
    real(dp), allocatable :: values(:), level(:), area(:), volume(:), gw_storage(:)
    integer(i8) :: i, j
    logical :: is_case1
    is_case1 = self%exchange%config%processes%lake == 1_i4
    nc = NcDataset(self%restart_input_path, "r")
    if (.not.nc%hasVariable("lake_id") .or. .not.nc%hasVariable("lake_outflow")) then
      log_fatal(*) "mLM restart is missing lake IDs or lake outflow: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_id")
    call var%getData(ids)
    var = nc%getVariable("lake_outflow")
    call var%getData(values)
    if (is_case1) then
      if (.not.nc%hasVariable("lake_level") .or. .not.nc%hasVariable("lake_surface_area") .or. &
          .not.nc%hasVariable("lake_volume") .or. .not.nc%hasVariable("lake_gw_storage")) then
        log_fatal(*) "mLM process 1 restart is missing one or more lake states: ", self%restart_input_path
        error stop 1
      end if
      var = nc%getVariable("lake_level")
      call var%getData(level)
      var = nc%getVariable("lake_surface_area")
      call var%getData(area)
      var = nc%getVariable("lake_volume")
      call var%getData(volume)
      var = nc%getVariable("lake_gw_storage")
      call var%getData(gw_storage)
    end if
    call nc%close()
    if (size(ids, kind=i8) /= size(self%lake_ids, kind=i8) .or. size(values, kind=i8) /= size(ids, kind=i8)) then
      log_fatal(*) "mLM restart lake count does not match configured lakes."
      error stop 1
    end if
    if (is_case1) then
      if (size(level, kind=i8) /= size(ids, kind=i8) .or. size(area, kind=i8) /= size(ids, kind=i8) .or. &
          size(volume, kind=i8) /= size(ids, kind=i8) .or. size(gw_storage, kind=i8) /= size(ids, kind=i8)) then
        log_fatal(*) "mLM restart lake states do not align with stable lake IDs."
        error stop 1
      end if
    end if
    do i = 1_i8, size(ids, kind=i8) - 1_i8
      if (any(ids(i + 1_i8:) == ids(i))) then
        log_fatal(*) "mLM restart contains duplicate stable lake IDs."
        error stop 1
      end if
    end do
    do i = 1_i8, size(self%lake_ids, kind=i8)
      do j = 1_i8, size(ids, kind=i8)
        if (ids(j) == self%lake_ids(i)) exit
      end do
      if (j > size(ids, kind=i8)) then
        log_fatal(*) "mLM restart is missing stable lake ID."
        error stop 1
      end if
      self%outflow(i) = values(j)
      if (is_case1) then
        self%lake_level(i) = level(j)
        self%lake_surface_area(i) = area(j)
        self%lake_volume(i) = volume(j)
        self%lake_gw_storage(i) = gw_storage(j)
      end if
    end do
  end subroutine mlm_read_restart_state

end module mo_mlm_container
