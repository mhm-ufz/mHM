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
!!     * Evaporation: Priestley-Taylor (1972) evaporation as a second method (evaporation_method selector).
!!     * Evaporation: Varying water albedo on projected grids via geographic lake coordinates (level-0 latlon aux).
!!     * Evaporation: Calibration-free Priestley-Taylor (Su & Singh, 2023) evaporation as a third method.
!!     * Percolatoin: Physical lakebed percolation: Darcy leakance Ks/b, linear in head (replaced the lumped
!!       percolation recession rate + depth-normaliser and the H^2 nonlinearity); self-consistent
!!       groundwater cold-start (initial baseflow = initial percolation; retired longterm_baseflow).
!!     * Parameter bounds: lakes_1 parameters registered in the process parameter registry (bounds-checked + optimizer-visible).
!!     * Percolation: Detached lake MPR (lakebed Ks): per-lake saturated hydraulic conductivity regionalized from the
!!       lake-footprint soil texture via the Cosby (1984) pedotransfer function (mo_mlm_mpr); borrows the
!!       land soil-moisture pedotransfer parameters when that process is active, else uses lakes_1's own.
!!     * Baseflow: Detached lake MPR (baseflow recession coefficient): per-lake k_bf regionalized from the
!!       lake-footprint geology (per-geological-unit lookup + footprint mean, mo_mlm_mpr); borrows the
!!       land baseflow parameters when that process is active, else uses lakes_1's own per-unit array.
!!     * Lake banks: reconcile the fixed DEM-footprint collection area against the dynamic water-balance
!!       area each step; the net atmospheric flux over the ring (A_footprint - A_water) enters as a
!!       same-step inflow correction. Exposed bank (A_bank > 0): adds (P - Es) with soil evaporation
!!       Es = min(c_soil * E, P) - the bank has no soil-moisture store, so Es consumes only same-step
!!       P and (P - Es) >= 0. Over-flooded ring (A_bank < 0): removes (P - E) with open-water E, the
!!       flux already counted on the mHM land side. Closes the budget over the true footprint; the
!!       DEM-footprint area is retained and restart-persisted (lake_footprint_area).
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
  use mo_lake_forcing, only: lake_grids_match
  use mo_lake_topology, only: lake_derive_spill_reference_surface_area_dem, lake_read_topology, lake_write_topology
  use mo_lake_bathymetry, only: generate_lake_bathymetry, lake_hav_t
  use mo_lake_spill, only: lake_spill
  use mo_lake_water_balance, only: lake_balance, interpolate_lake_states
  use mo_lake_groundwaterstorage, only: lake_seepage, lake_gws_baseflow
  use mo_input_lut, only: read_geoformation_lut
  use mo_mlm_mpr, only: lake_derive_ks, lake_derive_kbf, lake_derive_spillway_length, lake_check_geo_units
  use mo_lake_evaporation, only: lake_evaporation
  use mo_mlm_constants, only: bank_soil_evap_ratio
  use mo_meteo_temporal_tools, only: temporal_disagg_flux_daynight
  use mo_grid_io, only: var, daily, monthly, yearly, no_time, time_units_delta
  use mo_netcdf_utils, only: add_var
  use mo_points, only: points_t, cartesian, spherical
  use mo_points_io, only: points_output_dataset
  use mo_utils, only: is_close
  use mo_string_utils, only: n2s => num2str
  use mo_parameter_types, only: parameter_t
  use nml_config_mlm, only: nml_config_mlm_t
  use nml_output_mlm, only: nml_output_mlm_t
  use nml_helper, only: NML_OK
  implicit none
  private

  character(*), parameter :: s = "mlm" !< logging scope

  !> \class mlm_definition_t
  !> \brief Domain-owned lake definition data.
  type :: mlm_definition_t
    integer(i8), allocatable :: ids(:) !< stable IDs in point-set order
    integer(i8), allocatable :: map(:) !< stable IDs on packed level-0 lake cells
    real(dp), allocatable :: spill_reference_surface_area_dem_derived(:) !< maximum lake surface area from DEM footprint [m2]
    type(grid_t) :: grid !< lake footprint geometry retained for restart
    type(points_t) :: points !< static point metadata retained for restart output
    real(dp), allocatable :: spill_reference_level(:) !< maximum levels retained for restart output
    real(dp), allocatable :: spill_reference_surface_area(:) !< reported maximum surface areas [m2] retained for the lake points
    real(dp), allocatable :: spill_reference_volume(:) !< reported maximum volumes [m3] retained for the lake points
    real(dp), allocatable :: spillway_length(:) !< effective spillway crest length [m] retained for the lake points
    logical, allocatable :: endorheic(:) !< per-lake terminal-sink flag: no downstream river, so spill/percolation/baseflow are all suppressed (derived from topology; restart-persisted)
    real(dp), allocatable :: initial_level(:) !< initial lake level [m] retained for the lake points
    real(dp), allocatable :: ks(:) !< per-lake regionalized lakebed Ks [cm d-1] (lake MPR; restart-persisted)
    real(dp), allocatable :: kbf(:) !< per-lake regionalized baseflow recession coefficient [d] (lake MPR; restart-persisted)
    real(dp), allocatable :: footprint_area(:) !< per-lake DEM-footprint surface area [m2] (lake banks; restart-persisted)
    real(dp), allocatable :: latitude(:) !< geographic latitude [degree] retained for the lake points
    real(dp), allocatable :: longitude(:) !< geographic longitude [degree] retained for the lake points
  end type mlm_definition_t

  !> \class mlm_regionalization_t
  !> \brief Domain-owned lake regionalization data.
  type :: mlm_regionalization_t
    real(dp), allocatable :: sand(:) !< top-horizon sand content [%] on packed level-0 lake cells (lake MPR)
    real(dp), allocatable :: clay(:) !< top-horizon clay content [%] on packed level-0 lake cells (lake MPR)
    integer(i4), allocatable :: geo(:) !< geological unit ID on packed level-0 lake cells (lake MPR)
    integer(i4), allocatable :: geo_units(:) !< geological units read during connection
    integer(i4), allocatable :: geo_indices(:) !< geological parameter mapping
    integer(i8), allocatable :: lake_index(:) !< cached point-set slot per footprint cell
    integer(i4), allocatable :: geo_parameter(:) !< cached geological parameter slot per cell
    real(dp), allocatable :: catchment(:) !< lake catchment areas [m2]
  end type mlm_regionalization_t

  !> \class mlm_runtime_t
  !> \brief Domain-owned lake runtime data.
  type :: mlm_runtime_t
    real(dp) :: cd = 0.0_dp !< spill discharge coefficient
    real(dp) :: exponent = 0.0_dp !< spill exponent
    real(dp) :: multiplier = 0.0_dp !< catchment width multiplier
    real(dp) :: b = 0.0_dp !< lakebed thickness [m]
    real(dp) :: ptf(4) = 0.0_dp !< Cosby coefficients and unit conversion
    real(dp), allocatable :: recession(:) !< current effective geological coefficients [d]
  end type mlm_runtime_t

  !> \class mlm_state_t
  !> \brief Domain-owned lake state data.
  type :: mlm_state_t
    type(lake_hav_t), allocatable :: hav(:) !< per-lake height-area-volume curves built at initialize
    real(dp), allocatable :: level(:) !< current lake surface level [masl]
    real(dp), allocatable :: surface_area(:) !< current lake surface area [m2]
    real(dp), allocatable :: volume(:) !< current lake storage volume [m3]
    real(dp), allocatable :: gw_storage(:) !< current lake groundwater storage [m3]
    real(dp), allocatable :: level_previous(:) !< previous-step lake surface level [masl]
    real(dp), allocatable :: surface_area_previous(:) !< previous-step lake surface area [m2]
    real(dp), allocatable :: volume_previous(:) !< previous-step lake storage volume [m3]
    real(dp), allocatable :: gw_storage_previous(:) !< previous-step lake groundwater storage [m3]
    real(dp), allocatable :: outflow(:) !< current lake outflow [m3 s-1]
    real(dp), allocatable :: evaporation(:) !< current lake evaporation [mm per step] (diagnostic, for output)
    real(dp), allocatable :: seepage(:) !< current lake seepage/percolation [mm per step] (diagnostic, for output)
    real(dp), allocatable :: baseflow(:) !< current lake baseflow [m3 s-1] (diagnostic, for output)
    real(dp), allocatable :: spill(:) !< current lake crest spill [m3 s-1] (diagnostic, for output; outflow = spill + baseflow)
  end type mlm_state_t

  !> \class mlm_restart_t
  !> \brief Domain-owned lake restart data.
  type :: mlm_restart_t
    real(dp), allocatable :: generating_ks(:) !< restored effective Cosby coefficients
    real(dp), allocatable :: generating_kbf(:) !< restored effective geological coefficients
    real(dp), allocatable :: generating_length(:) !< restored effective spillway multiplier
    type(points_t) :: points !< point set restored from restart metadata
    integer(i8), allocatable :: ids(:) !< stable IDs in restart-file order
    real(dp), allocatable :: spill_reference_level(:) !< maximum levels restored from restart metadata
    real(dp), allocatable :: spill_reference_surface_area(:) !< reported maximum surface areas [m2] restored from restart metadata
    real(dp), allocatable :: spill_reference_volume(:) !< reported maximum volumes [m3] restored from restart metadata
    real(dp), allocatable :: spillway_length(:) !< effective spillway crest length [m] restored from restart metadata
    logical, allocatable :: endorheic(:) !< per-lake terminal-sink flag restored from restart metadata
    real(dp), allocatable :: initial_level(:) !< initial lake level [m] restored from restart metadata
    real(dp), allocatable :: ks(:) !< per-lake lakebed Ks [cm d-1] restored from restart metadata
    real(dp), allocatable :: kbf(:) !< per-lake baseflow recession coefficient [d] restored from restart metadata
    real(dp), allocatable :: footprint_area(:) !< per-lake DEM-footprint surface area [m2] restored from restart metadata
    real(dp), allocatable :: latitude(:) !< geographic latitude [degree] restored from restart metadata
    real(dp), allocatable :: longitude(:) !< geographic longitude [degree] restored from restart metadata
  end type mlm_restart_t

  !> \class mlm_t
  !> \brief Domain-owned natural lake definition, state and lifecycle.
  type, public :: mlm_t
    type(nml_config_mlm_t) :: config !< mLM configuration
    type(nml_output_mlm_t) :: output_config !< lake-point output configuration
    type(exchange_t), pointer :: exchange => null() !< owning domain exchange
    logical :: active = .false. !< whether process -1 is selected
    logical :: read_restart = .false. !< read mLM restart
    logical :: restart_save_topology = .false. !< include optional level-0 geometry
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
    type(mlm_definition_t) :: definition !< cohesive definition storage
    type(mlm_regionalization_t) :: regionalization !< cohesive regionalization storage
    type(mlm_runtime_t) :: runtime !< cohesive runtime storage
    type(mlm_state_t) :: state !< cohesive state storage
    type(mlm_restart_t) :: restart !< cohesive restart storage

  contains
    procedure :: set_dims => mlm_set_dims
    procedure :: configure => mlm_configure
    procedure :: prepare_restart => mlm_prepare_restart
    procedure :: connect => mlm_connect
    procedure :: initialize => mlm_initialize
    procedure :: update => mlm_update
    procedure :: finalize => mlm_finalize
    procedure :: destroy => mlm_destroy
    procedure, private :: cache_regionalization => mlm_cache_regionalization
    procedure, private :: fetch_parameters => mlm_fetch_parameters
    procedure, private :: create_restart => mlm_create_restart
    procedure, private :: read_restart_metadata => mlm_read_restart_metadata
    procedure, private :: derive_lake_ks => mlm_derive_lake_ks
    procedure, private :: derive_lake_kbf => mlm_derive_lake_kbf
    procedure, private :: derive_lake_spillway_length => mlm_derive_lake_spillway_length
    procedure, private :: derive_lake_endorheic => mlm_derive_lake_endorheic
    procedure, private :: read_restart_topology => mlm_read_restart_topology
    procedure, private :: read_restart_state => mlm_read_restart_state
    procedure, private :: validate_restart_metadata => mlm_validate_restart_metadata
    procedure, private :: create_output => mlm_create_output
    procedure, private :: validate_output_timing => mlm_validate_output_timing
    procedure, private :: at_output_boundary => mlm_at_output_boundary
    procedure, private :: update_output => mlm_update_output
  end type mlm_t

contains

  !> \brief Cache parameter-independent regionalization inputs once during connection.
  subroutine mlm_cache_regionalization(self)
    use mo_river, only: river_t
    class(mlm_t), target, intent(inout) :: self
    type(river_t), pointer :: river
    integer(i4) :: n_classes, class_slot
    integer(i4), allocatable :: geo_karstic(:)
    integer(i8) :: n, cell
    character(:), allocatable :: path
    if (allocated(self%definition%map)) then
      allocate(self%regionalization%lake_index(size(self%definition%map)))
      do cell = 1_i8,size(self%definition%map,kind=i8)
        self%regionalization%lake_index(cell) = findloc(self%definition%ids,self%definition%map(cell),dim=1,kind=i8)
      end do
      if (any(self%regionalization%lake_index < 1_i8)) then
        log_fatal(*) "mLM: lake map has unknown stable IDs."
        error stop 1
      end if
    end if
    if (.not.allocated(self%regionalization%geo)) then
      log_fatal(*) "mLM process 1 baseflow requires lake geology and geo_lut_path."
      error stop 1
    end if
    if (self%exchange%config%lut%is_set("geo_lut_path", idx=[self%exchange%nml_domain_id]) /= NML_OK) then
      log_fatal(*) "mLM process 1 baseflow requires geo_lut_path."
      error stop 1
    end if
    path = self%exchange%get_path(self%exchange%config%lut%geo_lut_path(self%exchange%nml_domain_id))
    call read_geoformation_lut(path, nGeo=n_classes, geo_unit=self%regionalization%geo_units, &
      geo_karstic=geo_karstic, geo_param_index=self%regionalization%geo_indices)
    call lake_check_geo_units(self%regionalization%geo, self%regionalization%geo_units)
    allocate(self%regionalization%geo_parameter(size(self%regionalization%geo)))
    do cell = 1_i8,size(self%regionalization%geo,kind=i8)
      if (count(self%regionalization%geo_units == self%regionalization%geo(cell)) /= 1) then
        log_fatal(*) "mLM: ambiguous geological unit in lake LUT."
        error stop 1
      end if
      class_slot = findloc(self%regionalization%geo_units,self%regionalization%geo(cell),dim=1)
      self%regionalization%geo_parameter(cell) = self%regionalization%geo_indices(class_slot)
    end do
    if (.not.associated(self%exchange%river_l0)) then
      log_fatal(*) "mLM process 1 spillway length requires the level-0 river."
      error stop 1
    end if
    river => self%exchange%river_l0
    if (.not.allocated(river%lake_outlet_nodes)) then
      log_fatal(*) "mLM process 1 spillway length requires lake outlet nodes."
      error stop 1
    end if
    n = size(self%definition%ids, kind=i8)
    if (size(river%lake_outlet_nodes, kind=i8) /= n) then
      log_fatal(*) "mLM: lake outlet count does not match stable IDs."
      error stop 1
    end if
    if (.not.allocated(river%upstream_area)) call river%calc_upstream_area()
    self%regionalization%catchment = river%upstream_area(river%lake_outlet_nodes)
  end subroutine mlm_cache_regionalization

  !> \brief Fetch the conditional registry layout on every initialization.
  subroutine mlm_fetch_parameters(self)
    class(mlm_t), target, intent(inout) :: self
    real(dp), allocatable :: values(:), soil(:)
    integer(i4) :: offset
    values = self%exchange%parameters%get_process("lake")
    if (size(values) < 4_i4) then
      log_fatal(*) "mLM: incomplete lake runtime parameter vector."
      error stop 1
    end if
    self%runtime%cd = values(1)
    self%runtime%exponent = values(2)
    self%runtime%multiplier = values(3)
    self%runtime%b = values(4)
    offset = 4_i4
    if (self%exchange%config%processes%soil_moisture /= 0_i4) then
      soil = self%exchange%parameters%get_process("soil_moisture")
      if (size(soil) < 12_i4) then
        log_fatal(*) "mLM: soil_moisture parameters must include Cosby coefficients at positions 10-12."
        error stop 1
      end if
      self%runtime%ptf(1:3) = soil(10:12)
    else
      if (size(values) < offset + 3_i4) then
        log_fatal(*) "mLM: missing lake-owned Cosby coefficients."
        error stop 1
      end if
      self%runtime%ptf(1:3) = values(offset+1:offset+3)
      offset = offset + 3_i4
    end if
    self%runtime%ptf(4) = 60.96_dp ! inch h-1 -> cm d-1; same conversion as land Cosby
    if (self%exchange%config%processes%baseflow /= 0_i4) then
      self%runtime%recession = self%exchange%parameters%get_process("baseflow")
    else
      self%runtime%recession = values(offset+1:)
    end if
    if (.not.self%read_restart) return
    call mlm_check_generating_parameters(self%restart%generating_ks, self%runtime%ptf(1:3))
    call mlm_check_generating_parameters(self%restart%generating_kbf, self%runtime%recession)
    call mlm_check_generating_parameters(self%restart%generating_length, [self%runtime%multiplier])
  end subroutine mlm_fetch_parameters

  !> \brief Reject changed or malformed generating parameter records.
  subroutine mlm_check_generating_parameters(saved, current)
    real(dp), intent(in) :: saved(:) !< effective generating values in restart
    real(dp), intent(in) :: current(:) !< current effective runtime values
    if (size(saved) /= size(current)) call mlm_changed_restart_parameters()
    if (.not.all(ieee_is_finite(saved))) call mlm_changed_restart_parameters()
    if (any(abs(saved-current) > 1.0e-12_dp * max(1.0_dp, abs(current)))) &
      call mlm_changed_restart_parameters()
  end subroutine mlm_check_generating_parameters

  !> \brief Diagnose an incompatible restart parameter override.
  subroutine mlm_changed_restart_parameters()
    log_fatal(*) "mLM: parameters generating restored Ks, recession or spillway length cannot be changed."
    error stop 1
  end subroutine mlm_changed_restart_parameters

  !> \brief Load restart-owned lake definition before field connection.
  subroutine mlm_prepare_restart(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i4) :: lake_case

    if (.not.self%read_restart) return
    call self%read_restart_metadata()
    ! Publish restart-owned definition before meteo connects. Keep the restored targets
    ! separate from the copies mLM uses internally; Input remains the owner when selected.
    if (self%exchange%lake_ids%provider == "mLM") then
      self%exchange%lake_points => self%restart%points
      call self%exchange%lake_ids%prepare_data("mLM", no_time)
      self%exchange%lake_ids%data => self%restart%ids
      call self%exchange%lake_spill_reference_level%prepare_data("mLM", no_time)
      self%exchange%lake_spill_reference_level%data => self%restart%spill_reference_level
      self%owns_restart_lake_points = .true.
      self%owns_restart_lake_metadata = .true.
      lake_case = self%exchange%config%processes%lake
      if (lake_case == 1_i4 .or. lake_case == -2_i4) then
        call self%read_restart_topology()
        if (allocated(self%definition%map)) then
          if (associated(self%exchange%level0_lake)) then
            if (.not.lake_grids_match(self%definition%grid,self%exchange%level0_lake)) then
              log_fatal(*) "mLM: restored lake topology conflicts with the supplied level-0 lake grid."
              error stop 1
            end if
          else
            self%exchange%level0_lake => self%definition%grid
          end if
          call self%exchange%lake_map%prepare_data("mLM", no_time)
          self%exchange%lake_map%data => self%definition%map
          self%owns_restart_lake_grid = .true.
        end if
      end if
    end if
  end subroutine mlm_prepare_restart

  !> \brief Set generated lake-configuration dimensions.
  subroutine mlm_set_dims(self)
    class(mlm_t), target, intent(inout) :: self
    character(1024) :: errmsg
    integer :: status
    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, max_lakes=self%exchange%max_lakes, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM: error setting config_mlm dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mlm_set_dims

  !> \brief Configure lake processes, restart paths, and optional lake-point output.
  subroutine mlm_configure(self, file, out_file)
    class(mlm_t), target, intent(inout) :: self
    character(*), optional, intent(in) :: file !< component configuration namelist path
    character(*), optional, intent(in) :: out_file !< lake output configuration namelist path
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
    self%restart_save_topology = self%config%restart_save_topology(id(1))
    if (self%read_restart) then
      if (self%config%is_set("restart_input_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart input path required when read_restart is true."
        error stop 1
      end if
      if (len_trim(self%config%restart_input_path(id(1))) == 0) then
        log_fatal(*) "mLM: restart input path must not be empty."
        error stop 1
      end if
      self%restart_input_path = self%exchange%get_path(self%config%restart_input_path(id(1)))
    end if
    if (self%write_restart) then
      if (self%config%is_set("restart_output_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart output path required when write_restart is true."
        error stop 1
      end if
      if (len_trim(self%config%restart_output_path(id(1))) == 0) then
        log_fatal(*) "mLM: restart output path must not be empty."
        error stop 1
      end if
      self%restart_output_path = self%exchange%get_path(self%config%restart_output_path(id(1)))
    end if
    if (self%read_restart .and. .not.self%exchange%lake_ids%provided) then
      call self%exchange%lake_ids%provide("mLM")
      call self%exchange%lake_spill_reference_level%provide("mLM")
      if (lake_case == -2_i4 .or. lake_case == 1_i4) call self%exchange%lake_map%provide("mLM")
    end if
    if (lake_case == -2_i4) call self%exchange%lake_spill_reference_surface_area_dem_derived%provide("mLM")
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
        if (.not.parameter_config%lakes_1%is_configured .and. allocated(self%exchange%parameter_file)) then
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
        ! Register the lake parameters with the optimizer. Two parameter sub-sets are conditional on the
        ! land-surface processes the lake borrows from, so a lake copy is not duplicated in the
        ! calibration space:
        !  - the Cosby pedotransfer parameters (ptf_ks_*) are registered ONLY when the land-surface
        !    soil-moisture process is inactive (otherwise the lakebed Ks borrows the land pedotransfer
        !    parameters, see mlm_derive_lake_ks);
        !  - the per-geological-unit baseflow recession array (baseflow_recession) is registered ONLY
        !    when the land-surface baseflow (geoparameter) process is inactive (otherwise the lake
        !    baseflow recession coefficient borrows the land baseflow parameters, see mlm_derive_lake_kbf).
        ! The spill parameters and the lakebed thickness b are always lake-owned and always registered.
        block
          type(parameter_t), allocatable :: lake_params(:)
          character(64), allocatable :: lake_names(:)
          ! always-registered lake-owned parameters
          lake_params = [ &
            parameter_config%lakes_1%spill_coefficient_of_discharge, &
            parameter_config%lakes_1%spill_exponent, &
            parameter_config%lakes_1%spillway_length_multiplier, &
            parameter_config%lakes_1%b]
          lake_names = [character(64) :: "spill_coefficient_of_discharge", "spill_exponent", &
            "spillway_length_multiplier", "b"]
          ! lakebed Ks pedotransfer parameters - only when soil_moisture is inactive
          if (self%exchange%config%processes%soil_moisture == 0_i4) then
            lake_params = [lake_params, &
              parameter_config%lakes_1%ptf_ks_constant, &
              parameter_config%lakes_1%ptf_ks_sand, &
              parameter_config%lakes_1%ptf_ks_clay]
            lake_names = [lake_names, [character(64) :: "ptf_ks_constant", "ptf_ks_sand", "ptf_ks_clay"]]
          end if
          ! baseflow recession per geological unit - only when the land baseflow process is inactive
          if (self%exchange%config%processes%baseflow == 0_i4) then
            block
              integer(i4) :: ig, n_old, n_new
              character(64), allocatable :: tmp_names(:)
              lake_params = [lake_params, parameter_config%lakes_1%baseflow_recession]
              ! grow lake_names via explicit temporary to avoid gfortran issues with
              ! length-mismatched array constructors and self-referential reallocation
              n_old = size(lake_names, kind=i4)
              n_new = size(parameter_config%lakes_1%baseflow_recession, kind=i4)
              allocate(tmp_names(n_old + n_new))
              tmp_names(1:n_old) = lake_names
              do ig = 1_i4, n_new
                tmp_names(n_old + ig) = "baseflow_recession(" // trim(adjustl(n2s(ig))) // ")"
              end do
              call move_alloc(tmp_names, lake_names)
            end block
          end if
          call self%exchange%parameters%add_process("lake", lake_params, lake_names, group="lakes_1")
          deallocate(lake_params, lake_names)
        end block
      end associate
    end if
    call self%exchange%lake_outflow%provide("mLM")
  end subroutine mlm_configure

  !> \brief Resolve static lake metadata and publish hourly outflow.
  subroutine mlm_connect(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: n_lakes
    integer(i8) :: k  !< target slot in point-set (lake_ids) order
    integer(i8), allocatable :: restart_perm(:)  !< restart-order index for each point-set slot
    integer(i4) :: lake_case
    logical :: has_points, has_ids, has_levels, has_area, has_volume, has_metadata

    lake_case = self%exchange%config%processes%lake
    has_points = associated(self%exchange%lake_points)
    has_ids = self%exchange%lake_ids%provided .and. self%exchange%lake_ids%has_data()
    has_levels = self%exchange%lake_spill_reference_level%provided .and. self%exchange%lake_spill_reference_level%has_data()
    has_area = self%exchange%lake_spill_reference_surface_area%provided .and. self%exchange%lake_spill_reference_surface_area%has_data()
    has_volume = self%exchange%lake_spill_reference_volume%provided .and. self%exchange%lake_spill_reference_volume%has_data()
    has_metadata = has_points .and. has_ids .and. has_levels
    if ((has_points .neqv. has_ids) .or. (has_points .neqv. has_levels)) then
      log_fatal(*) "mLM: lake points, stable IDs, and maximum levels must be published together."
      error stop 1
    end if

    if (has_metadata) then
      n_lakes = self%exchange%lake_points%n_points
      call self%exchange%check_data(self%exchange%lake_ids, "mLM")
      call self%exchange%check_data(self%exchange%lake_spill_reference_level, "mLM")
      if (self%read_restart) call self%validate_restart_metadata()
      self%definition%ids = self%exchange%lake_ids%data
      call self%definition%points%init( &
        self%exchange%lake_points%x, self%exchange%lake_points%y, coordsys=self%exchange%lake_points%coordsys)
      self%definition%spill_reference_level = self%exchange%lake_spill_reference_level%data
      if (has_area) self%definition%spill_reference_surface_area = self%exchange%lake_spill_reference_surface_area%data
      if (has_volume) self%definition%spill_reference_volume = self%exchange%lake_spill_reference_volume%data
      ! Restore the complete static restart definition in the supplied point-set order.
      ! prepare_restart keeps the published restart targets separate from these working copies.
      if (self%read_restart) then
        allocate(restart_perm(n_lakes))
        do k = 1_i8, n_lakes
          restart_perm(k) = findloc(self%restart%ids, self%definition%ids(k), dim=1, kind=i8)
        end do
        if (lake_case == 1_i4 .or. lake_case == -2_i4) &
          self%definition%footprint_area = self%restart%footprint_area(restart_perm)
        if (lake_case == 1_i4) then
          self%definition%ks = self%restart%ks(restart_perm)
          self%definition%kbf = self%restart%kbf(restart_perm)
          self%definition%spillway_length = self%restart%spillway_length(restart_perm)
          self%definition%endorheic = self%restart%endorheic(restart_perm)
          self%definition%initial_level = self%restart%initial_level(restart_perm)
          self%definition%spill_reference_surface_area = self%restart%spill_reference_surface_area(restart_perm)
          self%definition%spill_reference_volume = self%restart%spill_reference_volume(restart_perm)
          if (allocated(self%restart%latitude)) self%definition%latitude = self%restart%latitude(restart_perm)
          if (allocated(self%restart%longitude)) self%definition%longitude = self%restart%longitude(restart_perm)
        end if
        deallocate(restart_perm)
      end if
    else
      log_fatal(*) "mLM: lake process requires lake points, stable IDs and spill reference levels."
      error stop 1
    end if
    if (n_lakes < 1_i8) then
      log_fatal(*) "mLM: lake process requires at least one lake."
      error stop 1
    end if
    if (lake_case == -2_i4) then
      if (.not.self%owns_restart_lake_grid) then
        if (associated(self%exchange%level0_lake) .and. self%exchange%lake_map%has_data()) then
          call self%exchange%check_data(self%exchange%lake_map, "mLM")
          self%definition%grid = self%exchange%level0_lake
          self%definition%map = self%exchange%lake_map%data
        end if
      end if
      if (.not.self%read_restart .and. .not.allocated(self%definition%map)) then
        log_fatal(*) "mLM process -2 requires level-0 lake topology for a cold start."
        error stop 1
      end if
      if (allocated(self%definition%map)) then
        block
          real(dp), allocatable :: footprint(:)
          call lake_derive_spill_reference_surface_area_dem(self%definition%grid, self%definition%map, &
            self%definition%ids, footprint)
          if (allocated(self%definition%footprint_area)) then
            if (any(abs(footprint-self%definition%footprint_area) > 1.0e-10_dp*max(1.0_dp,footprint))) then
              log_fatal(*) "mLM: supplied footprint area conflicts with restart definition."
              error stop 1
            end if
          else
            call move_alloc(footprint, self%definition%footprint_area)
          end if
        end block
      end if
      self%definition%spill_reference_surface_area_dem_derived = self%definition%footprint_area
      call self%exchange%lake_spill_reference_surface_area_dem_derived%prepare_data("mLM", no_time)
      self%exchange%lake_spill_reference_surface_area_dem_derived%data => self%definition%spill_reference_surface_area_dem_derived
      call self%exchange%lake_pre%check_provided("mLM")
      call self%exchange%lake_pet%check_provided("mLM")
    end if
    if (lake_case == 1_i4) then
      ! The minimal natural lake needs the reported crest surface area and volume to build
      ! its height-area-volume (hAV) curves at initialize.
      if (.not.allocated(self%definition%spill_reference_surface_area) .or. &
          .not.allocated(self%definition%spill_reference_volume)) then
        log_fatal(*) "mLM process 1 requires reported maximum lake surface area and volume."
        error stop 1
      end if
      ! Capture the level-0 lake footprint, its stable lake map and its top-horizon soil class.
      ! The footprint grid + map drive the detached lake-MPR regionalization of the lakebed Ks at
      ! initialize (soil texture over the footprint -> Cosby Ks -> area-weighted mean per lake),
      ! and are retained for restart. On a restart-only run Input provides none of these; the
      ! lakebed Ks is then restored from the restart metadata instead (see the restart branch above).
      ! Skip the Input-data binding when prepare_restart already restored the topology from the restart
      ! file (restart continued without the original lake inputs). Otherwise bind the Input-published
      ! footprint grid, lake map and lake soil/geology for the detached lake-MPR regionalization.
      if (.not.self%owns_restart_lake_grid .and. &
          associated(self%exchange%level0_lake) .and. self%exchange%lake_map%provided) then
        call self%exchange%check_data(self%exchange%lake_map, "mLM")
        self%definition%grid = self%exchange%level0_lake
        self%definition%map = self%exchange%lake_map%data
        if (self%exchange%lake_sand%provided .and. self%exchange%lake_sand%has_data() .and. &
            self%exchange%lake_clay%provided .and. self%exchange%lake_clay%has_data()) then
          call self%exchange%check_data(self%exchange%lake_sand, "mLM")
          call self%exchange%check_data(self%exchange%lake_clay, "mLM")
          self%regionalization%sand = self%exchange%lake_sand%data
          self%regionalization%clay = self%exchange%lake_clay%data
        end if
        ! geological unit over the lake footprint for the baseflow recession regionalization
        if (self%exchange%lake_geo%provided .and. self%exchange%lake_geo%has_data()) then
          call self%exchange%check_data(self%exchange%lake_geo, "mLM")
          self%regionalization%geo = self%exchange%lake_geo%data
        end if
        ! The DEM-footprint surface area per lake (level-0 lake-cell area summed over the footprint)
        ! is the fixed area over which mRM collects the lake's inflow. It is retained here (and
        ! restart-persisted) so the per-step lake banks term can reconcile it against the dynamic
        ! water-balance area A_water (see mlm_update). Also logged against the reported area to record
        ! how far the level-0 delineation lies from the reported crest area.
        block
          real(dp), allocatable :: footprint(:)
          call lake_derive_spill_reference_surface_area_dem(self%definition%grid, self%definition%map, &
            self%definition%ids, footprint)
          if (allocated(self%definition%footprint_area)) then
            if (any(abs(footprint-self%definition%footprint_area) > 1.0e-10_dp*max(1.0_dp,footprint))) then
              log_fatal(*) "mLM: supplied footprint area conflicts with restart definition."
              error stop 1
            end if
          else
            call move_alloc(footprint, self%definition%footprint_area)
          end if
        end block
        block
          integer(i8) :: kk
          do kk = 1_i8, size(self%definition%ids, kind=i8)
            log_info(*) "mLM lake ", n2s(self%definition%ids(kk)), ": reported area [m2] = ", &
              n2s(self%definition%spill_reference_surface_area(kk)), " , DEM-footprint area [m2] = ", &
              n2s(self%definition%footprint_area(kk)), " ( ratio DEM/reported = ", &
              n2s(self%definition%footprint_area(kk) / self%definition%spill_reference_surface_area(kk)), " )"
          end do
        end block
      end if
      ! Allocate the lake states carried by the component (built at initialize, restart-persisted in c4).
      allocate(self%state%level(n_lakes), source=0.0_dp)
      allocate(self%state%surface_area(n_lakes), source=0.0_dp)
      allocate(self%state%volume(n_lakes), source=0.0_dp)
      allocate(self%state%gw_storage(n_lakes), source=0.0_dp)
      call self%exchange%lake_level%prepare_data("mLM", 1_i4)
      self%exchange%lake_level%data => self%state%level
      call self%exchange%lake_surface_area%prepare_data("mLM", 1_i4)
      self%exchange%lake_surface_area%data => self%state%surface_area
      call self%exchange%lake_volume%prepare_data("mLM", 1_i4)
      self%exchange%lake_volume%data => self%state%volume
      call self%exchange%lake_gw_storage%prepare_data("mLM", 1_i4)
      self%exchange%lake_gw_storage%data => self%state%gw_storage
      ! The minimal natural lake also needs the spillway length (crest spill) and an
      ! initial level (cold-start). The spillway length is regionalized per lake from the lake
      ! catchment (drainage) area (see derive_lake_spillway_length, called below alongside the Ks and
      ! baseflow-recession regionalizations); on a restart-only run (no Input) it was restored from the
      ! restart metadata at connect (the restart branch above). Nothing is read here: if it was not
      ! restart-restored it is derived below, so leave it unallocated for the derivation to fill.
      ! The initial level is optional: if it is absent (or carries the missing-value sentinel
      ! nodata_dp = -9999), the lake is cold-started by type (natural: just above the spill reference level;
      ! reservoir: 50% of its spill-reference volume) - see mlm_initialize.
      if (self%exchange%lake_initial_level%provided .and. self%exchange%lake_initial_level%has_data()) then
        self%definition%initial_level = self%exchange%lake_initial_level%data
      else if (.not.allocated(self%definition%initial_level)) then
        allocate(self%definition%initial_level(n_lakes), source=nodata_dp)
      end if
      ! Geographic latitude/longitude of the lake points, used by the varying-albedo solar geometry.
      ! The Input container derives them (from the point coordinates when spherical, or the level-0
      ! latlon auxiliary coordinates when projected); on a restart-only run they come from the restart
      ! metadata. Left as nodata when unavailable; varying_albedo then errors in the guard below.
      if (self%exchange%lake_latitude%provided .and. self%exchange%lake_latitude%has_data()) then
        self%definition%latitude = self%exchange%lake_latitude%data
      else if (.not.allocated(self%definition%latitude)) then
        allocate(self%definition%latitude(n_lakes), source=nodata_dp)
      end if
      if (self%exchange%lake_longitude%provided .and. self%exchange%lake_longitude%has_data()) then
        self%definition%longitude = self%exchange%lake_longitude%data
      else if (.not.allocated(self%definition%longitude)) then
        allocate(self%definition%longitude(n_lakes), source=nodata_dp)
      end if
      ! Validate the configured capacity before indexing per-lake options.
      if (n_lakes > int(self%config%max_lakes,i8)) then
        log_fatal(*) "mLM: lake count exceeds config_project max_lakes."
        error stop 1
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
          if (.not.allocated(self%definition%latitude) .or. .not.allocated(self%definition%longitude)) then
            log_fatal(*) "mLM process 1 varying_albedo requires lake geographic coordinates (set latlon_path for a projected grid)."
            error stop 1
          end if
          if (any(self%definition%latitude < -90.0_dp) .or. any(self%definition%latitude > 90.0_dp) .or. &
              any(self%definition%longitude < -360.0_dp) .or. any(self%definition%longitude > 360.0_dp)) then
            log_fatal(*) "mLM process 1 varying_albedo: lake geographic coordinates are missing or out of range (set latlon_path for a projected grid)."
            error stop 1
          end if
        end if
      end if
      ! Previous-step state storage for the water-balance integration.
      allocate(self%state%level_previous(n_lakes), source=0.0_dp)
      allocate(self%state%surface_area_previous(n_lakes), source=0.0_dp)
      allocate(self%state%volume_previous(n_lakes), source=0.0_dp)
      allocate(self%state%gw_storage_previous(n_lakes), source=0.0_dp)
      ! Diagnostic per-lake fluxes retained for the lake-point output.
      allocate(self%state%evaporation(n_lakes), source=0.0_dp)
      allocate(self%state%seepage(n_lakes), source=0.0_dp)
      allocate(self%state%baseflow(n_lakes), source=0.0_dp)
      allocate(self%state%spill(n_lakes), source=0.0_dp)
    end if
    if (lake_case == 1_i4 .and. .not.self%read_restart) then
      call self%cache_regionalization()
      call self%derive_lake_endorheic()
    end if
    if (self%write_restart .and. self%restart_save_topology) then
      if (.not.allocated(self%definition%map)) then
        log_fatal(*) "mLM: restart_save_topology requires available level-0 lake geometry."
        error stop 1
      end if
    end if
    allocate(self%state%outflow(n_lakes), source=0.0_dp)
    call self%exchange%lake_outflow%prepare_data("mLM", 1_i4)
    self%exchange%lake_outflow%data => self%state%outflow
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
    ! Cold-start head above the spill reference level for a NATURAL lake when no initial level is
    ! provided: a small positive head so the lake discharges over its outlet sill from the first step
    ! (natural lakes sit at/above the sill and spill continuously) rather than filling up to it first.
    real(dp), parameter :: natural_coldstart_head = 0.1_dp  ! [m]
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
      call self%exchange%check_data(self%exchange%lake_spill_reference_surface_area_dem_derived, "mLM")
    end if
    if (self%exchange%config%processes%lake == 1_i4) then
      call self%fetch_parameters()
      self%state%evaporation(:) = 0.0_dp
      self%state%seepage(:) = 0.0_dp
      self%state%baseflow(:) = 0.0_dp
      self%state%spill(:) = 0.0_dp
      ! Build the per-lake height-area-volume (hAV) curves from a built-in analytical shape.
      n_lakes = size(self%definition%ids, kind=i8)
      domain = self%exchange%nml_domain_id
      tsHydroSecs = real(self%exchange%step_hours, dp) * 3600.0_dp
      if (allocated(self%state%hav)) deallocate(self%state%hav)
      allocate(self%state%hav(n_lakes))
      do k = 1_i8, n_lakes
        bathymetry_type = self%config%bathymetry_type(self%exchange%nml_domain_id, k)
        if (bathymetry_type /= 1_i4 .and. bathymetry_type /= 2_i4) then
          log_fatal(*) "mLM process 1: bathymetry_type ", n2s(bathymetry_type), &
            " for lake ", n2s(k), " is not supported (0 = external bathymetry file is not yet implemented)."
          error stop 1
        end if
        call generate_lake_bathymetry(bathymetry_type, &
          self%definition%spill_reference_level(k), self%definition%spill_reference_surface_area(k), &
          self%definition%spill_reference_volume(k), &
          self%state%hav(k)%level, self%state%hav(k)%surface_area, self%state%hav(k)%volume)
      end do
      ! Rebuild parameter-dependent regionalization only for cold initialization.
      ! Restart definition and its effective generating parameters were validated before binding.
      if (.not.self%read_restart) then
        call self%derive_lake_ks()
        call self%derive_lake_kbf()
        call self%derive_lake_spillway_length()
      end if
      ! Cold-start the lake states from the reported initial level; the matching surface
      ! area and volume are read off the per-lake hAV curve. A restart overrides these below.
      if (.not. self%read_restart) then
        do k = 1_i8, n_lakes
          ! When the initial level was not provided (missing-value sentinel), cold-start by lake type:
          !  - a NATURAL lake sits at its outlet sill essentially all the time and discharges
          !    continuously, so cold-start it just above the spill reference level (a small head
          !    natural_coldstart_head above the sill) so it spills from the first step rather than
          !    filling up to the sill first;
          !  - a RESERVOIR is normally held below its full-supply level by controlled releases, so
          !    cold-start it at the level holding 50% of its spill-reference (full-supply) volume
          !    (mid-pool, the v5.13 behaviour): invert the hAV curve (volume -> level) at half that volume.
          if (self%definition%initial_level(k) < -9000.0_dp) then
            if (self%config%natural(self%exchange%nml_domain_id, int(k, i4))) then
              self%definition%initial_level(k) = self%definition%spill_reference_level(k) + natural_coldstart_head
              log_info(*) "mLM lake ", n2s(self%definition%ids(k)), &
                ": initial level not provided; natural lake cold-started at spill reference level + ", &
                n2s(natural_coldstart_head), " m -> level ", n2s(self%definition%initial_level(k)), " masl"
            else
              block
                real(dp) :: target_volume
                target_volume = 0.5_dp * self%definition%spill_reference_volume(k)
                call interpolate_lake_states(.false., self%state%hav(k)%volume, self%state%hav(k)%level, &
                  self%definition%spill_reference_volume(k), self%definition%spill_reference_level(k), &
                  target_volume, self%definition%initial_level(k))
              end block
              log_info(*) "mLM lake ", n2s(self%definition%ids(k)), &
                ": initial level not provided; reservoir cold-started at 50% of spill-reference volume -> level ", &
                n2s(self%definition%initial_level(k)), " masl"
            end if
          end if
          self%state%level(k) = self%definition%initial_level(k)
          level_query = self%definition%initial_level(k)
          call interpolate_lake_states(.false., self%state%hav(k)%level, self%state%hav(k)%surface_area, &
            self%definition%spill_reference_level(k), self%state%hav(k)%surface_area(size(self%state%hav(k)%surface_area)), &
            level_query, self%state%surface_area(k))
          level_query = self%definition%initial_level(k)
          call interpolate_lake_states(.false., self%state%hav(k)%level, self%state%hav(k)%volume, &
            self%definition%spill_reference_level(k), self%state%hav(k)%volume(size(self%state%hav(k)%volume)), &
            level_query, self%state%volume(k))
          ! Cold-start the lake groundwater storage self-consistently: the initial baseflow is
          ! set equal to the initial percolation, so the groundwater store starts at its
          ! own equilibrium and there is no spin-up drift. Evaluate the Darcy seepage at the initial
          ! level -> perc0 [mm/TS]; the per-step recession kbf0 = (tsHydroSecs/86400) /
          ! static_lake_kbf (the per-lake regionalized recession coefficient); then gw0 = (perc0/1000 * A0) / kbf0 gives
          ! baseflow0 = kbf0 * gw0 = perc0 volume. leakance = Ks/b, with Ks the per-lake
          ! regionalized lakebed conductivity (lake MPR) and b the global lakebed thickness.
          ! Endorheic lakes have no percolation (and no baseflow), so their groundwater store
          ! stays empty; skip the self-consistent cold-start (and its kbf0 divide) for them.
          if (self%definition%endorheic(k)) then
            self%state%gw_storage(k) = 0.0_dp
          else
            leakance = (self%definition%ks(k) * 0.01_dp &
              * tsHydroSecs / 86400.0_dp) / self%runtime%b
            call lake_seepage(.true., self%state%hav(k)%level, self%state%hav(k)%surface_area, &
              self%definition%spill_reference_level(k), &
              self%state%hav(k)%surface_area(size(self%state%hav(k)%surface_area)), &
              self%definition%initial_level(k), leakance, perc0)
            kbf0 = (tsHydroSecs / 86400.0_dp) / self%definition%kbf(k)
            self%state%gw_storage(k) = (perc0 / 1000.0_dp * self%state%surface_area(k)) / kbf0
          end if
        end do
      end if
    end if
    self%state%outflow(:) = 0.0_dp
    if (self%read_restart) call self%read_restart_state()
    if (self%exchange%config%processes%lake == 1_i4) then
      ! Seed the previous-step states for the water-balance integration.
      self%state%level_previous = self%state%level
      self%state%surface_area_previous = self%state%surface_area
      self%state%volume_previous = self%state%volume
      self%state%gw_storage_previous = self%state%gw_storage
    end if
    call self%validate_output_timing()
    call self%create_output()
  end subroutine mlm_initialize

  !> \brief Regionalize the per-lake lakebed saturated hydraulic conductivity (detached lake MPR).
  !> \details Applies current effective Cosby coefficients to cached topsoil texture,
  !! then forms a cell-count arithmetic mean over each lake footprint.
  subroutine mlm_derive_lake_ks(self)
    class(mlm_t), target, intent(inout) :: self
    if (.not.allocated(self%regionalization%sand) .or. .not.allocated(self%regionalization%clay)) then
      log_fatal(*) "mLM process 1 lakebed conductivity requires the lake footprint soil texture."
      error stop 1
    end if
    call lake_derive_ks(self%definition%grid, self%definition%map, self%definition%ids, &
      self%regionalization%sand, self%regionalization%clay, self%runtime%ptf, self%definition%ks, self%regionalization%lake_index)

    block
      integer(i8) :: kk
      do kk = 1_i8, size(self%definition%ids, kind=i8)
        log_info(*) "mLM lake ", n2s(self%definition%ids(kk)), ": regionalized lakebed Ks [cm d-1] = ", &
          n2s(self%definition%ks(kk))
      end do
    end block

  end subroutine mlm_derive_lake_ks

  !> \brief Regionalize the per-lake baseflow recession coefficient (detached lake MPR).
  !> \details Runs the geological baseflow-parameter lookup over the lake footprint geological units
  !!          (geo_unit -> geology LUT -> global baseflow-parameter slot -> per-geological-unit baseflow
  !!          recession value) and aggregates it to one recession coefficient per lake by an arithmetic
  !!          (cell-count) mean, the same path the v5 mpr_lake used for the lake baseflow coefficient and
  !!          matching the land-surface mpr_bridge_baseflow_param cell lookup. Self-contained and
  !!          mHM-independent: the lake reads its OWN geology LUT (config_lut geo_lut_path) and the
  !!          lake-footprint geological unit is surfaced by the Input container, so the regionalization
  !!          does not depend on the land-surface geology class definition or MPR. The per-geological-unit
  !!          baseflow recession values come from one of two sources: when the land-surface baseflow
  !!          (geoparameter) process is active the lake BORROWS its baseflow_recession parameters (so the
  !!          calibration space is not duplicated); when it is inactive the lake uses its own lakes_1
  !!          baseflow_recession parameters. Every geological unit on the lake footprint is checked
  !!          against the geology LUT first (lake-only units are not validated by the land-surface check).
  subroutine mlm_derive_lake_kbf(self)
    class(mlm_t), target, intent(inout) :: self
    if (.not.allocated(self%regionalization%geo_indices)) then
      log_fatal(*) "mLM process 1 baseflow requires cached lake geology."
      error stop 1
    end if
    if (any(self%regionalization%geo_indices > size(self%runtime%recession))) then
      log_fatal(*) "mLM lake geology references an unavailable baseflow parameter."
      error stop 1
    end if
    call lake_derive_kbf(self%definition%map, self%definition%ids, self%regionalization%geo, self%regionalization%geo_units, &
      self%regionalization%geo_indices, self%runtime%recession, self%definition%kbf, &
      self%regionalization%lake_index, self%regionalization%geo_parameter)

    block
      integer(i8) :: kk
      do kk = 1_i8, size(self%definition%ids, kind=i8)
        log_info(*) "mLM lake ", n2s(self%definition%ids(kk)), &
          ": regionalized baseflow recession coefficient [d] = ", n2s(self%definition%kbf(kk))
      end do
    end block

  end subroutine mlm_derive_lake_kbf

  !> \brief Regionalize the per-lake effective spillway crest length (detached lake MPR).
  !> \details Derives the spillway crest length Ls [m] per lake from the lake catchment (drainage)
  !!          area CA by the downstream-hydraulic-geometry power law Ls = m * 2.7 * CA^0.35. CA is the
  !!          routing upstream / contributing area at the lake outlet node (river%upstream_area at
  !!          river%lake_outlet_nodes, [m2], converted to km2 inside lake_derive_spillway_length), and
  !!          m is the global (pan-domain) spillway_length_multiplier (lakes_1). The outlet nodes are in
  !!          lake-point order (set together with lake_ids by river%label_lakes), so no footprint
  !!          packing is needed - this is one upstream-area value read per lake. The upstream area is
  !!          computed on demand if the routing output did not already request it.
  subroutine mlm_derive_lake_spillway_length(self)
    class(mlm_t), target, intent(inout) :: self
    if (.not.allocated(self%regionalization%catchment)) then
      log_fatal(*) "mLM process 1 spillway length requires cached lake catchment areas."
      error stop 1
    end if
    call lake_derive_spillway_length(self%regionalization%catchment, self%runtime%multiplier, &
      self%definition%spillway_length)
  end subroutine mlm_derive_lake_spillway_length

  !> \brief Flag endorheic (terminal, no-outlet) lakes from the routing topology.
  !!
  !! \details A lake whose snapped level-0 outlet node is a sink (its downstream node is zero) has
  !!          no downstream river: it is endorheic. Such a lake can shed water only through the
  !!          atmosphere, so its crest spill, lakebed percolation and groundwater baseflow are all
  !!          suppressed downstream (see mlm_update) and its balance closes by evaporation alone. The
  !!          outlet nodes are in lake-point order (set together with lake_ids by river%label_lakes),
  !!          the same handle used by the spillway-length regionalization, so no footprint packing is
  !!          needed. On a restart-only run (no river) the flag is restored from restart metadata at
  !!          connect and this routine is not reached.
  subroutine mlm_derive_lake_endorheic(self)
    use mo_river, only: river_t
    class(mlm_t), target, intent(inout) :: self
    type(river_t), pointer :: river
    integer(i8) :: i, k, n_lakes

    if (.not.associated(self%exchange%river_l0)) then
      log_fatal(*) "mLM endorheic flag requires the level-0 river (active routing)."
      error stop 1
    end if
    river => self%exchange%river_l0
    if (.not.allocated(river%lake_map) .or. .not.allocated(river%is_inland_sink)) then
      log_fatal(*) "mLM endorheic flag requires the level-0 river lake map and inland-sink flags."
      error stop 1
    end if

    ! A lake is endorheic iff its footprint contains an inland depression: a cell whose input flow
    ! direction is a true pit (not a valid D8 direction, e.g. the -1 depression marker), as opposed to
    ! a cell that drains off the domain edge. river%is_inland_sink carries exactly that distinction
    ! (captured from the input fdir before the edge-forcing in river_from_fdir), so pairing it with
    ! the lake footprint mask (river%lake_map) flags endorheic lakes without the domain-edge
    ! false-positive that a bare sink test (down == 0) would give.
    n_lakes = size(self%definition%ids, kind=i8)
    if (allocated(self%definition%endorheic)) deallocate(self%definition%endorheic)
    allocate(self%definition%endorheic(n_lakes), source=.false.)
    do i = 1_i8, size(river%lake_map, kind=i8)
      if (river%lake_map(i) <= 0_i8 .or. .not.river%is_inland_sink(i)) cycle
      k = findloc(self%definition%ids, river%lake_map(i), dim=1, kind=i8)
      if (k >= 1_i8) self%definition%endorheic(k) = .true.
    end do
    do k = 1_i8, n_lakes
      log_info(*) "mLM lake ", n2s(self%definition%ids(k)), ": endorheic = ", &
        merge("T", "F", self%definition%endorheic(k))
    end do

  end subroutine mlm_derive_lake_endorheic

  !> \brief Advance lake processes using lagged inflow on the shared hourly interval.
  subroutine mlm_update(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: k, n_lakes
    integer(i4) :: domain
    integer(i4) :: hour, doy
    real(dp) :: spill, volume_adjustment, tsHydroSecs, perc, baseflow, leakance, kbf
    real(dp) :: ntimesteps_day, evap_daily, evap_step, pr_step, qin_d, qout_d, pre_d, lat, lon
    real(dp) :: wind_k, tdew_k
    real(dp) :: a_bank, es_step, q_bank
    real(dp) :: deficit, spill_vol, perc_vol, evap_vol, avg_area, cut
    real(dp) :: evap_components(3)
    logical :: isday
    real(dp) :: input_windspeed_height, climatological_windspeed
    integer(i4) :: albedo_method
    logical :: varying_albedo, varying_lambda, varying_psychro
    select case (self%exchange%config%processes%lake)
      case(-1_i4)
        self%state%outflow(:) = self%exchange%lake_inflow%data
      case(-2_i4)
      ! if (self%exchange%config%processes%lake == -2_i4) then
        self%state%outflow(:) = lake_balance_outflow(self%exchange%lake_inflow%data, self%definition%spill_reference_surface_area_dem_derived, self%exchange%lake_pre%data, &
          self%exchange%lake_pet%data, 3600.0_dp)
      case(1_i4)
        ! Evaluate losses from the previous state, reconcile the lake water balance, then
        ! update groundwater. Routed outflow is corrected crest spill plus groundwater baseflow.
        n_lakes = size(self%definition%ids, kind=i8)
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
        !$omp           evap_daily, evap_step, pr_step, qin_d, qout_d, pre_d, lat, lon, wind_k, tdew_k, evap_components, &
        !$omp           a_bank, es_step, q_bank, deficit, spill_vol, perc_vol, evap_vol, avg_area, cut) &
        !$omp   schedule(dynamic)
        do k = 1_i8, n_lakes
          ! Darcy lakebed leakance Ks/b [TS-1 per unit head] and per-step baseflow recession.
          ! leakance replaces the former kp (recession rate) + exponent_for_percolation
          ! (depth-normaliser). Ks is the per-lake regionalized lakebed conductivity (lake MPR) and b
          ! is the global pan-domain lakebed thickness; Ks is cm/d, converted to m/TS then divided by b [m].
          leakance = (self%definition%ks(k) * 0.01_dp &
            * tsHydroSecs / 86400.0_dp) / self%runtime%b
          ! per-step baseflow recession from the per-lake regionalized recession coefficient (lake MPR)
          kbf = (tsHydroSecs / 86400.0_dp) / self%definition%kbf(k)
          ! Lake seepage from the previous-step level -> percolation to the groundwater store.
          ! Endorheic (terminal) lakes have no subsurface outlet, so percolation is suppressed: the
          ! only loss is evaporation (the physical closure of a terminal lake).
          if (self%definition%endorheic(k)) then
            perc = 0.0_dp
          else
            call lake_seepage(.true., self%state%hav(k)%level, self%state%hav(k)%surface_area, &
              self%definition%spill_reference_level(k), &
              self%state%hav(k)%surface_area(size(self%state%hav(k)%surface_area)), &
              self%state%level_previous(k), leakance, &
              perc)
          end if
          ! Crest-controlled spill from the previous-step level. Endorheic lakes never spill (no
          ! downstream river); the level floats above the sill on the water balance.
          if (self%definition%endorheic(k)) then
            spill = 0.0_dp
          else
            call lake_spill(self%state%level_previous(k), self%definition%spill_reference_level(k), &
              self%runtime%cd, &
              self%definition%spillway_length(k), &
              self%runtime%exponent, &
              spill)
          end if
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
            if (self%state%surface_area_previous(k) > 0.0_dp) then
              qin_d  = self%exchange%lake_inflow%data(k) / self%state%surface_area_previous(k) * 1000.0_dp * 86400.0_dp
              qout_d = self%state%outflow(k) / self%state%surface_area_previous(k) * 1000.0_dp * 86400.0_dp
            else
              qin_d  = 0.0_dp
              qout_d = 0.0_dp
            end if
            pre_d = pr_step * ntimesteps_day
            ! Geographic coordinates for the solar geometry (used only when varying_albedo); these are
            ! true lat/lon even on a projected grid (from the level-0 latlon auxiliary coordinates),
            ! unlike static_lake_points%x/y which hold projected eastings/northings on a cartesian grid.
            lat = self%definition%latitude(k)
            lon = self%definition%longitude(k)
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
              self%state%level_previous(k), qin_d, qout_d, pre_d, lat, lon, &
              input_windspeed_height, climatological_windspeed, varying_albedo, albedo_method, varying_lambda, varying_psychro, &
              evap_daily, evap_components)
            call temporal_disagg_flux_daynight(isday, ntimesteps_day, evap_daily, 0.5_dp, 0.5_dp, evap_step)
          end if
          ! Lake banks (Shrestha & Mueller, 2026): reconcile the fixed DEM-footprint collection area
          ! against the dynamic water-balance area. The atmospheric term of the water balance applies
          ! P and E over the dynamic water area A_water, while mRM collects inflow over the fixed
          ! DEM footprint A_footprint. The ring between them (A_footprint - A_water) is therefore
          ! unaccounted. Account its net atmospheric flux as a same-step inflow correction so the lake
          ! budget closes over the true footprint. The correction is symmetric:
          !   A_bank > 0 (water inside the footprint): the exposed-bank ring receives P and loses soil
          !     evaporation Es = c_soil * E (moist bank); its net (P - Es) is routed to the lake, on the
          !     basis that bank runoff and subsurface seepage both ultimately report to the lake. The
          !     bank has NO soil-moisture store, so Es can only consume water that fell this step:
          !     Es = min(c_soil * E, P). Hence (P - Es) >= 0 and the exposed bank never acts as a net
          !     sink (no dry-bank water removal without a store).
          !   A_bank < 0 (water beyond the footprint): the over-flooded ring's (P - E) was already
          !     accounted on the mHM land side, so it is removed from the lake to avoid double counting.
          !     This side uses open-water E and stays signed (a genuine correction of counted flux).
          ! A_water is the previous-step surface area (non-circular, as for perc/baseflow).
          a_bank = self%definition%footprint_area(k) - self%state%surface_area_previous(k)
          if (a_bank >= 0.0_dp) then
            es_step = min(bank_soil_evap_ratio * evap_step, pr_step)  ! exposed bank: soil E capped by same-step P (no store)
          else
            es_step = evap_step                                  ! over-flooded: open-water evaporation
          end if
          q_bank = (pr_step - es_step) / 1000.0_dp * a_bank / tsHydroSecs  ! [m3 s-1]; exposed side >=0, flooded side signed
          ! Iterative water balance -> updates level, surface area and volume (consumes perc, evap, spill).
          call lake_balance(.true., pr_step, evap_step, self%exchange%lake_inflow%data(k) + q_bank, perc, 0.0_dp, spill, 0.0_dp, &
            self%state%hav(k)%level, self%state%hav(k)%surface_area, self%state%hav(k)%volume, &
            self%definition%spill_reference_level(k), &
            self%state%hav(k)%surface_area(size(self%state%hav(k)%surface_area)), &
            self%state%hav(k)%volume(size(self%state%hav(k)%volume)), &
            self%state%level_previous(k), self%state%surface_area_previous(k), self%state%volume_previous(k), &
            self%state%level(k), self%state%surface_area(k), self%state%volume(k), &
            volume_adjustment, tsHydroSecs)
          ! Reconcile the published losses with the water actually available. When the step's net
          ! demand (evaporation + seepage + spill, less inflow and precipitation) exceeds the lake
          ! storage, lake_balance clamps the lake to empty and returns the shortfall as a positive
          ! volume_adjustment [m3]. The fluxes here were computed from the PREVIOUS-step level, so
          ! without this reconciliation the container would still route the full spill downstream and
          ! the full seepage to groundwater -- creating water the lake never held. Absorb the deficit
          ! in physical priority: spill first (a lake drawn below its sill cannot spill), then seepage
          ! (its driving head vanishes as the lake empties), then evaporation last (the terminal
          ! surface loss, and the only loss of an endorheic lake). Fluxes are rates; convert each to a
          ! step volume, cut it, and convert back so spill, perc and evap stay mutually consistent and
          ! the lake->groundwater transfer (perc) still matches on both sides.
          if (volume_adjustment > 0.0_dp) then
            deficit   = volume_adjustment
            avg_area  = 0.5_dp * (self%state%surface_area_previous(k) + self%state%surface_area(k))
            ! 1) spill [m3 s-1] -> step volume
            spill_vol = spill * tsHydroSecs
            cut       = min(deficit, spill_vol)
            spill     = (spill_vol - cut) / tsHydroSecs
            deficit   = deficit - cut
            ! 2) seepage [mm TS-1] over the previous-step area (the conserved lake->GW transfer volume)
            if (deficit > 0.0_dp) then
              perc_vol = perc / 1000.0_dp * self%state%surface_area_previous(k)
              cut      = min(deficit, perc_vol)
              if (self%state%surface_area_previous(k) > 0.0_dp) &
                perc = (perc_vol - cut) * 1000.0_dp / self%state%surface_area_previous(k)
              deficit  = deficit - cut
            end if
            ! 3) evaporation [mm TS-1] over the step-average area (as the balance applied it)
            if (deficit > 0.0_dp) then
              if (avg_area > 0.0_dp) then
                evap_vol  = evap_step / 1000.0_dp * avg_area
                cut       = min(deficit, evap_vol)
                evap_step = (evap_vol - cut) * 1000.0_dp / avg_area
                deficit   = deficit - cut
              end if
            end if
          end if
          ! Lake groundwater storage update -> baseflow (uses the previous-step surface area and GW storage).
          call lake_gws_baseflow(perc, self%state%surface_area_previous(k), &
            self%state%gw_storage_previous(k), self%state%gw_storage(k), &
            kbf, tsHydroSecs, self%definition%endorheic(k), baseflow)
          ! Natural-lake outflow is the crest spill plus the groundwater baseflow.
          self%state%outflow(k) = spill + baseflow
          ! Retain the per-lake fluxes for the lake-point output.
          self%state%evaporation(k) = evap_step
          self%state%seepage(k) = perc
          self%state%baseflow(k) = baseflow
          self%state%spill(k) = spill
        end do
        !$omp end parallel do
        ! Roll the current states into the previous-step slots for the next update.
        self%state%level_previous = self%state%level
        self%state%surface_area_previous = self%state%surface_area
        self%state%volume_previous = self%state%volume
        self%state%gw_storage_previous = self%state%gw_storage
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
    if (allocated(self%definition%ids)) deallocate(self%definition%ids)
    if (allocated(self%definition%map)) deallocate(self%definition%map)
    if (allocated(self%regionalization%sand)) deallocate(self%regionalization%sand)
    if (allocated(self%regionalization%clay)) deallocate(self%regionalization%clay)
    if (allocated(self%regionalization%geo)) deallocate(self%regionalization%geo)
    if (allocated(self%definition%spill_reference_surface_area_dem_derived)) deallocate(self%definition%spill_reference_surface_area_dem_derived)
    call self%definition%points%destroy()
    if (allocated(self%definition%spill_reference_level)) deallocate(self%definition%spill_reference_level)
    if (allocated(self%definition%spill_reference_surface_area)) deallocate(self%definition%spill_reference_surface_area)
    if (allocated(self%definition%spill_reference_volume)) deallocate(self%definition%spill_reference_volume)
    if (allocated(self%definition%spillway_length)) deallocate(self%definition%spillway_length)
    if (allocated(self%definition%endorheic)) deallocate(self%definition%endorheic)
    if (allocated(self%definition%initial_level)) deallocate(self%definition%initial_level)
    if (allocated(self%definition%ks)) deallocate(self%definition%ks)
    if (allocated(self%definition%kbf)) deallocate(self%definition%kbf)
    if (allocated(self%definition%footprint_area)) deallocate(self%definition%footprint_area)
    if (allocated(self%definition%latitude)) deallocate(self%definition%latitude)
    if (allocated(self%definition%longitude)) deallocate(self%definition%longitude)
    if (allocated(self%state%level)) deallocate(self%state%level)
    if (allocated(self%state%surface_area)) deallocate(self%state%surface_area)
    if (allocated(self%state%volume)) deallocate(self%state%volume)
    if (allocated(self%state%gw_storage)) deallocate(self%state%gw_storage)
    if (allocated(self%state%level_previous)) deallocate(self%state%level_previous)
    if (allocated(self%state%surface_area_previous)) deallocate(self%state%surface_area_previous)
    if (allocated(self%state%volume_previous)) deallocate(self%state%volume_previous)
    if (allocated(self%state%gw_storage_previous)) deallocate(self%state%gw_storage_previous)
    if (allocated(self%state%evaporation)) deallocate(self%state%evaporation)
    if (allocated(self%state%seepage)) deallocate(self%state%seepage)
    if (allocated(self%state%baseflow)) deallocate(self%state%baseflow)
    if (allocated(self%state%spill)) deallocate(self%state%spill)
    if (allocated(self%state%hav)) deallocate(self%state%hav)
    call self%restart%points%destroy()
    if (allocated(self%restart%ids)) deallocate(self%restart%ids)
    if (allocated(self%restart%spill_reference_level)) deallocate(self%restart%spill_reference_level)
    if (allocated(self%restart%spill_reference_surface_area)) deallocate(self%restart%spill_reference_surface_area)
    if (allocated(self%restart%spill_reference_volume)) deallocate(self%restart%spill_reference_volume)
    if (allocated(self%restart%spillway_length)) deallocate(self%restart%spillway_length)
    if (allocated(self%restart%endorheic)) deallocate(self%restart%endorheic)
    if (allocated(self%restart%initial_level)) deallocate(self%restart%initial_level)
    if (allocated(self%restart%ks)) deallocate(self%restart%ks)
    if (allocated(self%restart%kbf)) deallocate(self%restart%kbf)
    if (allocated(self%restart%footprint_area)) deallocate(self%restart%footprint_area)
    if (allocated(self%restart%latitude)) deallocate(self%restart%latitude)
    if (allocated(self%restart%longitude)) deallocate(self%restart%longitude)
    if (allocated(self%state%outflow)) deallocate(self%state%outflow)
    if (allocated(self%runtime%recession)) deallocate(self%runtime%recession)
    if (allocated(self%restart%generating_ks)) deallocate(self%restart%generating_ks)
    if (allocated(self%restart%generating_kbf)) deallocate(self%restart%generating_kbf)
    if (allocated(self%restart%generating_length)) deallocate(self%restart%generating_length)
    if (allocated(self%regionalization%geo_units)) deallocate(self%regionalization%geo_units)
    if (allocated(self%regionalization%geo_indices)) deallocate(self%regionalization%geo_indices)
    if (allocated(self%regionalization%catchment)) deallocate(self%regionalization%catchment)
    if (allocated(self%regionalization%lake_index)) deallocate(self%regionalization%lake_index)
    if (allocated(self%regionalization%geo_parameter)) deallocate(self%regionalization%geo_parameter)
    call self%definition%grid%destroy()
    if (allocated(self%restart_input_path)) deallocate(self%restart_input_path)
    if (allocated(self%restart_output_path)) deallocate(self%restart_output_path)
    if (allocated(self%output_path)) deallocate(self%output_path)
    self%active = .false.
    self%read_restart = .false.
    self%write_restart = .false.
    self%restart_save_topology = .false.
    self%evaporation_active = .false.
    self%output_active = .false.
  end subroutine mlm_destroy

  !> \brief Create the lake-point output dataset in stable exchange order.
  subroutine mlm_create_output(self)
    class(mlm_t), target, intent(inout) :: self
    type(var), allocatable :: vars(:)
    character(:), allocatable :: delta, dtype
    integer(i4) :: timestamp

    if (allocated(self%ds_out%vars)) call self%ds_out%close()
    if (.not.self%output_config%is_configured) return
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
    call self%ds_out%update("lake_id", self%definition%ids)
    call self%ds_out%write_static()
  end subroutine mlm_create_output

  !> \brief Validate mLM output cadence and restart/output boundary alignment.
  subroutine mlm_validate_output_timing(self)
    class(mlm_t), target, intent(in) :: self
    integer(i4) :: frequency

    if (.not.self%output_config%is_configured) return
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
    type(datetime), intent(in) :: time !< time
    integer(i4) :: elapsed_hours

    if (.not.self%output_config%is_configured) then
      boundary = .true.
      return
    end if
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

    if (.not.self%output_config%is_configured) return
    if (.not.self%output_active) return
    if (self%exchange%config%processes%lake == 1_i4) then
      if (self%output_config%out_lake_level) call self%ds_out%update("lake_level", self%state%level)
      if (self%output_config%out_lake_surface_area) &
        call self%ds_out%update("lake_surface_area", self%state%surface_area * 1.0e-6_dp)
      if (self%output_config%out_lake_volume) &
        call self%ds_out%update("lake_volume", self%state%volume * 1.0e-6_dp)
      if (self%output_config%out_lake_gw_storage) &
        call self%ds_out%update("lake_gw_storage", self%state%gw_storage * 1.0e-6_dp)
      if (self%output_config%out_lake_evaporation) call self%ds_out%update("lake_evaporation", self%state%evaporation)
      if (self%output_config%out_lake_seepage) call self%ds_out%update("lake_seepage", self%state%seepage)
      if (self%output_config%out_lake_baseflow) call self%ds_out%update("lake_baseflow", self%state%baseflow)
      if (self%output_config%out_lake_spill) call self%ds_out%update("lake_spill", self%state%spill)
    end if
    if (self%output_config%out_lake_inflow) call self%ds_out%update("lake_inflow", self%exchange%lake_inflow%data)
    if (self%output_config%out_lake_precipitation .and. self%exchange%lake_pre%provided) &
      call self%ds_out%update("lake_precipitation", self%exchange%lake_pre%data)
    if (self%output_config%out_lake_outflow) call self%ds_out%update("lake_outflow", self%state%outflow)
    if (self%at_output_boundary(self%exchange%time)) call self%ds_out%write(self%exchange%time)
  end subroutine mlm_update_output

  !> \brief Persist static lake metadata and current outflow with stable IDs.
  subroutine mlm_create_restart(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcDimension) :: dims(0), lake_dim
    type(NcVariable) :: var
    integer(i4), allocatable :: endorheic_flag(:) !< endorheic logical stored as i4 0/1 (NetCDF has no logical type)

    if (self%definition%points%n_points /= size(self%definition%ids, kind=i8) .or. &
        .not.allocated(self%definition%spill_reference_level)) then
      log_fatal(*) "mLM: cannot write restart without static lake metadata."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == 1_i4 .and. &
        (.not.allocated(self%definition%spill_reference_surface_area) .or. &
         .not.allocated(self%definition%spill_reference_volume))) then
      log_fatal(*) "mLM: cannot write process 1 restart without reported maximum surface area and volume."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == 1_i4 .and. &
        (.not.allocated(self%definition%spillway_length) .or. &
         .not.allocated(self%definition%initial_level))) then
      log_fatal(*) "mLM: cannot write process 1 restart without spillway length and initial level."
      error stop 1
    end if
    nc = NcDataset(self%restart_output_path, "w")
    lake_dim = nc%setDimension("lake", int(size(self%definition%ids), i4))
    var = nc%setVariable("lake_id", "i64", [lake_dim])
    call var%setData(self%definition%ids)
    var = nc%setVariable("lake_outlet_x", "f64", [lake_dim])
    call var%setData(self%definition%points%x)
    var = nc%setVariable("lake_outlet_y", "f64", [lake_dim])
    call var%setData(self%definition%points%y)
    var = nc%setVariable("lake_spill_reference_level", "f64", [lake_dim])
    call var%setAttribute("units", "m")
    call var%setData(self%definition%spill_reference_level)
    if (allocated(self%definition%spill_reference_surface_area)) then
      var = nc%setVariable("lake_spill_reference_surface_area", "f64", [lake_dim])
      call var%setAttribute("units", "m2")
      call var%setData(self%definition%spill_reference_surface_area)
    end if
    if (allocated(self%definition%spill_reference_volume)) then
      var = nc%setVariable("lake_spill_reference_volume", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%definition%spill_reference_volume)
    end if
    if (allocated(self%definition%spillway_length)) then
      var = nc%setVariable("lake_spillway_length", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%definition%spillway_length)
    end if
    ! The endorheic flag is written (as i4 0/1, since NetCDF has no logical type) so a restart-only
    ! run (no river, hence no topology to re-derive it from) restores it instead.
    if (allocated(self%definition%endorheic)) then
      allocate(endorheic_flag(size(self%definition%endorheic)))
      endorheic_flag = merge(1_i4, 0_i4, self%definition%endorheic)
      var = nc%setVariable("lake_endorheic", "i32", [lake_dim])
      call var%setAttribute("long_name", "1 if the lake is endorheic (terminal, no downstream river)")
      call var%setData(endorheic_flag)
      deallocate(endorheic_flag)
    end if
    if (allocated(self%definition%initial_level)) then
      var = nc%setVariable("lake_initial_level", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%definition%initial_level)
    end if
    ! The regionalized lakebed Ks is written so a restart-only run (no Input, no soil database)
    ! can restore it instead of re-deriving it from soil.
    if (allocated(self%definition%ks)) then
      var = nc%setVariable("lake_ks", "f64", [lake_dim])
      call var%setAttribute("units", "cm d-1")
      call var%setData(self%definition%ks)
    end if
    if (allocated(self%definition%kbf)) then
      var = nc%setVariable("lake_kbf", "f64", [lake_dim])
      call var%setAttribute("units", "d")
      call var%setData(self%definition%kbf)
    end if
    ! The DEM-footprint area is written so a restart-only run (no Input) can restore it for the
    ! lake banks reconciliation instead of re-deriving it from the lake footprint.
    if (allocated(self%definition%footprint_area)) then
      var = nc%setVariable("lake_footprint_area", "f64", [lake_dim])
      call var%setAttribute("units", "m2")
      call var%setData(self%definition%footprint_area)
    end if
    ! Geographic coordinates are written only when they were derived (i.e. a latlon file was present
    ! for a projected grid, or the points were spherical); they are required at read time only when
    ! varying_albedo is configured (checked in the connect guard).
    if (allocated(self%definition%latitude)) then
      var = nc%setVariable("lake_latitude", "f64", [lake_dim])
      call var%setAttribute("units", "degree")
      call var%setData(self%definition%latitude)
    end if
    if (allocated(self%definition%longitude)) then
      var = nc%setVariable("lake_longitude", "f64", [lake_dim])
      call var%setAttribute("units", "degree")
      call var%setData(self%definition%longitude)
    end if
    var = nc%setVariable("lake_coordsys", "i32", dims(:0))
    call var%setData(self%definition%points%coordsys)
    var = nc%setVariable("lake_outflow", "f64", [lake_dim])
    call var%setAttribute("units", "m3 s-1")
    call var%setData(self%state%outflow)
    if (self%exchange%config%processes%lake == 1_i4) then
      var = nc%setVariable("lake_level", "f64", [lake_dim])
      call var%setAttribute("units", "m")
      call var%setData(self%state%level)
      var = nc%setVariable("lake_surface_area", "f64", [lake_dim])
      call var%setAttribute("units", "m2")
      call var%setData(self%state%surface_area)
      var = nc%setVariable("lake_volume", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%state%volume)
      var = nc%setVariable("lake_gw_storage", "f64", [lake_dim])
      call var%setAttribute("units", "m3")
      call var%setData(self%state%gw_storage)
    end if
    if (self%restart_save_topology) then
      call lake_write_topology(nc, self%definition%grid, self%definition%map)
    end if
    if (self%exchange%config%processes%lake == 1_i4) then
      block
        type(NcDimension) :: parameter_dim
        parameter_dim = nc%setDimension("lake_ks_parameter", 3_i4)
        var = nc%setVariable("lake_generating_ks", "f64", [parameter_dim])
        call var%setData(self%runtime%ptf(1:3))
        parameter_dim = nc%setDimension("lake_kbf_parameter", size(self%runtime%recession))
        var = nc%setVariable("lake_generating_kbf", "f64", [parameter_dim])
        call var%setData(self%runtime%recession)
        parameter_dim = nc%setDimension("lake_length_parameter", 1_i4)
        var = nc%setVariable("lake_generating_length", "f64", [parameter_dim])
        call var%setData([self%runtime%multiplier])
      end block
    end if
    var = nc%setVariable("mlm_meta", "i32", dims(:0))
    call var%setData(0_i4)
    call var%setAttribute("format_version", 1_i4)
    call var%setAttribute("lake_case", self%exchange%config%processes%lake)
    call var%setAttribute("time_stamp", self%exchange%time%str())
    call nc%close()
  end subroutine mlm_create_restart

  !> \brief Read a complete current-format definition before resolving the exchange provider.
  subroutine mlm_read_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: var
    real(dp), allocatable :: outlet_x(:), outlet_y(:)
    integer(i4), allocatable :: endorheic_flag(:)
    integer(i4), allocatable :: coordsys_data
    integer(i4) :: version, lake_case
    integer(i8) :: i, n
    character(64) :: restart_time
    character(:), allocatable :: expected_time

    nc = NcDataset(self%restart_input_path, "r")
    call mlm_require_restart_fields(nc, [character(40) :: "mlm_meta"])
    var = nc%getVariable("mlm_meta")
    if (.not.var%hasAttribute("format_version")) then
      log_fatal(*) "mLM: unsupported lake restart format; regenerate the lake restart."
      error stop 1
    end if
    call var%getAttribute("format_version", version)
    if (version /= 1_i4) then
      log_fatal(*) "mLM: unsupported lake restart format; regenerate the lake restart."
      error stop 1
    end if
    if (.not.var%hasAttribute("lake_case") .or. .not.var%hasAttribute("time_stamp")) then
      log_fatal(*) "mLM: restart metadata requires lake_case and time_stamp."
      error stop 1
    end if
    call var%getAttribute("lake_case", lake_case)
    if (lake_case /= self%exchange%config%processes%lake) then
      log_fatal(*) "mLM: lake process case conflicts with restart."
      error stop 1
    end if
    call var%getAttribute("time_stamp", restart_time)
    expected_time = self%exchange%start_time%str()
    if (trim(restart_time) /= trim(expected_time)) then
      log_fatal(*) "mLM restart timestamp ", trim(restart_time), " does not match domain restart time ", trim(expected_time), "."
      error stop 1
    end if
    call mlm_require_restart_fields(nc, [character(40) :: "lake_id", "lake_outlet_x", "lake_outlet_y", &
      "lake_spill_reference_level", "lake_coordsys", "lake_outflow"])
    if (lake_case == 1_i4 .or. lake_case == -2_i4) &
      call mlm_require_restart_fields(nc, [character(40) :: "lake_footprint_area"])
    if (lake_case == 1_i4) then
      call mlm_require_restart_fields(nc, [character(40) :: "lake_spill_reference_surface_area", &
        "lake_spill_reference_volume", "lake_spillway_length", "lake_endorheic", "lake_initial_level", &
        "lake_ks", "lake_kbf", "lake_generating_ks", "lake_generating_kbf", "lake_generating_length", &
        "lake_level", "lake_surface_area", "lake_volume", "lake_gw_storage"])
    end if
    var = nc%getVariable("lake_id")
    call var%getData(self%restart%ids)
    n = size(self%restart%ids, kind=i8)
    if (n < 1_i8) then
      log_fatal(*) "mLM restart contains no lakes: ", self%restart_input_path
      error stop 1
    end if
    if (any(self%restart%ids <= 0_i8)) then
      log_fatal(*) "mLM restart contains a non-positive stable lake ID."
      error stop 1
    end if
    do i = 1_i8, n - 1_i8
      if (any(self%restart%ids(i + 1_i8:) == self%restart%ids(i))) then
        log_fatal(*) "mLM restart contains duplicate stable lake IDs."
        error stop 1
      end if
    end do
    call mlm_read_restart_field(nc, "lake_outlet_x", outlet_x, n)
    call mlm_read_restart_field(nc, "lake_outlet_y", outlet_y, n)
    call mlm_read_restart_field(nc, "lake_spill_reference_level", self%restart%spill_reference_level, n)
    if (lake_case == 1_i4 .or. lake_case == -2_i4) &
      call mlm_read_restart_field(nc, "lake_footprint_area", self%restart%footprint_area, n, positive=.true.)
    if (lake_case == 1_i4) then
      call mlm_read_restart_field(nc, "lake_spill_reference_surface_area", &
        self%restart%spill_reference_surface_area, n, positive=.true.)
      call mlm_read_restart_field(nc, "lake_spill_reference_volume", self%restart%spill_reference_volume, n, positive=.true.)
      call mlm_read_restart_field(nc, "lake_spillway_length", self%restart%spillway_length, n, positive=.true.)
      call mlm_read_restart_field(nc, "lake_initial_level", self%restart%initial_level, n)
      call mlm_read_restart_field(nc, "lake_ks", self%restart%ks, n, positive=.true.)
      call mlm_read_restart_field(nc, "lake_kbf", self%restart%kbf, n, positive=.true.)
      var = nc%getVariable("lake_endorheic")
      call var%getData(endorheic_flag)
      if (size(endorheic_flag, kind=i8) /= n) then
        log_fatal(*) "mLM restart endorheic flag does not align with stable lake IDs."
        error stop 1
      end if
      if (any(endorheic_flag /= 0_i4 .and. endorheic_flag /= 1_i4)) then
        log_fatal(*) "mLM restart: invalid endorheic flags."
        error stop 1
      end if
      self%restart%endorheic = endorheic_flag /= 0_i4
      var = nc%getVariable("lake_generating_ks")
      call var%getData(self%restart%generating_ks)
      var = nc%getVariable("lake_generating_kbf")
      call var%getData(self%restart%generating_kbf)
      var = nc%getVariable("lake_generating_length")
      call var%getData(self%restart%generating_length)
    end if
    ! Geographic coordinates remain optional for runs without varying albedo.
    if (nc%hasVariable("lake_latitude")) &
      call mlm_read_restart_field(nc, "lake_latitude", self%restart%latitude, n)
    if (nc%hasVariable("lake_longitude")) &
      call mlm_read_restart_field(nc, "lake_longitude", self%restart%longitude, n)
    var = nc%getVariable("lake_coordsys")
    call var%getData(coordsys_data)
    if (coordsys_data /= cartesian .and. coordsys_data /= spherical) then
      log_fatal(*) "mLM restart contains an unsupported lake coordinate system."
      error stop 1
    end if
    call self%restart%points%init(outlet_x, outlet_y, coordsys=coordsys_data)
    call nc%close()
  end subroutine mlm_read_restart_metadata

  !> \brief Reject incomplete restart definitions instead of reconstructing older formats.
  subroutine mlm_require_restart_fields(nc, names)
    type(NcDataset), intent(in) :: nc !< open current-format restart
    character(*), intent(in) :: names(:) !< fields required by the active lake case
    integer :: i
    do i = 1, size(names)
      if (.not.nc%hasVariable(trim(names(i)))) then
        log_fatal(*) "mLM: restart is missing required field ", trim(names(i)), "; regenerate the lake restart."
        error stop 1
      end if
    end do
  end subroutine mlm_require_restart_fields

  !> \brief Read and validate one per-lake definition field in stable-ID order.
  subroutine mlm_read_restart_field(nc, name, data, n, positive)
    type(NcDataset), intent(inout) :: nc !< open current-format restart
    character(*), intent(in) :: name !< required field name
    real(dp), allocatable, intent(out) :: data(:) !< owned definition storage, not yet published
    integer(i8), intent(in) :: n !< number of stable lake IDs
    logical, optional, intent(in) :: positive !< require strictly positive values
    type(NcVariable) :: var
    var = nc%getVariable(name)
    call var%getData(data)
    if (size(data, kind=i8) /= n) then
      log_fatal(*) "mLM restart field ", name, " does not align with stable lake IDs."
      error stop 1
    end if
    if (.not.all(ieee_is_finite(data))) then
      log_fatal(*) "mLM restart field ", name, " must be finite."
      error stop 1
    end if
    if (present(positive)) then
      if (positive) then
        if (any(data <= 0.0_dp)) then
          log_fatal(*) "mLM restart field ", name, " must be positive."
          error stop 1
        end if
      end if
    end if
  end subroutine mlm_read_restart_field

  !> \brief Restore optional static level-0 lake topology when present.
  subroutine mlm_read_restart_topology(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    nc = NcDataset(self%restart_input_path, "r")
    if (nc%hasVariable("lake_map")) call lake_read_topology(nc, self%definition%grid, self%definition%map)
    call nc%close()
  end subroutine mlm_read_restart_topology

  !> \brief Verify restart metadata against Input metadata without imposing restart-file ordering.
  subroutine mlm_validate_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: i, j, n_lakes

    n_lakes = self%exchange%lake_points%n_points
    if (size(self%restart%ids, kind=i8) /= n_lakes) then
      log_fatal(*) "mLM restart lake IDs do not match Input lake metadata."
      error stop 1
    end if
    if (self%restart%points%coordsys /= self%exchange%lake_points%coordsys) then
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
      do j = 1_i8, size(self%restart%ids, kind=i8)
        if (self%restart%ids(j) == self%exchange%lake_ids%data(i)) exit
      end do
      if (j > size(self%restart%ids, kind=i8)) then
        log_fatal(*) "mLM restart and Input lake metadata have different stable lake IDs."
        error stop 1
      end if
      if (.not.is_close(self%restart%points%x(j), self%exchange%lake_points%x(i)) .or. &
          .not.is_close(self%restart%points%y(j), self%exchange%lake_points%y(i)) .or. &
          .not.is_close(self%restart%spill_reference_level(j), self%exchange%lake_spill_reference_level%data(i))) then
        log_fatal(*) "mLM restart and Input lake metadata differ for stable lake ID."
        error stop 1
      end if
      if (allocated(self%restart%spill_reference_surface_area) .and. &
          self%exchange%lake_spill_reference_surface_area%has_data()) then
        if (.not.is_close(self%restart%spill_reference_surface_area(j), &
                          self%exchange%lake_spill_reference_surface_area%data(i))) then
          log_fatal(*) "mLM restart and Input reported maximum surface area differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart%spill_reference_volume) .and. &
          self%exchange%lake_spill_reference_volume%has_data()) then
        if (.not.is_close(self%restart%spill_reference_volume(j), &
                          self%exchange%lake_spill_reference_volume%data(i))) then
          log_fatal(*) "mLM restart and Input reported maximum volume differ for stable lake ID."
          error stop 1
        end if
      end if
      ! The initial level is a cold-start-only field: on a continuation the lake level comes from the
      ! restart, not from this input, and the restart carries the RESOLVED level (sentinel inputs are
      ! resolved to real levels at the first initialize) while the input may still hold the -9999
      ! sentinel. Comparing the two is meaningless and wrongly aborts a valid continuation, so it is
      ! not validated here.
      if (allocated(self%restart%latitude) .and. &
          self%exchange%lake_latitude%has_data()) then
        if (.not.is_close(self%restart%latitude(j), &
                          self%exchange%lake_latitude%data(i))) then
          log_fatal(*) "mLM restart and Input latitude differ for stable lake ID."
          error stop 1
        end if
      end if
      if (allocated(self%restart%longitude) .and. &
          self%exchange%lake_longitude%has_data()) then
        if (.not.is_close(self%restart%longitude(j), &
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
    if (size(ids, kind=i8) /= size(self%definition%ids, kind=i8) .or. size(values, kind=i8) /= size(ids, kind=i8)) then
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
    do i = 1_i8, size(self%definition%ids, kind=i8)
      do j = 1_i8, size(ids, kind=i8)
        if (ids(j) == self%definition%ids(i)) exit
      end do
      if (j > size(ids, kind=i8)) then
        log_fatal(*) "mLM restart is missing stable lake ID."
        error stop 1
      end if
      self%state%outflow(i) = values(j)
      if (is_case1) then
        self%state%level(i) = level(j)
        self%state%surface_area(i) = area(j)
        self%state%volume(i) = volume(j)
        self%state%gw_storage(i) = gw_storage(j)
      end if
    end do
  end subroutine mlm_read_restart_state

end module mo_mlm_container
