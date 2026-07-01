!> \file    mo_mhm_container.f90
!> \brief   \copybrief mo_mhm_container
!> \details \copydetails mo_mhm_container

!> \brief   Module for a mHM process container.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Aug 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_mhm_container
  use mo_logging
  use mo_canopy_interc, only: canopy_interc
  use mo_constants, only: nodata_dp
  use mo_grid_io, only: output_dataset, var, daily, monthly, yearly, no_time, time_units_delta
  use mo_runoff, only: L1_total_runoff, runoff_sat_zone, runoff_unsat_zone
  use mo_snow_accum_melt, only: snow_accum_melt
  use mo_soil_moisture, only: soil_moisture_direct_runoff, soil_moisture_pervious
  use mo_kind, only: i4, i8, dp
  use mo_constants, only: YearMonths
  use mo_common_constants, only: P1_InitStateFluxes, soilHorizonsVarName
  use mo_exchange_type, only: exchange_t
  use mo_datetime, only: one_hour
  use mo_grid, only: grid_t, cartesian
  use mo_mhm_constants, only: P2_InitStateFluxes, P3_InitStateFluxes, P4_InitStateFluxes
  use mo_neutrons, only: COSMIC, DesiletsN0, TabularIntegralAFast
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_string_utils, only: n2s => num2str
  use mo_utils, only: is_close
  use nml_config_mhm, only: nml_config_mhm_t, NML_OK
  use nml_output_mhm, only: nml_output_mhm_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  !> \class   mhm_canopy_state_t
  !> \brief   Grouped canopy process fields.
  type :: mhm_canopy_state_t
    real(dp), allocatable :: interception(:) !< canopy interception storage
    real(dp), allocatable :: throughfall(:) !< throughfall flux
    real(dp), allocatable :: aet(:) !< actual evapotranspiration from canopy
  end type mhm_canopy_state_t

  !> \class   mhm_snow_state_t
  !> \brief   Grouped snow process fields.
  type :: mhm_snow_state_t
    real(dp), allocatable :: snowpack(:) !< snowpack storage
    real(dp), allocatable :: melt(:) !< melting snow depth
    real(dp), allocatable :: pre_effect(:) !< effective precipitation depth
    real(dp), allocatable :: rain(:) !< rain precipitation depth
    real(dp), allocatable :: snow(:) !< snow precipitation depth
    real(dp), allocatable :: degday(:) !< degree-day factor
  end type mhm_snow_state_t

  !> \class   mhm_direct_runoff_state_t
  !> \brief   Grouped sealed-area/direct-runoff process fields.
  type :: mhm_direct_runoff_state_t
    real(dp), allocatable :: storage(:) !< retention storage of impervious areas
    real(dp), allocatable :: aet(:) !< actual evapotranspiration from free-water surfaces
    real(dp), allocatable :: runoff(:) !< direct runoff from impervious areas
  end type mhm_direct_runoff_state_t

  !> \class   mhm_soil_state_t
  !> \brief   Grouped soil-moisture process fields.
  type :: mhm_soil_state_t
    real(dp), allocatable :: horizon_bounds(:) !< soil-horizon boundary depths copied from MPR metadata
    real(dp), allocatable :: moisture(:, :) !< soil moisture by horizon
    real(dp), allocatable :: infiltration(:, :) !< infiltration intensity by horizon
    real(dp), allocatable :: aet(:, :) !< actual ET by soil horizon
  end type mhm_soil_state_t

  !> \class   mhm_runoff_state_t
  !> \brief   Grouped runoff/baseflow process fields.
  type :: mhm_runoff_state_t
    real(dp), allocatable :: unsat_storage(:) !< upper soil storage
    real(dp), allocatable :: sat_storage(:) !< groundwater storage
    real(dp), allocatable :: percolation(:) !< percolation flux
    real(dp), allocatable :: fast_interflow(:) !< fast runoff component
    real(dp), allocatable :: slow_interflow(:) !< slow runoff component
    real(dp), allocatable :: baseflow(:) !< baseflow
    real(dp), allocatable :: total_runoff(:) !< generated runoff
  end type mhm_runoff_state_t

  !> \class   mhm_neutron_state_t
  !> \brief   Grouped neutron process fields and static support data.
  type :: mhm_neutron_state_t
    real(dp), allocatable :: counts(:) !< simulated ground albedo neutrons
    real(dp), allocatable :: integral_afast(:) !< pre-calculated integrand for neutron flux approximation
  end type mhm_neutron_state_t

  !> \class   mhm_forcing_state_t
  !> \brief   Grouped current-step forcing caches and monthly evaporation coefficients.
  type :: mhm_forcing_state_t
    real(dp), allocatable :: pet(:) !< current-step potential evapotranspiration
    real(dp), allocatable :: temp(:) !< current-step temperature
    real(dp), allocatable :: prec(:) !< current-step precipitation
    real(dp) :: evap_coeff(int(YearMonths, i4)) = 0.0_dp !< monthly evaporation coefficients for free-water surfaces
  end type mhm_forcing_state_t

  !> \class   mhm_contract_state_t
  !> \brief   Internal ownership flags for mHM subprocess exchange outputs.
  type :: mhm_contract_state_t
    logical :: own_interception = .false. !< whether mHM publishes canopy interception storage
    logical :: own_throughfall = .false. !< whether mHM publishes throughfall
    logical :: own_aet_canopy = .false. !< whether mHM publishes canopy actual evapotranspiration
    logical :: own_snowpack = .false. !< whether mHM publishes snowpack storage
    logical :: own_melt = .false. !< whether mHM publishes snow melt
    logical :: own_pre_effect = .false. !< whether mHM publishes effective precipitation
    logical :: own_rain = .false. !< whether mHM publishes liquid precipitation
    logical :: own_snow = .false. !< whether mHM publishes solid precipitation
    logical :: own_degday = .false. !< whether mHM publishes degree-day factor
    logical :: own_sealed_storage = .false. !< whether mHM publishes sealed-area storage
    logical :: own_aet_sealed = .false. !< whether mHM publishes sealed-area actual evapotranspiration
    logical :: own_runoff_sealed = .false. !< whether mHM publishes direct runoff
    logical :: own_soil_moisture = .false. !< whether mHM publishes soil moisture
    logical :: own_infiltration = .false. !< whether mHM publishes soil infiltration
    logical :: own_aet_soil = .false. !< whether mHM publishes soil actual evapotranspiration
    logical :: own_unsat_storage = .false. !< whether mHM publishes unsaturated-zone storage
    logical :: own_sat_storage = .false. !< whether mHM publishes saturated-zone storage
    logical :: own_interflow_fast = .false. !< whether mHM publishes fast interflow
    logical :: own_interflow_slow = .false. !< whether mHM publishes slow interflow
    logical :: own_percolation = .false. !< whether mHM publishes percolation
    logical :: own_baseflow = .false. !< whether mHM publishes baseflow
    logical :: own_total_runoff = .false. !< whether mHM publishes total runoff
    logical :: own_neutrons = .false. !< whether mHM publishes neutron counts
  end type mhm_contract_state_t

  !> \class   mhm_io_state_t
  !> \brief   Grouped restart/output bookkeeping owned by the mHM container.
  type :: mhm_io_state_t
    logical :: read_restart = .false. !< whether to read restart data
    logical :: write_restart = .false. !< whether to write restart data
    logical :: output_active = .false. !< whether gridded mHM output is enabled
    logical :: calc_f_not_sealed = .false. !< whether output needs the derived non-sealed fraction
    character(:), allocatable :: restart_input_path !< resolved restart input path
    character(:), allocatable :: restart_output_path !< resolved restart output path
    character(:), allocatable :: output_path !< resolved output path
  end type mhm_io_state_t

  !> \class   mhm_runtime_state_t
  !> \brief   Grouped scalar runtime bookkeeping for the mHM container.
  type :: mhm_runtime_state_t
    real(dp) :: c2TSTu = 0.0_dp !< conversion factor from model time step to legacy day-based units
  end type mhm_runtime_state_t

  !> \class   mhm_t
  !> \brief   Class for a single mHM process container.
  type, public :: mhm_t
    type(nml_config_mhm_t) :: config !< configuration of the mHM process container
    type(nml_output_mhm_t) :: output_config !< output configuration of the mHM process container
    type(exchange_t), pointer :: exchange => null() !< exchange container of the domain
    type(output_dataset) :: ds_out !< output dataset for gridded mHM output
    type(mhm_canopy_state_t) :: canopy !< canopy process fields
    type(mhm_snow_state_t) :: snow !< snow process fields
    type(mhm_direct_runoff_state_t) :: direct_runoff !< direct-runoff/sealed-area process fields
    type(mhm_soil_state_t) :: soil !< soil-moisture process fields
    type(mhm_runoff_state_t) :: runoff !< runoff/baseflow process fields
    type(mhm_neutron_state_t) :: neutrons !< neutron process fields and static support data
    type(mhm_forcing_state_t) :: forcing !< forcing caches and monthly evap coefficients
    type(mhm_contract_state_t) :: contract !< internal ownership of couplable subprocess outputs
    type(mhm_io_state_t) :: io !< restart/output bookkeeping
    type(mhm_runtime_state_t) :: runtime !< scalar runtime bookkeeping
  contains
    procedure :: set_dims => mhm_set_dims
    procedure :: configure => mhm_configure
    procedure :: connect => mhm_connect
    procedure :: initialize => mhm_initialize
    procedure :: update => mhm_update
    procedure :: finalize => mhm_finalize
    procedure :: update_interception => mhm_update_interception
    procedure :: update_snow => mhm_update_snow
    procedure :: update_direct_runoff => mhm_update_direct_runoff
    procedure :: update_soil_moisture => mhm_update_soil_moisture
    procedure :: update_interflow => mhm_update_interflow
    procedure :: update_baseflow => mhm_update_baseflow
    procedure :: update_total_runoff => mhm_update_total_runoff
    procedure :: update_neutrons => mhm_update_neutrons
    procedure, private :: reset_fields => mhm_reset_fields
    procedure, private :: filter_output => mhm_filter_output
    procedure, private :: publish_exchange => mhm_publish_exchange
    procedure, private :: clear_exchange => mhm_clear_exchange
    procedure, private :: copy_horizon_bounds => mhm_copy_horizon_bounds
    procedure, private :: create_output => mhm_create_output
    procedure, private :: update_output => mhm_update_output
    procedure, private :: create_restart => mhm_create_restart
    procedure, private :: write_restart_data => mhm_write_restart_data
    procedure, private :: read_restart_data => mhm_read_restart_data
    procedure, private :: validate_restart_grid => mhm_validate_restart_grid
    procedure, private :: validate_restart_horizon_bounds => mhm_validate_restart_horizon_bounds
    procedure, private :: validate_restart_process_cases => mhm_validate_restart_process_cases
    procedure, private :: write_restart_field_2d => mhm_write_restart_field_2d
    procedure, private :: write_restart_field_3d => mhm_write_restart_field_3d
    procedure, private :: read_restart_field_2d => mhm_read_restart_field_2d
    procedure, private :: read_restart_field_3d => mhm_read_restart_field_3d
  end type mhm_t

contains

  !> \brief Set runtime dimensions for generated mHM namelists.
  subroutine mhm_set_dims(self)
    class(mhm_t), intent(inout), target :: self
    character(1024) :: errmsg
    integer :: status

    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Error setting mHM config dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mhm_set_dims

  !> \brief Configure the mHM process container.
  subroutine mhm_configure(self, file, out_file)
    class(mhm_t), intent(inout), target :: self
    character(*), intent(in), optional :: file !< file containing the namelists
    character(*), intent(in), optional :: out_file !< file containing the output namelists
    integer(i4) :: id(1)
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer :: status

    log_info(*) "Configure mhm"
    if (present(file)) then
      path = self%exchange%get_path(file)
      log_info(*) "Read mhm config: ", path
      status = self%config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading mHM config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%is_configured) then
      log_fatal(*) "mHM configuration not set."
      error stop 1
    end if
    status = self%config%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mHM config not valid: ", trim(errmsg)
      error stop 1
    end if

    id(1) = self%exchange%nml_domain_id

    self%io%read_restart = self%config%read_restart(id(1))
    self%io%write_restart = self%config%write_restart(id(1))

    self%io%output_active = .true.
    if (present(out_file)) then
      path = self%exchange%get_path(out_file, root=.true.)
      log_info(*) "Read mHM output config: ", path
      status = self%output_config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        self%io%output_active = .false.
        log_warn(*) "mHM output disabled, config not found: ", trim(errmsg)
      end if
    end if
    if (self%output_config%is_configured) then
      status = self%output_config%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mHM output config invalid: ", trim(errmsg)
        error stop 1
      end if
    else
      self%io%output_active = .false.
      log_warn(*) "mHM output disabled, config not set."
    end if

    if (self%io%output_active) then
      status = self%config%is_set("output_path", idx=id, errmsg=errmsg)
      self%io%output_active = self%io%output_active .and. (status == NML_OK)
      if (status /= NML_OK) then
        log_warn(*) "mHM output disabled, path not set for domain ", n2s(id(1)), ": ", trim(errmsg)
      else
        self%io%output_path = self%exchange%get_path(self%config%output_path(id(1)))
      end if
    end if
  end subroutine mhm_configure

  !> \brief Connect the mHM process container with other components.
  subroutine mhm_connect(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: id(1)
    integer(i4) :: n_cells
    integer(i4) :: n_horizons
    integer(i8) :: expected_shape_1d(1)
    integer(i8) :: expected_shape_2d(2)
    integer(i4) :: interception_case
    integer(i4) :: snow_case
    integer(i4) :: soil_case
    integer(i4) :: direct_runoff_case
    integer(i4) :: interflow_case
    integer(i4) :: percolation_case
    integer(i4) :: routing_case
    integer(i4) :: baseflow_case
    integer(i4) :: neutron_case
    logical :: allocate_throughfall
    logical :: allocate_pre_effect
    logical :: allocate_rain
    logical :: need_f_sealed
    integer :: status
    character(1024) :: errmsg

    log_info(*) "Connect mhm"

    id(1) = self%exchange%nml_domain_id
    self%io%read_restart = self%config%read_restart(id(1))
    self%io%write_restart = self%config%write_restart(id(1))
    if (self%io%read_restart) then
      status = self%config%is_set("restart_input_path", idx=id, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mHM restart input path not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg)
        error stop 1
      end if
      self%io%restart_input_path = self%exchange%get_path(self%config%restart_input_path(id(1)))
    end if
    if (self%io%write_restart) then
      status = self%config%is_set("restart_output_path", idx=id, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mHM restart output path not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg)
        error stop 1
      end if
      self%io%restart_output_path = self%exchange%get_path(self%config%restart_output_path(id(1)))
    end if

    if (.not.associated(self%exchange%level1)) then
      log_fatal(*) "mHM: level1 grid not available (check MPR/meteo setup)."
      error stop 1
    end if
    if (.not.associated(self%exchange%soil_horizon_bounds)) then
      log_fatal(*) "mHM: soil_horizon_bounds not provided (check MPR setup)."
      error stop 1
    end if

    call self%copy_horizon_bounds(self%exchange%soil_horizon_bounds)
    n_cells = int(self%exchange%level1%ncells, i4)
    n_horizons = size(self%soil%horizon_bounds) - 1_i4
    expected_shape_1d = [self%exchange%level1%ncells]
    expected_shape_2d = [self%exchange%level1%ncells, int(n_horizons, i8)]
    if (n_horizons < 1_i4) then
      log_fatal(*) "mHM: soil_horizon_bounds must contain at least two entries."
      error stop 1
    end if

    interception_case = self%exchange%parameters%process_matrix(1, 1)
    snow_case = self%exchange%parameters%process_matrix(2, 1)
    soil_case = self%exchange%parameters%process_matrix(3, 1)
    direct_runoff_case = self%exchange%parameters%process_matrix(4, 1)
    interflow_case = self%exchange%parameters%process_matrix(6, 1)
    percolation_case = self%exchange%parameters%process_matrix(7, 1)
    routing_case = self%exchange%parameters%process_matrix(8, 1)
    baseflow_case = self%exchange%parameters%process_matrix(9, 1)
    neutron_case = self%exchange%parameters%process_matrix(10, 1)

    if ((interflow_case == 0_i4) .neqv. (percolation_case == 0_i4)) then
      log_fatal(*) "mHM: interflow and percolation must currently be activated together on the exchange path."
      error stop 1
    end if

    self%contract%own_interception = interception_case /= 0_i4
    self%contract%own_throughfall = interception_case /= 0_i4
    self%contract%own_aet_canopy = interception_case /= 0_i4
    self%contract%own_snowpack = snow_case /= 0_i4
    self%contract%own_melt = snow_case /= 0_i4
    self%contract%own_pre_effect = snow_case /= 0_i4
    self%contract%own_rain = snow_case /= 0_i4
    self%contract%own_snow = snow_case /= 0_i4
    self%contract%own_degday = snow_case /= 0_i4
    self%contract%own_sealed_storage = direct_runoff_case /= 0_i4
    self%contract%own_aet_sealed = direct_runoff_case /= 0_i4
    self%contract%own_runoff_sealed = direct_runoff_case /= 0_i4
    self%contract%own_soil_moisture = soil_case /= 0_i4
    self%contract%own_infiltration = soil_case /= 0_i4
    self%contract%own_aet_soil = soil_case /= 0_i4
    self%contract%own_unsat_storage = interflow_case /= 0_i4
    self%contract%own_sat_storage = interflow_case /= 0_i4 .or. baseflow_case /= 0_i4
    self%contract%own_interflow_fast = interflow_case /= 0_i4
    self%contract%own_interflow_slow = interflow_case /= 0_i4
    self%contract%own_percolation = percolation_case /= 0_i4
    self%contract%own_baseflow = baseflow_case /= 0_i4
    self%contract%own_total_runoff = direct_runoff_case /= 0_i4 .or. interflow_case /= 0_i4 .or. baseflow_case /= 0_i4
    self%contract%own_neutrons = neutron_case /= 0_i4
    allocate_throughfall = interception_case == 1_i4
    allocate_pre_effect = snow_case == 1_i4
    allocate_rain = snow_case == 1_i4

    call self%filter_output()

    need_f_sealed = &
           (soil_case /= 0_i4) &
      .or. (direct_runoff_case /= 0_i4) &
      .or. (interflow_case /= 0_i4) &
      .or. (baseflow_case /= 0_i4) &
      .or. self%io%calc_f_not_sealed &
      .or. self%output_config%out_qd &
      .or. (self%output_config%out_aet_all .and. self%exchange%aet_sealed%available(owned=self%contract%own_aet_sealed))

    call self%exchange%pre%require("mHM", (interception_case /= 0_i4) .or. (snow_case == 1_i4), expected_shape_1d)
    call self%exchange%temp%require("mHM", snow_case == 1_i4, expected_shape_1d)
    call self%exchange%pet%require("mHM", &
      (interception_case == 1_i4) .or. (soil_case /= 0_i4) .or. (direct_runoff_case /= 0_i4) .or. &
      self%output_config%out_pet, expected_shape_1d)

    call self%exchange%max_interception%require("mHM", interception_case == 1_i4, expected_shape_1d)
    call self%exchange%thresh_temp%require("mHM", snow_case == 1_i4, expected_shape_1d)
    call self%exchange%degday_dry%require("mHM", snow_case == 1_i4, expected_shape_1d)
    call self%exchange%degday_inc%require("mHM", snow_case == 1_i4, expected_shape_1d)
    call self%exchange%degday_max%require("mHM", snow_case == 1_i4, expected_shape_1d)

    call self%exchange%f_sealed%require("mHM", need_f_sealed ,expected_shape_1d)
    call self%exchange%f_roots%require("mHM", soil_case /= 0_i4, expected_shape_2d)
    call self%exchange%sm_saturation%require("mHM", &
      (soil_case /= 0_i4) .or. self%output_config%out_sm .or. self%output_config%out_sm_all, expected_shape_2d)
    call self%exchange%sm_exponent%require("mHM", soil_case /= 0_i4, expected_shape_2d)
    call self%exchange%sm_field_capacity%require("mHM", soil_case /= 0_i4, expected_shape_2d)
    call self%exchange%wilting_point%require("mHM", soil_case /= 0_i4, expected_shape_2d)
    call self%exchange%thresh_jarvis%require("mHM", (soil_case == 2_i4) .or. (soil_case == 3_i4), expected_shape_1d)
    call self%exchange%thresh_sealed%require("mHM", direct_runoff_case /= 0_i4, expected_shape_1d)

    call self%exchange%alpha%require("mHM", interflow_case /= 0_i4, expected_shape_1d)
    call self%exchange%k_fastflow%require("mHM", interflow_case /= 0_i4, expected_shape_1d)
    call self%exchange%k_slowflow%require("mHM", interflow_case /= 0_i4, expected_shape_1d)
    call self%exchange%k_percolation%require("mHM", percolation_case /= 0_i4, expected_shape_1d)
    call self%exchange%f_karst_loss%require("mHM", percolation_case /= 0_i4, expected_shape_1d)
    call self%exchange%thresh_unsat%require("mHM", interflow_case /= 0_i4, expected_shape_1d)
    call self%exchange%k_baseflow%require("mHM", baseflow_case /= 0_i4, expected_shape_1d)

    call self%exchange%desilets_n0%require("mHM", neutron_case /= 0_i4, expected_shape_1d)
    call self%exchange%bulk_density%require("mHM", neutron_case /= 0_i4, expected_shape_2d)
    call self%exchange%lattice_water%require("mHM", neutron_case /= 0_i4, expected_shape_2d)
    call self%exchange%cosmic_l3%require("mHM", neutron_case == 2_i4, expected_shape_2d)

    call allocate_1d(self%canopy%interception, n_cells, self%contract%own_interception)
    call allocate_1d(self%canopy%throughfall, n_cells, allocate_throughfall)
    call allocate_1d(self%canopy%aet, n_cells, self%contract%own_aet_canopy)

    call allocate_1d(self%snow%snowpack, n_cells, self%contract%own_snowpack)
    call allocate_1d(self%snow%melt, n_cells, self%contract%own_melt)
    call allocate_1d(self%snow%pre_effect, n_cells, allocate_pre_effect)
    call allocate_1d(self%snow%rain, n_cells, allocate_rain)
    call allocate_1d(self%snow%snow, n_cells, self%contract%own_snow)
    call allocate_1d(self%snow%degday, n_cells, self%contract%own_degday)

    call allocate_1d(self%direct_runoff%storage, n_cells, self%contract%own_sealed_storage)
    call allocate_1d(self%direct_runoff%aet, n_cells, self%contract%own_aet_sealed)
    call allocate_1d(self%direct_runoff%runoff, n_cells, self%contract%own_runoff_sealed)

    call allocate_2d(self%soil%moisture, n_cells, n_horizons, self%contract%own_soil_moisture)
    call allocate_2d(self%soil%infiltration, n_cells, n_horizons, self%contract%own_infiltration)
    call allocate_2d(self%soil%aet, n_cells, n_horizons, self%contract%own_aet_soil)

    call allocate_1d(self%runoff%unsat_storage, n_cells, self%contract%own_unsat_storage)
    call allocate_1d(self%runoff%sat_storage, n_cells, self%contract%own_sat_storage)
    call allocate_1d(self%runoff%percolation, n_cells, self%contract%own_percolation)
    call allocate_1d(self%runoff%fast_interflow, n_cells, self%contract%own_interflow_fast)
    call allocate_1d(self%runoff%slow_interflow, n_cells, self%contract%own_interflow_slow)
    call allocate_1d(self%runoff%baseflow, n_cells, self%contract%own_baseflow)
    call allocate_1d(self%runoff%total_runoff, n_cells, self%contract%own_total_runoff)

    call allocate_1d(self%neutrons%counts, n_cells, self%contract%own_neutrons)

    call self%publish_exchange()
    call self%exchange%interception%expect_handoff("mHM", self%output_config%out_interception .or. (neutron_case /= 0_i4), expected_shape_1d)
    call self%exchange%throughfall%expect_handoff("mHM", snow_case /= 0_i4, expected_shape_1d)
    call self%exchange%aet_canopy%expect_handoff("mHM", (direct_runoff_case /= 0_i4) .or. (soil_case /= 0_i4), expected_shape_1d)
    call self%exchange%snowpack%expect_handoff("mHM", self%output_config%out_snowpack .or. (neutron_case /= 0_i4), expected_shape_1d)
    call self%exchange%melt%expect_handoff("mHM", self%output_config%out_qsm, expected_shape_1d)
    call self%exchange%pre_eff%expect_handoff("mHM", (direct_runoff_case /= 0_i4) .or. (soil_case /= 0_i4) .or. self%output_config%out_preeffect, expected_shape_1d)
    call self%exchange%sealed_storage%expect_handoff("mHM", self%output_config%out_sealedstw, expected_shape_1d)
    call self%exchange%infiltration%expect_handoff("mHM", (interflow_case /= 0_i4) .or. self%output_config%out_soil_infil, expected_shape_2d)
    call self%exchange%soil_moisture%expect_handoff("mHM", self%output_config%out_swc .or. self%output_config%out_sm .or. self%output_config%out_sm_all .or. (neutron_case /= 0_i4), expected_shape_2d)
    call self%exchange%aet_soil%expect_handoff("mHM", self%output_config%out_aet_layer, expected_shape_2d)
    call self%exchange%runoff_sealed%expect_handoff("mHM", self%output_config%out_qd, expected_shape_1d)
    call self%exchange%unsat_storage%expect_handoff("mHM", self%output_config%out_unsatstw, expected_shape_1d)
    call self%exchange%sat_storage%expect_handoff("mHM", self%output_config%out_satstw, expected_shape_1d)
    call self%exchange%interflow_fast%expect_handoff("mHM", self%output_config%out_qif, expected_shape_1d)
    call self%exchange%interflow_slow%expect_handoff("mHM", self%output_config%out_qis, expected_shape_1d)
    call self%exchange%percolation%expect_handoff("mHM", self%output_config%out_recharge, expected_shape_1d)
    call self%exchange%baseflow%expect_handoff("mHM", self%output_config%out_qb, expected_shape_1d)
    call self%exchange%runoff_total%expect_handoff("mHM", (routing_case /= 0_i4) .or. self%output_config%out_q, expected_shape_1d)
    call self%exchange%neutrons%expect_handoff("mHM", self%output_config%out_neutrons, expected_shape_1d)
  end subroutine mhm_connect

  !> \brief Initialize the mHM process container for the simulation.
  subroutine mhm_initialize(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: id(1)
    integer(i4) :: coeff_domain
    integer(i4) :: direct_runoff_case
    integer(i4) :: neutron_case

    log_info(*) "Initialize mhm"

    if (.not.allocated(self%soil%horizon_bounds)) then
      log_fatal(*) "mHM: soil horizon metadata not initialized before initialize."
      error stop 1
    end if

    coeff_domain = self%exchange%nml_domain_id
    id(1) = self%exchange%nml_domain_id
    if (self%config%share_evap_coeff) coeff_domain = 1_i4
    direct_runoff_case = self%exchange%parameters%process_matrix(4, 1)
    if (direct_runoff_case /= 0_i4) then
      self%forcing%evap_coeff = self%config%evap_coeff(:, coeff_domain)
      if (any(.not.ieee_is_finite(self%forcing%evap_coeff))) then
        log_fatal(*) "mHM: evap_coeff contains invalid values for domain ", n2s(coeff_domain), "."
        error stop 1
      end if
    else
      self%forcing%evap_coeff = 1.0_dp
    end if

    self%runtime%c2TSTu = real(int(self%exchange%step / one_hour(), i4), dp) / 24.0_dp
    if (.not.ieee_is_finite(self%runtime%c2TSTu) .or. self%runtime%c2TSTu <= 0.0_dp) then
      log_fatal(*) "mHM: invalid time-step conversion factor c2TSTu."
      error stop 1
    end if

    neutron_case = self%exchange%parameters%process_matrix(10, 1)
    if (allocated(self%neutrons%integral_afast)) deallocate(self%neutrons%integral_afast)
    if (neutron_case == 2_i4) then
      allocate(self%neutrons%integral_afast(10000 + 2))
      call TabularIntegralAFast(self%neutrons%integral_afast, 20.0_dp)
    else
      allocate(self%neutrons%integral_afast(1))
      self%neutrons%integral_afast = 0.0_dp
    end if

    call self%reset_fields()
    if (self%io%read_restart) call self%read_restart_data()

    call self%create_output()
  end subroutine mhm_initialize

  !> \brief Update the mHM process container for the current time step.
  subroutine mhm_update(self)
    class(mhm_t), intent(inout), target :: self
    log_trace(*) "Update mhm at step ", n2s(self%exchange%step_count)
    !$omp parallel default(shared)
    call self%update_interception()
    call self%update_snow()
    call self%update_direct_runoff()
    call self%update_soil_moisture()
    call self%update_interflow()
    call self%update_baseflow()
    call self%update_total_runoff()
    call self%update_neutrons()
    !$omp end parallel
    call self%update_output()
  end subroutine mhm_update

  !> \brief Create the gridded mHM output dataset from the configured output flags.
  subroutine mhm_create_output(self)
    class(mhm_t), intent(inout), target :: self
    type(var), allocatable :: vars(:)
    character(:), allocatable :: dtype
    character(:), allocatable :: flux_unit
    character(:), allocatable :: delta
    integer(i4) :: timestamp
    integer(i4) :: horizon

    if (.not.self%io%output_active) return

    timestamp = self%output_config%output_time_reference
    delta = time_units_delta(self%output_config%output_frequency, timestamp)
    dtype = "f64"
    if (.not.self%output_config%output_double_precision) dtype = "f32"
    flux_unit = mhm_flux_units(self)
    allocate(vars(0))

    if (self%output_config%out_interception) then
      vars = [vars, var(name="interception", long_name="canopy interception storage", units="mm", dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_snowpack) then
      vars = [vars, var(name="snowpack", long_name="depth of snowpack", units="mm", dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_swc) then
      do horizon = 1_i4, size(self%exchange%soil_moisture%data, 2)
        vars = [vars, var(name="SWC_L" // trim(n2s(horizon, '(i2.2)')), &
          long_name="soil water content of soil layer" // trim(n2s(horizon)), units="mm", dtype=dtype, avg=.true.)]
      end do
    end if
    if (self%output_config%out_sm) then
      do horizon = 1_i4, size(self%exchange%soil_moisture%data, 2)
        vars = [vars, var(name="SM_L" // trim(n2s(horizon, '(i2.2)')), &
          long_name="volumetric soil moisture of soil layer" // trim(n2s(horizon)), units="mm mm-1", dtype=dtype, avg=.true.)]
      end do
    end if
    if (self%output_config%out_sm_all) then
      vars = [vars, var(name="SM_Lall", long_name="average soil moisture over all layers", units="mm mm-1", &
        dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_sealedstw) then
      vars = [vars, var(name="sealedSTW", long_name="reservoir of sealed areas (sealedSTW)", units="mm", dtype=dtype, &
        avg=.true.)]
    end if
    if (self%output_config%out_unsatstw) then
      vars = [vars, var(name="unsatSTW", long_name="reservoir of unsaturated zone", units="mm", dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_satstw) then
      vars = [vars, var(name="satSTW", long_name="water level in groundwater reservoir", units="mm", dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_neutrons) then
      vars = [vars, var(name="neutrons", long_name="ground albedo neutrons", units="cph", dtype=dtype, avg=.true.)]
    end if
    if (self%output_config%out_pet) then
      vars = [vars, var(name="PET", long_name="potential Evapotranspiration", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_aet_all) then
      vars = [vars, var(name="aET", long_name="actual Evapotranspiration", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_q) then
      vars = [vars, var(name="Q", long_name="total runoff generated by every cell", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_qd) then
      vars = [vars, var(name="QD", long_name="direct runoff generated by every cell (runoffSeal)", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_qif) then
      vars = [vars, var(name="QIf", long_name="fast interflow generated by every cell (fastRunoff)", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_qis) then
      vars = [vars, var(name="QIs", long_name="slow interflow generated by every cell (slowRunoff)", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_qb) then
      vars = [vars, var(name="QB", long_name="baseflow generated by every cell", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_recharge) then
      vars = [vars, var(name="recharge", long_name="groundwater recharge", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_soil_infil) then
      do horizon = 1_i4, size(self%exchange%infiltration%data, 2)
        vars = [vars, var(name="soil_infil_L" // trim(n2s(horizon, '(i2.2)')), &
          long_name="infiltration flux from soil layer" // trim(n2s(horizon)), units=flux_unit, dtype=dtype)]
      end do
    end if
    if (self%output_config%out_aet_layer) then
      do horizon = 1_i4, size(self%exchange%aet_soil%data, 2)
        vars = [vars, var(name="aET_L" // trim(n2s(horizon, '(i2.2)')), &
          long_name="actual Evapotranspiration from soil layer" // trim(n2s(horizon)), units="mm " // trim(flux_unit), &
          dtype=dtype)]
      end do
    end if
    if (self%output_config%out_preeffect) then
      vars = [vars, var(name="preEffect", long_name="effective precipitation", units=flux_unit, dtype=dtype)]
    end if
    if (self%output_config%out_qsm) then
      vars = [vars, var(name="Qsm", &
        long_name="Average liquid water generated from solid to liquid phase change in the snow", units=flux_unit, dtype=dtype)]
    end if

    log_info(*) "Create mHM output file: ", self%io%output_path
    call self%ds_out%init( &
      path          = self%io%output_path, &
      grid          = self%exchange%level1, &
      vars          = vars, &
      start_time    = self%exchange%start_time, &
      delta         = delta, &
      timestamp     = timestamp, &
      deflate_level = self%output_config%output_deflate_level)
  end subroutine mhm_create_output

  !> \brief Update buffered mHM outputs and write them when the configured stamp is reached.
  subroutine mhm_update_output(self)
    class(mhm_t), intent(inout), target :: self
    real(dp), allocatable :: f_not_sealed(:)
    real(dp), allocatable :: tmp(:)
    logical :: write_stamp
    integer(i4) :: horizon

    if (.not.self%io%output_active) return

    if (self%io%calc_f_not_sealed) then
      allocate(f_not_sealed(size(self%exchange%f_sealed%data)))
      f_not_sealed = 1.0_dp - self%exchange%f_sealed%data
    end if

    if (self%output_config%out_interception) call self%ds_out%update("interception", self%exchange%interception%data)
    if (self%output_config%out_snowpack) call self%ds_out%update("snowpack", self%exchange%snowpack%data)
    if (self%output_config%out_swc) then
      do horizon = 1_i4, size(self%exchange%soil_moisture%data, 2)
        call self%ds_out%update("SWC_L" // trim(n2s(horizon, '(i2.2)')), self%exchange%soil_moisture%data(:, horizon))
      end do
    end if
    if (self%output_config%out_sm) then
      do horizon = 1_i4, size(self%exchange%soil_moisture%data, 2)
        call self%ds_out%update("SM_L" // trim(n2s(horizon, '(i2.2)')), &
          self%exchange%soil_moisture%data(:, horizon) / self%exchange%sm_saturation%data(:, horizon))
      end do
    end if
    if (self%output_config%out_sm_all) then
      allocate(tmp(size(self%exchange%soil_moisture%data, 1)))
      tmp = sum(self%exchange%soil_moisture%data, dim=2) / sum(self%exchange%sm_saturation%data, dim=2)
      call self%ds_out%update("SM_Lall", tmp)
      deallocate(tmp)
    end if
    if (self%output_config%out_sealedstw) call self%ds_out%update("sealedSTW", self%exchange%sealed_storage%data)
    if (self%output_config%out_unsatstw) call self%ds_out%update("unsatSTW", self%exchange%unsat_storage%data)
    if (self%output_config%out_satstw) call self%ds_out%update("satSTW", self%exchange%sat_storage%data)
    if (self%output_config%out_neutrons) call self%ds_out%update("neutrons", self%exchange%neutrons%data)
    if (self%output_config%out_pet) call self%ds_out%update("PET", self%exchange%pet%data)
    if (self%output_config%out_aet_all) then
      allocate(tmp(int(self%exchange%level1%ncells, i4)))
      tmp = 0.0_dp
      if (self%exchange%aet_soil%provided) tmp = tmp + sum(self%exchange%aet_soil%data, dim=2) * f_not_sealed
      if (self%exchange%aet_canopy%provided) tmp = tmp + self%exchange%aet_canopy%data
      if (self%exchange%aet_sealed%provided) tmp = tmp + self%exchange%aet_sealed%data * self%exchange%f_sealed%data
      call self%ds_out%update("aET", tmp)
      deallocate(tmp)
    end if
    if (self%output_config%out_q) call self%ds_out%update("Q", self%exchange%runoff_total%data)
    if (self%output_config%out_qd) call self%ds_out%update("QD", self%exchange%runoff_sealed%data * self%exchange%f_sealed%data)
    if (self%output_config%out_qif) call self%ds_out%update("QIf", self%exchange%interflow_fast%data * f_not_sealed)
    if (self%output_config%out_qis) call self%ds_out%update("QIs", self%exchange%interflow_slow%data * f_not_sealed)
    if (self%output_config%out_qb) call self%ds_out%update("QB", self%exchange%baseflow%data * f_not_sealed)
    if (self%output_config%out_recharge) call self%ds_out%update("recharge", self%exchange%percolation%data * f_not_sealed)
    if (self%output_config%out_soil_infil) then
      do horizon = 1_i4, size(self%exchange%infiltration%data, 2)
        call self%ds_out%update("soil_infil_L" // trim(n2s(horizon, '(i2.2)')), &
          self%exchange%infiltration%data(:, horizon) * f_not_sealed)
      end do
    end if
    if (self%output_config%out_aet_layer) then
      do horizon = 1_i4, size(self%exchange%aet_soil%data, 2)
        call self%ds_out%update("aET_L" // trim(n2s(horizon, '(i2.2)')), self%exchange%aet_soil%data(:, horizon) * f_not_sealed)
      end do
    end if
    if (self%output_config%out_preeffect) call self%ds_out%update("preEffect", self%exchange%pre_eff%data)
    if (self%output_config%out_qsm) call self%ds_out%update("Qsm", self%exchange%melt%data)

    write_stamp = .false.
    select case (self%output_config%output_frequency)
    case (daily)
      if (self%exchange%time%is_new_day()) write_stamp = .true.
    case (monthly)
      if (self%exchange%time%is_new_month()) write_stamp = .true.
    case (yearly)
      if (self%exchange%time%is_new_year()) write_stamp = .true.
    case (no_time)
      if (self%exchange%time == self%exchange%end_time) write_stamp = .true.
    case default
      if (mod(self%exchange%step_count, self%output_config%output_frequency) == 0_i4) write_stamp = .true.
    end select
    if (write_stamp) call self%ds_out%write(self%exchange%time)

    if (allocated(f_not_sealed)) deallocate(f_not_sealed)
  end subroutine mhm_update_output

  !> \brief Update canopy interception and throughfall for the current step.
  recursive subroutine mhm_update_interception(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: interception_case
    integer(i4) :: k

    interception_case = self%exchange%parameters%process_matrix(1, 1)
    select case (interception_case)
    case (-1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%aet_canopy%data)
        self%exchange%interception%data(k) = 0.0_dp
        self%exchange%aet_canopy%data(k) = 0.0_dp
      end do
      !$omp end do
    case (0_i4)
      continue
    case (1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%interception%data)
        call canopy_interc(self%exchange%pet%data(k), self%exchange%max_interception%data(k), self%exchange%pre%data(k), &
          self%exchange%interception%data(k), self%exchange%throughfall%data(k), self%exchange%aet_canopy%data(k))
      end do
      !$omp end do
    case default
      log_fatal(*) "mHM: unsupported interception case ", n2s(interception_case), "."
      error stop 1
    end select
  end subroutine mhm_update_interception

  !> \brief Update snow accumulation, melt, and effective precipitation.
  recursive subroutine mhm_update_snow(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: snow_case
    integer(i4) :: k

    snow_case = self%exchange%parameters%process_matrix(2, 1)
    select case (snow_case)
    case (-1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%melt%data)
        self%exchange%snowpack%data(k) = 0.0_dp
        self%exchange%melt%data(k) = 0.0_dp
        self%exchange%snow%data(k) = 0.0_dp
        self%exchange%degday%data(k) = 0.0_dp
      end do
      !$omp end do
    case (0_i4)
      continue
    case (1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%snowpack%data)
        call snow_accum_melt(self%exchange%degday_inc%data(k), self%exchange%degday_max%data(k) * self%runtime%c2TSTu, &
          self%exchange%degday_dry%data(k) * self%runtime%c2TSTu, self%exchange%pre%data(k), self%exchange%temp%data(k), &
          self%exchange%thresh_temp%data(k), self%exchange%throughfall%data(k), self%exchange%snowpack%data(k), &
          self%exchange%degday%data(k), self%exchange%melt%data(k), self%exchange%pre_eff%data(k), self%exchange%rain%data(k), &
          self%exchange%snow%data(k))
      end do
      !$omp end do
    case default
      log_fatal(*) "mHM: unsupported snow case ", n2s(snow_case), "."
      error stop 1
    end select
  end subroutine mhm_update_snow

  !> \brief Update direct runoff and sealed-area storage independently from pervious soil water.
  recursive subroutine mhm_update_direct_runoff(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: direct_runoff_case
    integer(i4) :: month
    integer(i4) :: k

    direct_runoff_case = self%exchange%parameters%process_matrix(4, 1)
    select case (direct_runoff_case)
    case (0_i4)
      continue
    case (1_i4)
      month = self%exchange%time%month
      if (month < 1_i4 .or. month > int(YearMonths, i4)) then
        log_fatal(*) "mHM: invalid month for evaporation coefficients: ", n2s(month), "."
        error stop 1
      end if
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%sealed_storage%data)
        call soil_moisture_direct_runoff(self%exchange%f_sealed%data(k), self%exchange%thresh_sealed%data(k), &
          self%exchange%pet%data(k), self%forcing%evap_coeff(month), self%exchange%aet_canopy%data(k), &
          self%exchange%pre_eff%data(k), &
          self%exchange%runoff_sealed%data(k), self%exchange%sealed_storage%data(k), self%exchange%aet_sealed%data(k))
      end do
      !$omp end do
    case default
      log_fatal(*) "mHM: unsupported direct runoff case ", n2s(direct_runoff_case), "."
      error stop 1
    end select
  end subroutine mhm_update_direct_runoff

  !> \brief Update pervious-soil infiltration, evapotranspiration, and soil moisture.
  recursive subroutine mhm_update_soil_moisture(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: soil_case
    integer(i4) :: k
    real(dp) :: jarvis_thresh
    real(dp), allocatable :: tmp_infiltration(:)
    real(dp), allocatable :: tmp_soil_moisture(:)
    real(dp), allocatable :: tmp_aet(:)

    soil_case = self%exchange%parameters%process_matrix(3, 1)
    if (soil_case == 0_i4) then
      return
    end if
    if (soil_case < 0_i4 .or. soil_case > 4_i4) then
      log_fatal(*) "mHM: unsupported soil moisture case ", n2s(soil_case), "."
      error stop 1
    end if

    allocate(tmp_infiltration(size(self%exchange%infiltration%data, 2)))
    allocate(tmp_soil_moisture(size(self%exchange%soil_moisture%data, 2)))
    allocate(tmp_aet(size(self%exchange%aet_soil%data, 2)))
    !$omp do schedule(static)
    do k = 1_i4, size(self%exchange%soil_moisture%data, 1)
      tmp_infiltration = self%exchange%infiltration%data(k, :)
      tmp_soil_moisture = self%exchange%soil_moisture%data(k, :)
      jarvis_thresh = 0.0_dp
      if (soil_case == 2_i4 .or. soil_case == 3_i4) jarvis_thresh = self%exchange%thresh_jarvis%data(k)

      call soil_moisture_pervious(soil_case, self%exchange%pet%data(k), self%exchange%sm_saturation%data(k, :), &
        self%exchange%f_roots%data(k, :), self%exchange%sm_field_capacity%data(k, :), self%exchange%wilting_point%data(k, :), &
        self%exchange%sm_exponent%data(k, :), jarvis_thresh, self%exchange%aet_canopy%data(k), self%exchange%pre_eff%data(k), &
        tmp_infiltration, tmp_soil_moisture, tmp_aet)

      self%exchange%infiltration%data(k, :) = tmp_infiltration
      self%exchange%soil_moisture%data(k, :) = tmp_soil_moisture
      self%exchange%aet_soil%data(k, :) = tmp_aet
    end do
    !$omp end do
    deallocate(tmp_infiltration, tmp_soil_moisture, tmp_aet)
  end subroutine mhm_update_soil_moisture

  !> \brief Update the unsaturated-zone runoff, interflow, and percolation state.
  recursive subroutine mhm_update_interflow(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: interflow_case
    integer(i4) :: percolation_case
    integer(i4) :: n_horizons
    integer(i4) :: k

    interflow_case = self%exchange%parameters%process_matrix(6, 1)
    percolation_case = self%exchange%parameters%process_matrix(7, 1)
    if (interflow_case == 0_i4 .and. percolation_case == 0_i4) then
      return
    end if
    if (interflow_case /= 1_i4 .or. percolation_case /= 1_i4) then
      log_fatal(*) "mHM: unsupported interflow/percolation cases ", n2s(interflow_case), "/", n2s(percolation_case), "."
      error stop 1
    end if

    n_horizons = size(self%exchange%infiltration%data, 2)
    !$omp do schedule(static)
    do k = 1_i4, size(self%exchange%unsat_storage%data)
      call runoff_unsat_zone(self%runtime%c2TSTu / self%exchange%k_slowflow%data(k), &
        self%runtime%c2TSTu / self%exchange%k_percolation%data(k), self%runtime%c2TSTu / self%exchange%k_fastflow%data(k), &
        self%exchange%alpha%data(k), self%exchange%f_karst_loss%data(k), self%exchange%infiltration%data(k, n_horizons), &
        self%exchange%thresh_unsat%data(k), self%exchange%sat_storage%data(k), self%exchange%unsat_storage%data(k), &
        self%exchange%interflow_slow%data(k), self%exchange%interflow_fast%data(k), self%exchange%percolation%data(k))
      end do
    !$omp end do
  end subroutine mhm_update_interflow

  !> \brief Update the saturated-zone storage and baseflow.
  recursive subroutine mhm_update_baseflow(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: baseflow_case
    integer(i4) :: k

    baseflow_case = self%exchange%parameters%process_matrix(9, 1)
    select case (baseflow_case)
    case (0_i4)
      continue
    case (1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%sat_storage%data)
        call runoff_sat_zone(self%runtime%c2TSTu / self%exchange%k_baseflow%data(k), self%exchange%sat_storage%data(k), &
          self%exchange%baseflow%data(k))
      end do
      !$omp end do
    case default
      log_fatal(*) "mHM: unsupported baseflow case ", n2s(baseflow_case), "."
      error stop 1
    end select
  end subroutine mhm_update_baseflow

  !> \brief Accumulate total runoff from the generated runoff components.
  recursive subroutine mhm_update_total_runoff(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: k
    real(dp) :: baseflow
    real(dp) :: direct_runoff
    real(dp) :: fast_interflow
    real(dp) :: slow_interflow

    if (.not.self%contract%own_total_runoff) return

    !$omp do schedule(static) private(baseflow, direct_runoff, fast_interflow, slow_interflow)
    do k = 1_i4, size(self%exchange%runoff_total%data)
      direct_runoff = 0.0_dp
      fast_interflow = 0.0_dp
      slow_interflow = 0.0_dp
      baseflow = 0.0_dp
      if (self%exchange%runoff_sealed%provided) direct_runoff = self%exchange%runoff_sealed%data(k)
      if (self%exchange%interflow_fast%provided) fast_interflow = self%exchange%interflow_fast%data(k)
      if (self%exchange%interflow_slow%provided) slow_interflow = self%exchange%interflow_slow%data(k)
      if (self%exchange%baseflow%provided) baseflow = self%exchange%baseflow%data(k)
      call L1_total_runoff(self%exchange%f_sealed%data(k), fast_interflow, slow_interflow, baseflow, direct_runoff, &
        self%exchange%runoff_total%data(k))
    end do
    !$omp end do
  end subroutine mhm_update_total_runoff

  !> \brief Update neutron counts from the current soil-water and surface states.
  recursive subroutine mhm_update_neutrons(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: neutron_case
    integer(i4) :: n_layers
    integer(i4) :: k

    neutron_case = self%exchange%parameters%process_matrix(10, 1)
    if (neutron_case == 0_i4) then
      return
    end if

    n_layers = size(self%exchange%soil_moisture%data, 2) - 1_i4
    if (n_layers < 1_i4) then
      log_fatal(*) "mHM: neutron calculations require at least two soil horizons."
      error stop 1
    end if

    select case (neutron_case)
    case (1_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%neutrons%data)
        call DesiletsN0(self%exchange%soil_moisture%data(k, 1:n_layers), self%soil%horizon_bounds(2:n_layers + 1_i4), &
          self%exchange%bulk_density%data(k, 1:n_layers), self%exchange%lattice_water%data(k, 1:n_layers), &
          self%exchange%desilets_n0%data(k), self%exchange%neutrons%data(k))
      end do
      !$omp end do
    case (2_i4)
      !$omp do schedule(static)
      do k = 1_i4, size(self%exchange%neutrons%data)
        call COSMIC(self%exchange%soil_moisture%data(k, 1:n_layers), self%soil%horizon_bounds(2:n_layers + 1_i4), &
          self%neutrons%integral_afast, self%exchange%interception%data(k), self%exchange%snowpack%data(k), &
          self%exchange%desilets_n0%data(k), &
          self%exchange%bulk_density%data(k, 1:n_layers), self%exchange%lattice_water%data(k, 1:n_layers), &
          self%exchange%cosmic_l3%data(k, 1:n_layers), self%exchange%neutrons%data(k))
      end do
      !$omp end do
    case default
      log_fatal(*) "mHM: unsupported neutrons case ", n2s(neutron_case), "."
      error stop 1
    end select
  end subroutine mhm_update_neutrons

  !> \brief Create an mHM restart file containing active physical state fields.
  subroutine mhm_create_restart(self)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset) :: nc
    type(NcVariable) :: nc_var
    type(NcDimension) :: dims0(0)

    if (.not.allocated(self%io%restart_output_path)) then
      log_fatal(*) "mHM: restart output path is not configured."
      error stop 1
    end if
    if (.not.associated(self%exchange%level1)) then
      log_fatal(*) "mHM: cannot write restart without a connected level1 grid."
      error stop 1
    end if

    log_info(*) "Write mHM restart to file: ", self%io%restart_output_path
    nc = NcDataset(self%io%restart_output_path, "w")
    call self%exchange%level1%to_restart(nc)
    call self%write_restart_data(nc)

    nc_var = nc%setVariable("mhm_meta", "i32", dims0(:0))
    call nc_var%setAttribute("time_stamp", self%exchange%time%str())
    call nc_var%setAttribute("domain", self%exchange%nml_domain_id)
    call nc_var%setAttribute("interception_case", self%exchange%parameters%process_matrix(1, 1))
    call nc_var%setAttribute("snow_case", self%exchange%parameters%process_matrix(2, 1))
    call nc_var%setAttribute("soil_moisture_case", self%exchange%parameters%process_matrix(3, 1))
    call nc_var%setAttribute("direct_runoff_case", self%exchange%parameters%process_matrix(4, 1))
    call nc_var%setAttribute("interflow_case", self%exchange%parameters%process_matrix(6, 1))
    call nc_var%setAttribute("percolation_case", self%exchange%parameters%process_matrix(7, 1))
    call nc_var%setAttribute("baseflow_case", self%exchange%parameters%process_matrix(9, 1))
    call nc%close()
  end subroutine mhm_create_restart

  !> \brief Write all mHM restart physical state fields in unpacked level1 form.
  subroutine mhm_write_restart_data(self, nc)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(inout) :: nc
    type(NcDimension) :: dims_xy(2)
    type(NcDimension) :: soil_dim
    real(dp), allocatable :: soil_bounds(:)
    integer(i4) :: interception_case
    integer(i4) :: snow_case
    integer(i4) :: soil_case
    integer(i4) :: direct_runoff_case
    integer(i4) :: interflow_case
    integer(i4) :: baseflow_case

    if (.not.associated(self%exchange%level1)) then
      log_fatal(*) "mHM: level1 grid not connected while writing restart."
      error stop 1
    end if
    if (.not.allocated(self%soil%horizon_bounds)) then
      log_fatal(*) "mHM: soil horizon bounds not available while writing restart."
      error stop 1
    end if

    if (self%exchange%level1%coordsys == cartesian) then
      dims_xy(1) = nc%getDimension("x")
      dims_xy(2) = nc%getDimension("y")
    else
      dims_xy(1) = nc%getDimension("lon")
      dims_xy(2) = nc%getDimension("lat")
    end if
    allocate(soil_bounds(size(self%soil%horizon_bounds)))
    soil_bounds = self%soil%horizon_bounds
    soil_dim = nc%setCoordinate(trim(soilHorizonsVarName), size(self%soil%horizon_bounds) - 1_i4, &
      bounds=soil_bounds, reference=2_i4)
    deallocate(soil_bounds)

    interception_case = self%exchange%parameters%process_matrix(1, 1)
    snow_case = self%exchange%parameters%process_matrix(2, 1)
    soil_case = self%exchange%parameters%process_matrix(3, 1)
    direct_runoff_case = self%exchange%parameters%process_matrix(4, 1)
    interflow_case = self%exchange%parameters%process_matrix(6, 1)
    baseflow_case = self%exchange%parameters%process_matrix(9, 1)

    if (interception_case == 1_i4) then
      call self%write_restart_field_2d(nc, dims_xy, "L1_Inter", "Interception storage at level 1", &
        self%exchange%interception%data)
    end if
    if (snow_case == 1_i4) then
      call self%write_restart_field_2d(nc, dims_xy, "L1_snowPack", "Snowpack at level 1", self%exchange%snowpack%data)
    end if
    if (direct_runoff_case == 1_i4) then
      call self%write_restart_field_2d(nc, dims_xy, "L1_sealSTW", &
        "Retention storage of impervious areas at level 1", self%exchange%sealed_storage%data)
    end if
    if (soil_case > 0_i4) then
      call self%write_restart_field_3d(nc, dims_xy, soil_dim, "L1_soilMoist", "soil moisture at level 1", &
        self%exchange%soil_moisture%data)
    end if
    if (interflow_case == 1_i4) then
      call self%write_restart_field_2d(nc, dims_xy, "L1_unsatSTW", "upper soil storage at level 1", &
        self%exchange%unsat_storage%data)
    end if
    if (interflow_case == 1_i4 .or. baseflow_case == 1_i4) then
      call self%write_restart_field_2d(nc, dims_xy, "L1_satSTW", "groundwater storage at level 1", &
        self%exchange%sat_storage%data)
    end if
  end subroutine mhm_write_restart_data

  !> \brief Restore mHM physical-process states from a restart file.
  subroutine mhm_read_restart_data(self)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset) :: nc
    type(grid_t) :: restart_grid
    integer(i4) :: interception_case
    integer(i4) :: snow_case
    integer(i4) :: soil_case
    integer(i4) :: direct_runoff_case
    integer(i4) :: interflow_case
    integer(i4) :: baseflow_case

    if (.not.allocated(self%io%restart_input_path)) then
      log_fatal(*) "mHM: restart input path is not configured."
      error stop 1
    end if
    if (.not.associated(self%exchange%level1)) then
      log_fatal(*) "mHM: level1 grid not connected before restart read."
      error stop 1
    end if

    log_info(*) "Read mHM restart from file: ", self%io%restart_input_path
    nc = NcDataset(self%io%restart_input_path, "r")
    call restart_grid%from_restart(nc)
    call self%validate_restart_grid(restart_grid)
    call self%validate_restart_horizon_bounds(nc)
    call self%validate_restart_process_cases(nc)

    interception_case = self%exchange%parameters%process_matrix(1, 1)
    snow_case = self%exchange%parameters%process_matrix(2, 1)
    soil_case = self%exchange%parameters%process_matrix(3, 1)
    direct_runoff_case = self%exchange%parameters%process_matrix(4, 1)
    interflow_case = self%exchange%parameters%process_matrix(6, 1)
    baseflow_case = self%exchange%parameters%process_matrix(9, 1)

    if (interception_case == 1_i4) call self%read_restart_field_2d(nc, "L1_Inter", self%canopy%interception)
    if (snow_case == 1_i4) call self%read_restart_field_2d(nc, "L1_snowPack", self%snow%snowpack)
    if (direct_runoff_case == 1_i4) call self%read_restart_field_2d(nc, "L1_sealSTW", self%direct_runoff%storage)
    if (soil_case > 0_i4) call self%read_restart_field_3d(nc, "L1_soilMoist", self%soil%moisture)
    if (interflow_case == 1_i4) call self%read_restart_field_2d(nc, "L1_unsatSTW", self%runoff%unsat_storage)
    if (interflow_case == 1_i4 .or. baseflow_case == 1_i4) then
      call self%read_restart_field_2d(nc, "L1_satSTW", self%runoff%sat_storage)
    end if

    call nc%close()
  end subroutine mhm_read_restart_data

  !> \brief Validate that the restart level1 grid matches the active level1 grid exactly.
  subroutine mhm_validate_restart_grid(self, restart_grid)
    class(mhm_t), intent(inout), target :: self
    type(grid_t), intent(in), target :: restart_grid

    if (restart_grid%coordsys /= self%exchange%level1%coordsys) then
      log_fatal(*) "mHM restart: restart grid coordinate system does not match current level1 grid."
      error stop 1
    end if
    if (restart_grid%nx /= self%exchange%level1%nx .or. restart_grid%ny /= self%exchange%level1%ny) then
      log_fatal(*) "mHM restart: restart grid dimensions do not match current level1 grid."
      error stop 1
    end if
    if (.not.is_close(restart_grid%cellsize, self%exchange%level1%cellsize) .or. &
      .not.is_close(restart_grid%xllcorner, self%exchange%level1%xllcorner) .or. &
      .not.is_close(restart_grid%yllcorner, self%exchange%level1%yllcorner)) then
      log_fatal(*) "mHM restart: restart grid geometry does not match current level1 grid."
      error stop 1
    end if
    if (restart_grid%y_direction /= self%exchange%level1%y_direction) then
      log_fatal(*) "mHM restart: restart grid y-direction does not match current level1 grid."
      error stop 1
    end if
    if (.not.allocated(restart_grid%mask) .or. .not.allocated(self%exchange%level1%mask)) then
      log_fatal(*) "mHM restart: mask information missing during restart-grid validation."
      error stop 1
    end if
    if (any(restart_grid%mask .neqv. self%exchange%level1%mask)) then
      log_fatal(*) "mHM restart: restart grid mask does not match current level1 grid."
      error stop 1
    end if
    if (self%exchange%level1%has_aux_coords()) then
      if (.not.restart_grid%has_aux_coords()) then
        log_fatal(*) "mHM restart: current level1 grid has auxiliary coordinates, but restart grid does not."
        error stop 1
      end if
      if (any(.not.is_close(restart_grid%lon, self%exchange%level1%lon)) .or. &
        any(.not.is_close(restart_grid%lat, self%exchange%level1%lat))) then
        log_fatal(*) "mHM restart: restart grid auxiliary coordinates do not match current level1 grid."
        error stop 1
      end if
    end if
  end subroutine mhm_validate_restart_grid

  !> \brief Validate the restart soil-horizon coordinate against the active MPR-provided soil axis.
  subroutine mhm_validate_restart_horizon_bounds(self, nc)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(inout) :: nc
    type(NcDimension) :: soil_dim
    type(NcVariable) :: soil_var
    type(NcVariable) :: bounds_var
    integer(i4), allocatable :: var_shape(:)
    real(dp), allocatable :: bounds_in(:, :)
    real(dp), allocatable :: bounds_2d(:, :)
    real(dp), allocatable :: restart_bounds(:)
    character(256) :: bounds_name
    integer(i4) :: n_horizons

    if (.not.allocated(self%soil%horizon_bounds)) then
      log_fatal(*) "mHM restart: soil_horizon_bounds not initialized before restart validation."
      error stop 1
    end if
    if (.not.nc%hasDimension(trim(soilHorizonsVarName))) then
      log_fatal(*) "mHM restart: required soil-horizon dimension missing: ", trim(soilHorizonsVarName)
      error stop 1
    end if
    soil_dim = nc%getDimension(trim(soilHorizonsVarName))
    n_horizons = soil_dim%getLength()
    if (n_horizons /= size(self%soil%horizon_bounds) - 1_i4) then
      log_fatal(*) "mHM restart: soil-horizon count ", n2s(n_horizons), &
        " does not match current mHM horizon count ", n2s(size(self%soil%horizon_bounds) - 1_i4), "."
      error stop 1
    end if
    if (.not.nc%hasVariable(trim(soilHorizonsVarName))) then
      log_fatal(*) "mHM restart: required soil-horizon coordinate missing: ", trim(soilHorizonsVarName)
      error stop 1
    end if
    soil_var = nc%getVariable(trim(soilHorizonsVarName))
    if (.not.soil_var%hasAttribute("bounds")) then
      log_fatal(*) "mHM restart: soil-horizon coordinate is missing bounds metadata."
      error stop 1
    end if
    call soil_var%getAttribute("bounds", bounds_name)
    bounds_name = trim(bounds_name)
    if (len_trim(bounds_name) < 1 .or. .not.nc%hasVariable(bounds_name)) then
      log_fatal(*) "mHM restart: required soil-horizon bounds variable missing: ", trim(bounds_name)
      error stop 1
    end if

    bounds_var = nc%getVariable(bounds_name)
    var_shape = bounds_var%getShape()
    if (size(var_shape) /= 2_i4) then
      log_fatal(*) "mHM restart: soil-horizon bounds variable has rank ", n2s(size(var_shape)), ", expected 2."
      error stop 1
    end if
    call bounds_var%getData(bounds_in)
    if (size(bounds_in, 1) == 2_i4 .and. size(bounds_in, 2) == n_horizons) then
      allocate(bounds_2d(2, n_horizons))
      bounds_2d = bounds_in
    else if (size(bounds_in, 2) == 2_i4 .and. size(bounds_in, 1) == n_horizons) then
      allocate(bounds_2d(2, n_horizons))
      bounds_2d = transpose(bounds_in)
    else
      log_fatal(*) "mHM restart: soil-horizon bounds variable has incompatible shape."
      error stop 1
    end if
    allocate(restart_bounds(n_horizons + 1_i4))
    restart_bounds(1) = bounds_2d(1, 1)
    restart_bounds(2:) = bounds_2d(2, :)
    if (any(.not.is_close(restart_bounds, self%soil%horizon_bounds))) then
      log_fatal(*) "mHM restart: soil-horizon bounds do not match the current MPR-provided soil axis."
      error stop 1
    end if
  end subroutine mhm_validate_restart_horizon_bounds

  !> \brief Validate that restart process cases are compatible with the current mHM configuration.
  subroutine mhm_validate_restart_process_cases(self, nc)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(inout) :: nc
    type(NcVariable) :: meta_var

    if (.not.nc%hasVariable("mhm_meta")) then
      log_fatal(*) "mHM restart: required metadata variable mhm_meta is missing. Regenerate the restart file."
      error stop 1
    end if

    meta_var = nc%getVariable("mhm_meta")
    call mhm_validate_restart_case(meta_var, "interception_case", "interception", &
      self%exchange%parameters%process_matrix(1, 1))
    call mhm_validate_restart_case(meta_var, "snow_case", "snow", self%exchange%parameters%process_matrix(2, 1))
    call mhm_validate_restart_case(meta_var, "soil_moisture_case", "soil_moisture", &
      self%exchange%parameters%process_matrix(3, 1))
    call mhm_validate_restart_case(meta_var, "direct_runoff_case", "direct_runoff", &
      self%exchange%parameters%process_matrix(4, 1))
    call mhm_validate_restart_case(meta_var, "interflow_case", "interflow", self%exchange%parameters%process_matrix(6, 1))
    call mhm_validate_restart_case(meta_var, "baseflow_case", "baseflow", self%exchange%parameters%process_matrix(9, 1))
  end subroutine mhm_validate_restart_process_cases

  !> \brief Validate one restart process-case metadata attribute.
  subroutine mhm_validate_restart_case(meta_var, attr_name, process_name, current_case)
    type(NcVariable), intent(in) :: meta_var
    character(*), intent(in) :: attr_name
    character(*), intent(in) :: process_name
    integer(i4), intent(in) :: current_case
    integer(i4) :: restart_case

    if (.not.meta_var%hasAttribute(trim(attr_name))) then
      log_fatal(*) "mHM restart: required metadata attribute ", trim(attr_name), " is missing. Regenerate the restart file."
      error stop 1
    end if

    call meta_var%getAttribute(trim(attr_name), restart_case)
    if ((restart_case > 0_i4 .or. current_case > 0_i4) .and. restart_case /= current_case) then
      log_fatal(*) "mHM restart: process case mismatch for ", trim(process_name), &
        ". Restart value: ", n2s(restart_case), ", current value: ", n2s(current_case), "."
      error stop 1
    end if
  end subroutine mhm_validate_restart_case

  !> \brief Write a packed L1 scalar field as unpacked 2D restart data.
  subroutine mhm_write_restart_field_2d(self, nc, dims_xy, var_name, long_name, data_packed)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(inout) :: nc
    type(NcDimension), intent(in) :: dims_xy(2)
    character(*), intent(in) :: var_name
    character(*), intent(in) :: long_name
    real(dp), intent(in) :: data_packed(:)
    type(NcVariable) :: nc_var
    real(dp), allocatable :: data_2d(:, :)

    allocate(data_2d(self%exchange%level1%nx, self%exchange%level1%ny))
    call self%exchange%level1%unpack_into(data_packed, data_2d)
    nc_var = nc%setVariable(trim(var_name), "f64", dims_xy)
    call nc_var%setFillValue(nodata_dp)
    call nc_var%setAttribute("missing_value", nodata_dp)
    call nc_var%setAttribute("long_name", trim(long_name))
    if (self%exchange%level1%has_aux_coords()) call nc_var%setAttribute("coordinates", "lon lat")
    call nc_var%setData(data_2d)
    deallocate(data_2d)
  end subroutine mhm_write_restart_field_2d

  !> \brief Write a packed L1 horizon-dependent field as unpacked 3D restart data.
  subroutine mhm_write_restart_field_3d(self, nc, dims_xy, dim3, var_name, long_name, data_packed)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(inout) :: nc
    type(NcDimension), intent(in) :: dims_xy(2)
    type(NcDimension), intent(in) :: dim3
    character(*), intent(in) :: var_name
    character(*), intent(in) :: long_name
    real(dp), intent(in) :: data_packed(:, :)
    type(NcDimension) :: dims(3)
    type(NcVariable) :: nc_var
    real(dp), allocatable :: data_3d(:, :, :)
    integer(i4) :: idx

    dims(1:2) = dims_xy
    dims(3) = dim3
    allocate(data_3d(self%exchange%level1%nx, self%exchange%level1%ny, size(data_packed, 2)))
    do idx = 1_i4, size(data_packed, 2)
      call self%exchange%level1%unpack_into(data_packed(:, idx), data_3d(:, :, idx))
    end do
    nc_var = nc%setVariable(trim(var_name), "f64", dims)
    call nc_var%setFillValue(nodata_dp)
    call nc_var%setAttribute("missing_value", nodata_dp)
    call nc_var%setAttribute("long_name", trim(long_name))
    if (self%exchange%level1%has_aux_coords()) call nc_var%setAttribute("coordinates", "lon lat")
    call nc_var%setData(data_3d)
    deallocate(data_3d)
  end subroutine mhm_write_restart_field_3d

  !> \brief Read a packed L1 scalar field from unpacked 2D restart data.
  subroutine mhm_read_restart_field_2d(self, nc, var_name, data_packed)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(in) :: nc
    character(*), intent(in) :: var_name
    real(dp), intent(out) :: data_packed(:)
    type(NcVariable) :: nc_var
    integer(i4), allocatable :: var_shape(:)
    real(dp), allocatable :: data_2d(:, :)

    if (.not.nc%hasVariable(trim(var_name))) then
      log_fatal(*) "mHM restart: required variable missing: ", trim(var_name)
      error stop 1
    end if
    nc_var = nc%getVariable(trim(var_name))
    var_shape = nc_var%getShape()
    if (size(var_shape) /= 2_i4) then
      log_fatal(*) "mHM restart: variable ", trim(var_name), " has rank ", n2s(size(var_shape)), ", expected 2."
      error stop 1
    end if
    if (any(var_shape /= [self%exchange%level1%nx, self%exchange%level1%ny])) then
      log_fatal(*) "mHM restart: variable ", trim(var_name), " has incompatible x/y shape."
      error stop 1
    end if
    if (size(data_packed) /= self%exchange%level1%ncells) then
      log_fatal(*) "mHM restart: target array for ", trim(var_name), " has unexpected packed level1 size."
      error stop 1
    end if
    call nc_var%getData(data_2d)
    call self%exchange%level1%pack_into(data_2d, data_packed)
  end subroutine mhm_read_restart_field_2d

  !> \brief Read a packed L1 horizon-dependent field from unpacked 3D restart data.
  subroutine mhm_read_restart_field_3d(self, nc, var_name, data_packed)
    class(mhm_t), intent(inout), target :: self
    type(NcDataset), intent(in) :: nc
    character(*), intent(in) :: var_name
    real(dp), intent(out) :: data_packed(:, :)
    type(NcVariable) :: nc_var
    integer(i4), allocatable :: var_shape(:)
    real(dp), allocatable :: data_3d(:, :, :)
    integer(i4) :: idx

    if (.not.nc%hasVariable(trim(var_name))) then
      log_fatal(*) "mHM restart: required variable missing: ", trim(var_name)
      error stop 1
    end if
    nc_var = nc%getVariable(trim(var_name))
    var_shape = nc_var%getShape()
    if (size(var_shape) /= 3_i4) then
      log_fatal(*) "mHM restart: variable ", trim(var_name), " has rank ", n2s(size(var_shape)), ", expected 3."
      error stop 1
    end if
    if (any(var_shape(1:2) /= [self%exchange%level1%nx, self%exchange%level1%ny])) then
      log_fatal(*) "mHM restart: variable ", trim(var_name), " has incompatible x/y shape."
      error stop 1
    end if
    if (var_shape(3) /= size(data_packed, 2)) then
      log_fatal(*) "mHM restart: variable ", trim(var_name), " has incompatible soil-horizon shape."
      error stop 1
    end if
    if (size(data_packed, 1) /= self%exchange%level1%ncells) then
      log_fatal(*) "mHM restart: target array for ", trim(var_name), " has unexpected packed level1 size."
      error stop 1
    end if
    call nc_var%getData(data_3d)
    do idx = 1_i4, var_shape(3)
      call self%exchange%level1%pack_into(data_3d(:, :, idx), data_packed(:, idx))
    end do
  end subroutine mhm_read_restart_field_3d

  !> \brief Finalize the mHM process container after the simulation.
  subroutine mhm_finalize(self)
    class(mhm_t), intent(inout), target :: self
    if (self%io%write_restart) call self%create_restart()
    if (self%io%output_active) then
      call self%ds_out%close()
      log_info(*) "Close mHM output file: ", self%io%output_path
    else
      log_info(*) "No mHM output file will be written"
    end if
    call self%clear_exchange()
    if (allocated(self%canopy%interception)) deallocate(self%canopy%interception)
    if (allocated(self%canopy%throughfall)) deallocate(self%canopy%throughfall)
    if (allocated(self%canopy%aet)) deallocate(self%canopy%aet)
    if (allocated(self%snow%snowpack)) deallocate(self%snow%snowpack)
    if (allocated(self%snow%melt)) deallocate(self%snow%melt)
    if (allocated(self%snow%pre_effect)) deallocate(self%snow%pre_effect)
    if (allocated(self%snow%rain)) deallocate(self%snow%rain)
    if (allocated(self%snow%snow)) deallocate(self%snow%snow)
    if (allocated(self%snow%degday)) deallocate(self%snow%degday)
    if (allocated(self%direct_runoff%storage)) deallocate(self%direct_runoff%storage)
    if (allocated(self%direct_runoff%aet)) deallocate(self%direct_runoff%aet)
    if (allocated(self%direct_runoff%runoff)) deallocate(self%direct_runoff%runoff)
    if (allocated(self%soil%horizon_bounds)) deallocate(self%soil%horizon_bounds)
    if (allocated(self%soil%moisture)) deallocate(self%soil%moisture)
    if (allocated(self%soil%infiltration)) deallocate(self%soil%infiltration)
    if (allocated(self%soil%aet)) deallocate(self%soil%aet)
    if (allocated(self%runoff%unsat_storage)) deallocate(self%runoff%unsat_storage)
    if (allocated(self%runoff%sat_storage)) deallocate(self%runoff%sat_storage)
    if (allocated(self%runoff%percolation)) deallocate(self%runoff%percolation)
    if (allocated(self%runoff%fast_interflow)) deallocate(self%runoff%fast_interflow)
    if (allocated(self%runoff%slow_interflow)) deallocate(self%runoff%slow_interflow)
    if (allocated(self%runoff%baseflow)) deallocate(self%runoff%baseflow)
    if (allocated(self%runoff%total_runoff)) deallocate(self%runoff%total_runoff)
    if (allocated(self%neutrons%counts)) deallocate(self%neutrons%counts)
    if (allocated(self%neutrons%integral_afast)) deallocate(self%neutrons%integral_afast)
    self%runtime%c2TSTu = 0.0_dp
    log_info(*) "Finalize mhm"
  end subroutine mhm_finalize

  !> \brief Build the legacy flux-unit string for mHM gridded output metadata.
  function mhm_flux_units(self) result(units)
    class(mhm_t), intent(in), target :: self
    character(:), allocatable :: units
    integer(i4) :: step_hours
    integer(i4) :: total_steps

    step_hours = max(1_i4, int(self%exchange%step / one_hour(), i4))
    select case (self%output_config%output_frequency)
    case (daily)
      units = "mm d-1"
    case (monthly)
      units = "mm month-1"
    case (yearly)
      units = "mm a-1"
    case (no_time)
      total_steps = max(1_i4, nint((self%exchange%end_time - self%exchange%start_time) / self%exchange%step))
      units = "mm " // trim(n2s(total_steps)) // "h-1"
    case default
      if (step_hours * self%output_config%output_frequency == 1_i4) then
        units = "mm h-1"
      else
        units = "mm " // trim(n2s(step_hours * self%output_config%output_frequency)) // "h-1"
      end if
    end select
  end function mhm_flux_units

  !> \brief Disable an unavailable output flag while keeping the run configuration valid.
  subroutine disable_unavailable(flag, name, available)
    logical, intent(inout) :: flag
    character(*), intent(in) :: name
    logical, intent(in) :: available

    if (flag .and. .not.available) then
      log_warn(*) "mHM output ", trim(name), " requested but the required exchange field is not provided; output disabled."
      flag = .false.
    end if
  end subroutine disable_unavailable

  !> \brief Drop optional output requests that are not backed by provided exchange fields.
  subroutine mhm_filter_output(self)
    class(mhm_t), intent(inout), target :: self
    logical :: aet_available

    if (.not.self%io%output_active) return

    call disable_unavailable(self%output_config%out_interception, "interception", &
      self%exchange%interception%available(owned=self%contract%own_interception))
    call disable_unavailable(self%output_config%out_snowpack, "snowpack", &
      self%exchange%snowpack%available(owned=self%contract%own_snowpack))
    call disable_unavailable(self%output_config%out_swc, "SWC_L", &
      self%exchange%soil_moisture%available(owned=self%contract%own_soil_moisture))
    call disable_unavailable(self%output_config%out_sm, "SM_L", &
      self%exchange%soil_moisture%available(owned=self%contract%own_soil_moisture))
    call disable_unavailable(self%output_config%out_sm_all, "SM_Lall", &
      self%exchange%soil_moisture%available(owned=self%contract%own_soil_moisture))
    call disable_unavailable(self%output_config%out_sealedstw, "sealedSTW", &
      self%exchange%sealed_storage%available(owned=self%contract%own_sealed_storage))
    call disable_unavailable(self%output_config%out_unsatstw, "unsatSTW", &
      self%exchange%unsat_storage%available(owned=self%contract%own_unsat_storage))
    call disable_unavailable(self%output_config%out_satstw, "satSTW", &
      self%exchange%sat_storage%available(owned=self%contract%own_sat_storage))
    call disable_unavailable(self%output_config%out_neutrons, "neutrons", &
      self%exchange%neutrons%available(owned=self%contract%own_neutrons))
    call disable_unavailable(self%output_config%out_pet, "PET", &
      self%exchange%pet%available())

    aet_available = self%exchange%aet_canopy%available(owned=self%contract%own_aet_canopy) .or. &
      self%exchange%aet_sealed%available(owned=self%contract%own_aet_sealed) .or. &
      self%exchange%aet_soil%available(owned=self%contract%own_aet_soil)
    call disable_unavailable(self%output_config%out_aet_all, "aET", aet_available)

    call disable_unavailable(self%output_config%out_q, "Q", &
      self%exchange%runoff_total%available(owned=self%contract%own_total_runoff))
    call disable_unavailable(self%output_config%out_qd, "QD", &
      self%exchange%runoff_sealed%available(owned=self%contract%own_runoff_sealed))
    call disable_unavailable(self%output_config%out_qif, "QIf", &
      self%exchange%interflow_fast%available(owned=self%contract%own_interflow_fast))
    call disable_unavailable(self%output_config%out_qis, "QIs", &
      self%exchange%interflow_slow%available(owned=self%contract%own_interflow_slow))
    call disable_unavailable(self%output_config%out_qb, "QB", &
      self%exchange%baseflow%available(owned=self%contract%own_baseflow))
    call disable_unavailable(self%output_config%out_recharge, "recharge", &
      self%exchange%percolation%available(owned=self%contract%own_percolation))
    call disable_unavailable(self%output_config%out_soil_infil, "soil_infil_L", &
      self%exchange%infiltration%available(owned=self%contract%own_infiltration))
    call disable_unavailable(self%output_config%out_aet_layer, "aET_L", &
      self%exchange%aet_soil%available(owned=self%contract%own_aet_soil))
    call disable_unavailable(self%output_config%out_preeffect, "preEffect", &
      self%exchange%pre_eff%available(owned=self%contract%own_pre_effect))
    call disable_unavailable(self%output_config%out_qsm, "Qsm", &
      self%exchange%melt%available(owned=self%contract%own_melt))

    self%io%calc_f_not_sealed = self%output_config%out_qif .or. self%output_config%out_qis .or. &
      self%output_config%out_qb .or. self%output_config%out_recharge .or. self%output_config%out_soil_infil .or. &
      self%output_config%out_aet_layer .or. &
      (self%output_config%out_aet_all .and. self%exchange%aet_soil%available(owned=self%contract%own_aet_soil))
  end subroutine mhm_filter_output

  !> \brief Allocate, resize, or release a 1D field array.
  subroutine allocate_1d(data, n_cells, required)
    real(dp), allocatable, intent(inout) :: data(:)
    integer(i4), intent(in) :: n_cells
    logical, intent(in) :: required

    if (.not.required) then
      if (allocated(data)) deallocate(data)
      return
    end if
    if (allocated(data)) then
      if (size(data) /= n_cells) deallocate(data)
    end if
    if (.not.allocated(data)) allocate(data(n_cells))
  end subroutine allocate_1d

  !> \brief Allocate, resize, or release a 2D field array.
  subroutine allocate_2d(data, n_cells, n_horizons, required)
    real(dp), allocatable, intent(inout) :: data(:, :)
    integer(i4), intent(in) :: n_cells
    integer(i4), intent(in) :: n_horizons
    logical, intent(in) :: required

    if (.not.required) then
      if (allocated(data)) deallocate(data)
      return
    end if
    if (allocated(data)) then
      if (size(data, 1) /= n_cells .or. size(data, 2) /= n_horizons) deallocate(data)
    end if
    if (.not.allocated(data)) allocate(data(n_cells, n_horizons))
  end subroutine allocate_2d

  !> \brief Copy soil-horizon metadata from the exchange contract into container-owned state.
  subroutine mhm_copy_horizon_bounds(self, bounds)
    class(mhm_t), intent(inout), target :: self
    real(dp), dimension(:), pointer, intent(in) :: bounds
    integer(i4) :: i

    if (size(bounds) < 2_i4) then
      log_fatal(*) "mHM: soil_horizon_bounds must contain at least two entries."
      error stop 1
    end if
    if (allocated(self%soil%horizon_bounds)) then
      if (size(self%soil%horizon_bounds) /= size(bounds)) deallocate(self%soil%horizon_bounds)
    end if
    if (.not.allocated(self%soil%horizon_bounds)) allocate(self%soil%horizon_bounds(size(bounds)))
    self%soil%horizon_bounds = bounds
    do i = 2_i4, size(self%soil%horizon_bounds)
      if (.not.ieee_is_finite(self%soil%horizon_bounds(i)) .or. &
        self%soil%horizon_bounds(i) <= self%soil%horizon_bounds(i - 1_i4)) then
        log_fatal(*) "mHM: soil_horizon_bounds must be finite and strictly increasing."
        error stop 1
      end if
    end do
  end subroutine mhm_copy_horizon_bounds

  !> \brief Publish mHM-owned fields through the exchange contract.
  subroutine mhm_publish_exchange(self)
    class(mhm_t), intent(inout), target :: self
    integer(i4) :: interception_case
    integer(i4) :: snow_case

    interception_case = self%exchange%parameters%process_matrix(1, 1)
    snow_case = self%exchange%parameters%process_matrix(2, 1)

    select case (interception_case)
    case (1_i4)
      call self%exchange%interception%publish_local("mHM", self%canopy%interception)
      call self%exchange%throughfall%publish_local("mHM", self%canopy%throughfall)
      call self%exchange%aet_canopy%publish_local("mHM", self%canopy%aet)
    case (-1_i4)
      call self%exchange%interception%publish_local("mHM", self%canopy%interception)
      call self%exchange%throughfall%publish_alias("mHM", self%exchange%pre)
      call self%exchange%aet_canopy%publish_local("mHM", self%canopy%aet)
    case (0_i4)
    end select

    select case (snow_case)
    case (1_i4)
      call self%exchange%snowpack%publish_local("mHM", self%snow%snowpack)
      call self%exchange%melt%publish_local("mHM", self%snow%melt)
      call self%exchange%pre_eff%publish_local("mHM", self%snow%pre_effect)
      call self%exchange%rain%publish_local("mHM", self%snow%rain)
      call self%exchange%snow%publish_local("mHM", self%snow%snow)
      call self%exchange%degday%publish_local("mHM", self%snow%degday)
    case (-1_i4)
      call self%exchange%snowpack%publish_local("mHM", self%snow%snowpack)
      call self%exchange%melt%publish_local("mHM", self%snow%melt)
      call self%exchange%pre_eff%publish_alias("mHM", self%exchange%throughfall)
      call self%exchange%rain%publish_alias("mHM", self%exchange%throughfall)
      call self%exchange%snow%publish_local("mHM", self%snow%snow)
      call self%exchange%degday%publish_local("mHM", self%snow%degday)
    case (0_i4)
    end select

    if (self%contract%own_sealed_storage) call self%exchange%sealed_storage%publish_local("mHM", self%direct_runoff%storage)
    if (self%contract%own_aet_sealed) call self%exchange%aet_sealed%publish_local("mHM", self%direct_runoff%aet)
    if (self%contract%own_runoff_sealed) call self%exchange%runoff_sealed%publish_local("mHM", self%direct_runoff%runoff)

    if (self%contract%own_soil_moisture) call self%exchange%soil_moisture%publish_local("mHM", self%soil%moisture)
    if (self%contract%own_infiltration) call self%exchange%infiltration%publish_local("mHM", self%soil%infiltration)
    if (self%contract%own_aet_soil) call self%exchange%aet_soil%publish_local("mHM", self%soil%aet)

    if (self%contract%own_unsat_storage) call self%exchange%unsat_storage%publish_local("mHM", self%runoff%unsat_storage)
    if (self%contract%own_sat_storage) call self%exchange%sat_storage%publish_local("mHM", self%runoff%sat_storage)
    if (self%contract%own_percolation) call self%exchange%percolation%publish_local("mHM", self%runoff%percolation)
    if (self%contract%own_interflow_fast) call self%exchange%interflow_fast%publish_local("mHM", self%runoff%fast_interflow)
    if (self%contract%own_interflow_slow) call self%exchange%interflow_slow%publish_local("mHM", self%runoff%slow_interflow)
    if (self%contract%own_baseflow) call self%exchange%baseflow%publish_local("mHM", self%runoff%baseflow)
    if (self%contract%own_total_runoff) call self%exchange%runoff_total%publish_local("mHM", self%runoff%total_runoff)

    if (self%contract%own_neutrons) call self%exchange%neutrons%publish_local("mHM", self%neutrons%counts)
  end subroutine mhm_publish_exchange

  !> \brief Reset all mHM fields to default values.
  subroutine mhm_reset_fields(self)
    class(mhm_t), intent(inout), target :: self

    if (allocated(self%canopy%interception)) self%canopy%interception = P1_InitStateFluxes

    if (allocated(self%snow%snowpack)) self%snow%snowpack = P2_InitStateFluxes

    if (allocated(self%direct_runoff%storage)) self%direct_runoff%storage = P1_InitStateFluxes

    if (allocated(self%soil%moisture)) then
      if (.not.associated(self%exchange%sm_field_capacity%data)) then
        log_fatal(*) "mHM: sm_field_capacity data not connected for soil-moisture field initialization."
        error stop 1
      end if
      self%soil%moisture = 0.5_dp * self%exchange%sm_field_capacity%data
    end if

    if (allocated(self%runoff%unsat_storage)) self%runoff%unsat_storage = P3_InitStateFluxes
    if (allocated(self%runoff%sat_storage)) self%runoff%sat_storage = P4_InitStateFluxes

    if (allocated(self%canopy%throughfall)) self%canopy%throughfall = P1_InitStateFluxes
    if (allocated(self%canopy%aet)) self%canopy%aet = P1_InitStateFluxes

    if (allocated(self%snow%melt)) self%snow%melt = P1_InitStateFluxes
    if (allocated(self%snow%pre_effect)) self%snow%pre_effect = P1_InitStateFluxes
    if (allocated(self%snow%rain)) self%snow%rain = P1_InitStateFluxes
    if (allocated(self%snow%snow)) self%snow%snow = P1_InitStateFluxes
    if (allocated(self%snow%degday)) self%snow%degday = P1_InitStateFluxes

    if (allocated(self%direct_runoff%aet)) self%direct_runoff%aet = P1_InitStateFluxes
    if (allocated(self%direct_runoff%runoff)) self%direct_runoff%runoff = P1_InitStateFluxes

    if (allocated(self%soil%infiltration)) self%soil%infiltration = P1_InitStateFluxes
    if (allocated(self%soil%aet)) self%soil%aet = P1_InitStateFluxes

    if (allocated(self%runoff%percolation)) self%runoff%percolation = P1_InitStateFluxes
    if (allocated(self%runoff%fast_interflow)) self%runoff%fast_interflow = P1_InitStateFluxes
    if (allocated(self%runoff%slow_interflow)) self%runoff%slow_interflow = P1_InitStateFluxes
    if (allocated(self%runoff%baseflow)) self%runoff%baseflow = P1_InitStateFluxes
    if (allocated(self%runoff%total_runoff)) self%runoff%total_runoff = P1_InitStateFluxes

    if (allocated(self%neutrons%counts)) self%neutrons%counts = P1_InitStateFluxes
  end subroutine mhm_reset_fields

  !> \brief Clear mHM-owned exchange publications.
  subroutine mhm_clear_exchange(self)
    class(mhm_t), intent(inout), target :: self

    call self%exchange%interception%clear(self%contract%own_interception)
    call self%exchange%throughfall%clear(self%contract%own_throughfall)
    call self%exchange%aet_canopy%clear(self%contract%own_aet_canopy)
    call self%exchange%snowpack%clear(self%contract%own_snowpack)
    call self%exchange%melt%clear(self%contract%own_melt)
    call self%exchange%pre_eff%clear(self%contract%own_pre_effect)
    call self%exchange%rain%clear(self%contract%own_rain)
    call self%exchange%snow%clear(self%contract%own_snow)
    call self%exchange%degday%clear(self%contract%own_degday)
    call self%exchange%sealed_storage%clear(self%contract%own_sealed_storage)
    call self%exchange%aet_sealed%clear(self%contract%own_aet_sealed)
    call self%exchange%runoff_sealed%clear(self%contract%own_runoff_sealed)
    call self%exchange%soil_moisture%clear(self%contract%own_soil_moisture)
    call self%exchange%infiltration%clear(self%contract%own_infiltration)
    call self%exchange%aet_soil%clear(self%contract%own_aet_soil)
    call self%exchange%unsat_storage%clear(self%contract%own_unsat_storage)
    call self%exchange%sat_storage%clear(self%contract%own_sat_storage)
    call self%exchange%percolation%clear(self%contract%own_percolation)
    call self%exchange%interflow_fast%clear(self%contract%own_interflow_fast)
    call self%exchange%interflow_slow%clear(self%contract%own_interflow_slow)
    call self%exchange%baseflow%clear(self%contract%own_baseflow)
    call self%exchange%runoff_total%clear(self%contract%own_total_runoff)
    call self%exchange%neutrons%clear(self%contract%own_neutrons)
  end subroutine mhm_clear_exchange

end module mo_mhm_container
