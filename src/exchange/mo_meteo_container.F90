!> \file    mo_meteo_container.f90
!> \copydoc mo_meteo_container

!> \brief   Container-side meteorological processing for the exchange runtime.
!> \details This module consumes raw meteorological forcings published on the
!! exchange, remaps them from level2 to level1, applies temporal
!! disaggregation, and computes PET corrections or PET estimates for the active
!! PET process case.
!> \version 0.1
!> \changelog
!! - Matthias Cuntz (2012): temporal disaggregation of daily meteorological forcings.
!! - Rohini Kumar (2013): meteorological forcing preparation and spatial remapping.
!! - Matthias Zink, Christoph Schneider, Matthias Cuntz (2014): PET process formulations.
!! - Stephan Thober (2014-2022): chunked forcing input, prescribed weights, and hourly forcing support.
!! - Sebastian Mueller (2023): object-oriented meteorological handler.
!! - Sebastian Mueller (2026): exchange-side meteorology container rewrite.
!> \authors Matthias Cuntz
!> \authors Rohini Kumar
!> \authors Matthias Zink
!> \authors Christoph Schneider
!> \authors Stephan Thober
!> \authors Sebastian Mueller
!> \date    2012 - 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_meteo_container
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_constants, only: T0_dp
  use mo_datetime, only: datetime
  use mo_exchange_type, only: exchange_t, var_dp, l1
  use mo_grid, only: grid_t, spherical
  use mo_grid_io, only: no_time, daily, monthly, yearly
  use mo_grid_scaler, only: scaler_t
  use mo_lake_forcing, only: lake_support_t, lake_support_build
  use mo_kind, only: i4, i8, dp
  use mo_meteo_temporal_tools, only: temporal_disagg_meteo_weights, temporal_disagg_flux_daynight, temporal_disagg_state_daynight
  use mo_mhm_constants, only: HarSamConst
  use mo_pet, only: pet_hargreaves, pet_penman, pet_priestly
  use mo_read_nc, only: read_weights_nc
  use mo_string_utils, only: n2s => num2str
  use mo_utils, only: is_close
  use nml_config_meteo, only: nml_config_meteo_t, NML_OK, NML_ERR_NML_NOT_FOUND

  character(len=*), parameter :: s = "meteo" !< module scope for logging
  public :: meteo_is_day_step, meteo_supports_model_step, meteo_supports_forcing_step

  !> \class   meteo_weight_state_t
  !> \brief   Cached hourly disaggregation weights on level1.
  type :: meteo_weight_state_t
    real(dp), allocatable :: pre(:, :, :) !< precipitation weights (nCells1, 12, 24)
    real(dp), allocatable :: pet(:, :, :) !< PET weights (nCells1, 12, 24)
    real(dp), allocatable :: temp(:, :, :) !< temperature weights (nCells1, 12, 24)
    real(dp), allocatable :: ssrd(:, :, :) !< short-wave radiation weights (nCells1, 12, 24)
    real(dp), allocatable :: strd(:, :, :) !< long-wave radiation weights (nCells1, 12, 24)
  end type meteo_weight_state_t

  !> \class   meteo_output_state_t
  !> \brief   Current processed meteorological outputs on level1.
  type :: meteo_output_state_t
    real(dp), allocatable :: pre(:) !< current precipitation on level1
    real(dp), allocatable :: temp(:) !< current temperature on level1
    real(dp), allocatable :: pet(:) !< current PET on level1
    real(dp), allocatable :: ssrd(:) !< current short-wave radiation on level1
    real(dp), allocatable :: strd(:) !< current long-wave radiation on level1
    real(dp), allocatable :: tann(:) !< current annual mean temperature on level1
  end type meteo_output_state_t

  !> \class   meteo_scratch_state_t
  !> \brief   Reusable remapped raw forcings and PET work arrays on level1.
  type :: meteo_scratch_state_t
    real(dp), allocatable :: pre(:) !< remapped raw precipitation on level1
    real(dp), allocatable :: temp(:) !< remapped raw temperature on level1
    real(dp), allocatable :: pet(:) !< remapped raw or computed PET on level1
    real(dp), allocatable :: tann(:) !< remapped raw annual mean temperature on level1
    real(dp), allocatable :: tmin(:) !< remapped raw minimum temperature on level1
    real(dp), allocatable :: tmax(:) !< remapped raw maximum temperature on level1
    real(dp), allocatable :: ssrd(:) !< remapped raw short-wave radiation on level1
    real(dp), allocatable :: strd(:) !< remapped raw long-wave radiation on level1
    real(dp), allocatable :: netrad(:) !< remapped raw net radiation on level1
    real(dp), allocatable :: eabs(:) !< remapped raw vapor pressure on level1
    real(dp), allocatable :: wind(:) !< remapped raw wind speed on level1
    real(dp), allocatable :: latitude(:) !< packed level1 latitude for Hargreaves PET
  end type meteo_scratch_state_t

  !> \class   meteo_lake_state_t
  !> \brief   Lake-specific forcing support and current forcing caches.
  type :: meteo_lake_state_t
    type(lake_support_t), allocatable :: forcing(:) !< sparse level2 support per lake
    real(dp), allocatable :: pre(:) !< current precipitation on lakes
    real(dp), allocatable :: pet(:) !< current PET on lakes
    real(dp), allocatable :: pre_weights(:, :, :) !< precipitation weights (lake,12,24)
    real(dp), allocatable :: pet_weights(:, :, :) !< PET weights (lake,12,24)
    logical :: active = .false. !< whether lake forcing is part of the selected process case
  end type meteo_lake_state_t

  !> \class   meteo_t
  !> \brief   Class for a single meteorology process container.
  !> \authors Sebastian Mueller
  type, public :: meteo_t
    type(nml_config_meteo_t) :: config !< configuration of the meteorology process container
    type(exchange_t), pointer :: exchange => null() !< exchange container of the domain
    type(grid_t) :: tgt_level1_land !< internal level1 grid derived from level0 when needed
    type(scaler_t) :: regrid !< level2-to-level1 remapper for packed fields
    type(meteo_weight_state_t) :: weights !< cached disaggregation weights
    type(meteo_output_state_t) :: out !< processed meteo outputs
    type(meteo_scratch_state_t) :: scratch !< reusable remapped raw forcings
    type(meteo_lake_state_t) :: lake !< lake-specific forcing state
    logical :: active = .false. !< whether meteorological processing participates in the configured domain
  contains
    procedure :: set_dims => meteo_set_dims
    procedure :: configure => meteo_configure
    procedure :: prepare_restart => meteo_prepare_restart
    procedure :: connect => meteo_connect
    procedure :: initialize => meteo_initialize
    procedure :: update => meteo_update
    procedure :: finalize => meteo_finalize
    procedure :: destroy => meteo_destroy
    procedure, private :: steps_per_day => meteo_steps_per_day
    procedure, private :: weight_mode_active => meteo_weight_mode_active
    procedure, private :: fraction_domain => meteo_fraction_domain
    procedure, private :: ensure_level1_grid => meteo_ensure_level1_grid
    procedure, private :: ensure_size => meteo_ensure_size
    procedure, private :: remap_raw => meteo_remap_raw
    procedure, private :: load_weight_cache => meteo_load_weight_cache
    procedure, private :: validate_step => meteo_validate_step
    procedure, private :: require_fraction => meteo_require_fraction
    procedure, private :: warn_unused_raw => meteo_warn_unused_raw
    procedure, private :: load_level1_latitude => meteo_load_level1_latitude
    procedure, private :: update_pre => meteo_update_pre
    procedure, private :: update_temp => meteo_update_temp
    procedure, private :: update_pet => meteo_update_pet
    procedure, private :: update_ssrd => meteo_update_ssrd
    procedure, private :: update_strd => meteo_update_strd
    procedure, private :: update_tann => meteo_update_tann
    procedure, private :: setup_lake_forcing => meteo_setup_lake_forcing
    procedure, private :: update_lake_forcing => meteo_update_lake_forcing
    procedure, private :: aggregate_lake => meteo_aggregate_lake
    procedure, private :: load_lake_weight_cache => meteo_load_lake_weight_cache
  end type meteo_t

contains

  !> \brief Meteorology owns no restart definition; retained for the uniform lifecycle.
  subroutine meteo_prepare_restart(self)
    class(meteo_t), intent(inout), target :: self
  end subroutine meteo_prepare_restart

  !> \brief Classify an interval with the literal legacy predicate using its start hour.
  pure logical function meteo_is_day_step(hour) result(is_day)
    integer(i4), intent(in) :: hour !< hour at the start of the represented interval

    is_day = (hour >= 7_i4) .and. (hour < 19_i4)
  end function meteo_is_day_step

  !> \brief Report whether meteo can process the configured global model cadence without resampling.
  pure logical function meteo_supports_model_step(step_hours) result(supported)
    integer(i4), intent(in) :: step_hours !< global model cadence in hours

    supported = step_hours == 1_i4 .or. step_hours == 24_i4
  end function meteo_supports_model_step

  !> \brief Report whether forcing support is compatible without temporal resampling.
  pure logical function meteo_supports_forcing_step(model_step_hours, forcing_stepping) result(supported)
    integer(i4), intent(in) :: model_step_hours !< global model cadence in hours
    integer(i4), intent(in) :: forcing_stepping !< forcing support encoding

    supported = forcing_stepping <= 0_i4 .or. forcing_stepping == model_step_hours
  end function meteo_supports_forcing_step

  !> \brief Set runtime dimensions for generated meteo namelists.
  subroutine meteo_set_dims(self)
    class(meteo_t), intent(inout), target :: self
    character(1024) :: errmsg
    integer :: status

    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Error setting meteo config dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine meteo_set_dims

  !> \brief Configure the meteorology process container.
  subroutine meteo_configure(self, file)
    class(meteo_t), intent(inout), target :: self
    character(*), intent(in), optional :: file !< file containing the namelists
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer :: status

    log_info(*) "Configure meteo"
    self%active = any([ &
      self%exchange%config%processes%interception, self%exchange%config%processes%snow, &
      self%exchange%config%processes%soil_moisture, self%exchange%config%processes%direct_runoff, &
      self%exchange%config%processes%pet, self%exchange%config%processes%interflow, &
      self%exchange%config%processes%percolation, self%exchange%config%processes%baseflow, &
      self%exchange%config%processes%neutrons, self%exchange%config%processes%temperature_routing, &
      self%exchange%config%processes%lake] /= 0_i4)
    if (.not.self%active) return
    if (.not.meteo_supports_model_step(self%exchange%step_hours)) then
      log_fatal(*) "Meteo supports only 1-hour and 24-hour model steps; temporal aggregation/disaggregation for intermediate steps is not implemented."
      error stop 1
    end if
    if (present(file)) then
      path = self%exchange%get_path(file)
      log_info(*) "Read meteo config: ", path
      status = self%config%from_file(file=path, errmsg=errmsg)
      if (status == NML_ERR_NML_NOT_FOUND) then
        status = self%config%set(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error setting default meteo config: ", trim(errmsg)
          error stop 1
        end if
      else if (status /= NML_OK) then
        log_fatal(*) "Error reading meteo config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%is_configured) then
      log_fatal(*) "Meteo config not set."
      error stop 1
    end if
    status = self%config%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Meteo config not valid: ", trim(errmsg)
      error stop 1
    end if
    self%lake%active = self%exchange%config%processes%lake == -2_i4
    if (any([self%exchange%config%processes%interception, self%exchange%config%processes%snow, &
        self%exchange%config%processes%soil_moisture, self%exchange%config%processes%direct_runoff, &
        self%exchange%config%processes%pet, self%exchange%config%processes%interflow, &
        self%exchange%config%processes%percolation, self%exchange%config%processes%baseflow, &
        self%exchange%config%processes%neutrons, self%exchange%config%processes%temperature_routing] /= 0_i4)) then
      call self%exchange%pre%provide("Meteo")
    end if
    if (self%exchange%config%processes%snow == 1_i4 .or. any(self%exchange%config%processes%pet == [1_i4, 2_i4, 3_i4])) call self%exchange%temp%provide("Meteo")
    if (self%exchange%config%processes%pet /= 0_i4) call self%exchange%pet%provide("Meteo")
    if (self%exchange%config%processes%temperature_routing > 0_i4) then
      call self%exchange%ssrd%provide("Meteo")
      call self%exchange%strd%provide("Meteo")
      call self%exchange%tann%provide("Meteo")
    end if
    if (self%lake%active .and. self%exchange%raw_pre%provided) call self%exchange%lake_pre%provide("Meteo")
    if (self%lake%active .and. self%exchange%raw_pet%provided) call self%exchange%lake_pet%provide("Meteo")
  end subroutine meteo_configure

  !> \brief Connect the meteorology process container with other components.
  subroutine meteo_connect(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: pet_process
    integer(i4) :: snow_process
    integer(i4) :: riv_temp_process
    integer(i4) :: lake_process
    integer(i4) :: steps_day
    integer(i4) :: step_hours
    integer(i4) :: frac_domain_id
    integer(i4) :: id(1)
    integer :: status
    character(1024) :: errmsg
    character(:), allocatable :: path
    logical :: need_raw_pre
    logical :: need_raw_temp
    logical :: need_raw_pet
    logical :: need_raw_tmin
    logical :: need_raw_tmax
    logical :: need_raw_netrad
    logical :: need_raw_eabs
    logical :: need_raw_wind
    logical :: need_raw_ssrd
    logical :: need_raw_strd
    logical :: need_raw_tann
    logical :: need_raw_pre_land
    logical :: need_raw_temp_land
    logical :: need_raw_pet_land
    logical :: need_raw_tmin_land
    logical :: need_raw_tmax_land
    logical :: need_raw_netrad_land
    logical :: need_raw_eabs_land
    logical :: need_raw_wind_land
    logical :: need_raw_ssrd_land
    logical :: need_raw_strd_land
    logical :: need_raw_tann_land
    logical :: need_raw_pre_lake
    logical :: need_raw_temp_lake
    logical :: need_raw_pet_lake
    logical :: need_raw_tmin_lake
    logical :: need_raw_tmax_lake
    logical :: need_raw_netrad_lake
    logical :: need_raw_eabs_lake
    logical :: need_raw_wind_lake
    logical :: need_raw_ssrd_lake
    logical :: need_raw_strd_lake
    logical :: need_raw_tann_lake
    logical :: land_active, lake_active

    log_info(*) "Connect meteo"

    lake_active = self%lake%active
    land_active = any([                                   &
      self%exchange%config%processes%interception,        &
      self%exchange%config%processes%snow,                &
      self%exchange%config%processes%soil_moisture,       &
      self%exchange%config%processes%direct_runoff,       &
      self%exchange%config%processes%pet,                 &
      self%exchange%config%processes%interflow,           &
      self%exchange%config%processes%percolation,         &
      self%exchange%config%processes%baseflow,            &
      self%exchange%config%processes%neutrons,            &
      self%exchange%config%processes%temperature_routing] /= 0_i4)

    pet_process = self%exchange%config%processes%pet
    snow_process = self%exchange%config%processes%snow
    riv_temp_process = self%exchange%config%processes%temperature_routing
    lake_process = self%exchange%config%processes%lake

    ! Determine every raw forcing need before connecting any field. Raw PET is
    ! needed on land only by the pass-through PET cases; calculated PET uses
    ! its own temperature and radiation inputs. The lake process case -2
    ! currently consumes precipitation and PET only.
    need_raw_pre_land = land_active
    need_raw_temp_land = snow_process == 1_i4 .or. any(pet_process == [1_i4, 2_i4, 3_i4])
    need_raw_pet_land = any(pet_process == [-2_i4, -1_i4])
    need_raw_tmin_land = pet_process == 1_i4
    need_raw_tmax_land = pet_process == 1_i4
    need_raw_netrad_land = any(pet_process == [2_i4, 3_i4])
    need_raw_eabs_land = pet_process == 3_i4
    need_raw_wind_land = pet_process == 3_i4
    need_raw_ssrd_land = riv_temp_process > 0_i4
    need_raw_strd_land = riv_temp_process > 0_i4
    need_raw_tann_land = riv_temp_process > 0_i4
    need_raw_pre_lake = lake_process == -2_i4
    need_raw_temp_lake = .false.
    need_raw_pet_lake = lake_process == -2_i4
    need_raw_tmin_lake = .false.
    need_raw_tmax_lake = .false.
    need_raw_netrad_lake = .false.
    need_raw_eabs_lake = .false.
    need_raw_wind_lake = .false.
    need_raw_ssrd_lake = .false.
    need_raw_strd_lake = .false.
    need_raw_tann_lake = .false.
    need_raw_pre = need_raw_pre_land .or. need_raw_pre_lake
    need_raw_temp = need_raw_temp_land .or. need_raw_temp_lake
    need_raw_pet = need_raw_pet_land .or. need_raw_pet_lake
    need_raw_tmin = need_raw_tmin_land .or. need_raw_tmin_lake
    need_raw_tmax = need_raw_tmax_land .or. need_raw_tmax_lake
    need_raw_netrad = need_raw_netrad_land .or. need_raw_netrad_lake
    need_raw_eabs = need_raw_eabs_land .or. need_raw_eabs_lake
    need_raw_wind = need_raw_wind_land .or. need_raw_wind_lake
    need_raw_ssrd = need_raw_ssrd_land .or. need_raw_ssrd_lake
    need_raw_strd = need_raw_strd_land .or. need_raw_strd_lake
    need_raw_tann = need_raw_tann_land .or. need_raw_tann_lake

    call self%warn_unused_raw("raw_pre", self%exchange%raw_pre%provided, need_raw_pre)
    call self%warn_unused_raw("raw_temp", self%exchange%raw_temp%provided, need_raw_temp)
    call self%warn_unused_raw("raw_pet", self%exchange%raw_pet%provided, need_raw_pet)
    call self%warn_unused_raw("raw_tmin", self%exchange%raw_tmin%provided, need_raw_tmin)
    call self%warn_unused_raw("raw_tmax", self%exchange%raw_tmax%provided, need_raw_tmax)
    call self%warn_unused_raw("raw_netrad", self%exchange%raw_netrad%provided, need_raw_netrad)
    call self%warn_unused_raw("raw_eabs", self%exchange%raw_eabs%provided, need_raw_eabs)
    call self%warn_unused_raw("raw_wind", self%exchange%raw_wind%provided, need_raw_wind)
    call self%warn_unused_raw("raw_ssrd", self%exchange%raw_ssrd%provided, need_raw_ssrd)
    call self%warn_unused_raw("raw_strd", self%exchange%raw_strd%provided, need_raw_strd)
    call self%warn_unused_raw("raw_tann", self%exchange%raw_tann%provided, need_raw_tann)

    if (.not.land_active .and. .not.lake_active) return

    if (.not.associated(self%exchange%level2)) then
      log_fatal(*) "Meteo: level2 grid not connected."
      error stop 1
    end if

    steps_day = self%steps_per_day()
    frac_domain_id = self%fraction_domain()

    if (need_raw_pre) then
      call self%exchange%raw_pre%check_provided("Meteo")
      call self%validate_step("raw_pre", self%exchange%raw_pre%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_pre%stepping == daily) then
        call self%require_fraction("frac_night_pre", frac_domain_id)
      end if
    end if
    if (need_raw_pet) then
      call self%exchange%raw_pet%check_provided("Meteo")
      call self%validate_step("raw_pet", self%exchange%raw_pet%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_pet%stepping == daily) then
        call self%require_fraction("frac_night_pet", frac_domain_id)
      end if
    end if
    if (need_raw_temp) then
      call self%exchange%raw_temp%check_provided("Meteo")
      call self%validate_step("raw_temp", self%exchange%raw_temp%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_temp%stepping == daily) then
        call self%require_fraction("frac_night_temp", frac_domain_id)
      end if
    end if
    if (need_raw_ssrd) then
      call self%exchange%raw_ssrd%check_provided("Meteo")
      call self%validate_step("raw_ssrd", self%exchange%raw_ssrd%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_ssrd%stepping == daily) then
        call self%require_fraction("frac_night_ssrd", frac_domain_id)
      end if
    end if
    if (need_raw_strd) then
      call self%exchange%raw_strd%check_provided("Meteo")
      call self%validate_step("raw_strd", self%exchange%raw_strd%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_strd%stepping == daily) then
        call self%require_fraction("frac_night_strd", frac_domain_id)
      end if
    end if
    if (need_raw_tann) then
      call self%exchange%raw_tann%check_provided("Meteo")
      call self%validate_step("raw_tann", self%exchange%raw_tann%stepping, allow_static=.true., allow_daily=.true., &
        allow_monthly=.true., allow_yearly=.true., allow_hourly=.true.)
    end if

    if (need_raw_tmin) call self%exchange%raw_tmin%check_provided("Meteo")
    if (need_raw_tmax) call self%exchange%raw_tmax%check_provided("Meteo")
    if (need_raw_netrad) call self%exchange%raw_netrad%check_provided("Meteo")
    if (need_raw_eabs) call self%exchange%raw_eabs%check_provided("Meteo")
    if (need_raw_wind) call self%exchange%raw_wind%check_provided("Meteo")

    if (need_raw_pre_lake .or. need_raw_pet_lake) then
      call self%exchange%lake_ids%check_provided("Meteo")
      call self%exchange%lake_map%check_provided("Meteo")
      call self%setup_lake_forcing()
    end if

    if (.not.land_active) return

    call self%ensure_level1_grid()
    call self%regrid%init(self%exchange%level2, self%exchange%level1_land)
    domain_id = self%exchange%nml_domain_id
    id(1) = domain_id
    step_hours = self%exchange%step_hours

    if (need_raw_pre_land) then
      call self%ensure_size(self%out%pre, self%exchange%level1_land%ncells)
      call self%exchange%pre%prepare_data("Meteo", step_hours)
      self%exchange%pre%data => self%out%pre
    end if

    if (need_raw_temp_land) then
      call self%ensure_size(self%out%temp, self%exchange%level1_land%ncells)
      call self%exchange%temp%prepare_data("Meteo", step_hours)
      self%exchange%temp%data => self%out%temp
    end if

    if (pet_process /= 0_i4) then
      call self%ensure_size(self%out%pet, self%exchange%level1_land%ncells)
      call self%exchange%pet%prepare_data("Meteo", step_hours)
      self%exchange%pet%data => self%out%pet
    end if

    if (pet_process == 1_i4) then
      call self%validate_step("raw_temp", self%exchange%raw_temp%stepping, allow_daily=.true.)
      call self%validate_step("raw_tmin", self%exchange%raw_tmin%stepping, allow_daily=.true.)
      call self%validate_step("raw_tmax", self%exchange%raw_tmax%stepping, allow_daily=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4) then
        call self%require_fraction("frac_night_pet", frac_domain_id)
      end if
      call self%load_level1_latitude()
    else if (pet_process == 2_i4) then
      call self%validate_step("raw_temp", self%exchange%raw_temp%stepping, allow_daily=.true.)
      call self%validate_step("raw_netrad", self%exchange%raw_netrad%stepping, allow_daily=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4) then
        call self%require_fraction("frac_night_pet", frac_domain_id)
      end if
    else if (pet_process == 3_i4) then
      call self%validate_step("raw_temp", self%exchange%raw_temp%stepping, allow_daily=.true.)
      call self%validate_step("raw_netrad", self%exchange%raw_netrad%stepping, allow_daily=.true.)
      call self%validate_step("raw_eabs", self%exchange%raw_eabs%stepping, allow_daily=.true.)
      call self%validate_step("raw_wind", self%exchange%raw_wind%stepping, allow_daily=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4) then
        call self%require_fraction("frac_night_pet", frac_domain_id)
      end if
    end if

    if (need_raw_ssrd_land) then
      call self%ensure_size(self%out%ssrd, self%exchange%level1_land%ncells)
      call self%exchange%ssrd%prepare_data("Meteo", step_hours)
      self%exchange%ssrd%data => self%out%ssrd
    end if
    if (need_raw_strd_land) then
      call self%ensure_size(self%out%strd, self%exchange%level1_land%ncells)
      call self%exchange%strd%prepare_data("Meteo", step_hours)
      self%exchange%strd%data => self%out%strd
    end if
    if (need_raw_tann_land) then
      call self%ensure_size(self%out%tann, self%exchange%level1_land%ncells)
      call self%exchange%tann%prepare_data("Meteo", step_hours)
      self%exchange%tann%data => self%out%tann
    end if

    if (self%weight_mode_active() .and. steps_day > 1_i4) then
      if (need_raw_pre_land .and. self%exchange%raw_pre%stepping == daily) then
        status = self%config%is_set("pre_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: pre_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%pre_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%pre_weights_var(domain_id)), self%weights%pre)
      end if
      if (need_raw_temp_land .and. self%exchange%raw_temp%stepping == daily) then
        status = self%config%is_set("temp_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: temp_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%temp_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%temp_weights_var(domain_id)), self%weights%temp)
      end if
      if (pet_process /= 0_i4 .and. (pet_process > 0_i4 .or. self%exchange%raw_pet%stepping == daily)) then
        status = self%config%is_set("pet_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: pet_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%pet_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%pet_weights_var(domain_id)), self%weights%pet)
      end if
      if (need_raw_ssrd_land .and. self%exchange%raw_ssrd%stepping == daily) then
        status = self%config%is_set("ssrd_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: ssrd_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%ssrd_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%ssrd_weights_var(domain_id)), self%weights%ssrd)
      end if
      if (need_raw_strd_land .and. self%exchange%raw_strd%stepping == daily) then
        status = self%config%is_set("strd_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: strd_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%strd_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%strd_weights_var(domain_id)), self%weights%strd)
      end if
    end if
  end subroutine meteo_connect

  !> \brief Initialize the meteorology process container for the simulation.
  subroutine meteo_initialize(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: pet_process

    log_info(*) "Initialize meteo"

    ! A pass-through lake-only configuration activates meteo for lifecycle consistency,
    ! but owns no meteorological fields and does not establish a level1 grid.
    if (.not.associated(self%exchange%level1_land)) return

    pet_process = self%exchange%config%processes%pet
    select case (pet_process)
    case (1_i4)
      call self%exchange%check_data(self%exchange%pet_fac_aspect, "Meteo")
      call self%exchange%check_data(self%exchange%pet_coeff_hs, "Meteo")
    case (2_i4)
      call self%exchange%check_data(self%exchange%pet_coeff_pt, "Meteo")
    case (3_i4)
      call self%exchange%check_data(self%exchange%resist_aero, "Meteo")
      call self%exchange%check_data(self%exchange%resist_surf, "Meteo")
    case (-2_i4)
      call self%exchange%check_data(self%exchange%pet_fac_aspect, "Meteo")
    case (-1_i4)
      call self%exchange%check_data(self%exchange%pet_fac_lai, "Meteo")
    end select

    if (allocated(self%out%pre)) self%out%pre(:) = 0.0_dp
    if (allocated(self%out%temp)) self%out%temp(:) = 0.0_dp
    if (allocated(self%out%pet)) self%out%pet(:) = 0.0_dp
    if (allocated(self%out%ssrd)) self%out%ssrd(:) = 0.0_dp
    if (allocated(self%out%strd)) self%out%strd(:) = 0.0_dp
    if (allocated(self%out%tann)) self%out%tann(:) = 0.0_dp
  end subroutine meteo_initialize

  !> \brief Update the meteorology process container for the current time step.
  subroutine meteo_update(self)
    class(meteo_t), intent(inout), target :: self

    log_trace(*) "Update meteo"
    if (allocated(self%lake%forcing)) call self%update_lake_forcing()
    if (allocated(self%out%pre)) call self%update_pre()
    if (allocated(self%out%temp)) call self%update_temp()
    if (allocated(self%out%pet)) call self%update_pet()
    if (allocated(self%out%ssrd)) call self%update_ssrd()
    if (allocated(self%out%strd)) call self%update_strd()
    if (allocated(self%out%tann)) call self%update_tann()
  end subroutine meteo_update

  !> \brief Finalize the meteorology process container after the simulation.
  subroutine meteo_finalize(self)
    class(meteo_t), intent(inout), target :: self

    log_info(*) "Finalize meteo"
  end subroutine meteo_finalize

  !> \brief Release meteorological output and remapping caches.
  subroutine meteo_destroy(self)
    class(meteo_t), intent(inout), target :: self

    if (allocated(self%weights%pre)) deallocate(self%weights%pre)
    if (allocated(self%weights%pet)) deallocate(self%weights%pet)
    if (allocated(self%weights%temp)) deallocate(self%weights%temp)
    if (allocated(self%weights%ssrd)) deallocate(self%weights%ssrd)
    if (allocated(self%weights%strd)) deallocate(self%weights%strd)

    if (allocated(self%out%pre)) deallocate(self%out%pre)
    if (allocated(self%out%temp)) deallocate(self%out%temp)
    if (allocated(self%out%pet)) deallocate(self%out%pet)
    if (allocated(self%out%ssrd)) deallocate(self%out%ssrd)
    if (allocated(self%out%strd)) deallocate(self%out%strd)
    if (allocated(self%out%tann)) deallocate(self%out%tann)

    if (allocated(self%scratch%pre)) deallocate(self%scratch%pre)
    if (allocated(self%scratch%temp)) deallocate(self%scratch%temp)
    if (allocated(self%scratch%pet)) deallocate(self%scratch%pet)
    if (allocated(self%scratch%tann)) deallocate(self%scratch%tann)
    if (allocated(self%scratch%tmin)) deallocate(self%scratch%tmin)
    if (allocated(self%scratch%tmax)) deallocate(self%scratch%tmax)
    if (allocated(self%scratch%ssrd)) deallocate(self%scratch%ssrd)
    if (allocated(self%scratch%strd)) deallocate(self%scratch%strd)
    if (allocated(self%scratch%netrad)) deallocate(self%scratch%netrad)
    if (allocated(self%scratch%eabs)) deallocate(self%scratch%eabs)
    if (allocated(self%scratch%wind)) deallocate(self%scratch%wind)
    if (allocated(self%scratch%latitude)) deallocate(self%scratch%latitude)
    if (allocated(self%lake%forcing)) deallocate(self%lake%forcing)
    if (allocated(self%lake%pre)) deallocate(self%lake%pre)
    if (allocated(self%lake%pet)) deallocate(self%lake%pet)
    if (allocated(self%lake%pre_weights)) deallocate(self%lake%pre_weights)
    if (allocated(self%lake%pet_weights)) deallocate(self%lake%pet_weights)
    self%lake%active = .false.
    call self%regrid%destroy()
    call self%tgt_level1_land%destroy()
  end subroutine meteo_destroy

  !> \brief Return the number of model steps per day.
  integer(i4) function meteo_steps_per_day(self) result(steps_day)
    class(meteo_t), intent(in) :: self
    steps_day = 24_i4 / self%exchange%step_hours
    if (steps_day < 1_i4) then
      log_fatal(*) "Meteo: invalid model step size for temporal disaggregation."
      error stop 1
    end if
  end function meteo_steps_per_day

  !> \brief Check whether weight-based temporal disaggregation is active for this domain.
  logical function meteo_weight_mode_active(self) result(active)
    class(meteo_t), intent(in) :: self
    active = self%config%read_meteo_weights(self%exchange%nml_domain_id)
  end function meteo_weight_mode_active

  !> \brief Resolve which domain supplies the day/night fractions.
  integer(i4) function meteo_fraction_domain(self) result(domain_id)
    class(meteo_t), intent(in) :: self
    domain_id = merge(1_i4, self%exchange%nml_domain_id, self%config%share_frac)
  end function meteo_fraction_domain

  !> \brief Ensure level1 is available before remapping level2 forcings.
  subroutine meteo_ensure_level1_grid(self)
    class(meteo_t), intent(inout), target :: self
    real(dp) :: l1_res

    if (.not.associated(self%exchange%level0_land)) then
      log_fatal(*) "Meteo: level0 land grid not connected."
      error stop 1
    end if

    l1_res = self%exchange%level1_resolution
    if (associated(self%exchange%level1_land)) then
      if (ieee_is_finite(l1_res) .and. l1_res > 0.0_dp .and. &
          .not.is_close(self%exchange%level1_land%cellsize, l1_res)) then
        log_fatal(*) "Meteo: level1 grid cellsize (", n2s(self%exchange%level1_land%cellsize), &
          ") conflicts with configured level1_resolution (", n2s(l1_res), ")."
        error stop 1
      end if
      call self%exchange%level1_land%check_is_filled_by(self%exchange%level0_land, check_mask=.true.)
      call self%exchange%alias_full_grid_no_lakes(l1)
      return
    end if

    if (.not.ieee_is_finite(l1_res) .or. l1_res <= 0.0_dp) then
      log_fatal(*) "Meteo: level1 resolution not configured."
      error stop 1
    end if

    call self%exchange%level0_land%gen_grid(self%tgt_level1_land, target_resolution=l1_res)
    self%exchange%level1_land => self%tgt_level1_land
    call self%exchange%alias_full_grid_no_lakes(l1)
    log_info(*) "Meteo: derive level1 grid from level0 with resolution ", n2s(l1_res)
  end subroutine meteo_ensure_level1_grid

  !> \brief Ensure a packed level1 work array has the expected size.
  subroutine meteo_ensure_size(self, arr, n_cells)
    class(meteo_t), intent(inout) :: self
    real(dp), allocatable, intent(inout) :: arr(:)
    integer(i8), intent(in) :: n_cells

    if (.not.allocated(arr)) then
      allocate(arr(n_cells))
    else if (size(arr, kind=i8) /= n_cells) then
      log_fatal(*) "Meteo: work array size changed after allocation."
      error stop 1
    end if
  end subroutine meteo_ensure_size

  !> \brief Remap the current raw packed level2 field to packed level1.
  subroutine meteo_remap_raw(self, raw_var, l1_data, name)
    class(meteo_t), intent(inout), target :: self
    type(var_dp), intent(in), target :: raw_var
    real(dp), allocatable, intent(inout) :: l1_data(:)
    character(*), intent(in) :: name

    if (.not.associated(raw_var%data)) then
      log_fatal(*) "Meteo: raw field not connected for ", trim(name), "."
      error stop 1
    end if
    call self%ensure_size(l1_data, self%exchange%level1_land%ncells)
    call self%regrid%execute(raw_var%data, l1_data)
  end subroutine meteo_remap_raw

  !> \brief Build the static lake-to-level2 support from connected lake topology.
  subroutine meteo_setup_lake_forcing(self)
    class(meteo_t), target, intent(inout) :: self
    integer(i8) :: n_lakes
    integer(i4) :: domain_id
    character(:), allocatable :: path

    if (.not.associated(self%exchange%level0_lake)) then
      log_fatal(*) "Meteo lake forcing requires a connected level-0 lake grid."
      error stop 1
    end if
    n_lakes = size(self%exchange%lake_ids%data, kind=i8)
    call self%exchange%check_data(self%exchange%lake_map, "Meteo")
    call lake_support_build(self%exchange%level2, self%exchange%level0_lake, self%exchange%lake_map%data, &
      self%exchange%lake_ids%data, self%lake%forcing)

    if (self%exchange%lake_pre%provided) then
      allocate(self%lake%pre(n_lakes), source=0.0_dp)
      call self%exchange%lake_pre%prepare_data("Meteo", self%exchange%step_hours)
      self%exchange%lake_pre%data => self%lake%pre
      if (self%weight_mode_active() .and. self%exchange%raw_pre%stepping == daily) then
        domain_id = self%exchange%nml_domain_id
        path = self%exchange%get_path(self%config%pre_weights_path(domain_id))
        call self%load_lake_weight_cache(path, trim(self%config%pre_weights_var(domain_id)), self%lake%pre_weights)
      end if
    end if

    if (self%exchange%lake_pet%provided) then
      allocate(self%lake%pet(n_lakes), source=0.0_dp)
      call self%exchange%lake_pet%prepare_data("Meteo", self%exchange%step_hours)
      self%exchange%lake_pet%data => self%lake%pet
      if (self%weight_mode_active() .and. self%exchange%raw_pet%stepping == daily) then
        domain_id = self%exchange%nml_domain_id
        path = self%exchange%get_path(self%config%pet_weights_path(domain_id))
        call self%load_lake_weight_cache(path, trim(self%config%pet_weights_var(domain_id)), self%lake%pet_weights)
      end if
    end if

  end subroutine meteo_setup_lake_forcing

  !> \brief Area-average a level2 field on every lake support.
  subroutine meteo_aggregate_lake(self, field, result)
    class(meteo_t), intent(in) :: self
    real(dp), intent(in) :: field(:)
    real(dp), intent(out) :: result(:)
    integer(i8) :: lake
    do lake = 1_i8, size(self%lake%forcing, kind=i8)
      result(lake) = self%lake%forcing(lake)%aggregate(field)
    end do
  end subroutine meteo_aggregate_lake

  !> \brief Update current lake precipitation and PET with standard temporal disaggregation.
  subroutine meteo_update_lake_forcing(self)
    class(meteo_t), target, intent(inout) :: self
    integer(i4) :: month, hour, steps_day, domain_id
    integer(i8) :: i
    logical :: isday
    real(dp), allocatable :: daily_values(:)

    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    domain_id = self%fraction_domain()
    isday = meteo_is_day_step(hour)

    if (self%exchange%lake_pre%provided) then
      call self%aggregate_lake(self%exchange%raw_pre%data, self%lake%pre)
      if (self%exchange%raw_pre%stepping == daily .and. steps_day > 1_i4) then
        daily_values = self%lake%pre
        if (self%weight_mode_active()) then
          do i = 1_i8, size(self%lake%pre, kind=i8)
            call temporal_disagg_meteo_weights([daily_values(i)], self%lake%pre_weights(i, month, hour + 1_i4), self%lake%pre(i:i))
          end do
        else
          call temporal_disagg_flux_daynight(isday, real(steps_day, dp), daily_values, &
            1.0_dp-self%config%frac_night_pre(month, domain_id), self%config%frac_night_pre(month, domain_id), self%lake%pre)
        end if
      end if
    end if

    if (self%exchange%lake_pet%provided) then
      call self%aggregate_lake(self%exchange%raw_pet%data, self%lake%pet)
      if (self%exchange%raw_pet%stepping == daily .and. steps_day > 1_i4) then
        daily_values = self%lake%pet
        if (self%weight_mode_active()) then
          do i = 1_i8, size(self%lake%pet, kind=i8)
            call temporal_disagg_meteo_weights([daily_values(i)], self%lake%pet_weights(i, month, hour + 1_i4), self%lake%pet(i:i))
          end do
        else
          call temporal_disagg_flux_daynight(isday, real(steps_day, dp), daily_values, &
            1.0_dp-self%config%frac_night_pet(month, domain_id), self%config%frac_night_pet(month, domain_id), self%lake%pet)
        end if
      end if
    end if

  end subroutine meteo_update_lake_forcing

  !> \brief Read a level2 temporal weight cube and aggregate it onto lake supports.
  subroutine meteo_load_lake_weight_cache(self, path, var_name, cache)
    class(meteo_t), target, intent(inout) :: self
    character(*), intent(in) :: path, var_name
    real(dp), allocatable, intent(inout) :: cache(:, :, :)
    real(dp), allocatable :: l2_data(:, :, :, :), packed(:), tmp(:)
    character(len=256) :: path_fixed
    integer(i4) :: month, hour

    path_fixed = trim(path)
    call read_weights_nc("", self%exchange%level2%nx, self%exchange%level2%ny, trim(var_name), l2_data, &
      self%exchange%level2%mask, fileName=path_fixed)
    if (size(l2_data, 3) /= 12_i4 .or. size(l2_data, 4) /= 24_i4) then
      log_fatal(*) "Meteo: lake temporal weights must have dimensions (12,24)."
      error stop 1
    end if
    allocate(cache(size(self%lake%forcing), 12, 24), tmp(size(self%lake%forcing)), packed(self%exchange%level2%ncells))
    do month = 1_i4, 12_i4
      do hour = 1_i4, 24_i4
        call self%exchange%level2%pack_into(l2_data(:, :, month, hour), packed)
        call self%aggregate_lake(packed, tmp)
        cache(:, month, hour) = tmp
      end do
    end do
    deallocate(l2_data, packed, tmp)
  end subroutine meteo_load_lake_weight_cache

  !> \brief Read and regrid one weight cube to cached packed level1 weights.
  subroutine meteo_load_weight_cache(self, path, var_name, cache)
    class(meteo_t), intent(inout), target :: self
    character(*), intent(in) :: path
    character(*), intent(in) :: var_name
    character(len=256) :: path_fixed
    real(dp), allocatable, intent(inout) :: cache(:, :, :)
    real(dp), allocatable :: l2_data(:, :, :, :)
    real(dp), allocatable :: packed_l2(:)
    real(dp), allocatable :: packed_l1(:)
    integer(i4) :: n_months
    integer(i4) :: n_hours
    integer(i4) :: month
    integer(i4) :: hour

    path_fixed = trim(path)
    call read_weights_nc("", self%exchange%level2%nx, self%exchange%level2%ny, trim(var_name), l2_data, &
      self%exchange%level2%mask, fileName=path_fixed)

    n_months = size(l2_data, 3)
    n_hours = size(l2_data, 4)
    if (n_months /= 12_i4 .or. n_hours /= 24_i4) then
      log_fatal(*) "Meteo: weight cube for ", trim(var_name), " must have shape (12, 24), got (", &
        n2s(n_months), ", ", n2s(n_hours), ")."
      error stop 1
    end if

    if (allocated(cache)) deallocate(cache)
    allocate(cache(self%exchange%level1_land%ncells, n_months, n_hours))
    allocate(packed_l2(self%exchange%level2%ncells))
    allocate(packed_l1(self%exchange%level1_land%ncells))

    do month = 1_i4, n_months
      do hour = 1_i4, n_hours
        call self%exchange%level2%pack_into(l2_data(:, :, month, hour), packed_l2)
        call self%regrid%execute(packed_l2, packed_l1)
        cache(:, month, hour) = packed_l1
      end do
    end do

    deallocate(packed_l2)
    deallocate(packed_l1)
    deallocate(l2_data)
  end subroutine meteo_load_weight_cache

  !> \brief Warn when a supplied raw forcing is not used by the selected Meteo processes.
  subroutine meteo_warn_unused_raw(self, name, provided, needed)
    class(meteo_t), intent(in), target :: self
    character(*), intent(in) :: name
    logical, intent(in) :: provided
    logical, intent(in) :: needed

    if (provided .and. .not.needed) then
      log_warn(*) "Meteo: ", trim(name), " is provided but not used by the selected meteorological processes."
    end if
  end subroutine meteo_warn_unused_raw

  !> \brief Require that monthly day/night fractions are explicitly configured for one domain.
  subroutine meteo_require_fraction(self, name, domain_id)
    class(meteo_t), intent(inout), target :: self
    character(*), intent(in) :: name
    integer(i4), intent(in) :: domain_id
    integer(i4) :: month
    integer(i4) :: idx(2)
    integer :: status
    character(1024) :: errmsg

    do month = 1_i4, 12_i4
      idx = [month, domain_id]
      status = self%config%is_set(name, idx=idx, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Meteo: missing day/night fractions for ", trim(name), &
          " in domain ", n2s(domain_id), ". Error: ", trim(errmsg)
        error stop 1
      end if
    end do
  end subroutine meteo_require_fraction

  !> \brief Validate the allowed stepping contract for one raw forcing.
  subroutine meteo_validate_step(self, name, stepping, allow_static, allow_daily, allow_monthly, allow_yearly, allow_hourly)
    class(meteo_t), intent(in) :: self
    character(*), intent(in) :: name
    integer(i4), intent(in) :: stepping
    logical, intent(in), optional :: allow_static
    logical, intent(in), optional :: allow_daily
    logical, intent(in), optional :: allow_monthly
    logical, intent(in), optional :: allow_yearly
    logical, intent(in), optional :: allow_hourly
    logical :: valid

    valid = .false.
    if (present(allow_static)) valid = valid .or. (allow_static .and. stepping == no_time)
    if (present(allow_daily)) valid = valid .or. (allow_daily .and. stepping == daily)
    if (present(allow_monthly)) valid = valid .or. (allow_monthly .and. stepping == monthly)
    if (present(allow_yearly)) valid = valid .or. (allow_yearly .and. stepping == yearly)
    if (present(allow_hourly)) valid = valid .or. (allow_hourly .and. stepping > 0_i4)

    if (.not.valid) then
      log_fatal(*) "Meteo: unsupported stepping for ", trim(name), ": ", n2s(stepping)
      error stop 1
    end if
    if (.not.meteo_supports_forcing_step(self%exchange%step_hours, stepping)) then
      log_fatal(*) "Meteo: ", trim(name), " has fixed support of ", n2s(stepping), &
        "h, but the model step is ", n2s(self%exchange%step_hours), &
        "h; temporal forcing resampling is not implemented."
      error stop 1
    end if
  end subroutine meteo_validate_step

  !> \brief Pack latitude on level1 for PET formulations that need it.
  subroutine meteo_load_level1_latitude(self)
    class(meteo_t), intent(inout), target :: self
    real(dp), allocatable :: y_axis(:)
    integer(i8) :: k
    integer(i4) :: j

    call self%ensure_size(self%scratch%latitude, self%exchange%level1_land%ncells)
    if (self%exchange%level1_land%coordsys == spherical) then
      y_axis = self%exchange%level1_land%y_axis()
      do k = 1_i8, self%exchange%level1_land%ncells
        j = self%exchange%level1_land%cell_ij(k, 2)
        self%scratch%latitude(k) = y_axis(j)
      end do
      deallocate(y_axis)
    else if (self%exchange%level1_land%has_aux_coords()) then
      call self%exchange%level1_land%pack_into(self%exchange%level1_land%lat, self%scratch%latitude)
    else
      log_fatal(*) "Meteo: PET Hargreaves requires latitude on level1."
      error stop 1
    end if
  end subroutine meteo_load_level1_latitude

  !> \brief Update precipitation on level1 for the current model step.
  subroutine meteo_update_pre(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: month
    integer(i4) :: hour
    integer(i4) :: steps_day
    logical :: isday

    domain_id = self%fraction_domain()
    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    isday = meteo_is_day_step(hour)
    call self%remap_raw(self%exchange%raw_pre, self%scratch%pre, "raw_pre")
    select case (self%exchange%raw_pre%stepping)
      case (daily)
        if (steps_day == 1_i4) then
          self%out%pre(:) = self%scratch%pre
        else if (self%weight_mode_active()) then
          if (.not.allocated(self%weights%pre)) then
            log_fatal(*) "Meteo: precipitation weights not loaded."
            error stop 1
          end if
          call temporal_disagg_meteo_weights(self%scratch%pre, self%weights%pre(:, month, hour + 1_i4), self%out%pre)
        else
          call temporal_disagg_flux_daynight(isday, real(steps_day, dp), self%scratch%pre, &
            1.0_dp - self%config%frac_night_pre(month, domain_id), self%config%frac_night_pre(month, domain_id), self%out%pre)
        end if
      case default
        self%out%pre(:) = self%scratch%pre
    end select
  end subroutine meteo_update_pre

  !> \brief Update temperature on level1 for the current model step.
  subroutine meteo_update_temp(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: month
    integer(i4) :: hour
    integer(i4) :: steps_day
    logical :: isday

    domain_id = self%fraction_domain()
    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    isday = meteo_is_day_step(hour)

    call self%remap_raw(self%exchange%raw_temp, self%scratch%temp, "raw_temp")
    select case (self%exchange%raw_temp%stepping)
      case (daily)
        if (steps_day == 1_i4) then
          self%out%temp(:) = self%scratch%temp
        else if (self%weight_mode_active()) then
          if (.not.allocated(self%weights%temp)) then
            log_fatal(*) "Meteo: temperature weights not loaded."
            error stop 1
          end if
          call temporal_disagg_meteo_weights(self%scratch%temp, self%weights%temp(:, month, hour + 1_i4), self%out%temp, &
            weights_correction=T0_dp)
        else
          call temporal_disagg_state_daynight(isday, real(steps_day, dp), self%scratch%temp, &
            -1.0_dp * self%config%frac_night_temp(month, domain_id), self%config%frac_night_temp(month, domain_id), &
            self%out%temp, add_correction=.true.)
        end if
      case default
        self%out%temp(:) = self%scratch%temp
    end select
  end subroutine meteo_update_temp

  !> \brief Update PET on level1 for the current model step.
  subroutine meteo_update_pet(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: pet_process
    integer(i4) :: domain_id
    integer(i4) :: month
    integer(i4) :: hour
    integer(i4) :: steps_day
    integer(i4) :: pet_stepping
    logical :: isday

    pet_process = self%exchange%config%processes%pet
    domain_id = self%fraction_domain()
    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    isday = meteo_is_day_step(hour)

    select case (pet_process)
      case (-2_i4)
        call self%remap_raw(self%exchange%raw_pet, self%scratch%pet, "raw_pet")
        self%scratch%pet = self%exchange%pet_fac_aspect%data * self%scratch%pet
        pet_stepping = self%exchange%raw_pet%stepping
      case (-1_i4)
        call self%remap_raw(self%exchange%raw_pet, self%scratch%pet, "raw_pet")
        self%scratch%pet = self%exchange%pet_fac_lai%data * self%scratch%pet
        pet_stepping = self%exchange%raw_pet%stepping
      case (1_i4)
        call self%remap_raw(self%exchange%raw_temp, self%scratch%temp, "raw_temp")
        call self%remap_raw(self%exchange%raw_tmin, self%scratch%tmin, "raw_tmin")
        call self%remap_raw(self%exchange%raw_tmax, self%scratch%tmax, "raw_tmax")
        if (any(self%scratch%tmax < self%scratch%tmin)) then
          log_warn(*) "Meteo: tmax smaller than tmin for at least one cell at ", self%exchange%time%str()
        end if
        call self%ensure_size(self%scratch%pet, self%exchange%level1_land%ncells)
        self%scratch%pet = self%exchange%pet_fac_aspect%data * pet_hargreaves( &
          HarSamCoeff=self%exchange%pet_coeff_hs%data, &
          HarSamConst=HarSamConst, &
          tavg=self%scratch%temp, &
          tmax=self%scratch%tmax, &
          tmin=self%scratch%tmin, &
          latitude=self%scratch%latitude, &
          doy=self%exchange%time_step_start%doy())
        pet_stepping = daily
      case (2_i4)
        call self%remap_raw(self%exchange%raw_temp, self%scratch%temp, "raw_temp")
        call self%remap_raw(self%exchange%raw_netrad, self%scratch%netrad, "raw_netrad")
        call self%ensure_size(self%scratch%pet, self%exchange%level1_land%ncells)
        self%scratch%pet = pet_priestly(PrieTayParam=self%exchange%pet_coeff_pt%data, &
          Rn=max(self%scratch%netrad, 0.0_dp), tavg=self%scratch%temp)
        pet_stepping = daily
      case (3_i4)
        call self%remap_raw(self%exchange%raw_temp, self%scratch%temp, "raw_temp")
        call self%remap_raw(self%exchange%raw_netrad, self%scratch%netrad, "raw_netrad")
        call self%remap_raw(self%exchange%raw_eabs, self%scratch%eabs, "raw_eabs")
        call self%remap_raw(self%exchange%raw_wind, self%scratch%wind, "raw_wind")
        call self%ensure_size(self%scratch%pet, self%exchange%level1_land%ncells)
        self%scratch%pet = pet_penman( &
          net_rad=max(self%scratch%netrad, 0.0_dp), &
          tavg=self%scratch%temp, &
          act_vap_pressure=self%scratch%eabs / 1000.0_dp, &
          aerodyn_resistance=self%exchange%resist_aero%data / self%scratch%wind, &
          bulksurface_resistance=self%exchange%resist_surf%data, &
          a_s=1.0_dp, &
          a_sh=1.0_dp)
        pet_stepping = daily
      case default
        log_fatal(*) "Meteo: unsupported PET process case ", n2s(pet_process), "."
        error stop 1
    end select

    select case (pet_stepping)
      case (daily)
        if (steps_day == 1_i4) then
          self%out%pet(:) = self%scratch%pet
        else if (self%weight_mode_active()) then
          if (.not.allocated(self%weights%pet)) then
            log_fatal(*) "Meteo: PET weights not loaded."
            error stop 1
          end if
          call temporal_disagg_meteo_weights(self%scratch%pet, self%weights%pet(:, month, hour + 1_i4), self%out%pet)
        else
          call temporal_disagg_flux_daynight(isday, real(steps_day, dp), self%scratch%pet, &
            1.0_dp - self%config%frac_night_pet(month, domain_id), self%config%frac_night_pet(month, domain_id), self%out%pet)
        end if
      case default
        self%out%pet(:) = self%scratch%pet
    end select
  end subroutine meteo_update_pet

  !> \brief Update short-wave radiation on level1 for the current model step.
  subroutine meteo_update_ssrd(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: month
    integer(i4) :: hour
    integer(i4) :: steps_day
    logical :: isday

    domain_id = self%fraction_domain()
    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    isday = meteo_is_day_step(hour)

    call self%remap_raw(self%exchange%raw_ssrd, self%scratch%ssrd, "raw_ssrd")
    select case (self%exchange%raw_ssrd%stepping)
      case (daily)
        if (steps_day == 1_i4) then
          self%out%ssrd(:) = self%scratch%ssrd
        else if (self%weight_mode_active()) then
          if (.not.allocated(self%weights%ssrd)) then
            log_fatal(*) "Meteo: short-wave radiation weights not loaded."
            error stop 1
          end if
          call temporal_disagg_meteo_weights(self%scratch%ssrd, self%weights%ssrd(:, month, hour + 1_i4), self%out%ssrd)
        else
          call temporal_disagg_state_daynight(isday, real(steps_day, dp), self%scratch%ssrd, &
            1.0_dp - self%config%frac_night_ssrd(month, domain_id), self%config%frac_night_ssrd(month, domain_id), self%out%ssrd)
        end if
      case default
        self%out%ssrd(:) = self%scratch%ssrd
    end select
  end subroutine meteo_update_ssrd

  !> \brief Update long-wave radiation on level1 for the current model step.
  subroutine meteo_update_strd(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: month
    integer(i4) :: hour
    integer(i4) :: steps_day
    logical :: isday

    domain_id = self%fraction_domain()
    month = self%exchange%time_step_start%month
    hour = self%exchange%time_step_start%hour
    steps_day = self%steps_per_day()
    isday = meteo_is_day_step(hour)

    call self%remap_raw(self%exchange%raw_strd, self%scratch%strd, "raw_strd")
    select case (self%exchange%raw_strd%stepping)
      case (daily)
        if (steps_day == 1_i4) then
          self%out%strd(:) = self%scratch%strd
        else if (self%weight_mode_active()) then
          if (.not.allocated(self%weights%strd)) then
            log_fatal(*) "Meteo: long-wave radiation weights not loaded."
            error stop 1
          end if
          call temporal_disagg_meteo_weights(self%scratch%strd, self%weights%strd(:, month, hour + 1_i4), self%out%strd)
        else
          call temporal_disagg_state_daynight(isday, real(steps_day, dp), self%scratch%strd, &
            1.0_dp - self%config%frac_night_strd(month, domain_id), self%config%frac_night_strd(month, domain_id), self%out%strd)
        end if
      case default
        self%out%strd(:) = self%scratch%strd
    end select
  end subroutine meteo_update_strd

  !> \brief Update annual mean temperature on level1 for the current model step.
  subroutine meteo_update_tann(self)
    class(meteo_t), intent(inout), target :: self
    call self%remap_raw(self%exchange%raw_tann, self%scratch%tann, "raw_tann")
    self%out%tann(:) = self%scratch%tann
  end subroutine meteo_update_tann

end module mo_meteo_container
