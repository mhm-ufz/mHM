!> \file    mo_meteo_container.f90
!> \copydoc mo_meteo_container

!> \brief   Container-side meteorological processing for the exchange runtime.
!> \details This module consumes raw meteorological forcings published on the
!! exchange, remaps them from level2 to level1, applies temporal
!! disaggregation, and computes PET corrections or PET estimates for the active
!! PET process case.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Mar 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_meteo_container
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_constants, only: T0_dp
  use mo_datetime, only: datetime
  use mo_exchange_type, only: exchange_t, var_dp
  use mo_grid, only: grid_t, spherical
  use mo_grid_io, only: no_time, daily, monthly, yearly
  use mo_grid_scaler, only: scaler_t
  use mo_kind, only: i4, i8, dp
  use mo_meteo_temporal_tools, only: temporal_disagg_meteo_weights, temporal_disagg_flux_daynight, temporal_disagg_state_daynight
  use mo_mhm_constants, only: HarSamConst
  use mo_pet, only: pet_hargreaves, pet_penman, pet_priestly
  use mo_read_nc, only: read_weights_nc
  use mo_string_utils, only: n2s => num2str
  use mo_utils, only: is_close
  use nml_config_meteo, only: nml_config_meteo_t, NML_OK

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

  !> \class   meteo_t
  !> \brief   Class for a single meteorology process container.
  type, public :: meteo_t
    type(nml_config_meteo_t) :: config !< configuration of the meteorology process container
    type(exchange_t), pointer :: exchange => null() !< exchange container of the domain
    type(grid_t) :: tgt_level1 !< internal level1 grid derived from level0 when needed
    type(scaler_t) :: regrid !< level2-to-level1 remapper for packed fields
    type(meteo_weight_state_t) :: weights !< cached disaggregation weights
    type(meteo_output_state_t) :: out !< processed meteo outputs
    type(meteo_scratch_state_t) :: scratch !< reusable remapped raw forcings
    logical :: active = .false. !< whether meteorological processing participates in the configured domain
  contains
    procedure :: set_dims => meteo_set_dims
    procedure :: configure => meteo_configure
    procedure :: connect => meteo_connect
    procedure :: initialize => meteo_initialize
    procedure :: update => meteo_update
    procedure :: finalize => meteo_finalize
    procedure, private :: steps_per_day => meteo_steps_per_day
    procedure, private :: weight_mode_active => meteo_weight_mode_active
    procedure, private :: fraction_domain => meteo_fraction_domain
    procedure, private :: ensure_level1_grid => meteo_ensure_level1_grid
    procedure, private :: ensure_size => meteo_ensure_size
    procedure, private :: remap_raw => meteo_remap_raw
    procedure, private :: load_weight_cache => meteo_load_weight_cache
    procedure, private :: validate_step => meteo_validate_step
    procedure, private :: require_fraction => meteo_require_fraction
    procedure, private :: load_level1_latitude => meteo_load_level1_latitude
    procedure, private :: update_pre => meteo_update_pre
    procedure, private :: update_temp => meteo_update_temp
    procedure, private :: update_pet => meteo_update_pet
    procedure, private :: update_ssrd => meteo_update_ssrd
    procedure, private :: update_strd => meteo_update_strd
    procedure, private :: update_tann => meteo_update_tann
  end type meteo_t

contains

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
      self%exchange%config%processes%neutrons, self%exchange%config%processes%temperature_routing] /= 0_i4)
    if (.not.self%active) return
    if (.not.meteo_supports_model_step(self%exchange%step_hours)) then
      log_fatal(*) "Meteo supports only 1-hour and 24-hour model steps; temporal aggregation/disaggregation for intermediate steps is not implemented."
      error stop 1
    end if
    if (present(file)) then
      path = self%exchange%get_path(file)
      log_info(*) "Read meteo config: ", path
      status = self%config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
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
  end subroutine meteo_configure

  !> \brief Connect the meteorology process container with other components.
  subroutine meteo_connect(self)
    class(meteo_t), intent(inout), target :: self
    integer(i4) :: domain_id
    integer(i4) :: pet_process
    integer(i4) :: snow_process
    integer(i4) :: riv_temp_process
    integer(i4) :: steps_day
    integer(i4) :: step_hours
    integer(i4) :: frac_domain_id
    integer(i4) :: id(1)
    integer :: status
    character(1024) :: errmsg
    character(:), allocatable :: path
    logical :: need_pre
    logical :: need_temp

    log_info(*) "Connect meteo"

    if (.not.associated(self%exchange%level2)) then
      log_fatal(*) "Meteo: level2 grid not connected."
      error stop 1
    end if
    call self%ensure_level1_grid()
    call self%regrid%init(self%exchange%level2, self%exchange%level1)
    domain_id = self%exchange%nml_domain_id
    id(1) = domain_id
    pet_process = self%exchange%config%processes%pet
    snow_process = self%exchange%config%processes%snow
    riv_temp_process = self%exchange%config%processes%temperature_routing
    steps_day = self%steps_per_day()
    step_hours = self%exchange%step_hours
    frac_domain_id = self%fraction_domain()

    need_pre = self%active
    need_temp = (snow_process == 1_i4) .or. any(pet_process == [1_i4, 2_i4, 3_i4])

    call self%exchange%raw_pre%require("Meteo", need_pre, check_data=.false.)
    call self%exchange%raw_temp%require("Meteo", need_temp, check_data=.false.)
    call self%exchange%raw_pet%require("Meteo", any(pet_process == [-2_i4, -1_i4]), check_data=.false.)
    call self%exchange%raw_tann%require("Meteo", riv_temp_process > 0_i4, check_data=.false.)
    call self%exchange%raw_tmin%require("Meteo", pet_process == 1_i4, check_data=.false.)
    call self%exchange%raw_tmax%require("Meteo", pet_process == 1_i4, check_data=.false.)
    call self%exchange%raw_ssrd%require("Meteo", riv_temp_process > 0_i4, check_data=.false.)
    call self%exchange%raw_strd%require("Meteo", riv_temp_process > 0_i4, check_data=.false.)
    call self%exchange%raw_netrad%require("Meteo", any(pet_process == [2_i4, 3_i4]), check_data=.false.)
    call self%exchange%raw_eabs%require("Meteo", pet_process == 3_i4, check_data=.false.)
    call self%exchange%raw_wind%require("Meteo", pet_process == 3_i4, check_data=.false.)

    if (need_pre) then
      call self%validate_step("raw_pre", self%exchange%raw_pre%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_pre%stepping == daily) then
        call self%require_fraction("frac_night_pre", frac_domain_id)
      end if
      call self%ensure_size(self%out%pre, self%exchange%level1%ncells)
      call self%exchange%pre%publish_local("Meteo", self%out%pre, step_hours)
    end if

    if (need_temp) then
      call self%validate_step("raw_temp", self%exchange%raw_temp%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_temp%stepping == daily) then
        call self%require_fraction("frac_night_temp", frac_domain_id)
      end if
      call self%ensure_size(self%out%temp, self%exchange%level1%ncells)
      call self%exchange%temp%publish_local("Meteo", self%out%temp, step_hours)
    end if

    if (pet_process /= 0_i4) then
      call self%ensure_size(self%out%pet, self%exchange%level1%ncells)
      call self%exchange%pet%publish_local("Meteo", self%out%pet, step_hours)
    end if

    if (any(pet_process == [-2_i4, -1_i4])) then
      call self%validate_step("raw_pet", self%exchange%raw_pet%stepping, allow_daily=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4 .and. self%exchange%raw_pet%stepping == daily) then
        call self%require_fraction("frac_night_pet", frac_domain_id)
      end if
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

    if (riv_temp_process > 0_i4) then
      call self%validate_step("raw_ssrd", self%exchange%raw_ssrd%stepping, allow_daily=.true., allow_hourly=.true.)
      call self%validate_step("raw_strd", self%exchange%raw_strd%stepping, allow_daily=.true., allow_hourly=.true.)
      call self%validate_step("raw_tann", self%exchange%raw_tann%stepping, allow_static=.true., allow_daily=.true., &
        allow_monthly=.true., allow_yearly=.true., allow_hourly=.true.)
      if (.not.self%weight_mode_active() .and. steps_day > 1_i4) then
        if (self%exchange%raw_ssrd%stepping == daily) call self%require_fraction("frac_night_ssrd", frac_domain_id)
        if (self%exchange%raw_strd%stepping == daily) call self%require_fraction("frac_night_strd", frac_domain_id)
      end if
      call self%ensure_size(self%out%ssrd, self%exchange%level1%ncells)
      call self%ensure_size(self%out%strd, self%exchange%level1%ncells)
      call self%ensure_size(self%out%tann, self%exchange%level1%ncells)
      call self%exchange%ssrd%publish_local("Meteo", self%out%ssrd, step_hours)
      call self%exchange%strd%publish_local("Meteo", self%out%strd, step_hours)
      call self%exchange%tann%publish_local("Meteo", self%out%tann, step_hours)
    end if

    if (self%weight_mode_active() .and. steps_day > 1_i4) then
      if (need_pre .and. self%exchange%raw_pre%stepping == daily) then
        status = self%config%is_set("pre_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: pre_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%pre_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%pre_weights_var(domain_id)), self%weights%pre)
      end if
      if (need_temp .and. self%exchange%raw_temp%stepping == daily) then
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
      if (riv_temp_process > 0_i4 .and. self%exchange%raw_ssrd%stepping == daily) then
        status = self%config%is_set("ssrd_weights_path", idx=id, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Meteo: ssrd_weights_path not set for domain ", n2s(domain_id), ". Error: ", trim(errmsg)
          error stop 1
        end if
        path = self%exchange%get_path(self%config%ssrd_weights_path(domain_id))
        call self%load_weight_cache(path, trim(self%config%ssrd_weights_var(domain_id)), self%weights%ssrd)
      end if
      if (riv_temp_process > 0_i4 .and. self%exchange%raw_strd%stepping == daily) then
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
    integer(i8) :: n_l1

    log_info(*) "Initialize meteo"

    pet_process = self%exchange%config%processes%pet
    n_l1 = self%exchange%level1%ncells
    select case (pet_process)
    case (1_i4)
      call self%exchange%pet_fac_aspect%require("Meteo", .true., [n_l1])
      call self%exchange%pet_coeff_hs%require("Meteo", .true., [n_l1])
    case (2_i4)
      call self%exchange%pet_coeff_pt%require("Meteo", .true., [n_l1])
    case (3_i4)
      call self%exchange%resist_aero%require("Meteo", .true., [n_l1])
      call self%exchange%resist_surf%require("Meteo", .true., [n_l1])
    case (-2_i4)
      call self%exchange%pet_fac_aspect%require("Meteo", .true., [n_l1])
    case (-1_i4)
      call self%exchange%pet_fac_lai%require("Meteo", .true., [n_l1])
    end select

    if (allocated(self%out%pre)) self%out%pre = 0.0_dp
    if (allocated(self%out%temp)) self%out%temp = 0.0_dp
    if (allocated(self%out%pet)) self%out%pet = 0.0_dp
    if (allocated(self%out%ssrd)) self%out%ssrd = 0.0_dp
    if (allocated(self%out%strd)) self%out%strd = 0.0_dp
    if (allocated(self%out%tann)) self%out%tann = 0.0_dp
  end subroutine meteo_initialize

  !> \brief Update the meteorology process container for the current time step.
  subroutine meteo_update(self)
    class(meteo_t), intent(inout), target :: self

    log_trace(*) "Update meteo"
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

    call self%exchange%pre%clear(owned=.true.)
    call self%exchange%temp%clear(owned=.true.)
    call self%exchange%pet%clear(owned=.true.)
    call self%exchange%ssrd%clear(owned=.true.)
    call self%exchange%strd%clear(owned=.true.)
    call self%exchange%tann%clear(owned=.true.)

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
  end subroutine meteo_finalize

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

    if (.not.associated(self%exchange%level0)) then
      log_fatal(*) "Meteo: level0 grid not connected."
      error stop 1
    end if

    l1_res = self%exchange%level1_resolution
    if (associated(self%exchange%level1)) then
      if (ieee_is_finite(l1_res) .and. l1_res > 0.0_dp .and. &
          .not.is_close(self%exchange%level1%cellsize, l1_res)) then
        log_fatal(*) "Meteo: level1 grid cellsize (", n2s(self%exchange%level1%cellsize), &
          ") conflicts with configured level1_resolution (", n2s(l1_res), ")."
        error stop 1
      end if
      call self%exchange%level1%check_is_filled_by(self%exchange%level0, check_mask=.true.)
      return
    end if

    if (.not.ieee_is_finite(l1_res) .or. l1_res <= 0.0_dp) then
      log_fatal(*) "Meteo: level1 resolution not configured."
      error stop 1
    end if

    call self%exchange%level0%gen_grid(self%tgt_level1, target_resolution=l1_res)
    self%exchange%level1 => self%tgt_level1
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
      deallocate(arr)
      allocate(arr(n_cells))
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
    call self%ensure_size(l1_data, self%exchange%level1%ncells)
    call self%regrid%execute(raw_var%data, l1_data)
  end subroutine meteo_remap_raw

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
    allocate(cache(self%exchange%level1%ncells, n_months, n_hours))
    allocate(packed_l2(self%exchange%level2%ncells))
    allocate(packed_l1(self%exchange%level1%ncells))

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

    call self%ensure_size(self%scratch%latitude, self%exchange%level1%ncells)
    if (self%exchange%level1%coordsys == spherical) then
      y_axis = self%exchange%level1%y_axis()
      do k = 1_i8, self%exchange%level1%ncells
        j = self%exchange%level1%cell_ij(k, 2)
        self%scratch%latitude(k) = y_axis(j)
      end do
      deallocate(y_axis)
    else if (self%exchange%level1%has_aux_coords()) then
      call self%exchange%level1%pack_into(self%exchange%level1%lat, self%scratch%latitude)
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
          self%out%pre = self%scratch%pre
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
        self%out%pre = self%scratch%pre
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
          self%out%temp = self%scratch%temp
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
        self%out%temp = self%scratch%temp
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
        call self%ensure_size(self%scratch%pet, self%exchange%level1%ncells)
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
        call self%ensure_size(self%scratch%pet, self%exchange%level1%ncells)
        self%scratch%pet = pet_priestly(PrieTayParam=self%exchange%pet_coeff_pt%data, &
          Rn=max(self%scratch%netrad, 0.0_dp), tavg=self%scratch%temp)
        pet_stepping = daily
      case (3_i4)
        call self%remap_raw(self%exchange%raw_temp, self%scratch%temp, "raw_temp")
        call self%remap_raw(self%exchange%raw_netrad, self%scratch%netrad, "raw_netrad")
        call self%remap_raw(self%exchange%raw_eabs, self%scratch%eabs, "raw_eabs")
        call self%remap_raw(self%exchange%raw_wind, self%scratch%wind, "raw_wind")
        call self%ensure_size(self%scratch%pet, self%exchange%level1%ncells)
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
          self%out%pet = self%scratch%pet
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
        self%out%pet = self%scratch%pet
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
          self%out%ssrd = self%scratch%ssrd
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
        self%out%ssrd = self%scratch%ssrd
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
          self%out%strd = self%scratch%strd
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
        self%out%strd = self%scratch%strd
    end select
  end subroutine meteo_update_strd

  !> \brief Update annual mean temperature on level1 for the current model step.
  subroutine meteo_update_tann(self)
    class(meteo_t), intent(inout), target :: self
    call self%remap_raw(self%exchange%raw_tann, self%scratch%tann, "raw_tann")
    self%out%tann = self%scratch%tann
  end subroutine meteo_update_tann

end module mo_meteo_container
