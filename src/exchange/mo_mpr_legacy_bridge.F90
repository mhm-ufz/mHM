!> \file    mo_mpr_legacy_bridge.f90
!> \copydoc mo_mpr_legacy_bridge

!> \brief   Bridge helpers between the exchange MPR container and legacy MPR routines.
!> \details This module isolates the compatibility boundary between the new
!! exchange-side MPR container and the legacy MPR science routines. It reshapes
!! grids and packed fields into the conventions expected by the older code,
!! reuses the legacy PET/soil/runoff/neutron implementations, and keeps the
!! remaining legacy-global interaction scoped to this bridge layer.
!> \changelog
!! - Juliane Mai          Dec 2012
!!   - implemented the soil-database and LUT readers that now sit behind the bridge setup helpers
!! - Stephan Thober, Rohini Kumar
!!   Dec 2012 / Jan 2013
!!   - established the core MPR PET, soil-moisture, runoff, and upscaling helper routines reused here
!! - Luis Samaniego       Sep 2015
!!   - refined PET aspect handling, including southern-hemisphere support
!! - Rohini Kumar         Mar 2016
!!   - added support for multiple soil database options
!! - Matthias Zink        Jun 2017
!!   - moved PET helper functionality into the dedicated MPR PET module
!! - Mehmet Cuneyd Demirel, Simon Stisen
!!   Apr 2017 / Jun 2020 / Feb 2021
!!   - extended PET and soil-moisture helper routines with Jarvis/Feddes and FC/root variants
!! - Maren Kaluza         Dec 2017
!!   - added neutron helper functionality
!! - Robert Schweppe      Jun 2018
!!   - refactored and reformatted the donor MPR helper modules
!! - Sebastian Mueller    Mar 2026
!!   - wrapped the legacy helper routines behind the exchange-side bridge API used by mpr_t
!> \authors Juliane Mai
!> \authors Stephan Thober
!> \authors Rohini Kumar
!> \authors Luis Samaniego
!> \authors Matthias Zink
!> \authors Mehmet Cuneyd Demirel
!> \authors Simon Stisen
!> \authors Maren Kaluza
!> \authors Robert Schweppe
!> \authors Sebastian Mueller
!> \date    2012 - 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_mpr_legacy_bridge
  use mo_logging
  use mo_common_constants, only: nProcesses
  use mo_constants, only: nodata_dp, nodata_i4
  use mo_grid, only: grid_t, cartesian
  use mo_kind, only: i4, i8, dp
  use mo_grid_scaler, only: scaler_t, up_fraction, up_h_mean, up_a_mean
  use mo_mpr_global_variables, only: soilDB, HorizonDepth_mHM, iFlag_soilDB, nSoilHorizons_mHM, nSoilTypes, tillageDepth
  use mo_mpr_neutrons, only: mpr_neutrons
  use mo_mpr_pet, only: pet_correctbyASP, priestley_taylor_alpha, bulksurface_resistance
  use mo_mpr_runoff, only: mpr_runoff
  use mo_mpr_smhorizons, only: mpr_SMhorizons
  use mo_mpr_soilmoist, only: mpr_sm
  use mo_multi_param_reg, only: aerodynamical_resistance
  use mo_soil_database, only: read_soil_LUT, generate_soil_database
  use mo_string_utils, only: n2s => num2str
  use nml_config_mpr, only: nml_config_mpr_t

  implicit none
  private

  real(dp), parameter :: inch_per_hour_to_centimetre_per_day = 60.96_dp

  public :: mpr_bridge_land_cover_fraction
  public :: mpr_bridge_snow_param
  public :: mpr_bridge_pet_lai
  public :: mpr_bridge_pet_aspect
  public :: mpr_bridge_pet_hargreaves
  public :: mpr_bridge_pet_priestley_taylor
  public :: mpr_bridge_pet_penman_monteith
  public :: mpr_bridge_setup_soil_database
  public :: mpr_bridge_soil_moisture
  public :: mpr_bridge_runoff_param
  public :: mpr_bridge_karstic_param
  public :: mpr_bridge_baseflow_param
  public :: mpr_bridge_sealed_threshold

contains

  !> \brief Build packed L0 ids and packed L1 coarse bounds in the legacy x/y index convention.
  subroutine mpr_bridge_build_l0_to_l1_indices(level0, upscaler, id0, upper_bound1, lower_bound1, left_bound1, right_bound1)
    type(grid_t), intent(in), target :: level0
    class(scaler_t), intent(inout), target :: upscaler
    integer(i4), allocatable, intent(out) :: id0(:)
    integer(i4), allocatable, intent(out) :: upper_bound1(:)
    integer(i4), allocatable, intent(out) :: lower_bound1(:)
    integer(i4), allocatable, intent(out) :: left_bound1(:)
    integer(i4), allocatable, intent(out) :: right_bound1(:)
    integer(i4) :: i
    integer(i4) :: x_lb
    integer(i4) :: x_ub
    integer(i4) :: y_lb
    integer(i4) :: y_ub

    allocate(id0(level0%ncells))
    do i = 1_i4, level0%ncells
      id0(i) = i
    end do

    allocate(upper_bound1(upscaler%target_grid%ncells))
    allocate(lower_bound1(upscaler%target_grid%ncells))
    allocate(left_bound1(upscaler%target_grid%ncells))
    allocate(right_bound1(upscaler%target_grid%ncells))
    do i = 1_i4, upscaler%target_grid%ncells
      call upscaler%coarse_bounds(int(i, i8), x_lb, x_ub, y_lb, y_ub)
      ! Reorder the new coarse bounds into the legacy PET routine convention.
      upper_bound1(i) = x_lb
      lower_bound1(i) = x_ub
      left_bound1(i) = y_lb
      right_bound1(i) = y_ub
    end do
  end subroutine mpr_bridge_build_l0_to_l1_indices

  !> \brief Resolve packed level0 latitude values for legacy PET helpers.
  subroutine mpr_bridge_get_level0_latitude(level0, latitude_l0)
    type(grid_t), intent(in), target :: level0
    real(dp), allocatable, intent(out) :: latitude_l0(:)
    real(dp), allocatable :: y_axis(:)
    integer(i4) :: i
    integer(i4) :: ix
    integer(i4) :: iy

    allocate(latitude_l0(level0%ncells))
    if (level0%coordsys == cartesian) then
      if (.not.level0%has_aux_coords()) then
        log_fatal(*) "MPR bridge: PET aspect correction on a cartesian grid requires auxiliary latitude coordinates."
        error stop 1
      end if
      ! PET aspect correction needs physical latitude even on projected input grids.
      do i = 1_i4, level0%ncells
        ix = level0%cell_ij(i, 1)
        iy = level0%cell_ij(i, 2)
        latitude_l0(i) = level0%lat(ix, iy)
      end do
    else
      y_axis = level0%y_axis()
      do i = 1_i4, level0%ncells
        latitude_l0(i) = y_axis(level0%cell_ij(i, 2))
      end do
    end if
  end subroutine mpr_bridge_get_level0_latitude

  !> \brief Clear the legacy-global soil database state before rebuilding it from exchange inputs.
  subroutine mpr_bridge_reset_soil_database()
    if (allocated(HorizonDepth_mHM)) deallocate(HorizonDepth_mHM)
    if (allocated(soilDB%id)) deallocate(soilDB%id)
    if (allocated(soilDB%nHorizons)) deallocate(soilDB%nHorizons)
    if (allocated(soilDB%is_present)) deallocate(soilDB%is_present)
    if (allocated(soilDB%UD)) deallocate(soilDB%UD)
    if (allocated(soilDB%LD)) deallocate(soilDB%LD)
    if (allocated(soilDB%clay)) deallocate(soilDB%clay)
    if (allocated(soilDB%sand)) deallocate(soilDB%sand)
    if (allocated(soilDB%DbM)) deallocate(soilDB%DbM)
    if (allocated(soilDB%depth)) deallocate(soilDB%depth)
    if (allocated(soilDB%RZdepth)) deallocate(soilDB%RZdepth)
    if (allocated(soilDB%Wd)) deallocate(soilDB%Wd)
    if (allocated(soilDB%nTillHorizons)) deallocate(soilDB%nTillHorizons)
    if (allocated(soilDB%thetaS_till)) deallocate(soilDB%thetaS_till)
    if (allocated(soilDB%thetaS)) deallocate(soilDB%thetaS)
    if (allocated(soilDB%Db)) deallocate(soilDB%Db)
    if (allocated(soilDB%thetaFC_till)) deallocate(soilDB%thetaFC_till)
    if (allocated(soilDB%thetaFC)) deallocate(soilDB%thetaFC)
    if (allocated(soilDB%thetaPW_till)) deallocate(soilDB%thetaPW_till)
    if (allocated(soilDB%thetaPW)) deallocate(soilDB%thetaPW)
    if (allocated(soilDB%Ks)) deallocate(soilDB%Ks)
    nSoilTypes = 0_i4
    iFlag_soilDB = 0_i4
    nSoilHorizons_mHM = 0_i4
    tillageDepth = 0.0_dp
  end subroutine mpr_bridge_reset_soil_database

  !> \brief Initialize the legacy soil database from the new MPR configuration and packed L0 soil ids.
  subroutine mpr_bridge_setup_soil_database(config, domain, soil_lut_path, soil_id_l0, horizon_bounds)
    type(nml_config_mpr_t), intent(in) :: config
    integer(i4), intent(in) :: domain
    character(*), intent(in) :: soil_lut_path
    integer(i4), dimension(:, :), intent(in) :: soil_id_l0
    real(dp), allocatable, intent(out), optional :: horizon_bounds(:)
    integer(i4) :: soil_layers
    integer(i4) :: required_depths
    integer(i4) :: i
    integer(i4) :: j
    integer(i4) :: k

    if (domain < 1_i4) then
      log_fatal(*) "MPR bridge: invalid domain for soil database setup: ", n2s(domain), "."
      error stop 1
    end if

    call mpr_bridge_reset_soil_database()

    ! Translate the new config into the legacy-global soil database controls before reading the LUT.
    iFlag_soilDB = config%soil_db_mode(domain)
    nSoilHorizons_mHM = config%n_layers(domain)
    tillageDepth = real(config%tillage_depth(domain), dp)
    if (nSoilHorizons_mHM < 1_i4) then
      log_fatal(*) "MPR bridge: n_layers must be >= 1 for soil database setup."
      error stop 1
    end if

    allocate(HorizonDepth_mHM(nSoilHorizons_mHM))
    HorizonDepth_mHM = 0.0_dp
    select case (iFlag_soilDB)
      case (0_i4)
        required_depths = max(0_i4, nSoilHorizons_mHM - 1_i4)
      case (1_i4)
        required_depths = nSoilHorizons_mHM
      case default
        log_fatal(*) "MPR bridge: unsupported soil_db_mode=", n2s(iFlag_soilDB), "."
        error stop 1
    end select
    do i = 1_i4, required_depths
      HorizonDepth_mHM(i) = real(config%soil_depth(i, domain), dp)
    end do

    call read_soil_LUT(trim(soil_lut_path))
    call generate_soil_database()

    ! Mark only the soil classes that are actually present in the packed L0 input.
    if (allocated(soilDB%is_present)) deallocate(soilDB%is_present)
    allocate(soilDB%is_present(nSoilTypes))
    soilDB%is_present = 0_i4

    soil_layers = 1_i4
    if (iFlag_soilDB == 1_i4) soil_layers = nSoilHorizons_mHM
    if (size(soil_id_l0, 2) < soil_layers) then
      log_fatal(*) "MPR bridge: soil_id input provides ", n2s(size(soil_id_l0, 2)), &
        " layers, expected at least ", n2s(soil_layers), "."
      error stop 1
    end if

    do i = 1_i4, soil_layers
      do k = 1_i4, size(soil_id_l0, 1)
        j = soil_id_l0(k, i)
        if (j == nodata_i4) cycle
        if (j < 1_i4 .or. j > nSoilTypes) then
          log_fatal(*) "MPR bridge: soil ID ", n2s(j), " at packed L0 cell ", n2s(k), &
            " and layer ", n2s(i), " is outside the soil LUT range 1..", n2s(nSoilTypes), "."
          error stop 1
        end if
        soilDB%is_present(j) = 1_i4
      end do
    end do

    if (present(horizon_bounds)) then
      allocate(horizon_bounds(nSoilHorizons_mHM + 1_i4))
      horizon_bounds(1) = 0.0_dp
      horizon_bounds(2:) = HorizonDepth_mHM(:)
    end if
  end subroutine mpr_bridge_setup_soil_database

  !> \brief Upscale land-cover fractions from packed L0 ids to packed L1 fractions.
  subroutine mpr_bridge_land_cover_fraction(upscaler, land_cover_l0, frac_sealed_cityarea, forest_l1, sealed_l1, pervious_l1)
    class(scaler_t), intent(inout), target :: upscaler
    integer(i4), dimension(:), intent(in) :: land_cover_l0
    real(dp), intent(in) :: frac_sealed_cityarea
    real(dp), dimension(:), intent(out) :: forest_l1
    real(dp), dimension(:), intent(out) :: sealed_l1
    real(dp), dimension(:), intent(out) :: pervious_l1

    call upscaler%execute(land_cover_l0, forest_l1, upscaling_operator=up_fraction, class_id=1_i4)
    call upscaler%execute(land_cover_l0, sealed_l1, upscaling_operator=up_fraction, class_id=2_i4)

    sealed_l1 = frac_sealed_cityarea * sealed_l1
    pervious_l1 = 1.0_dp - forest_l1 - sealed_l1
  end subroutine mpr_bridge_land_cover_fraction

  !> \brief Bridge snow parameter generation through the legacy degree-day routine.
  subroutine mpr_bridge_snow_param(param, forest_l1, sealed_l1, pervious_l1, thresh_temp_l1, degday_dry_l1, degday_inc_l1, degday_max_l1)
    real(dp), dimension(:), intent(in) :: param
    real(dp), dimension(:), intent(in) :: forest_l1
    real(dp), dimension(:), intent(in) :: sealed_l1
    real(dp), dimension(:), intent(in) :: pervious_l1
    real(dp), dimension(:), intent(out) :: thresh_temp_l1
    real(dp), dimension(:), intent(out) :: degday_dry_l1
    real(dp), dimension(:), intent(out) :: degday_inc_l1
    real(dp), dimension(:), intent(out) :: degday_max_l1
    real(dp) :: degday_dry_forest
    real(dp) :: degday_dry_sealed
    real(dp) :: degday_dry_pervious
    real(dp) :: degday_max_forest
    real(dp) :: degday_max_sealed
    real(dp) :: degday_max_pervious

    if (size(param) < 8_i4) then
      log_fatal(*) "MPR bridge: snow parameter set must contain 8 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if

    degday_dry_forest = param(2)
    degday_dry_sealed = param(2) + param(4) + param(3)
    degday_dry_pervious = param(2) + param(4)
    degday_max_forest = param(2) + param(6)
    degday_max_sealed = param(2) + param(4) + param(3) + param(7)
    degday_max_pervious = param(2) + param(4) + param(8)

    thresh_temp_l1 = param(1)
    degday_inc_l1 = param(5)
    degday_dry_l1 = degday_dry_forest * forest_l1 + degday_dry_sealed * sealed_l1 + degday_dry_pervious * pervious_l1
    degday_max_l1 = degday_max_forest * forest_l1 + degday_max_sealed * sealed_l1 + degday_max_pervious * pervious_l1
  end subroutine mpr_bridge_snow_param

  !> \brief Bridge PET-LAI correction with legacy-compatible formulas on top of the new scaler.
  subroutine mpr_bridge_pet_lai(param, land_cover_l0, lai_l0, upscaler, pet_fac_lai_l1)
    real(dp), dimension(:), intent(in) :: param
    integer(i4), dimension(:), intent(in) :: land_cover_l0
    real(dp), dimension(:), intent(in) :: lai_l0
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: pet_fac_lai_l1
    real(dp), allocatable :: pet_fac_lai_l0(:)
    integer(i4) :: i
    integer(i4) :: land_cover_id

    if (size(param) < 5_i4) then
      log_fatal(*) "MPR bridge: PET-LAI parameter set must contain 5 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(land_cover_l0) /= size(lai_l0)) then
      log_fatal(*) "MPR bridge: PET-LAI inputs must have matching packed L0 sizes."
      error stop 1
    end if

    allocate(pet_fac_lai_l0(size(land_cover_l0)))
    do i = 1_i4, size(land_cover_l0)
      land_cover_id = land_cover_l0(i)
      select case (land_cover_id)
        case (1_i4)
          pet_fac_lai_l0(i) = param(1) + (param(4) * (1.0_dp - exp(param(5) * lai_l0(i))))
        case (2_i4)
          pet_fac_lai_l0(i) = param(2) + (param(4) * (1.0_dp - exp(param(5) * lai_l0(i))))
        case (3_i4)
          pet_fac_lai_l0(i) = param(3) + (param(4) * (1.0_dp - exp(param(5) * lai_l0(i))))
        case default
          log_fatal(*) "MPR bridge: unsupported land-cover ID ", n2s(land_cover_id), &
            " for PET-LAI correction at packed L0 cell ", n2s(i), "."
          error stop 1
      end select
    end do

    call upscaler%execute(pet_fac_lai_l0, pet_fac_lai_l1, upscaling_operator=up_h_mean)
    deallocate(pet_fac_lai_l0)
  end subroutine mpr_bridge_pet_lai

  !> \brief Bridge PET aspect correction through the legacy aspect routine.
  subroutine mpr_bridge_pet_aspect(level0, aspect_l0, param, upscaler, pet_fac_aspect_l1)
    type(grid_t), intent(in), target :: level0
    real(dp), dimension(:), intent(in) :: aspect_l0
    real(dp), dimension(:), intent(in) :: param
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: pet_fac_aspect_l1
    integer(i4), allocatable :: id0(:)
    real(dp), allocatable :: latitude_l0(:)
    real(dp), allocatable :: pet_fac_aspect_l0(:)
    integer(i4) :: i

    if (size(param) < 3_i4) then
      log_fatal(*) "MPR bridge: PET aspect parameter set must contain 3 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(aspect_l0) /= level0%ncells) then
      log_fatal(*) "MPR bridge: PET aspect input size does not match packed level0 cells."
      error stop 1
    end if
    if (size(pet_fac_aspect_l1) /= upscaler%target_grid%ncells) then
      log_fatal(*) "MPR bridge: PET aspect output size does not match packed level1 cells."
      error stop 1
    end if

    allocate(id0(level0%ncells))
    id0 = [(i, i = 1_i4, level0%ncells)]
    call mpr_bridge_get_level0_latitude(level0, latitude_l0)
    allocate(pet_fac_aspect_l0(level0%ncells))

    call pet_correctbyASP(id0, latitude_l0, max(aspect_l0, 1.0_dp), param(1:3), nodata_dp, pet_fac_aspect_l0)
    call upscaler%execute(pet_fac_aspect_l0, pet_fac_aspect_l1, upscaling_operator=up_a_mean)

    deallocate(id0)
    deallocate(latitude_l0)
    deallocate(pet_fac_aspect_l0)
  end subroutine mpr_bridge_pet_aspect

  !> \brief Bridge Hargreaves-Samani PET coefficients.
  subroutine mpr_bridge_pet_hargreaves(level0, aspect_l0, param, upscaler, pet_fac_aspect_l1, pet_coeff_hs_l1)
    type(grid_t), intent(in), target :: level0
    real(dp), dimension(:), intent(in) :: aspect_l0
    real(dp), dimension(:), intent(in) :: param
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: pet_fac_aspect_l1
    real(dp), dimension(:), intent(out) :: pet_coeff_hs_l1

    if (size(param) < 4_i4) then
      log_fatal(*) "MPR bridge: PET Hargreaves parameter set must contain 4 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(pet_coeff_hs_l1) /= size(pet_fac_aspect_l1)) then
      log_fatal(*) "MPR bridge: PET Hargreaves outputs must have matching packed L1 sizes."
      error stop 1
    end if

    call mpr_bridge_pet_aspect(level0, aspect_l0, param(1:3), upscaler, pet_fac_aspect_l1)
    pet_coeff_hs_l1 = param(4)
  end subroutine mpr_bridge_pet_hargreaves

  !> \brief Bridge Priestley-Taylor PET coefficients.
  subroutine mpr_bridge_pet_priestley_taylor(level0, lai_l0, param, upscaler, pet_coeff_pt_l1)
    type(grid_t), intent(in), target :: level0
    real(dp), dimension(:, :), intent(in) :: lai_l0
    real(dp), dimension(:), intent(in) :: param
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:, :), intent(out) :: pet_coeff_pt_l1
    integer(i4), allocatable :: id0(:)
    integer(i4), allocatable :: upper_bound1(:)
    integer(i4), allocatable :: lower_bound1(:)
    integer(i4), allocatable :: left_bound1(:)
    integer(i4), allocatable :: right_bound1(:)

    if (size(param) < 2_i4) then
      log_fatal(*) "MPR bridge: PET Priestley-Taylor parameter set must contain 2 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(lai_l0, 1) /= level0%ncells) then
      log_fatal(*) "MPR bridge: PET Priestley-Taylor LAI input size does not match packed level0 cells."
      error stop 1
    end if
    if (size(pet_coeff_pt_l1, 1) /= upscaler%target_grid%ncells) then
      log_fatal(*) "MPR bridge: PET Priestley-Taylor output size does not match packed level1 cells."
      error stop 1
    end if
    if (size(pet_coeff_pt_l1, 2) /= size(lai_l0, 2)) then
      log_fatal(*) "MPR bridge: PET Priestley-Taylor output LAI dimension does not match LAI cache."
      error stop 1
    end if

    call mpr_bridge_build_l0_to_l1_indices(level0, upscaler, id0, upper_bound1, lower_bound1, left_bound1, right_bound1)
    call priestley_taylor_alpha(lai_l0, param(1:2), level0%mask, nodata_dp, id0, upscaler%n_subcells, &
      upper_bound1, lower_bound1, left_bound1, right_bound1, pet_coeff_pt_l1)

    deallocate(id0)
    deallocate(upper_bound1)
    deallocate(lower_bound1)
    deallocate(left_bound1)
    deallocate(right_bound1)
  end subroutine mpr_bridge_pet_priestley_taylor

  !> \brief Bridge Penman-Monteith PET resistances.
  subroutine mpr_bridge_pet_penman_monteith(level0, land_cover_l0, lai_l0, param, upscaler, resist_aero_l1, resist_surf_l1)
    type(grid_t), intent(in), target :: level0
    integer(i4), dimension(:), intent(in) :: land_cover_l0
    real(dp), dimension(:, :), intent(in) :: lai_l0
    real(dp), dimension(:), intent(in) :: param
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:, :), intent(out) :: resist_aero_l1
    real(dp), dimension(:, :), intent(out) :: resist_surf_l1
    integer(i4), allocatable :: id0(:)
    integer(i4), allocatable :: upper_bound1(:)
    integer(i4), allocatable :: lower_bound1(:)
    integer(i4), allocatable :: left_bound1(:)
    integer(i4), allocatable :: right_bound1(:)
    if (size(param) < 7_i4) then
      log_fatal(*) "MPR bridge: PET Penman-Monteith parameter set must contain 7 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(land_cover_l0) /= level0%ncells .or. size(lai_l0, 1) /= level0%ncells) then
      log_fatal(*) "MPR bridge: PET Penman-Monteith L0 inputs must match packed level0 cells."
      error stop 1
    end if
    if (size(resist_aero_l1, 1) /= upscaler%target_grid%ncells .or. size(resist_surf_l1, 1) /= upscaler%target_grid%ncells) then
      log_fatal(*) "MPR bridge: PET Penman-Monteith outputs must match packed level1 cells."
      error stop 1
    end if
    if (size(resist_aero_l1, 2) /= size(lai_l0, 2) .or. size(resist_surf_l1, 2) /= size(lai_l0, 2)) then
      log_fatal(*) "MPR bridge: PET Penman-Monteith output LAI dimensions must match the LAI cache."
      error stop 1
    end if

    call mpr_bridge_build_l0_to_l1_indices(level0, upscaler, id0, upper_bound1, lower_bound1, left_bound1, right_bound1)
    call aerodynamical_resistance(lai_l0, land_cover_l0, param(1:6), level0%mask, id0, upscaler%n_subcells, &
      upper_bound1, lower_bound1, left_bound1, right_bound1, resist_aero_l1)
    call bulksurface_resistance(lai_l0, param(7), level0%mask, nodata_dp, id0, upscaler%n_subcells, &
      upper_bound1, lower_bound1, left_bound1, right_bound1, resist_surf_l1)

    deallocate(id0)
    deallocate(upper_bound1)
    deallocate(lower_bound1)
    deallocate(left_bound1)
    deallocate(right_bound1)
  end subroutine mpr_bridge_pet_penman_monteith

  !> \brief Bridge soil-moisture parameter generation through the legacy soil routines.
  subroutine mpr_bridge_soil_moisture(soil_case, neutron_case, param, land_cover_l0, soil_id_l0, level0, upscaler, &
      thresh_jarvis_l1, sm_exponent_l1, sm_saturation_l1, sm_field_capacity_l1, wilting_point_l1, f_roots_l1, &
      sm_deficit_fc_l0, ks_var_h_l0, ks_var_v_l0, bulk_density_l1, lattice_water_l1, cosmic_l3_l1, neutron_param)
    integer(i4), intent(in) :: soil_case
    integer(i4), intent(in) :: neutron_case
    real(dp), dimension(:), intent(in) :: param
    integer(i4), dimension(:), intent(in) :: land_cover_l0
    integer(i4), dimension(:, :), intent(in) :: soil_id_l0
    type(grid_t), intent(in), target :: level0
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: thresh_jarvis_l1
    real(dp), dimension(:, :), intent(out) :: sm_exponent_l1
    real(dp), dimension(:, :), intent(out) :: sm_saturation_l1
    real(dp), dimension(:, :), intent(out) :: sm_field_capacity_l1
    real(dp), dimension(:, :), intent(out) :: wilting_point_l1
    real(dp), dimension(:, :), intent(out) :: f_roots_l1
    real(dp), dimension(:), intent(out), optional :: sm_deficit_fc_l0
    real(dp), dimension(:), intent(out), optional :: ks_var_h_l0
    real(dp), dimension(:), intent(out), optional :: ks_var_v_l0
    real(dp), dimension(:, :), intent(out), optional :: bulk_density_l1
    real(dp), dimension(:, :), intent(out), optional :: lattice_water_l1
    real(dp), dimension(:, :), intent(out), optional :: cosmic_l3_l1
    real(dp), dimension(:), intent(in), optional :: neutron_param
    integer(i4), allocatable :: id0(:)
    integer(i4) :: i
    integer(i4), dimension(nProcesses, 3) :: legacy_process_matrix
    integer(i4) :: msoil
    integer(i4) :: m_lc
    integer(i4) :: mtill
    integer(i4) :: m_hor
    integer(i4) :: n_cells0
    integer(i4) :: n_cells1
    integer(i4) :: soil_param_end
    integer(i4) :: horizon_param_start
    integer(i4) :: horizon_param_end
    integer(i4) :: x_lb
    integer(i4) :: x_ub
    integer(i4) :: y_lb
    integer(i4) :: y_ub
    real(dp), allocatable :: thetaS_till(:, :, :)
    real(dp), allocatable :: thetaFC_till(:, :, :)
    real(dp), allocatable :: thetaPW_till(:, :, :)
    real(dp), allocatable :: thetaS(:, :)
    real(dp), allocatable :: thetaFC(:, :)
    real(dp), allocatable :: thetaPW(:, :)
    real(dp), allocatable :: Ks(:, :, :)
    real(dp), allocatable :: Db(:, :, :)
    real(dp), allocatable :: KsVar_H0(:)
    real(dp), allocatable :: KsVar_V0(:)
    real(dp), allocatable :: SMs_FC0(:)
    real(dp), allocatable :: latWat_till(:, :, :)
    real(dp), allocatable :: COSMIC_L3_till(:, :, :)
    real(dp), allocatable :: latWat(:, :)
    real(dp), allocatable :: COSMIC_L3(:, :)
    real(dp), allocatable :: bulk_density_tmp(:, :)
    real(dp), allocatable :: lattice_water_tmp(:, :)
    real(dp), allocatable :: cosmic_l3_tmp(:, :)
    real(dp), allocatable :: legacy_param(:)
    integer(i4), allocatable :: upper_bound1(:)
    integer(i4), allocatable :: lower_bound1(:)
    integer(i4), allocatable :: left_bound1(:)
    integer(i4), allocatable :: right_bound1(:)

    if (.not.allocated(soilDB%is_present)) then
      log_fatal(*) "MPR bridge: soil database not initialized before soil-moisture generation."
      error stop 1
    end if

    legacy_process_matrix = 0_i4
    legacy_process_matrix(3, 1) = soil_case
    legacy_process_matrix(10, 1) = neutron_case

    select case (soil_case)
      case (1_i4)
        if (size(param) < 16_i4) then
          log_fatal(*) "MPR bridge: soil moisture case 1 expects 16 parameters, got ", n2s(size(param, kind=i4)), "."
          error stop 1
        end if
        soil_param_end = 13_i4
        horizon_param_start = 14_i4
        horizon_param_end = 17_i4
        thresh_jarvis_l1 = nodata_dp
      case (2_i4)
        if (size(param) < 17_i4) then
          log_fatal(*) "MPR bridge: soil moisture case 2 expects 17 parameters, got ", n2s(size(param, kind=i4)), "."
          error stop 1
        end if
        soil_param_end = 13_i4
        horizon_param_start = 14_i4
        horizon_param_end = 17_i4
        thresh_jarvis_l1 = param(17)
      case (3_i4)
        if (size(param) < 21_i4) then
          log_fatal(*) "MPR bridge: soil moisture case 3 expects 21 parameters, got ", n2s(size(param, kind=i4)), "."
          error stop 1
        end if
        soil_param_end = 13_i4
        horizon_param_start = 14_i4
        horizon_param_end = 21_i4
        thresh_jarvis_l1 = param(21)
      case (4_i4)
        if (size(param) < 20_i4) then
          log_fatal(*) "MPR bridge: soil moisture case 4 expects 20 parameters, got ", n2s(size(param, kind=i4)), "."
          error stop 1
        end if
        soil_param_end = 13_i4
        horizon_param_start = 14_i4
        horizon_param_end = 21_i4
        thresh_jarvis_l1 = nodata_dp
      case default
        log_fatal(*) "MPR bridge: unsupported soil moisture process case ", n2s(soil_case), "."
        error stop 1
    end select

    ! The legacy PTF API expects the unit conversion as its thirteenth parameter.
    allocate(legacy_param(size(param) + 1_i4))
    legacy_param(:12) = param(:12)
    legacy_param(13) = inch_per_hour_to_centimetre_per_day
    legacy_param(14:) = param(13:)

    n_cells0 = size(land_cover_l0)
    if (size(soil_id_l0, 1) /= n_cells0) then
      log_fatal(*) "MPR bridge: soil_id and land_cover inputs must have matching packed L0 sizes."
      error stop 1
    end if
    n_cells1 = size(sm_exponent_l1, 1)
    if (size(sm_exponent_l1, 2) /= nSoilHorizons_mHM) then
      log_fatal(*) "MPR bridge: soil cache horizon dimension ", n2s(size(sm_exponent_l1, 2)), &
        " does not match n_layers=", n2s(nSoilHorizons_mHM), "."
      error stop 1
    end if

    msoil = size(soilDB%is_present)
    m_lc = maxval(land_cover_l0, land_cover_l0 /= nodata_i4)
    if (m_lc < 1_i4) then
      log_fatal(*) "MPR bridge: land-cover map contains no valid IDs for soil-moisture generation."
      error stop 1
    end if
    if (iFlag_soilDB == 0_i4) then
      mtill = maxval(soilDB%nTillHorizons, soilDB%nTillHorizons /= nodata_i4)
      m_hor = maxval(soilDB%nHorizons, soilDB%nHorizons /= nodata_i4)
    else if (iFlag_soilDB == 1_i4) then
      mtill = 1_i4
      m_hor = 1_i4
    else
      log_fatal(*) "MPR bridge: unsupported legacy soil_db_mode=", n2s(iFlag_soilDB), "."
      error stop 1
    end if

    allocate(id0(n_cells0))
    do i = 1_i4, n_cells0
      id0(i) = i
    end do

    allocate(thetaS_till(msoil, mtill, m_lc))
    allocate(thetaFC_till(msoil, mtill, m_lc))
    allocate(thetaPW_till(msoil, mtill, m_lc))
    allocate(thetaS(msoil, m_hor))
    allocate(thetaFC(msoil, m_hor))
    allocate(thetaPW(msoil, m_hor))
    allocate(Ks(msoil, m_hor, m_lc))
    allocate(Db(msoil, m_hor, m_lc))
    allocate(KsVar_H0(n_cells0))
    allocate(KsVar_V0(n_cells0))
    allocate(SMs_FC0(n_cells0))
    allocate(latWat_till(msoil, mtill, m_lc))
    allocate(COSMIC_L3_till(msoil, mtill, m_lc))
    allocate(latWat(msoil, m_hor))
    allocate(COSMIC_L3(msoil, m_hor))
    allocate(bulk_density_tmp(n_cells1, size(sm_exponent_l1, 2)))
    allocate(lattice_water_tmp(n_cells1, size(sm_exponent_l1, 2)))
    allocate(cosmic_l3_tmp(n_cells1, size(sm_exponent_l1, 2)))
    allocate(upper_bound1(n_cells1))
    allocate(lower_bound1(n_cells1))
    allocate(left_bound1(n_cells1))
    allocate(right_bound1(n_cells1))

    latWat_till = 1.0e-6_dp
    COSMIC_L3_till = 1.0e-6_dp
    latWat = 1.0e-6_dp
    COSMIC_L3 = 1.0e-6_dp
    bulk_density_tmp = nodata_dp
    lattice_water_tmp = nodata_dp
    cosmic_l3_tmp = nodata_dp
    do i = 1_i4, n_cells1
      call upscaler%coarse_bounds(int(i, i8), x_lb, x_ub, y_lb, y_ub)
      ! Legacy MPR upscaling routines slice arrays as (dim1, dim2). The new forces grid stores data as (x, y),
      ! so feed x bounds into the first legacy index pair and y bounds into the second one.
      upper_bound1(i) = x_lb
      lower_bound1(i) = x_ub
      left_bound1(i) = y_lb
      right_bound1(i) = y_ub
    end do

    call mpr_sm(legacy_param(:soil_param_end), legacy_process_matrix, &
      soilDB%is_present, soilDB%nHorizons, soilDB%nTillHorizons, &
      soilDB%sand, soilDB%clay, soilDB%DbM, id0, soil_id_l0, land_cover_l0, &
      thetaS_till, thetaFC_till, thetaPW_till, thetaS, thetaFC, thetaPW, Ks, Db, KsVar_H0, KsVar_V0, SMs_FC0)

    if (present(sm_deficit_fc_l0)) then
      if (size(sm_deficit_fc_l0) /= size(SMs_FC0)) then
        log_fatal(*) "MPR bridge: sm_deficit_fc_l0 output has wrong size."
        error stop 1
      end if
      sm_deficit_fc_l0 = SMs_FC0
    end if
    if (present(ks_var_h_l0)) then
      if (size(ks_var_h_l0) /= size(KsVar_H0)) then
        log_fatal(*) "MPR bridge: ks_var_h_l0 output has wrong size."
        error stop 1
      end if
      ks_var_h_l0 = KsVar_H0
    end if
    if (present(ks_var_v_l0)) then
      if (size(ks_var_v_l0) /= size(KsVar_V0)) then
        log_fatal(*) "MPR bridge: ks_var_v_l0 output has wrong size."
        error stop 1
      end if
      ks_var_v_l0 = KsVar_V0
    end if

    if (neutron_case > 0_i4) then
      if (.not.present(neutron_param)) then
        log_fatal(*) "MPR bridge: neutron parameters must be provided when the neutrons process is active."
        error stop 1
      end if
      call mpr_neutrons( &
        neutron_case, &
        neutron_param, &
        soilDB%is_present, soilDB%nHorizons, soilDB%nTillHorizons, land_cover_l0, &
        soilDB%clay, soilDB%DbM, Db, COSMIC_L3_till, latWat_till, COSMIC_L3, latWat)
    end if

    call mpr_SMhorizons(legacy_param(horizon_param_start:horizon_param_end), legacy_process_matrix, &
      iFlag_soilDB, nSoilHorizons_mHM, HorizonDepth_mHM, land_cover_l0, soil_id_l0, &
      soilDB%nHorizons, soilDB%nTillHorizons, thetaS_till, thetaFC_till, thetaPW_till, thetaS, thetaFC, thetaPW, &
      soilDB%Wd, Db, soilDB%DbM, soilDB%RZdepth, level0%mask, id0, &
      upper_bound1, lower_bound1, left_bound1, right_bound1, upscaler%n_subcells, &
      sm_exponent_l1, sm_saturation_l1, sm_field_capacity_l1, wilting_point_l1, f_roots_l1, &
      latWat_till, COSMIC_L3_till, latWat, COSMIC_L3, bulk_density_tmp, lattice_water_tmp, cosmic_l3_tmp)

    if (present(bulk_density_l1)) then
      if (any(shape(bulk_density_l1) /= shape(bulk_density_tmp))) then
        log_fatal(*) "MPR bridge: bulk_density_l1 output has wrong shape."
        error stop 1
      end if
      bulk_density_l1 = bulk_density_tmp
    end if
    if (present(lattice_water_l1)) then
      if (any(shape(lattice_water_l1) /= shape(lattice_water_tmp))) then
        log_fatal(*) "MPR bridge: lattice_water_l1 output has wrong shape."
        error stop 1
      end if
      lattice_water_l1 = lattice_water_tmp
    end if
    if (present(cosmic_l3_l1)) then
      if (any(shape(cosmic_l3_l1) /= shape(cosmic_l3_tmp))) then
        log_fatal(*) "MPR bridge: cosmic_l3_l1 output has wrong shape."
        error stop 1
      end if
      cosmic_l3_l1 = cosmic_l3_tmp
    end if

    deallocate(id0)
    deallocate(thetaS_till)
    deallocate(thetaFC_till)
    deallocate(thetaPW_till)
    deallocate(thetaS)
    deallocate(thetaFC)
    deallocate(thetaPW)
    deallocate(Ks)
    deallocate(Db)
    deallocate(KsVar_H0)
    deallocate(KsVar_V0)
    deallocate(SMs_FC0)
    deallocate(latWat_till)
    deallocate(COSMIC_L3_till)
    deallocate(latWat)
    deallocate(COSMIC_L3)
    deallocate(bulk_density_tmp)
    deallocate(lattice_water_tmp)
    deallocate(cosmic_l3_tmp)
    deallocate(legacy_param)
    deallocate(upper_bound1)
    deallocate(lower_bound1)
    deallocate(left_bound1)
    deallocate(right_bound1)
  end subroutine mpr_bridge_soil_moisture

  !> \brief Bridge runoff/interflow parameter generation with legacy mpr_runoff.
  subroutine mpr_bridge_runoff_param(land_cover_l0, slope_emp_l0, sm_deficit_fc_l0, ks_var_h_l0, level0, upscaler, &
      param, thresh_unsat_l1, k_fastflow_l1, k_slowflow_l1, alpha_l1)
    integer(i4), dimension(:), intent(in) :: land_cover_l0
    real(dp), dimension(:), intent(in) :: slope_emp_l0
    real(dp), dimension(:), intent(in) :: sm_deficit_fc_l0
    real(dp), dimension(:), intent(in) :: ks_var_h_l0
    type(grid_t), intent(in), target :: level0
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(in) :: param
    real(dp), dimension(:), intent(out) :: thresh_unsat_l1
    real(dp), dimension(:), intent(out) :: k_fastflow_l1
    real(dp), dimension(:), intent(out) :: k_slowflow_l1
    real(dp), dimension(:), intent(out) :: alpha_l1
    integer(i4), allocatable :: id0(:)
    integer(i4), allocatable :: upper_bound1(:)
    integer(i4), allocatable :: lower_bound1(:)
    integer(i4), allocatable :: left_bound1(:)
    integer(i4), allocatable :: right_bound1(:)
    integer(i4) :: i
    integer(i4) :: x_lb
    integer(i4) :: x_ub
    integer(i4) :: y_lb
    integer(i4) :: y_ub
    integer(i4) :: n_cells0
    integer(i4) :: n_cells1

    if (size(param) < 5_i4) then
      log_fatal(*) "MPR bridge: interflow parameter set must contain 5 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if

    n_cells0 = size(land_cover_l0)
    if (size(slope_emp_l0) /= n_cells0 .or. size(sm_deficit_fc_l0) /= n_cells0 .or. size(ks_var_h_l0) /= n_cells0) then
      log_fatal(*) "MPR bridge: runoff L0 inputs must have matching packed sizes."
      error stop 1
    end if
    n_cells1 = size(alpha_l1)
    if (size(thresh_unsat_l1) /= n_cells1 .or. size(k_fastflow_l1) /= n_cells1 .or. size(k_slowflow_l1) /= n_cells1) then
      log_fatal(*) "MPR bridge: runoff L1 outputs must have matching sizes."
      error stop 1
    end if

    allocate(id0(n_cells0))
    do i = 1_i4, n_cells0
      id0(i) = i
    end do
    allocate(upper_bound1(n_cells1), lower_bound1(n_cells1), left_bound1(n_cells1), right_bound1(n_cells1))
    do i = 1_i4, n_cells1
      call upscaler%coarse_bounds(int(i, i8), x_lb, x_ub, y_lb, y_ub)
      upper_bound1(i) = x_lb
      lower_bound1(i) = x_ub
      left_bound1(i) = y_lb
      right_bound1(i) = y_ub
    end do

    call mpr_runoff( &
      land_cover_l0, level0%mask, sm_deficit_fc_l0, slope_emp_l0, ks_var_h_l0, param(:5), id0, &
      upper_bound1, lower_bound1, left_bound1, right_bound1, upscaler%n_subcells, &
      thresh_unsat_l1, k_fastflow_l1, k_slowflow_l1, alpha_l1)

    deallocate(id0)
    deallocate(upper_bound1)
    deallocate(lower_bound1)
    deallocate(left_bound1)
    deallocate(right_bound1)
  end subroutine mpr_bridge_runoff_param

  !> \brief Bridge percolation and karst-loss parameter generation with legacy-compatible formulas.
  subroutine mpr_bridge_karstic_param(param, geo_unit_l0, geo_unit_list, geo_karstic, sm_deficit_fc_l0, ks_var_v_l0, &
      upscaler, f_karst_loss_l1, k_percolation_l1)
    real(dp), dimension(:), intent(in) :: param
    integer(i4), dimension(:), intent(in) :: geo_unit_l0
    integer(i4), dimension(:), intent(in) :: geo_unit_list
    integer(i4), dimension(:), intent(in) :: geo_karstic
    real(dp), dimension(:), intent(in) :: sm_deficit_fc_l0
    real(dp), dimension(:), intent(in) :: ks_var_v_l0
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: f_karst_loss_l1
    real(dp), dimension(:), intent(out) :: k_percolation_l1
    real(dp), allocatable :: percolation_l0(:)
    real(dp), allocatable :: karst_fraction_l1(:)
    real(dp), allocatable :: karst_class_fraction_l1(:)
    integer(i4) :: i

    if (size(param) < 2_i4) then
      log_fatal(*) "MPR bridge: percolation parameter set must contain 2 values, got ", n2s(size(param, kind=i4)), "."
      error stop 1
    end if
    if (size(geo_unit_l0) /= size(sm_deficit_fc_l0) .or. size(geo_unit_l0) /= size(ks_var_v_l0)) then
      log_fatal(*) "MPR bridge: karst/percolation L0 inputs must have matching packed sizes."
      error stop 1
    end if
    if (size(geo_unit_list) /= size(geo_karstic)) then
      log_fatal(*) "MPR bridge: geo_unit and geo_karstic LUT sizes do not match."
      error stop 1
    end if

    allocate(percolation_l0(size(geo_unit_l0)))
    allocate(karst_fraction_l1(size(f_karst_loss_l1)))
    allocate(karst_class_fraction_l1(size(f_karst_loss_l1)))

    percolation_l0 = param(1) * (1.0_dp + sm_deficit_fc_l0) / (1.0_dp + ks_var_v_l0)
    call upscaler%execute(percolation_l0, k_percolation_l1, upscaling_operator=up_a_mean)
    k_percolation_l1 = merge(2.0_dp, k_percolation_l1, k_percolation_l1 < 2.0_dp)

    karst_fraction_l1 = 0.0_dp
    do i = 1_i4, size(geo_unit_list)
      if (geo_karstic(i) == 0_i4) cycle
      call upscaler%execute(geo_unit_l0, karst_class_fraction_l1, upscaling_operator=up_fraction, class_id=geo_unit_list(i))
      karst_fraction_l1 = karst_fraction_l1 + karst_class_fraction_l1
    end do
    f_karst_loss_l1 = 1.0_dp - karst_fraction_l1 * param(2)

    deallocate(percolation_l0)
    deallocate(karst_fraction_l1)
    deallocate(karst_class_fraction_l1)
  end subroutine mpr_bridge_karstic_param

  !> \brief Bridge the sealed-surface runoff threshold parameter.
  subroutine mpr_bridge_sealed_threshold(param, thresh_sealed_l1)
    real(dp), dimension(:), intent(in) :: param
    real(dp), dimension(:), intent(out) :: thresh_sealed_l1

    if (size(param) < 1_i4) then
      log_fatal(*) "MPR bridge: direct-runoff parameter set must contain at least 1 value."
      error stop 1
    end if
    thresh_sealed_l1 = param(1)
  end subroutine mpr_bridge_sealed_threshold

  !> \brief Bridge geological baseflow parameter generation with legacy-compatible class mapping.
  subroutine mpr_bridge_baseflow_param(param, geo_unit_l0, geo_unit_list, geo_param_index, upscaler, k_baseflow_l1)
    real(dp), dimension(:), intent(in) :: param
    integer(i4), dimension(:), intent(in) :: geo_unit_l0
    integer(i4), dimension(:), intent(in) :: geo_unit_list
    integer(i4), dimension(:), intent(in) :: geo_param_index
    class(scaler_t), intent(inout), target :: upscaler
    real(dp), dimension(:), intent(out) :: k_baseflow_l1
    real(dp), allocatable :: baseflow_l0(:)
    integer(i4) :: i
    integer(i4) :: j
    integer(i4) :: class_pos

    if (size(geo_param_index) /= size(geo_unit_list)) then
      log_fatal(*) "MPR bridge: geo_param_index and geo_unit LUT sizes do not match."
      error stop 1
    end if
    if (any(geo_param_index < 1_i4) .or. any(geo_param_index > size(param))) then
      log_fatal(*) "MPR bridge: geology LUT GeoParam index is outside available baseflow parameters."
      error stop 1
    end if

    allocate(baseflow_l0(size(geo_unit_l0)))
    do i = 1_i4, size(geo_unit_l0)
      class_pos = 0_i4
      do j = 1_i4, size(geo_unit_list)
        if (geo_unit_list(j) == geo_unit_l0(i)) then
          class_pos = j
          exit
        end if
      end do
      if (class_pos < 1_i4) then
        log_fatal(*) "MPR bridge: geological unit ", n2s(geo_unit_l0(i)), " missing in geology LUT."
        error stop 1
      end if
      baseflow_l0(i) = param(geo_param_index(class_pos))
    end do

    call upscaler%execute(baseflow_l0, k_baseflow_l1, upscaling_operator=up_a_mean)
    deallocate(baseflow_l0)
  end subroutine mpr_bridge_baseflow_param

end module mo_mpr_legacy_bridge
