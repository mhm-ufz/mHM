!> \file mo_meteo_constants.f90
!> \copydoc mo_meteo_constants

!> \brief Provides meteorological constants
!> \details Provides coefficients for PET, extraterrestrial radiation, and vapour pressure.
! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting
!> \authors Matthias Cuntz
!> \date Nov 2011
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_meteo
module mo_meteo_constants
  use mo_kind, only: dp
  implicit none
  private

  !> Hargreaves-Samani ref. ET formula [deg C]
  real(dp), public, parameter :: HarSamConst = 17.800_dp

  ! Duffie formula for computing extraterrestrial radiation
  real(dp), public, parameter :: DuffieDr = 0.0330_dp
  real(dp), public, parameter :: DuffieDelta1 = 0.4090_dp
  real(dp), public, parameter :: DuffieDelta2 = 1.3900_dp

  !> Tetens's formula to calculate saturated vapour pressure
  real(dp), public, parameter :: tetens_c1 = 0.6108_dp
  real(dp), public, parameter :: tetens_c2 = 17.270_dp
  real(dp), public, parameter :: tetens_c3 = 237.30_dp
  !> calculation of the slope of the saturation vapour pressure curve following Tetens
  real(dp), public, parameter :: satpressureslope1 = 4098.0_dp

  ! Physical
  !> Refractive index of water [-]
  real(dp), public, parameter :: RefracIndex_w = 1.33_dp
  !> Emissivity of water [-]
  real(dp), public, parameter :: Emissivity_w = 0.96_dp

  !> Latent heat of vaporisation of water
  real(dp), public, parameter :: lat_heat_vap_c1 = 2.501_dp
  real(dp), public, parameter :: lat_heat_vap_c2 = 0.002361_dp

  !> Atmospheric pressure
  real(dp), public, parameter :: atm_pressure_c1 = 101.3_dp
  real(dp), public, parameter :: atm_pressure_c2 = 293.0_dp
  real(dp), public, parameter :: atm_pressure_c3 = 0.0065_dp
  real(dp), public, parameter :: atm_pressure_c4 = 5.256_dp

  !> Psychrometric constant
  real(dp), public, parameter :: psychro_constant_c1 = 0.0016286_dp

  ! Lake evaporation formulae
  !> Penman (1948) formula to calculate lake evaporation
  real(dp), public, parameter :: penman1948_1  = 6.43_dp
  real(dp), public, parameter :: penman1948_2  = 0.536_dp

  !> Priestley-Taylor (1972) coefficient
  real(dp), public, parameter :: ptcoeff       = 1.26_dp

  !> Calibration-free Priestley-Taylor, after Su and Singh (2023) adapted to the
  !> open-water Penman formulation: the constant ptcoeff is replaced by an analytical coefficient
  !> alpha = 1 + (gamma/Delta) * penman1948_1 (1 + penman1948_2 u2) (es - e) / (Rn + Ah), evaluated
  !> with a climatological (not instantaneous) windspeed so that noisy per-step wind does not swing
  !> the evaporation, and applied to the live equilibrium term. alpha is bounded to [1, cfpt_coeff_max]:
  !> the lower bound of 1 matches the no-condensation design of the leaf routines (the es - e < 0
  !> condensation regime is already zeroed there), and the upper cap guards the A = Rn + Ah -> 0 limit,
  !> where the coefficient would otherwise diverge.
  !> Available-energy floor [J m-2 d-1] (about 1.2 W m-2 daily mean; small but strictly positive).
  real(dp), public, parameter :: cfpt_avail_floor = 1.0e5_dp
  !> Upper cap on the calibration-free Priestley-Taylor coefficient [-] (open-water alpha is physically
  !> well below this; the cap only guards the A -> 0 limit).
  real(dp), public, parameter :: cfpt_coeff_max   = 5.0_dp

  ! Meteorology
  !> Constant water albedo (Cogley)
  real(dp), public, parameter :: albedo_w_const = 0.08_dp

  !> Fresnel's formula for water albedo estimation
  real(dp), public, parameter :: fresnel_c    = 50._dp

  !> Still water albedo cutoff values for the Fresnel equation
  real(dp), public, parameter :: albedo_w_low = 0.025_dp
  real(dp), public, parameter :: albedo_w_high= 0.5_dp

  !> Still water albedo fitted Weibull parameters (Grishchenko 1959)
  real(dp), public, parameter :: grishchenko_k = 1.517_dp
  real(dp), public, parameter :: grishchenko_l = 1.236_dp
  real(dp), public, parameter :: grishchenko_s = 52.279_dp
  real(dp), public, parameter :: grishchenko_c = 4.493_dp

  !> Boltzmann equation exponent for outgoing longwave radiation
  real(dp), public, parameter :: boltzmann_e  = 4._dp

  !> Hellmann's exponential law to correct windspeed height
  real(dp), public, parameter :: std_windspeed_ht = 2.0_dp
  real(dp), public, parameter :: Hellmann_exp = 0.10_dp

                                ! Note on Hellmann_exp:
                                !
                                ! Unstable air above open water surface:  0.06
                                ! Neutral air above open water surface: 0.10
                                ! Stable air above open water surface:  0.27
                                !
                                !       Kaltschmitt, M., Streicher, W., Wiese, Q. (2007).
                                !       Renewable energy: technology, economics, and environment. page 55
                                !       Springer, https://10.1007/3-540-70949-5

  !> Equation of time (Meeus)
  real(dp), public, parameter :: meeus_1 = 229.18_dp
  real(dp), public, parameter :: meeus_2 = 0.000075_dp
  real(dp), public, parameter :: meeus_3 = 0.001868_dp
  real(dp), public, parameter :: meeus_4 = -0.032077_dp
  real(dp), public, parameter :: meeus_5 = -0.014615_dp
  real(dp), public, parameter :: meeus_6 = -0.040849_dp

end module mo_meteo_constants
