!>       \file mo_mlm_constants.f90

!>       \brief Provides mLM specific constants

!>       \details Provides mLM specific constants such as the empirical coefficients of
!>                the lake evaporation formulae.

!>       \authors Pallav Kumar Shrestha

!>       \date Jul 2018

! Modifications:
! Pallav Shrestha, Oct 2026 - v6 migration: recreated in the v6 tree with the minimal-scope
!     subset actually used by the ported lake routines (the Penman-1948 evaporation path).
!     Values are unchanged from v5; generic physical constants continue to come from
!     mo_constants (forces). Entries are added here as further leaf routines are ported.
! Pallav Shrestha, Oct 2026 - v6 migration (Phase 5): moved the lake water-balance solver
!     constants (error_lake_vol, nIterations_lake) here from mo_lake_water_balance.
! Pallav Shrestha, Oct 2026 - added the lake bank soil-evaporation ratio (bank_soil_evap_ratio)
!     used by the lake banks footprint reconciliation in mo_mlm_container.

module mo_mlm_constants

  use mo_kind, only : i4, dp

  implicit none

  ! Lake water-balance solver
  !> Allowed lake volume closure error [m3]. In v5 this was 1e-6 mcm (= 1 m3), carried in
  !> mcm; with the v6 SI volume (m3) the same physical tolerance is 1 m3.
  real(dp), public, parameter :: error_lake_vol = 1.0_dp
  !> Maximum number of secant iterations in the lake water-balance solver [-]
  integer(i4), public, parameter :: nIterations_lake = 20_i4

  ! Lake banks (footprint reconciliation)
  !> Bank soil evaporation as a fraction of open-water lake evaporation (Es = c_soil * E), used to
  !> account the net atmospheric flux over the exposed DEM-footprint ring (A_footprint - A_water)
  !> as same-step lake inflow. Literature on saturated bare soil versus open water (same lysimeter
  !> group) finds the POTENTIAL evaporation of a saturated bare-soil surface is slightly HIGHER
  !> than open water, more so in spring/summer, because of the lower soil albedo and the surface
  !> energy balance (Li et al., 2020, J. Hydrol. 590:125140; Li et al., 2022, J. Hydrol. 613:128194,
  !> "role of soil texture"; Li et al., 2023, Hydrol. Process. 37(6):e14894). So a freshly exposed,
  !> near-saturated bank evaporates at close to - or just above - the open-water rate (c_soil ~ 1).
  !> A real exposed bank also dries down between wettings, which pulls actual evaporation below the
  !> saturated-potential rate. The single constant 0.8 is a moist-bank value between the dry limit
  !> and the saturated (>=1) limit; make it level/moisture-dependent later for endorheic beds. [-]
  real(dp), public, parameter :: bank_soil_evap_ratio = 0.8_dp

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

end module mo_mlm_constants
