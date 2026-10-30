!> \file mo_mlm_constants.f90
!> \copydoc mo_mlm_constants

!>       \brief Provides mLM specific constants

!>       \details Provides mLM specific constants such as the water-balance solver tolerances
!>                and the empirical exposed-bank evaporation ratio.

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

!> \version 0.1
!> \changelog
!! - October 2026: v6 component ownership and restart integration; retained donor arithmetic.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
module mo_mlm_constants

  use mo_kind, only : i4, dp

  implicit none
  private

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

end module mo_mlm_constants
