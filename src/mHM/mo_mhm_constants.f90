!> \file mo_mhm_constants.f90
!> \brief   \copybrief mo_mhm_constants
!> \details \copydetails mo_mhm_constants

!> \brief Provides mHM specific constants
!> \details Provides mHM specific constants such as flood plain elevation.
! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting
!> \authors Matthias Cuntz
!> \date Nov 2011
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_mhm_constants

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  ! natural
  real(dp), public, parameter :: H2Odens = 1000.0_dp ! Density of water (kg/m3)

  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P2_InitStateFluxes = 15.00_dp
  real(dp), public, parameter :: P3_InitStateFluxes = 10.00_dp
  real(dp), public, parameter :: P4_InitStateFluxes = 75.00_dp
  real(dp), public, parameter :: P5_InitStateFluxes = 1500.00_dp
  real(dp), public, parameter :: C1_InitStateSM = 0.25_dp

  ! maximum number of outputs (fluxes states) for mHM
  integer(i4), public, parameter :: nOutFlxState = 21_i4     ! max. number of outputs to write into a netcdf file

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

  !> Neutrons and moisture: N0 formula, Desilets et al. 2010
  real(dp), public, parameter :: Desilets_a0 = 0.0808_dp
  real(dp), public, parameter :: Desilets_a1 = 0.372_dp
  real(dp), public, parameter :: Desilets_a2 = 0.115_dp

  !> Neutrons and moisture: COSMIC, Shuttleworth et al. 2013
  real(dp), public, parameter :: COSMIC_N = 348.33_dp           ! High energy neutron flux (cph), original was 510.51737902_dp
  real(dp), public, parameter :: COSMIC_alpha = 0.2392421548_dp ! Ratio of Fast Neutron Creation Factor (Soil to Water)
  real(dp), public, parameter :: COSMIC_L1 = 161.98621864_dp    ! High Energy Soil Attenuation Length (g/cm2)
  real(dp), public, parameter :: COSMIC_L2 = 129.14558985_dp    ! High Energy Water Attenuation Length (g/cm2)
  real(dp), public, parameter :: COSMIC_L3 = 107.82204562_dp    ! Fast Neutron Soil Attenuation Length (g/cm2)
  real(dp), public, parameter :: COSMIC_L4 = 3.1627190566_dp    ! Fast Neutron Water Attenuation Length (g/cm2)

END MODULE mo_mhm_constants
