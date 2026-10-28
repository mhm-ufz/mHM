!> \file mo_mhm_constants.f90
!> \copydoc mo_mhm_constants

!> \brief Provides mHM specific constants
!> \details Provides mHM specific constants for hydrological states and neutron formulations.
! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting
!> \authors Matthias Cuntz
!> \date Nov 2011
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_mhm_constants

  USE mo_kind, ONLY : dp

  IMPLICIT NONE

  PRIVATE

  ! natural
  real(dp), public, parameter :: H2Odens = 1000.0_dp ! Density of water (kg/m3)

  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P2_InitStateFluxes = 15.00_dp
  real(dp), public, parameter :: P3_InitStateFluxes = 10.00_dp
  real(dp), public, parameter :: P4_InitStateFluxes = 75.00_dp

  !> Neutrons and moisture: N0 formula, Desilets et al. 2010
  real(dp), public, parameter :: Desilets_a0 = 0.0808_dp
  real(dp), public, parameter :: Desilets_a1 = 0.372_dp
  real(dp), public, parameter :: Desilets_a2 = 0.115_dp

  !> Neutrons and moisture: COSMIC, Shuttleworth et al. 2013
  real(dp), public, parameter :: COSMIC_alpha = 0.2392421548_dp ! Ratio of Fast Neutron Creation Factor (Soil to Water)
  real(dp), public, parameter :: COSMIC_L1 = 161.98621864_dp    ! High Energy Soil Attenuation Length (g/cm2)
  real(dp), public, parameter :: COSMIC_L2 = 129.14558985_dp    ! High Energy Water Attenuation Length (g/cm2)
  real(dp), public, parameter :: COSMIC_L4 = 3.1627190566_dp    ! Fast Neutron Water Attenuation Length (g/cm2)

END MODULE mo_mhm_constants
