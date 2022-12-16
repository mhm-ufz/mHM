!> \file mo_mpr_constants.f90
!> \brief \copybrief mo_mpr_constants
!> \details \copydetails mo_mpr_constants

!> \brief Provides MPR specific constants
!> \details Provides MPR specific constants such as flood plain elevation.
!> \authors Matthias Cuntz
!> \date Nov 2011
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_mpr_constants

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  ! hydrologic modeling
  integer(i4), public, parameter :: nLCover_class = 3_i4      ! [-]     Number of land cover class
  integer(i4), public, parameter :: maxGeoUnit = 25_i4     ! maximum number of allowed geological classes
  integer(i4), public, parameter :: maxNoSoilHorizons = 10_i4     ! maximum number of allowed soil layers

  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P2_InitStateFluxes = 15.00_dp
  real(dp), public, parameter :: P3_InitStateFluxes = 10.00_dp
  real(dp), public, parameter :: P4_InitStateFluxes = 75.00_dp
  real(dp), public, parameter :: P5_InitStateFluxes = 1500.00_dp
  real(dp), public, parameter :: C1_InitStateSM = 0.25_dp

  ! soil paramterization (mo_mpr_soilmoist)
  ! organic matter constant for calculation of mineral bulk density following RAWL
  real(dp), public, parameter :: BulkDens_OrgMatter = 0.224_dp     ! [g/cm3] from W.R. RAWLS
  ! constants for determinination of the field capacity following Twarakavi
  real(dp), public, parameter :: field_cap_c1 = -0.60_dp     ! field capacity constant 1
  real(dp), public, parameter :: field_cap_c2 = 2.0_dp      ! field capacity constant 2
  ! constants for determinination of the van Genuchten parameter n and sand treshold
  real(dp), public, parameter :: vGenuchten_sandtresh = 66.5_dp     ! van Genuchten snad treshold
  real(dp), public, parameter :: vGenuchtenN_c1 = 1.392_dp   ! constants for van Genuchten n
  real(dp), public, parameter :: vGenuchtenN_c2 = 0.418_dp
  real(dp), public, parameter :: vGenuchtenN_c3 = -0.024_dp
  real(dp), public, parameter :: vGenuchtenN_c4 = 1.212_dp
  real(dp), public, parameter :: vGenuchtenN_c5 = -0.704_dp
  real(dp), public, parameter :: vGenuchtenN_c6 = -0.648_dp
  real(dp), public, parameter :: vGenuchtenN_c7 = 0.023_dp
  real(dp), public, parameter :: vGenuchtenN_c8 = 0.044_dp
  real(dp), public, parameter :: vGenuchtenN_c9 = 3.168_dp
  real(dp), public, parameter :: vGenuchtenN_c10 = -2.562_dp
  real(dp), public, parameter :: vGenuchtenN_c11 = 7.0E-9_dp
  real(dp), public, parameter :: vGenuchtenN_c12 = 4.004_dp
  real(dp), public, parameter :: vGenuchtenN_c13 = 3.750_dp
  real(dp), public, parameter :: vGenuchtenN_c14 = -0.016_dp
  real(dp), public, parameter :: vGenuchtenN_c15 = -4.197_dp
  real(dp), public, parameter :: vGenuchtenN_c16 = 0.013_dp
  real(dp), public, parameter :: vGenuchtenN_c17 = 0.076_dp
  real(dp), public, parameter :: vGenuchtenN_c18 = 0.276_dp
  ! determinination Ks
  real(dp), public, parameter :: Ks_c = 10.0_dp
  ! permanent wiltung point (PWP)
  real(dp), public, parameter :: PWP_c = 1.0_dp
  real(dp), public, parameter :: PWP_matPot_ThetaR = 15000.0_dp ! [hPa] matrix potential of -1500 kPa, assumed as thetaR=0

  !> assumed meteorol. measurement hight for estimation of aeroResist and surfResist
  real(dp), public, parameter :: WindMeasHeight = 10.0_dp
  !> von karman constant
  real(dp), public, parameter :: karman = 0.41_dp

  !> LAI factor for bulk surface resistance formulation
  real(dp), public, parameter :: LAI_factor_surfResi = 0.3_dp
  !> LAI offset for bulk surface resistance formulation
  real(dp), public, parameter :: LAI_offset_surfResi = 1.2_dp
  !> maximum bulk surface resistance
  real(dp), public, parameter :: max_surfResist = 250.0_dp

END MODULE mo_mpr_constants
