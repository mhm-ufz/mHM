!> \file mo_common_constants.f90
!> \brief \copybrief mo_common_constants
!> \details \copydetails mo_common_constants

!> \brief Provides constants commonly used by mHM, mRM and MPR
!> \details Provides commonly used by mHM, mRM and MPR such as no_data values and eps
!> \changelog
!!  - Robert Schweppe Jun 2018
!!    - refactoring and reformatting
!> \authors Robert Schweppe
!> \date Dec 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_common_constants

  USE mo_kind, ONLY : i4, dp, sp

  IMPLICIT NONE

  PRIVATE

  ! Computational
  !> epsilon(1.0) in double precision
  REAL(dp), public, PARAMETER :: eps_dp = epsilon(1.0_dp)
  !> epsilon(1.0) in single precision
  REAL(sp), public, PARAMETER :: eps_sp = epsilon(1.0_sp)

  ! computational, these values need to be the same!!!
  integer(i4), public, parameter :: nodata_i4 = -9999_i4  ! [-]     global no data value
  real(dp), public, parameter :: nodata_dp = -9999._dp ! [-]     global no data value

  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P1_InitStateFluxes = 0.00_dp

  ! hydrologic modeling
  integer(i4), public, parameter :: nColPars = 5_i4      ! number of properties of the global variables
  integer(i4), public, parameter :: maxNoDomains = 50_i4     ! maximum number of allowed domains
  integer(i4), public, parameter :: maxNLcovers = 50_i4     ! maximum number of allowed LCover scenes

  character(64), public, parameter :: soilHorizonsVarName = "L1_SoilHorizons"
  character(64), public, parameter :: landCoverPeriodsVarName = "L1_LandCoverPeriods"
  character(64), public, parameter :: LAIVarName = "L1_LAITimesteps"


END MODULE mo_common_constants
