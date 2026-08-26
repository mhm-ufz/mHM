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

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  integer(i4), public, parameter :: nerror_model = 2    !< # possible parameters in error model

  ! PROCESSES description
  integer(i4), parameter, public :: nProcesses = 11 !< Number of possible processes to consider

  !> default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P1_InitStateFluxes = 0.00_dp

  ! hydrologic modeling
  integer(i4), public, parameter :: nColPars = 5_i4      !< number of properties of the global variables
  integer(i4), public, parameter :: maxNoDomains = 50_i4     !< maximum number of allowed domains
  integer(i4), public, parameter :: maxNLcovers = 50_i4     !< maximum number of allowed LCover scenes

  character(64), public, parameter :: soilHorizonsVarName = "L1_SoilHorizons"
  character(64), public, parameter :: landCoverPeriodsVarName = "L1_LandCoverPeriods"
  character(64), public, parameter :: LAIVarName = "L1_LAITimesteps"

  ! default tolerance for header or value comparisons
  real(dp), public, parameter:: defaultTolerance_dp = 1.0e-7_dp

END MODULE mo_common_constants
