!> \file mo_common_constants.f90
!> \copydoc mo_common_constants

!> \brief Shared run defaults and constants used by mHM and MPR.
!> \details Defines the default run paths, field reset value, and restart coordinate names.
!> \changelog
!!  - Robert Schweppe Jun 2018
!!    - refactoring and reformatting
!!  - Oct 2026: removed unused legacy constants and centralized run defaults.
!> \authors Robert Schweppe
!> \date Dec 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_constants

  use mo_kind, only: i4, dp

  implicit none

  private

  character(*), parameter, public :: default_main_file = "mhm.nml" !< Default main namelist path.
  character(*), parameter, public :: default_para_file = "parameters.nml" !< Default parameter namelist path.
  character(*), parameter, public :: default_out_file = "outputs.nml" !< Default output namelist path.
  character(*), parameter, public :: default_cwd = "." !< Default run root: the process working directory.

  integer(i4), parameter, public :: nProcesses = 11 !< Number of possible processes to consider

  real(dp), parameter, public :: P1_InitStateFluxes = 0.00_dp !< Default reset value for mHM states and fluxes.

  character(64), parameter, public :: soilHorizonsVarName = "L1_SoilHorizons" !< Soil-horizon restart coordinate name.
  character(64), parameter, public :: landCoverPeriodsVarName = "L1_LandCoverPeriods" !< Land-cover restart coordinate name.
  character(64), parameter, public :: LAIVarName = "L1_LAITimesteps" !< LAI restart coordinate name.

end module mo_common_constants
