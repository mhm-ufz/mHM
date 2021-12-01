!>       \file mo_common_constants.f90

!>       \brief Provides constants commonly used by mHM, mRM and MPR

!>       \details Provides commonly used by mHM, mRM and MPR such as no_data values and eps

!>       \authors Robert Schweppe

!>       \date Dec 2017

! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting

MODULE mo_common_constants

  USE mo_kind, ONLY : i4, dp, sp

  IMPLICIT NONE

  PRIVATE

  ! Computational
  !> epsilon(1.0) in double precision
  real(dp), public, parameter :: eps_dp = epsilon(1.0_dp)
  !> epsilon(1.0) in single precision
  real(sp), public, parameter :: eps_sp = epsilon(1.0_sp)
  !> the default precision used for comparing floats for equality (e.g. when read from netcdf files)
  real(dp), public, parameter :: floatComparisonPrecision = 1.e-6_dp
  !> the actual precision used for comparing floats for equality (e.g. when read from netcdf files)
  real(dp), public :: float_comparison_precision


  ! computational, these values need to be the same!!!
  integer(i4), public, parameter :: nodata_i4 = -9999_i4  ! [-]     global no data value
  real(dp), public, parameter :: nodata_dp = -9999._dp ! [-]     global no data value

  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P1_InitStateFluxes = 0.00_dp

  ! hydrologic modeling
  integer(i4), public, parameter :: nColPars = 5_i4      ! number of properties of the global variables
  integer(i4), public, parameter :: maxNoDomains = 50_i4     ! maximum number of allowed domains
  integer(i4), public, parameter :: maxNLcovers = 50_i4     ! maximum number of allowed LCover scenes
  integer(i4), public, parameter :: maxNLais = 366_i4     ! maximum number of allowed LAI periods
  integer(i4), public, parameter :: nLCover_class = 3_i4      ! [-]     Number of land cover class

  character(64), public, parameter :: soilHorizonsVarName = "L1_SoilHorizons"
  character(64), public, parameter :: landCoverPeriodsVarName = "L1_LandCoverPeriods"
  character(64), public, parameter :: LAIVarName = "L1_LAITimesteps"

  ! flags controlling the selection of input landcover and LAI periods based on simulation period
  logical, public, parameter :: keepUnneededPeriodsLandCover = .true.
  logical, public, parameter :: keepUnneededPeriodsLAI = .false.

END MODULE mo_common_constants
