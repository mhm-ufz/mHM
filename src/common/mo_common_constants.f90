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
  integer(i4), public, parameter :: maxNoBasins = 50_i4     ! maximum number of allowed basins
  integer(i4), public, parameter :: maxNLcovers = 50_i4     ! maximum number of allowed LCover scenes

  ! temporal
  real(dp), public, parameter :: DayHours = 24.0_dp  ! hours per day
  real(dp), public, parameter :: YearMonths = 12.0_dp  ! months per year
  integer(i4), public, parameter :: YearMonths_i4 = 12       ! months per year
  real(dp), public, parameter :: YearDays = 365.0_dp  ! days in a year
  real(dp), public, parameter :: DaySecs = 86400.0_dp  ! sec in a day
  real(dp), public, parameter :: HourSecs = 3600.0_dp  ! seconds per hour

END MODULE mo_common_constants
