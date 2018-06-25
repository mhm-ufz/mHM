!>       \file mo_common_file.f90

!>       \brief Provides file names and units for mRM

!>       \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.

!>       \authors Matthias Cuntz, Stephan Thober

!>       \date Aug 2015

! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting


MODULE mo_common_file

  IMPLICIT NONE
  !> DEM input data file
  CHARACTER(len=*), PARAMETER :: file_dem                = 'dem.asc'                     ! DEM
  !> Unit for  DEM input data file
  INTEGER,          PARAMETER :: udem                    = 53                            !
  !> Unit for  LCover input data file
  INTEGER,          PARAMETER :: ulcoverclass            = 61                            !

  !> file defining mHM's outputs
  CHARACTER(len=*), PARAMETER :: file_config             = 'ConfigFile.log'              ! configuration
  !> Unit for file defining mHM's outputs
  INTEGER,          PARAMETER :: uconfig                 = 68                            !



END MODULE mo_common_file
