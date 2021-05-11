!>       \file mo_common_file.f90

!>       \brief Provides file names and units for mRM

!>       \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.

!>       \authors Matthias Cuntz, Stephan Thober

!>       \date Aug 2015

! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting


MODULE mo_common_file

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: varNameDem                = 'dem'                         ! DEM
  CHARACTER(len=*), PARAMETER :: varNameLandCover          = 'land_cover'                  ! landcover

  !> file defining mHM's outputs
  CHARACTER(len=*), PARAMETER :: file_config             = 'ConfigFile.log'              ! configuration
  !> Unit for file defining mHM's outputs
  INTEGER,          PARAMETER :: uconfig                 = 68                            !



END MODULE mo_common_file
