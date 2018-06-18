!> \file mo_file.f90

!> \brief Provides file names and units for mHM

!> \details Provides all filenames as well as all units used for the hydrologic model mHM.

!> \author Matthias Cuntz
!> \date Jan 2012

MODULE mo_file

  IMPLICIT NONE

  !> Current mHM model version
  CHARACTER(len = *), PARAMETER :: version = '5.9'                         ! Version
  !> Time of current mHM model version release
  CHARACTER(len = *), PARAMETER :: version_date = 'June 2018'                   ! Release date
  !> Driver file
  CHARACTER(len = *), PARAMETER :: file_main = 'mhm_driver.f90'              ! Driver
  !> Namelist file name
  CHARACTER(len = *), PARAMETER :: file_namelist_mhm = 'mhm.nml'                     ! Namelist
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mhm = 30                            !
  !> Parameter namelists file name
  CHARACTER(len = *), PARAMETER :: file_namelist_mhm_param = 'mhm_parameter.nml'           ! Parameter namelists
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mhm_param = 31                            !

  !> file defining mHM's outputs
  CHARACTER(len = *), PARAMETER :: file_defOutput = 'mhm_outputs.nml'             ! output states and fluxes
  !> Unit for file defining mHM's outputs
  INTEGER, PARAMETER :: udefOutput = 67                            !

  !> unit for tws time series
  INTEGER, PARAMETER :: utws = 77                            !

END MODULE mo_file
