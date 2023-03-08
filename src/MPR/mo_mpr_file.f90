!> \file mo_mpr_file.f90
!> \brief \copybrief mo_mpr_file
!> \details \copydetails mo_mpr_file

!> \brief Provides file names and units for mRM
!> \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.
!> \authors Matthias Cuntz, Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_mpr_file

  IMPLICIT NONE

  !> Current mHM model version
  CHARACTER(len = *), PARAMETER :: version = '0.1'                         ! Version
  !> Time of current mHM model version release
  CHARACTER(len = *), PARAMETER :: version_date = 'Jun 2019'                    ! Release date
  !> Driver file
  CHARACTER(len = *), PARAMETER :: file_main = 'mpr_driver.f90'              ! Driver
  !> Namelist file name
  CHARACTER(len = *), PARAMETER :: file_namelist_mpr = 'mpr.nml'                     ! Namelist
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mpr = 80                            ! set different from mhm
  !> Parameter namelists file name
  CHARACTER(len = *), PARAMETER :: file_namelist_mpr_param = 'mpr_parameter.nml'           ! Parameter namelists
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mpr_param = 31                            !


  !> Soil database file (iFlag_soilDB = 0) = classical mHM format
  CHARACTER(len = *), PARAMETER :: file_soil_database = 'soil_classdefinition.txt'    ! Soil data base
  !> Soil database file (iFlag_soilDB = 1)
  CHARACTER(len = *), PARAMETER :: file_soil_database_1 = 'soil_classdefinition_iFlag_soilDB_1.txt'
  !> Unit for soil data base
  INTEGER, PARAMETER :: usoil_database = 52                            !
  !> slope input data file
  CHARACTER(len = *), PARAMETER :: file_slope = 'slope.asc'                   ! slope
  !> Unit for  slope input data file
  INTEGER, PARAMETER :: uslope = 54                            !
  !> aspect input data file
  CHARACTER(len = *), PARAMETER :: file_aspect = 'aspect.asc'                  ! aspect
  !> Unit for  aspect input data file
  INTEGER, PARAMETER :: uaspect = 55                            !
  !> hydrogeological classes input data file
  CHARACTER(len = *), PARAMETER :: file_hydrogeoclass = 'geology_class.asc'           ! hydrogeological classes
  !> Unit for  hydrogeological classes input data file
  INTEGER, PARAMETER :: uhydrogeoclass = 58                            !
  !> soil classes input data file
  CHARACTER(len = *), PARAMETER :: file_soilclass = 'soil_class.asc'              ! soil classes
  !> Unit for  soil classes input data file
  INTEGER, PARAMETER :: usoilclass = 59                            !
  !> LAI classes input data file
  CHARACTER(len = *), PARAMETER :: file_laiclass = 'LAI_class.asc'               ! LAI classes
  !> Unit for  LAI input data file
  INTEGER, PARAMETER :: ulaiclass = 60                            !

  !> geological formation lookup table file
  CHARACTER(len = *), PARAMETER :: file_geolut = 'geology_classdefinition.txt' ! geolog. formation lookup table
  !> Unit for geological formation lookup table file
  INTEGER, PARAMETER :: ugeolut = 64                            !

  !> LAI classes lookup table file
  CHARACTER(len = *), PARAMETER :: file_lailut = 'LAI_classdefinition.txt'     ! LAI classes lookup table
  !> Unit for LAI classes lookup table file
  INTEGER, PARAMETER :: ulailut = 65                            !

END MODULE mo_mpr_file
