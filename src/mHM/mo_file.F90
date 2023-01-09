#ifndef MHMVERSION
#define MHMVERSION "0.0.0-dev0"
#endif
#ifndef MHMDATE
#define MHMDATE __DATE__
#endif
#define set_version(x) CHARACTER(len = *), PARAMETER :: version = x
#define set_date(x) CHARACTER(len = *), PARAMETER :: version_date = x

!> \file mo_file.F90
!> \brief \copybrief mo_file
!> \details \copydetails mo_file

!> \brief Provides file names and units for mHM
!> \details Provides all filenames as well as all units used for the hydrologic model mHM.
!! The \c version parameter will be set during compilation to
!! \"\htmlinclude version.txt \latexinclude version.txt\".
!! The \c version_date parameter will be set during compilation to
!! \"\htmlinclude version_date.txt \latexinclude version_date.txt\",
!! if it is a release version, otherwise it will be the current date.
!> \changelog
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!! - Sebastian Mueller Sep 2020
!!   - setting version with pre-processor from version file
!> \authors Matthias Cuntz
!> \authors Sebastian Mueller
!> \date Jan 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_file

  IMPLICIT NONE

  set_version(MHMVERSION)
  !< Current mHM model version (will be set to \htmlinclude version.txt \latexinclude version.txt)

  set_date(MHMDATE)
  !< Time of current mHM model version release

  !> Driver file
  CHARACTER(len = *), PARAMETER :: file_main = 'mhm_driver.f90'
  !> Namelist file name
  character(:), allocatable :: file_namelist_mhm ! = 'mhm.nml'
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mhm = 30
  !> Parameter namelists file name
  character(:), allocatable :: file_namelist_mhm_param ! = 'mhm_parameter.nml'
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mhm_param = 31
  !> file defining mHM's outputs
  character(:), allocatable :: file_defOutput ! = 'mhm_outputs.nml'
  !> Unit for file defining mHM's outputs
  INTEGER, PARAMETER :: udefOutput = 67
  !> file containing mhm output
  character(len = *), PARAMETER :: file_mhm_output = 'mHM_Fluxes_States.nc'

END MODULE mo_file
