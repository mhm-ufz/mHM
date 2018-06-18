!>       \file mo_common_mHM_mRM_file.f90

!>       \brief Provides file names and units for mHM

!>       \details Provides all filenames as well as all units used for the hydrologic model mHM.

!>       \authors Matthias Cuntz

!>       \date Jan 2012

! Modifications:

MODULE mo_common_mhm_mrm_file

  IMPLICIT NONE

  !> file defining optimization outputs (objective and parameter set)
  CHARACTER(len = *), PARAMETER :: file_opti = 'FinalParam.out'              ! final parameters & objective
  !> Unit for file optimization outputs (objective and parameter set)
  INTEGER, PARAMETER :: uopti = 72                            !
  !> file defining optimization outputs in a namelist format (parameter set)
  CHARACTER(len = *), PARAMETER :: file_opti_nml = 'FinalParam.nml'              ! final parameters
  !> Unit for file optimization outputs in a namelist format (parameter set)
  INTEGER, PARAMETER :: uopti_nml = 73                            !

END MODULE mo_common_mhm_mrm_file
