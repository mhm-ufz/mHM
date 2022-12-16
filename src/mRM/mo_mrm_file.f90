!> \file mo_mrm_file.f90
!> \brief \copybrief mo_mrm_file
!> \details \copydetails mo_mrm_file

!> \brief Provides file names and units for mRM
!> \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.
!> \authors Matthias Cuntz, Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
MODULE mo_mrm_file

  IMPLICIT NONE

  !> Current mHM model version
  CHARACTER(len = *), PARAMETER :: version = '1.0'                         ! Version
  !> Time of current mHM model version release
  CHARACTER(len = *), PARAMETER :: version_date = 'May 2019'                    ! Release date
  !> Driver file
  CHARACTER(len = *), PARAMETER :: file_main = 'mrm_driver.f90'              ! Driver
  !> Namelist file name
  CHARACTER(len = *), PARAMETER :: file_namelist_mrm = 'mrm.nml'                     ! Namelist
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_mrm = 40                            ! set different from mhm
  !> Parameter namelists file name
  CHARACTER(len = *), PARAMETER :: file_namelist_param_mrm = 'mrm_parameter.nml'           ! Parameter namelists
  !> Unit for namelist
  INTEGER, PARAMETER :: unamelist_param_mrm = 41                            ! set different from mhm

  CHARACTER(len = *), PARAMETER :: file_facc = 'facc.asc'                    ! flow accumulation
  !> Unit for  flow accumulation input data file
  INTEGER, PARAMETER :: ufacc = 56                            !
  !> flow direction input data file
  CHARACTER(len = *), PARAMETER :: file_fdir = 'fdir.asc'                    ! flow direction
  !> Unit for  flow direction input data file
  INTEGER,          PARAMETER :: ufdir = 57                            !
  !> flow direction input data file
  CHARACTER(len=*), PARAMETER :: file_slope              = 'slope.asc'                    ! slope
  !> Unit for  flow direction input data file
  INTEGER,          PARAMETER :: uslope                  = 59                            !

  !> gauge location input data file
  CHARACTER(len = *), PARAMETER :: file_gaugeloc = 'idgauges.asc'                ! gauge location
  !> Unit for  gauge location input data file
  INTEGER, PARAMETER :: ugaugeloc = 62                            !

  !> unit for discharge time series
  INTEGER, PARAMETER :: udischarge = 66                            !

  !> file defining mRM's outputs
  character(:), allocatable :: file_defOutput ! = 'mrm_outputs.nml'             ! output states and fluxes
  !> Unit for file defining mRM's outputs
  INTEGER, PARAMETER :: udefOutput = 67                            !

  !> file defining mHM's outputs
  CHARACTER(len = *), PARAMETER :: file_config = 'ConfigFile.log'              ! configuration
  !> Unit for file defining mHM's outputs
  INTEGER, PARAMETER :: uconfig = 68                            !

  !> file defining optimazation outputs
  CHARACTER(len = *), PARAMETER :: file_daily_discharge = 'daily_discharge.out'         ! daily discharge file
  !> Unit for file optimazation outputs
  INTEGER, PARAMETER :: udaily_discharge = 74                            !
  !> file defining optimazation outputs
  CHARACTER(len = *), PARAMETER :: ncfile_discharge = 'discharge.nc'                ! discharge file as netcdf

  !> file containing mrm output
  character(len = *), PARAMETER :: file_mrm_output = 'mRM_Fluxes_States.nc'

!> file containing mrm output for groundwater coupling
  character(len = *), PARAMETER :: file_gw_output = 'mRM_gw_Fluxes_States.nc'

  !> file defining optimazation outputs
  CHARACTER(len = *), PARAMETER :: file_subdaily_discharge = 'subdaily_discharge.out'         ! input_timestep discharge file
  !> Unit for file optimazation outputs
  INTEGER, PARAMETER :: usubdaily_discharge = 75                            !
  !> file containing simulated discharge at observat time step
  CHARACTER(len = *), PARAMETER :: ncfile_subdaily_discharge = 'subdaily_discharge.nc'    ! discharge file as netcdf

END MODULE mo_mrm_file
