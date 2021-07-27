!>       \file mo_mrm_file.f90

!>       \brief Provides file names and units for mRM

!>       \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.

!>       \authors Matthias Cuntz, Stephan Thober

!>       \date Aug 2015

! Modifications:

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

  !> DEM input data file
  CHARACTER(len=*), PARAMETER :: file_dem = 'dem.nc'                     ! DEM

  !> land cover input data file
  CHARACTER(len=*), PARAMETER :: file_lcover = 'land_cover.nc'              ! land_cover

  CHARACTER(len = *), PARAMETER :: file_facc = 'facc.nc'                    ! flow accumulation
  !> flow direction input data file
  CHARACTER(len = *), PARAMETER :: file_fdir = 'fdir.nc'                    ! flow direction
  !> flow direction input data file
  CHARACTER(len=*), PARAMETER :: file_slope = 'slope.nc'                    ! slope

  !> gauge location input data file
  CHARACTER(len = *), PARAMETER :: file_gaugeloc = 'idgauges.nc'                ! gauge location
  !> unit for discharge time series 
  INTEGER, PARAMETER :: udischarge = 66                            !
  !> file defining mRM's outputs
  CHARACTER(len = *), PARAMETER :: file_defOutput = 'mrm_outputs.nml'             ! output states and fluxes
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

END MODULE mo_mrm_file
