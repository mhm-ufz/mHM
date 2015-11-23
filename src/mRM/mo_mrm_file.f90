!> \file mo_mrm_file.f90

!> \brief Provides file names and units for mRM

!> \details Provides all filenames as well as all units used for the multiscale Routing Model mRM.

!> \author Matthias Cuntz, Stephan Thober
!> \date Aug 2015

MODULE mo_mrm_file

  IMPLICIT NONE
  
  !> Current mHM model version
  CHARACTER(len=*), PARAMETER :: version                 = '0.9'                         ! Version
  !> Time of current mHM model version release
  CHARACTER(len=*), PARAMETER :: version_date            = 'Aug 2015'                    ! Release date
  !> Driver file
  CHARACTER(len=*), PARAMETER :: file_main               = 'mrm_driver.f90'              ! Driver
  !> Namelist file name
  CHARACTER(len=*), PARAMETER :: file_namelist_mrm       = 'mhm.nml'                     ! Namelist
  !> Unit for namelist
  INTEGER,          PARAMETER :: unamelist_mrm           = 40                            ! set different from mhm
  !> Parameter namelists file name
  CHARACTER(len=*), PARAMETER :: file_namelist_param_mrm = 'mhm_parameter.nml'           ! Parameter namelists
  !> Unit for namelist
  INTEGER,          PARAMETER :: unamelist_param         = 41                            ! set different from mhm

  !> DEM input data file
  CHARACTER(len=*), PARAMETER :: file_dem                = 'dem.asc'                     ! DEM
  !> Unit for  DEM input data file
  INTEGER,          PARAMETER :: udem                    = 53                            ! 
  !> flow accumulation input data file
  CHARACTER(len=*), PARAMETER :: file_facc               = 'facc.asc'                    ! flow accumulation
  !> Unit for  flow accumulation input data file
  INTEGER,          PARAMETER :: ufacc                   = 56                            ! 
  !> flow direction input data file
  CHARACTER(len=*), PARAMETER :: file_fdir               = 'fdir.asc'                    ! flow direction
  !> Unit for  flow direction input data file
  INTEGER,          PARAMETER :: ufdir                   = 57                            ! 
  !> Unit for  LCover input data file
  INTEGER,          PARAMETER :: ulcoverclass            = 61                            ! 

  !> gauge location input data file
  CHARACTER(len=*), PARAMETER :: file_gaugeloc           = 'idgauges.asc'                ! gauge location
  !> Unit for  gauge location input data file
  INTEGER,          PARAMETER :: ugaugeloc               = 62                            ! 

  !> unit for discharge time series 
  INTEGER,          PARAMETER :: udischarge              = 66                            ! 

  !> file defining mRM's outputs
  CHARACTER(len=*), PARAMETER :: file_defOutput         = 'mrm_outputs.nml'             ! output states and fluxes
  !> Unit for file defining mRM's outputs
  INTEGER,          PARAMETER :: udefOutput             = 67                            ! 

  !> file defining mHM's outputs
  CHARACTER(len=*), PARAMETER :: file_config             = 'ConfigFile.log'              ! configuration
  !> Unit for file defining mHM's outputs
  INTEGER,          PARAMETER :: uconfig                 = 68                            ! 

  !> file defining optimization outputs (objective and p arameter set)
  CHARACTER(len=*), PARAMETER :: file_opti               = 'FinalParam.out'              ! final parameters & objective
  !> Unit for file optimization outputs (objective and p arameter set)
  INTEGER,          PARAMETER :: uopti                   = 72                            ! 
  !> file defining optimization outputs in a namelist fo rmat (parameter set)
  CHARACTER(len=*), PARAMETER :: file_opti_nml           = 'FinalParam.nml'              ! final parameters
  !> Unit for file optimization outputs in a namelist fo rmat (parameter set)
  INTEGER,          PARAMETER :: uopti_nml               = 73                            ! 

  !> file defining optimazation outputs
  CHARACTER(len=*), PARAMETER :: file_daily_discharge    = 'daily_discharge.out'         ! daily discharge file
  !> Unit for file optimazation outputs
  INTEGER,          PARAMETER :: udaily_discharge        = 74                            ! 
  !> file defining optimazation outputs
  CHARACTER(len=*), PARAMETER :: ncfile_discharge        = 'discharge.nc'                ! discharge file as netcdf

  !> file containing mrm output
  character(len=*), PARAMETER :: file_mrm_output         = 'mRM_Fluxes_States.nc'


END MODULE mo_mrm_file
