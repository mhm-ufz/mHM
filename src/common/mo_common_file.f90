!>       \file mo_common_file.f90
!>       \brief provides file names and units for mrm
!>       \details provides all shared filenames of input and output of mHM and mRM

module mo_common_file

  implicit none
  !> variable name of lat lon for level 11 in netCDF file
  character(len = *), parameter :: varNameLatLon = 'latlon'
  !> variable name of digital elevation model in netCDF file
  character(len = *), parameter :: varNameDem = 'dem'
  !> variable name of land cover in netCDF file
  character(len = *), parameter :: varNameLandcover = 'land_cover'
  !> file defining mhm's outputs
  character(len = *), parameter :: file_config = 'ConfigFile.log'
  !> unit for file defining mhm's outputs
  integer, parameter :: uconfig = 68
  !> file defining optimization outputs (objective and parameter set)
  character(len = *), parameter :: file_opti = 'FinalParam.out'
  !> unit for file optimization outputs (objective and parameter set)
  integer, parameter :: uopti = 72                            !
  !> file defining optimization outputs in a namelist format (parameter set)
  character(len = *), parameter :: file_opti_nml = 'FinalParam.nml'              ! final parameters
  !> unit for file optimization outputs in a namelist format (parameter set)
  integer, parameter :: uopti_nml = 73                            !

end module mo_common_file
