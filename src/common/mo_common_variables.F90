!> \dir common
!> \brief \copybrief f_common
!> \details \copydetails f_common

!> \defgroup   f_common common - Fortran modules
!> \brief      Common modules used by mHM, mRM and MPR.
!> \details    This module provides different routines, constants and structures for all components of mHM.

!> \file mo_common_variables.f90
!> \brief \copybrief mo_common_variables
!> \details \copydetails mo_common_variables

!> \brief Provides structures needed by mHM, mRM and/or mpr.
!> \details Provides the global structure period that is used by both mHM and mRM.
!> \changelog
!! - Stephan Thober  Nov 2016
!!   - moved processdescription from mo_global_variables to here
!! - Robert Schweppe Dec 2017
!!   - merged more duplicated variables from mhm and mrm global variables
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!> \authors Stephan Thober
!> \date Sep 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_variables

  use mo_kind, only : i4, dp
  use mo_common_types, only: grid, gridremapper, domain_meta
#ifdef MPI
  USE mpi_f08
#endif
  implicit none

  integer(i4) :: itimer           ! Current timer number

  ! -------------------------------------------------------------------
  ! PROJECT DESCRIPTION for the NETCDF output file
  ! -------------------------------------------------------------------
  character(1024), public :: project_details            ! project including funding instituion., PI, etc.
  character(1024), public :: setup_description          ! any specific description of simulation
  character(1024), public :: simulation_type            ! e.g. seasonal forecast, climate projection, ...
  character(256), public :: Conventions                ! convention used for dataset
  character(1024), public :: contact                    ! contact details, incl. PI name
  character(1024), public :: mHM_details                ! developing institution, specific mHM revision
  character(1024), public :: history                    ! details on version/creation date

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of main part
  ! -------------------------------------------------------------------
  integer(i4), public :: iFlag_cordinate_sys        ! options model for the run cordinate system
  real(dp), dimension(:), allocatable, public :: resolutionHydrology        ! [m or degree] resolution of hydrology - Level 1
  integer(i4), dimension(:), allocatable, public :: L0_Domain
  logical, public :: write_restart              ! flag

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nDomains
  character(256), dimension(:), allocatable, public :: mhmFileRestartOut ! Directory where output of restart is written
  character(256), dimension(:), allocatable, public :: mrmFileRestartOut ! Directory where output of restart is written
  character(256), public :: dirConfigOut
  character(256), public :: dirCommonFiles ! directory where common input files should be located
  character(256), dimension(:), allocatable, public :: dirMorpho ! Directory where morphological files are located
  character(256), dimension(:), allocatable, public :: dirLCover ! Directory where land cover files are located
  character(256), dimension(:), allocatable, public :: dirOut ! Directory where output is written to
  character(256), dimension(:), allocatable, public :: fileLatLon ! Directory where the Lat Lon Files are located

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type(Grid), dimension(:), target, allocatable, public :: level0 ! grid information at morphological level (e.g., dem, fDir)
  type(Grid), dimension(:), target, allocatable, public :: level1 ! grid information at hydrologic level

  type(GridRemapper), dimension(:), allocatable, public :: l0_l1_remap  ! grid information at morphological level (e.g., dem, fDir)

  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain>
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  real(dp), public, dimension(:), allocatable :: L0_elev    ! [m]      Elevation (sinks removed)
  !          target variable for coupling to mRM
  integer(i4), public, dimension(:, :), allocatable :: L0_LCover      ! Classic mHM landcover class (upto 3 classes)
  !                                                                          ! dim1=number grid cells, dim2=Number of land cover scenes
  !                                                                          ! target variable for coupling to mRM

#ifdef MPI
  ! -------------------------------------------------------------------
  ! MPI variables
  type(MPI_Comm)      :: comm                ! MPI communicator
#endif
  ! -------------------------------------------------------------------
  !
  ! -------------------------------------------------------------------
  ! DOMAIN general description
  ! -------------------------------------------------------------------
  type(domain_meta), public :: domainMeta
  integer(i4), public :: nuniqueL0Domains ! Number of unique domains for L0

  ! -----------------------------------------------------------------
  ! LAND COVER DATA
  ! -----------------------------------------------------------------
  ! Land cover information
  integer(i4), public :: nLCoverScene        ! Number of land cover scene (lcs)
  character(256), dimension(:), allocatable, public :: LCfilename          ! file names for the different lcs
  integer(i4), dimension(:), allocatable, public :: LC_year_start       ! vector of start years for lcs
  integer(i4), dimension(:), allocatable, public :: LC_year_end         ! vector of end years for lcs

  ! -------------------------------------------------------------------
  ! PROCESSES description
  ! -------------------------------------------------------------------
  integer(i4), parameter, public :: nProcesses = 11 ! Number of possible processes to consider
  !                                                                !   process 1 :: interception
  !                                                                !   process 2 :: snow
  !                                                                !   process 3 :: soilmoisture
  !                                                                !   process 4 :: sealed area direct runoff
  !                                                                !   process 5 :: potential evapotranspiration
  !                                                                !   process 6 :: interflow
  !                                                                !   process 7 :: percolation
  !                                                                !   process 8 :: routing
  !                                                                !   process 9 :: baseflow
  !                                                                !   process 10:: neutrons
  !                                                                !   process 11:: river temperature routing
  integer(i4), dimension(nProcesses, 3), public :: processMatrix   ! Info about which process runs in which option and
  !                                                                ! number of parameters necessary for this option
  !                                                                !   col1: process_switch
  !                                                                !   col2: no. of parameters
  !                                                                !   col3: cum. no. of parameters

  ! -------------------------------------------------------------------
  ! PARAMETERS
  ! -------------------------------------------------------------------
  real(dp), dimension(:, :), allocatable, public, target :: global_parameters
  !                                                               ! Matrix of global parameters (former: gamma)
  !                                                               !   col1: min,  col2: max, col3: initial,
  !                                                               !   col4: flag, col5: scaling
  character(256), dimension(:), allocatable, public :: global_parameters_name
  !                                                               ! Matrix of global parameters (former: gamma)
  !                                                               !   col1: names
  ! -------------------------------------------------------------------
  ! ALMA convention
  ! -------------------------------------------------------------------
  ! TODO: this is currently used only be mRM, but could be useful for MPR and mHM also, ...
  ! so it is already in common_variables
  logical :: ALMA_convention ! flag for ALMA convention
  !                          ! see http://www.lmd.jussieu.fr/~polcher/ALMA/convention_3.html
  !                          ! .True.: ALMA convention is used for Input/Output
  !                          ! .False.: default mHM units are used
  !                          ! CAUTION: only Qall is considered at the moment

end module mo_common_variables
