!>       \file mo_common_variables.f90

!>       \brief Provides structures needed by mHM, mRM and/or mpr.

!>       \details Provides the global structure period that is used
!>       by both mHM and mRM.

!>       \authors Stephan Thober

!>       \date Sep 2015

! Modifications:
! Stephan Thober  Nov 2016 - moved processdescription from mo_global_variables to here
! Robert Schweppe Dec 2017 - merged more duplicated variables from mhm and mrm global variables
! Robert Schweppe Jun 2018 - refactoring and reformatting


module mo_common_variables

  use mo_kind, only : i4, i8, dp
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
  ! PERIOD description
  ! -------------------------------------------------------------------
  type period
    integer(i4) :: dStart      ! first day
    integer(i4) :: mStart      ! first month
    integer(i4) :: yStart      ! first year
    integer(i4) :: dEnd        ! last  day
    integer(i4) :: mEnd        ! last  month
    integer(i4) :: yEnd        ! last  year
    integer(i4) :: julStart    ! first julian day
    integer(i4) :: julEnd      ! last  julian day
    integer(i4) :: nObs        ! total number of observations
  end type period

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type Grid
    ! general domain information
    integer(i4) :: ncols     ! Number of columns
    integer(i4) :: nrows     ! Number of rows
    integer(i4) :: nCells    ! Number of cells in mask
    real(dp) :: xllcorner    ! x coordinate of the lowerleft corner
    real(dp) :: yllcorner    ! y coordinate of the lowerleft corner
    real(dp) :: cellsize     ! Cellsize x = cellsize y
    real(dp) :: nodata_value ! Code to define the mask
    real(dp), dimension(:, :), allocatable :: x  ! 2d longitude array (unmasked version is needed for output anyway)
    real(dp), dimension(:, :), allocatable :: y  ! 2d latitude  array (unmasked version is needed for output anyway)
    logical, dimension(:, :), allocatable :: mask  ! the mask for valid cells in the original grid (nrows*ncols)
    ! for referencing values in the nValidCells vector
    integer(i4) :: iStart          ! Starting cell index of a given domain
    integer(i4) :: iEnd            ! Ending cell index of a given domain
    ! dimension(nCells, (x,y) )
    integer(i4), dimension(:, :), allocatable :: CellCoor  ! this is only used for mRM
    real(dp), dimension(:), allocatable :: CellArea  ! area of the cell in sq m
    integer(i4), dimension(:), allocatable :: Id

  end type Grid

  type(Grid), dimension(:), target, allocatable, public :: level0 ! grid information at morphological level (e.g., dem, fDir)
  type(Grid), dimension(:), target, allocatable, public :: level1 ! grid information at hydrologic level

  type GridRemapper
    type(Grid), pointer :: high_res_grid
    type(Grid), pointer :: low_res_grid

    ! dimension nCells
    integer(i4), dimension(:), allocatable :: lower_bound  ! 1d index of lower side subgrid
    integer(i4), dimension(:), allocatable :: upper_bound  ! 1d index of upper side subgrid
    integer(i4), dimension(:), allocatable :: left_bound  ! 1d index of left side subgrid
    integer(i4), dimension(:), allocatable :: right_bound  ! 1d index of right side subgrid
    integer(i4), dimension(:), allocatable :: n_subcells   ! 1d numberof valid subgrid cells
    integer(i4), dimension(:, :), allocatable :: lowres_id_on_highres   ! 2d index array of lowres id

  end type GridRemapper

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
  type domain_meta
    integer(i4)                            :: nDomains
    integer(i4)                            :: overallNumberOfDomains  ! Number of domains for multi-domain optimization
    integer(i4), dimension(:), allocatable :: indices
    integer(i4), dimension(:), allocatable :: L0DataFrom
    ! optidata saves for each domain which optional data is assigned to it
    ! (0) default: the program decides. If you are confused, choose 0
    ! (1) runoff
    ! (2) sm
    ! (3) tws
    ! (4) neutons
    ! (5) et
    ! (6) et & tws
    integer(i4), dimension(:), allocatable :: optidata
    logical,     dimension(:), allocatable :: doRouting
#ifdef MPI
    logical                                :: isMasterInComLocal  ! true if the process is master proc in comLocal
    type(MPI_Comm)                         :: comMaster ! the communicater the domains are using to send messages to each other
                                                        ! here are all processes wich have rank 0 in comLocal
    type(MPI_Comm)                         :: comLocal  ! the communicater the domain internal communication takes place
#endif
  end type domain_meta

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
