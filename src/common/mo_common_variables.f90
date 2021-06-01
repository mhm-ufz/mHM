!>       \file mo_common_variables.f90
!>       \brief Provides structures needed by mHM, mRM and/or mpr.
!>       \details Provides the global variables needed by both mHM and mRM
!>       by both mHM and mRM.

module mo_common_variables

  use mo_kind, only : i4, i8, dp
  use mo_common_datetime_type, only : period
  use mo_grid, only : Grid, GridRemapper
#ifdef MPI
  USE mpi_f08
#endif
  implicit none

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
  real(dp), dimension(:), allocatable, public :: resolutionHydrology        ! [m or degree] resolution of hydrology - Level 1
  integer(i4), dimension(:), allocatable, public :: L0_Domain
  logical, public :: write_restart              ! flag

   integer(i4) :: mrm_coupling_mode !-1 = no mrm (mHM only)
  !                                ! 0 = stand alone (mRM only)
  !                                ! 1 = general coupling to a model (not used)
  !                                ! 2 = specific coupling to mHM (mHM and mRM)

  real(dp), public :: c2TSTu            !       Unit transformation = timeStep/24
  real(dp), dimension(:), allocatable, public :: resolutionRouting          ! [m or degree] resolution of routing - Level 11
  logical, public :: read_restart               ! flag
  logical, public :: mrm_read_river_network               ! flag

  type(period), dimension(:), allocatable, public :: warmPer     ! time period for warming
  type(period), dimension(:), allocatable, public :: evalPer     ! time period for model evaluation
  type(period), public :: readPer     ! start and end dates of read period
  integer(i4), dimension(:), allocatable, public :: warmingDays ! number of days for warm up period

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nDomains
  character(256), dimension(:), allocatable, public :: mhmFileRestartIn! Directory where input of restart is read from
  character(256), dimension(:), allocatable, public :: mrmFileRestartIn! Directory where input of restart is read from
  character(256), dimension(:), allocatable, public :: mhmFileRestartOut ! Directory where output of restart is written
  character(256), dimension(:), allocatable, public :: mrmFileRestartOut ! Directory where output of restart is written
  character(256), public :: dirConfigOut
  character(256), public :: dirCommonFiles ! directory where common input files should be located
  character(256), dimension(:), allocatable, public :: dirOut ! Directory where output is written to
  !TODO: MPR this will go, level0 is moved to mRM
  character(256), dimension(:), allocatable, public :: dirMorpho ! Directory where morphological files are located
  character(256), dimension(:), allocatable, public :: dirLCover ! Directory where land cover files are located
  character(256), dimension(:), allocatable, public :: fileLatLon ! Directory where the Lat Lon Files are located
  type(Grid), dimension(:), target, allocatable, public :: level0 ! grid information at morphological level (e.g., dem, fDir)
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

  type(Grid), dimension(:), target, allocatable, public :: level1 ! grid information at hydrologic level



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

  !TODO: MPR this will go
  ! -----------------------------------------------------------------
  ! LAND COVER DATA
  ! -----------------------------------------------------------------
  ! Land cover information
  integer(i4), public :: nLandCoverPeriods        ! Number of land cover scene (lcs)
  character(256), dimension(:), allocatable, public :: LCfilename          ! file names for the different lcs
  integer(i4), dimension(:), allocatable, public :: LC_year_start       ! vector of start years for lcs
  integer(i4), dimension(:), allocatable, public :: LC_year_end         ! vector of end years for lcs
  integer(i4), dimension(:, :), allocatable, public :: LCyearId  ! Mapping of landcover scenes (1, 2,..) for each basin
                                                                 ! dimensions ( year, iBasin) (year as integer!!!)

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
  character(64), dimension(:), allocatable, public :: global_parameters_name
  !                                                               ! Matrix of global parameters (former: gamma)
  !                                                               !   col1: names
  real(dp), dimension(:), allocatable, public:: dummy_global_parameters
  !                                                               ! Vector of global dummy parameter default values
  character(64), dimension(:), allocatable, public :: dummy_global_parameters_name
  !                                                               ! Vector of global dummy parameter names
  ! -------------------------------------------------------------------
  ! OPTIMIZATION
  ! -------------------------------------------------------------------
  integer(i4), public :: opti_method         ! Optimization algorithm:
  !                                                                       ! 1 - DDS
  !                                                                       ! 2 - Simulated Annealing
  !                                                                       ! 3 - SCE
  integer(i4), public :: opti_function       ! Objective function:
  !                                                                       ! 1 - 1.0-NSE
  !                                                                       ! 2 - 1.0-lnNSE
  !                                                                       ! 3 - 1.0-0.5*(NSE+lnNSE)
  logical, public :: optimize            ! Optimization   (.true. ) or
  !                                                                       ! Evaluation run (.false.)
  logical, public :: optimize_restart    ! Optimization will be restarted from
  !                                                                       ! mo_<opti_method>.restart file (.true.)
  ! settings for optimization algorithms:
  integer(i8), public :: seed                ! seed used for optimization
  !                                                                       ! default: -9 --> system time
  integer(i4), public :: nIterations         ! number of iterations for optimization
  real(dp), public :: dds_r               ! DDS: perturbation rate
  !                                                                       !      default: 0.2
  real(dp), public :: sa_temp             ! SA:  initial temperature
  !                                                                       !      default: -9.0 --> estimated
  integer(i4), public :: sce_ngs             ! SCE: # of complexes
  !                                                                       !      default: 2
  integer(i4), public :: sce_npg             ! SCE: # of points per complex
  !                                                                       !      default: -9 --> 2n+1
  integer(i4), public :: sce_nps             ! SCE: # of points per subcomplex
  !                                                                       !      default: -9 --> n+1
  logical, public :: mcmc_opti           ! MCMC: Optimization (.true. ) or
  !                                                                       !       Only parameter uncertainty (.false.)
  integer(i4), public, parameter :: nerror_model = 2    ! # possible parameters in error model
  !                                                                       !       e.g. for opti_function=8: 2
  real(dp), public, dimension(nerror_model) :: mcmc_error_params   !       Parameters of error model if mcmc_opti=.false.
  !                                                                       !       e.g. for opti_function=8: 0.01, 0.3

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
