!> \file mo_global_variables.f90

!> \brief Global variables ONLY used in reading, writing and startup.

!> \details 

!> \authors Luis Samaniego
!> \date Dec 2012

MODULE mo_global_variables

  ! This module provides 

  ! 
  ! Written   Luis Samaniego, Dec 2005
  ! Modified  Luis Samaniego, Feb 2013  - new variable names, new modules, units
  !           Rohini Kumar,   Jul 2013  - fraction of perfectly sealed area within city added
  !           Rohini Kumar,   Aug 2013  - name changed "inputFormat" to "inputFormat_meteo_forcings"
  !           Rohini Kumar,   Aug 2013  - name changed from "L0_LAI" to "L0_LCover_LAI"
  !           Rohini Kumar,   Aug 2013  - added dirSoil_LUT and dirGeology_LUT
  !           Luis Samaniego, Nov 2013  - documentation of dimensions
  !           Matthias Zink,  Nov 2013  - added "InflowGauge" and inflow gauge variabels in basin 
  !           Rohini Kumar,   May 2014  - added options for the model run cordinate system
  !           Stephan Thober, Jun 2014  - added timeStep_model_inputs and readPer
  !           Stephan Thober, Jun 2014  - added perform_mpr, updated restart flags
  !           Matthias Cuntz & Juliane Mai, Nov 2014 - LAI input from daily, monthly or yearly files

  USE mo_kind,          ONLY: i4, i8, dp
  USE mo_mhm_constants, ONLY: nOutFlxState, YearMonths, maxNoBasins

  IMPLICIT NONE

  ! Types have to be public
  PUBLIC :: period
  PUBLIC :: gaugingStation
  PUBLIC :: soilType
  PUBLIC :: gridGeoRef
  PUBLIC :: basinInfo

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mHM
  ! -------------------------------------------------------------------
  integer(i4),   public                              :: timeStep                   ! [h] simulation time step (= TS) in [h]
  integer(i4),   public                              :: timeStep_model_inputs      ! frequency for reading meteo input
  real(dp),      dimension(:), allocatable, public   :: resolutionHydrology        ! [m or °] resolution of hydrology - Level 1
  real(dp),      dimension(:), allocatable, public   :: resolutionRouting          ! [m or °] resolution of routing - Level 11
  integer(i4),   dimension(:), allocatable, public   :: L0_Basin
  logical,       public                              :: read_restart               ! flag 
  logical,       public                              :: write_restart              ! flag 
  logical,       public                              :: perform_mpr                ! switch for performing
                                                                                   ! multiscale parameter regionalization
  character(256),public                              :: inputFormat_meteo_forcings ! format of meteo input data(bin or nc)
  ! LAI information
  character(256), public                             :: inputFormat_gridded_LAI    ! format of gridded LAI data(bin or nc)
  integer(i4),    public                             :: timeStep_LAI_input         ! time step of gridded LAI input
  integer(i4),    public                             :: iFlag_cordinate_sys        ! options model for the run cordinate system
  ! -------------------------------------------------------------------
  ! OPTIMIZATION
  ! -------------------------------------------------------------------
  integer(i4), public                              :: opti_method         ! Optimization algorithm:
  !                                                                       ! 1 - DDS
  !                                                                       ! 2 - Simulated Annealing
  !                                                                       ! 3 - SCE
  integer(i4), public                              :: opti_function       ! Objective function:
  !                                                                       ! 1 - 1.0-NSE
  !                                                                       ! 2 - 1.0-lnNSE
  !                                                                       ! 3 - 1.0-0.5*(NSE+lnNSE)
  logical,     public                              :: optimize            ! Optimization   (.true. ) or
  !                                                                       ! Evaluation run (.false.)
  ! settings for optimization algorithms: 
  integer(i8), public                              :: seed                ! seed used for optimization
  !                                                                       ! default: -9 --> system time 
  integer(i4), public                              :: nIterations         ! number of iterations for optimization
  real(dp),    public                              :: dds_r               ! DDS: perturbation rate
  !                                                                       !      default: 0.2
  real(dp),    public                              :: sa_temp             ! SA:  initial temperature
  !                                                                       !      default: -9.0 --> estimated
  integer(i4), public                              :: sce_ngs             ! SCE: # of complexes
  !                                                                       !      default: 2
  integer(i4), public                              :: sce_npg             ! SCE: # of points per complex
  !                                                                       !      default: -9 --> 2n+1
  integer(i4), public                              :: sce_nps             ! SCE: # of points per subcomplex
  !                                                                       !      default: -9 --> n+1

  ! -------------------------------------------------------------------
  ! PROCESSES description
  ! -------------------------------------------------------------------
  integer(i4), parameter,                      public :: nProcesses = 9         ! Number of possible processes to consider
  !                                                                             !   process 1 :: interception
  !                                                                             !   process 2 :: snow
  !                                                                             !   process 3 :: soilmoisture
  !                                                                             !   process 4 :: sealed area direct runoff
  !                                                                             !   process 5 :: potential evapotranspiration
  !                                                                             !   process 6 :: interflow
  !                                                                             !   process 7 :: percolation
  !                                                                             !   process 8 :: routing
  !                                                                             !   process 9 :: baseflow  
  integer(i4),    dimension(nProcesses, 3),    public :: processMatrix          ! Info about which process runs in which option and
  !                                                                             ! number of parameters necessary for this option
  !                                                                             !   col1: process_switch 
  !                                                                             !   col2: no. of parameters
  !                                                                             !   col3: cum. no. of parameters
  real(dp),       dimension(:,:), allocatable, public :: global_parameters      ! Matrix of global parameters (former: gamma)
  !                                                                             !   col1: min,  col2: max, col3: initial, 
  !                                                                             !   col4: flag, col5: scaling
  character(256), dimension(:), allocatable,   public :: global_parameters_name ! Matrix of global parameters (former: gamma)
  !                                                                             !   col1: names

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nBasins
  character(256), dimension(:), allocatable, public :: dirMorpho          ! Directory where morphological files are located
  character(256), dimension(:), allocatable, public :: dirLCover          ! Directory where land cover files are located
  character(256), dimension(:), allocatable, public :: dirGauges          ! Directory where discharge files are located
  character(256), dimension(:), allocatable, public :: dirPrecipitation   ! Directory where precipitation files are located
  character(256), dimension(:), allocatable, public :: dirTemperature     ! Directory where temperature files are located
  character(256), dimension(:), allocatable, public :: dirMinTemperature  ! Directory where minimum temp. files are located
  character(256), dimension(:), allocatable, public :: dirMaxTemperature  ! Directory where maximum temp. files are located
  character(256), dimension(:), allocatable, public :: dirNetRadiation    ! Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirabsVapPressure  ! Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirwindspeed       ! Directory where windspeed files are located
  character(256), dimension(:), allocatable, public :: dirReferenceET     ! Directory where reference-ET files are located
  character(256), dimension(:), allocatable, public :: dirOut             ! Directory where output is written to
  character(256), dimension(:), allocatable, public :: dirRestartOut      ! Directory where output of restart is written to
  character(256), dimension(:), allocatable, public :: dirRestartIn       ! Directory where input of restart is read from
  character(256), dimension(:), allocatable, public :: dirgridded_LAI     ! directory where gridded LAI is located
                                                                          ! used when timeStep_LAI_input < 0
  character(256), dimension(:), allocatable, public :: dirLatLon          ! directory to lat lon files

  ! directory common to all basins 
  character(256),                            public :: dirConfigOut       ! Directory where config run output is written to
  character(256),                            public :: dirCommonFiles     ! directory where common input files should be located
  !                                                                       ! for all modeled basins
  ! ------------------------------------------------------------------
  ! NETCDF OUTPUT
  ! ------------------------------------------------------------------
  real(dp), dimension(:),   allocatable, target, public :: yCoor               ! GK4 (DHDN3-zone 4) easting
  real(dp), dimension(:),   allocatable, target, public :: xCoor               ! GK4 (DHDN3-zone 4) northing
  real(dp), dimension(:,:), allocatable, target, public :: lons                ! WGS84 lons 
  real(dp), dimension(:,:), allocatable, target, public :: lats                ! WGS84 lats

  ! ------------------------------------------------------------------
  ! CONSTANT 
  ! ------------------------------------------------------------------
  real(dp), public                                  :: c2TSTu             !       Unit transformation = timeStep/24 
  integer(i4), public                               :: nTstepDay          !       Number of time intervals per day
  !                                                                       !       (was previously NAGG)
  integer(i4), public, parameter                    :: routingStates = 2  ! [-]   Routing states (2=current, 1=past)

  ! ------------------------------------------------------------------
  ! SOIL DATA
  ! ------------------------------------------------------------------
  real(dp), public                                  :: tillageDepth       ! [mm]  Soil depth down to which organic
  !                                                                               matter is possible
  integer(i4), public                               :: nSoilHorizons_mHM  !       Number of horizons to model
  real(dp), dimension(:), allocatable, public       :: HorizonDepth_mHM   ! [mm]  Horizon depth from surface,
  !                                                                               positive downwards
  integer(i4), public                               :: nSoilTypes         !       Number of soil types

  type soilType
     ! dim1 =  nSoilType (e.g. i=1..72 for BUEK)
     ! dim2 =  the maximum of nHorizons
     ! dim3 =  land cover classes
     ! input data
     integer(i4), dimension(:), allocatable     :: id                 !            Soil Id
     integer(i4), dimension(:), allocatable     :: nHorizons          !            Number of horizons
     integer(i4), dimension(:), allocatable     :: is_present         !            Wether this soil type is present in 
     !                                                                !            this basin or not
     real(dp), dimension(:,:),   allocatable    :: UD                 ! [mm]       Upper Bound of depth      
     real(dp), dimension(:,:),   allocatable    :: LD                 ! [mm]       Lower Bound of depth       
     real(dp), dimension(:,:),   allocatable    :: clay               ! [%]        Clay content                  
     real(dp), dimension(:,:),   allocatable    :: sand               ! [%]        Sand content                 
     real(dp), dimension(:,:),   allocatable    :: DbM                ! [g/cm2]    Mineral Bulk density
     real(dp), dimension(:,:),   allocatable    :: depth              ! [mm]       Depth of the soil Horizon    
     real(dp), dimension(:),     allocatable    :: RZdepth            ! [mm]       Total soil depth       
     real(dp), dimension(:,:,:), allocatable    :: Wd                 ! [1]        Weights of mHM Horizons according to 
     !                                                                !            horizons provided in soil database  
     integer(i4), dimension(:),  allocatable    :: nTillHorizons      ! [1]        Number of tillage horizons 

     ! derived soil hydraulic properties 
     real(dp), dimension(:,:,:), allocatable    :: thetaS_Till        ! [1]        Saturated water content of soil horizons
     !                                                                !            tillage depth - f(OM, management)
     real(dp), dimension(:,:),   allocatable    :: thetaS             ! [1]        Saturated water content of soil horizons 
     !                                                                !            after tillage depth
     real(dp), dimension(:,:,:), allocatable    :: Db                 ! [g/cm2]    Bulk density, LUC dependent
     !                                                                !            = f( OM, management)
     real(dp), dimension(:,:,:), allocatable    :: thetaFC_Till       ! [1]        Field capacity of tillage layers; 
     !                                                                !            LUC dependent - f(OM, management)
     real(dp), dimension(:,:),   allocatable    :: thetaFC            ! [1]        Field capacity of deeper layers
     real(dp), dimension(:,:,:), allocatable    :: thetaPW_Till       ! [1]        Permament wilting point of tillage layers; 
     !                                                                !            LUC dependent - f(OM, management)
     real(dp), dimension(:,:),   allocatable    :: thetaPW            ! [1]        Permanent wilting point of deeper layers
     real(dp), dimension(:,:,:), allocatable    :: Ks                 ! [cm/d]     Saturated hydaulic conductivity
  end type soilType
  type(soilType), public                        :: soilDB             !            The soil database

  ! -----------------------------------------------------------------
  ! GAUGED station data
  ! -----------------------------------------------------------------
  integer(i4), public                            :: nGaugesTotal       ! Number of evaluation gauges for all basins 
  integer(i4), public                            :: nInflowGaugesTotal ! Number of evaluation gauges for all basins 

  integer(i4), public                            :: nMeasPerDay        ! Number of observations per day,
  !                                                                    ! e.g. 24 -> hourly discharge, 1 -> daily discharge
  type gaugingStation
     integer(i4),    dimension(:),   allocatable :: basinId            ! Basin Id
     integer(i4),    dimension(:),   allocatable :: gaugeId            ! Gauge Id (e.g. 0000444)
     character(256), dimension(:),   allocatable :: fname              ! Name runoff file
     real(dp),       dimension(:,:), allocatable :: Q                  ! [m3 s-1] observed daily mean discharge (simPer)
     !                                                                 ! dim1=number observations, dim2=number of gauges
  end type gaugingStation
  type(gaugingStation), public                   :: gauge              ! Gauging station information
  type(gaugingStation), public                   :: InflowGauge        ! inflow gauge information
  
  ! -----------------------------------------------------------------
  ! GEOLOGICAL FORMATION data
  ! -----------------------------------------------------------------
  integer(i4), public                              :: nGeoUnits   ! Number of geological formations
  integer(i4), dimension(:), allocatable, public   :: GeoUnitList ! List of ids of each geological formations
  integer(i4), dimension(:), allocatable, public   :: GeoUnitKar  ! Id of Karstic formation (0 == does not exist)

  ! -----------------------------------------------------------------
  ! Land cover, LAI LUT data
  ! -----------------------------------------------------------------
  ! Land cover information
  real(dp), public                                    :: fracSealed_cityArea ! fraction of area within city assumed to be
                                                                             ! perfectly sealed [0-1] 
  integer(i4),    public                              :: nLCover_scene       ! Number of land cover scene
  character(256), public, dimension(:), allocatable   :: LCfilename          ! file names for the different land cover scenes
  integer(i4),    public, dimension(:), allocatable   :: LCyearId            ! Mapping of landcover scenes (1, 2, ...)
                                                                             ! to the actual year(1960, 1961, ...)
  ! LAI data
  ! variables used when timeStep_LAI_input == 0
  integer(i4),    public                              :: nLAIclass         ! Number of LAI classes
  integer(i4),    public, dimension(:),   allocatable :: LAIUnitList       ! List of ids of each LAI class in LAILUT
  real(dp),       public, dimension(:,:), allocatable :: LAILUT            ! [m2/m2] Leaf area index for LAIUnit
  !                                                                        ! dim1=land cover class, dim2=month of year 

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type gridGeoRef
  !  dim1=basin id
     integer(i4), dimension(:), allocatable :: ncols        ! Number of columns
     integer(i4), dimension(:), allocatable :: nrows        ! Number of rows
     real(dp), dimension(:), allocatable    :: xllcorner    ! x coordinate of the lowerleft corner
     real(dp), dimension(:), allocatable    :: yllcorner    ! y coordinate of the lowerleft corner
     real(dp), dimension(:), allocatable    :: cellsize     ! Cellsize x = cellsize y
     real(dp), dimension(:), allocatable    :: nodata_value ! Code to define the mask
  end type gridGeoRef

  type(gridGeoRef), public                  :: level0       ! Reference of the input data files
  type(gridGeoRef), public                  :: level1       ! Reference of the hydrological variables
  type(gridGeoRef), public                  :: level11      ! Reference of the routing variables
  type(gridGeoRef), public                  :: level110     ! Reference of the routing variables at L0 scale (e.g., L0_floodPlain)
  type(gridGeoRef), public                  :: level2       ! Reference of the metereological variables

  real(dp), dimension(:), allocatable, public :: longitude  ! 1d longitude array
  real(dp), dimension(:), allocatable, public :: latitude   ! 1d latitude  array
  
  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type period
      integer(i4)      :: dStart      ! first day
      integer(i4)      :: mStart      ! first month
      integer(i4)      :: yStart      ! first year
      integer(i4)      :: dEnd        ! last  day
      integer(i4)      :: mEnd        ! last  month
      integer(i4)      :: yEnd        ! last  year
      integer(i4)      :: julStart    ! first julian day 
      integer(i4)      :: julEnd      ! last  julian day 
      integer(i4)      :: nObs        ! total number of observations
  end type period

  type(period), public :: warmPer     ! time period for warming
  type(period), public :: evalPer     ! time period for model evaluation
  type(period), public :: simPer      ! warmPer + evalPer
  type(period), public :: readPer     ! start and end dates of read period
  integer(i4), public  :: warmingDays ! number of days for warm up period

  ! -------------------------------------------------------------------
  ! BASIN general description
  ! -------------------------------------------------------------------
  integer(i4), public                           :: nBasins            ! Number of basins for multi-basin optimization
  type basinInfo
     ! dim1 = basin id
     ! dim2 = maximum number of gauges in a given basin
     ! discharge measurement gauges
     integer(i4), dimension(:),   allocatable   :: nGauges            ! Number of gauges within a basin
     integer(i4), dimension(:,:), allocatable   :: gaugeIdList        ! Gauge Id list (e.g. 0000444 0000445)
     integer(i4), dimension(:,:), allocatable   :: gaugeIndexList     ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
     integer(i4), dimension(:,:), allocatable   :: gaugeNodeList      ! Gauge node list at L11

     ! discharge inflow gauges (e.g if headwater bsins are missing)
     integer(i4), dimension(:),   allocatable   :: nInflowGauges      ! Number of gauges within a basin
     integer(i4), dimension(:,:), allocatable   :: InflowGaugeIdList  ! Gauge Id list (e.g. 0000444 0000445)
     integer(i4), dimension(:,:), allocatable   :: InflowGaugeIndexList ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
     integer(i4), dimension(:,:), allocatable   :: InflowGaugeNodeList ! Gauge node list at L11

     ! basin outlet
     integer(i4), dimension(:), allocatable     :: L0_rowOutlet       ! Outlet location in L0 
     integer(i4), dimension(:), allocatable     :: L0_colOutlet       ! Outlet location in L0 

     ! for remapping                                                   
     integer(i4), dimension(:), allocatable     :: L0_iStart          ! Starting cell index of a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iEnd            ! Ending cell index of a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iStartMask      ! Starting cell index of mask a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iEndMask        ! Ending cell index of mask a given basin at L0
     logical,     dimension(:), allocatable     :: L0_mask            ! Mask of level0 based on DEM

     integer(i4), dimension(:), allocatable     :: L1_iStart          ! Starting cell index of a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iEnd            ! Ending cell index of a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iStartMask      ! Starting cell index of mask a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iEndMask        ! Ending cell index of mask a given basin at L1
     logical,     dimension(:), allocatable     :: L1_mask            ! Mask of level1

     integer(i4), dimension(:), allocatable     :: L11_iStart         ! Sarting cell index of a given basin at L11 = node
     integer(i4), dimension(:), allocatable     :: L11_iEnd           ! Ending cell index of a given basin at L11   = node
     integer(i4), dimension(:), allocatable     :: L11_iStartMask     ! Starting cell index of mask a given basin at L11
     integer(i4), dimension(:), allocatable     :: L11_iEndMask       ! Ending cell index of mask a given basin at L11
     logical,     dimension(:), allocatable     :: L11_mask           ! Mask of level11

     integer(i4), dimension(:), allocatable     :: L110_iStart        ! Sarting cell index of L0_floodPlain 
     !                                                                ! at a given basin at L110 = node
     integer(i4), dimension(:), allocatable     :: L110_iEnd          ! Ending cell index of L0_floodPlain 
     !                                                                ! at a given basin at L110   = node

     Integer(i4), dimension(:), allocatable     :: L2_iStart          ! Starting cell index of a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iEnd            ! Ending cell index of a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iStartMask      ! Starting cell index of mask a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iEndMask        ! Ending cell index of mask a given basin at L2
     logical,     dimension(:), allocatable     :: L2_mask            ! Mask of level2

  end type basinInfo
  type(basinInfo), public                       :: basin              ! Basin structure

  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain> 
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  real(dp), public, dimension(:), allocatable      :: L0_elev       ! [m]       Elevation (sinks removed)   
  real(dp), public, dimension(:), allocatable      :: L0_slope      ! [%]       Slope 
  real(dp), public, dimension(:), allocatable      :: L0_asp        ! [degree]  Aspect degree
  integer(i4), public, dimension(:), allocatable   :: L0_fAcc       !           Flow accumulation
  integer(i4), public, dimension(:), allocatable   :: L0_fDir       !           Flow direction (standard ArcGIS)
  integer(i4), public, dimension(:), allocatable   :: L0_soilId     !           Soil id
  integer(i4), public, dimension(:), allocatable   :: L0_geoUnit    !           Geologic formation (unit)
  integer(i4), public, dimension(:), allocatable   :: L0_gaugeLoc   !           Location of gauges within the catchment
  integer(i4), public, dimension(:), allocatable   :: L0_InflowGaugeLoc         ! Location of inflow gauges within the catchment
 
  ! input data - land cover
  integer(i4), public, dimension(:), allocatable   :: L0_LCover_LAI ! Special landcover id (upto 10 classes) only for the LAI index
  integer(i4), public, dimension(:,:), allocatable :: L0_LCover     ! Normal  landcover id (upto 3 classes) 
  !                                                                 ! dim1=number grid cells, dim2=Number of land cover scenes

  ! mHM derived variables
  ! dim1 = number grid cells L0
  integer(i4), public                              :: L0_nCells     !      Number of valid cells 
  real(dp), public, dimension(:), allocatable      :: L0_areaCell   ! [m2] Area of a cell at level-0 
  integer(i4), public, dimension(:,:), allocatable :: L0_cellCoor   !      Cell coordinates (row,col) for each grid cell, dim2=2
  integer(i4), public, dimension(:), allocatable   :: L0_Id         !      Level-0 id
  integer(i4), public, dimension(:), allocatable   :: L0_L11_Id     !      Mapping of L11 Id on L0  
  !                                                                 !      (sub-cat. id. == cell Id L11)
  real(dp), public, dimension(:), allocatable      :: L0_slope_emp  !      Empirical quantiles of slope
  integer(i4), public, dimension(:), allocatable   :: L0_draSC      !      Index of draining cell of each sub catchment 
  !                                                                 !      i.e. a routing cell L11
  integer(i4), public, dimension(:), allocatable   :: L0_draCell    !      Draining cell id at L11 of ith cell of L0
  integer(i4), public, dimension(:), allocatable   :: L0_streamNet  !      Stream network
  integer(i4), public, dimension(:), allocatable   :: L0_floodPlain !      Floodplains of stream i
  !
  real(dp),    public, dimension(:,:), allocatable :: L0_gridded_LAI !      gridded LAI data used when timeStep_LAI_input<0
  !                                                                  !      dim1=number of grid cells, dim2=number of LAI time steps
  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  ! dim2 = 2
  integer(i4), public                              :: L1_nCells        !       No. of cells at this level 
  real(dp),    public, dimension(:), allocatable   :: L1_areaCell      ! [km2] Effective area of cell at this level
  integer(i4), public, dimension(:,:), allocatable :: L1_cellCoor      !       Cell coordinates (row,col)
  !                                                                    !       -> <only for the domain> L1 modelling
  integer(i4), public, dimension(:), allocatable   :: L1_Id            !       Level-1 id
  integer(i4), public, dimension(:), allocatable   :: L1_L11_Id        !       Mapping of L11 Id on L1

  integer(i4), public, dimension(:), allocatable   :: L1_upBound_L0    ! Row start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L1_downBound_L0  ! Row end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L1_leftBound_L0  ! Col start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L1_rightBound_L0 ! Col end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L1_nTCells_L0    ! Total number of valid L0 cells
  !                                                                    ! in a given L1 cell 
  ! Forcings
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  real(dp), public, dimension(:,:), allocatable    :: L1_pre           ! [mm]    Precipitation
  real(dp), public, dimension(:,:), allocatable    :: L1_temp          ! [degC]  Air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_pet           ! [mm]    Potential evapotranspiration
  real(dp), public, dimension(:,:), allocatable    :: L1_tmin          ! [degC]  minimum daily air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_tmax          ! [degC]  maximum daily air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_netrad        ! [W m2]  net radiation
  real(dp), public, dimension(:,:), allocatable    :: L1_absvappress   ! [hPa]   absolute vapour pressure
  real(dp), public, dimension(:,:), allocatable    :: L1_windspeed     ! [m s-1] windspeed

  ! Land cover
  ! dim1 = number grid cells L1
  real(dp), public, dimension(:), allocatable      :: L1_fSealed       ! [1]  Fraction of sealed area
  real(dp), public, dimension(:), allocatable      :: L1_fForest       ! [1]  Fraction of forest cover
  real(dp), public, dimension(:), allocatable      :: L1_fPerm         ! [1]  Fraction of permeable cover
  
  ! State variables
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable   :: L1_inter        ! [mm]  Interception               
  real(dp), public, dimension(:), allocatable   :: L1_snowPack     ! [mm]  Snowpack
  real(dp), public, dimension(:), allocatable   :: L1_sealSTW      ! [mm]  Retention storage of impervious areas
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoist    ! [mm]  Soil moisture of each horizon
  real(dp), public, dimension(:), allocatable   :: L1_unsatSTW     ! [mm]  upper soil storage
  real(dp), public, dimension(:), allocatable   :: L1_satSTW       ! [mm]  groundwater storage

  ! Fluxes
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:,:), allocatable :: L1_aETSoil      ! [mm]       Actual ET from soil layers
  real(dp), public, dimension(:), allocatable   :: L1_aETCanopy    ! [mm d-1]   Real evaporation intensity from canopy
  real(dp), public, dimension(:), allocatable   :: L1_aETSealed    ! [mm d-1]   Real evap. from free water surfaces
  real(dp), public, dimension(:), allocatable   :: L1_baseflow     ! [mm d-1]   Baseflow
  real(dp), public, dimension(:,:), allocatable :: L1_infilSoil    ! [mm d-1]   Infiltration intensity each soil horizon 
  real(dp), public, dimension(:), allocatable   :: L1_fastRunoff   ! [mm d-1]   Fast runoff component
  real(dp), public, dimension(:), allocatable   :: L1_melt         ! [mm d-1]   Melting snow depth.
  real(dp), public, dimension(:), allocatable   :: L1_percol       ! [mm d-1]   Percolation.
  real(dp), public, dimension(:), allocatable   :: L1_preEffect    ! [mm d-1]   Effective precip. depth (snow melt + rain)
  real(dp), public, dimension(:), allocatable   :: L1_rain         ! [mm]       Rain precipitation depth
  real(dp), public, dimension(:), allocatable   :: L1_runoffSeal   ! [mm d-1]   Direct runoff from impervious areas
  real(dp), public, dimension(:), allocatable   :: L1_slowRunoff   ! [mm d-1]   Slow runoff component
  real(dp), public, dimension(:), allocatable   :: L1_snow         ! [mm d-1]   Snow precipitation depth
  real(dp), public, dimension(:), allocatable   :: L1_Throughfall  ! [mm d-1]   Throughfall.             
  real(dp), public, dimension(:), allocatable   :: L1_total_runoff ! [m3 s-1]   Generated runoff

  ! Effective parameters
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable   :: L1_alpha        ! [1]            Exponent for the upper reservoir
  real(dp), public, dimension(:), allocatable   :: L1_degDayInc    ! [d-1 degC-1]   Increase of the Degree-day factor
  !                                                                !                per mm of increase in precipitation 
  real(dp), public, dimension(:), allocatable   :: L1_degDayMax    ! [mm-1 degC-1]  Maximum Degree-day factor 
  real(dp), public, dimension(:), allocatable   :: L1_degDayNoPre  ! [mm-1 degC-1]  Degree-day factor with no precipitation.
  real(dp), public, dimension(:), allocatable   :: L1_degDay       ! [mm d-1degC-1] Degree-day factor.
  real(dp), public, dimension(:), allocatable   :: L1_karstLoss    ! [1]    Karstic percolation loss
  real(dp), public, dimension(:), allocatable   :: L1_fAsp         ! [1]    PET correction for aspect
  real(dp), public, dimension(:), allocatable   :: L1_HarSamCoeff  ! [1]    Hargreaves Samani coeffiecient
  real(dp), public, dimension(:,:), allocatable :: L1_PrieTayAlpha ! [1]    Priestley Taylor coeffiecient
  real(dp), public, dimension(:,:), allocatable :: L1_aeroResist   ! [s m-1] aerodynamical resitance
  real(dp), public, dimension(:,:), allocatable :: L1_surfResist   ! [s m-1] bulk surface resitance
  !                                                                ! dim1 = No cells for basin, dim2 = No of Months in year
  real(dp), public, dimension(:,:), allocatable :: L1_fRoots       ! [1]    Fraction of roots in soil horizons   
  real(dp), public, dimension(:), allocatable   :: L1_maxInter     ! [mm]   Maximum interception 
  real(dp), public, dimension(:), allocatable   :: L1_kfastFlow    ! [d-1]  Fast interflow recession coefficient 
  real(dp), public, dimension(:), allocatable   :: L1_kSlowFlow    ! [d-1]  Slow interflow recession coefficient 
  real(dp), public, dimension(:), allocatable   :: L1_kBaseFlow    ! [d-1]  Baseflow recession coefficient 
  real(dp), public, dimension(:), allocatable   :: L1_kPerco       ! [d-1]  percolation coefficient
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistFC  ! [mm]   Soil moisture below which actual ET
  !                                                                !        is reduced linearly till PWP
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistSat ! [mm]   Saturation soil moisture for each horizon [mm]  
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistExp ! [1]    Exponential parameter to how non-linear
  !                                                                !        is the soil water retention
  real(dp), public, dimension(:), allocatable   :: L1_tempThresh   ! [degC]   Threshold temperature for snow/rain
  real(dp), public, dimension(:), allocatable   :: L1_unsatThresh  ! [mm]   Threshhold water depth controlling fast interflow
  real(dp), public, dimension(:), allocatable   :: L1_sealedThresh ! [mm]   Threshhold water depth for surface runoff
  !                                                                !        in sealed surfaces
  real(dp), public, dimension(:,:), allocatable :: L1_wiltingPoint ! [mm]   Permanent wilting point: below which neither 
  !                                                                !        plant can take water nor water can drain in

  ! -------------------------------------------------------------------
  ! L11 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L11
  ! dim2 = 2
  integer(i4), public                              :: L11_nCells        ! No. of routing cells  (= nSC = nNodes)
  integer(i4), public, dimension(:,:), allocatable :: L11_cellCoor      ! Cell coordinates (row,col)
  !                                                                     ! -> <only domain> Routing
  integer(i4), public, dimension(:), allocatable   :: L11_Id            ! Ids of grid at level-11           
  integer(i4), public, dimension(:), allocatable   :: L11_fDir          ! Flow direction (standard notation)

  ! Reference
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable   :: L11_upBound_L0    ! Row start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_downBound_L0  ! Row end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_leftBound_L0  ! Col start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_rightBound_L0 ! Col end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_upBound_L1    ! Row start at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_downBound_L1  ! Row end at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_leftBound_L1  ! Col start at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable   :: L11_rightBound_L1 ! Col end at finer level-1 scale 

  ! Constants
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable   :: L11_rowOut        ! Grid vertical location of the Outlet
  integer(i4), public, dimension(:), allocatable   :: L11_colOut        ! Grid horizontal location  of the Outlet

  ! -------------------------------------------------------------------
  ! L11 NETWORK description
  ! -------------------------------------------------------------------
  ! Fluxes
  ! dim1 = number grid cells L11
  ! dim2 = 2
  real(dp), public, dimension(:), allocatable     :: L11_Qmod        ! [m3 s-1] Simulated discharge
  real(dp), public, dimension(:), allocatable     :: L11_qOUT        ! [m3 s-1] Total outflow from cells L11 at time tt
  real(dp), public, dimension(:,:), allocatable   :: L11_qTIN        !          Total discharge inputs at t-1 and t
  real(dp), public, dimension(:,:), allocatable   :: L11_qTR         !          Routed outflow leaving a node

  ! Stream link description ( drainage network topology )
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable  :: L11_fromN       !         From node 
  integer(i4), public, dimension(:), allocatable  :: L11_toN         !         To node
  integer(i4), public, dimension(:), allocatable  :: L11_netPerm     !         Routing sequence (permutation of L11_rOrder)
  integer(i4), public, dimension(:), allocatable  :: L11_fRow        !         From row in L0 grid 
  integer(i4), public, dimension(:), allocatable  :: L11_fCol        !         From col in L0 grid
  integer(i4), public, dimension(:), allocatable  :: L11_tRow        !         To row in L0 grid 
  integer(i4), public, dimension(:), allocatable  :: L11_tCol        !         To col in L0 grid 
  integer(i4), public, dimension(:), allocatable  :: L11_rOrder      !         Network routing order  
  integer(i4), public, dimension(:), allocatable  :: L11_label       !         Label Id [0='', 1=HeadWater, 2=Sink]
  logical, public, dimension(:), allocatable      :: L11_sink        !         .true. if sink node reached
  real(dp), public, dimension(:), allocatable     :: L11_length      ! [m]     Total length of river link
  real(dp), public, dimension(:), allocatable     :: L11_aFloodPlain ! [m2]    Area of the flood plain
  real(dp), public, dimension(:), allocatable     :: L11_FracFPimp   ! [1]     Fraction of the flood plain with
  !                                                                  !         impervious cover
  real(dp), public, dimension(:), allocatable     :: L11_slope       ! [1]     Average slope of river link

  ! Parameters
  ! dim1 = number grid cells L11
  real(dp), public, dimension(:), allocatable     :: L11_K           ! [d]     kappa: Muskingum travel time parameter.
  real(dp), public, dimension(:), allocatable     :: L11_xi          ! [1]     xi:    Muskingum diffusion parameter
  !                                                                  !                (attenuation).
  real(dp), public, dimension(:), allocatable     :: L11_C1          ! [-]     Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
  real(dp), public, dimension(:), allocatable     :: L11_C2          ! [-]     Routing parameter C2 (")

  ! -------------------------------------------------------------------
  ! Monthly day/night variation of Meteorological variables
  ! for temporal disaggregation
  ! -------------------------------------------------------------------
  ! dim1 = number of months in a year
  real(dp), public, dimension(int(YearMonths,i4)) :: evap_coeff     ! [-] Evap. coef. for free-water surfaces
  real(dp), public, dimension(int(YearMonths,i4)) :: fday_prec      ! [-] Day ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths,i4)) :: fnight_prec    ! [-] Night ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths,i4)) :: fday_pet       ! [-] Day ratio PET  < 1
  real(dp), public, dimension(int(YearMonths,i4)) :: fnight_pet     ! [-] Night ratio PET  < 1
  real(dp), public, dimension(int(YearMonths,i4)) :: fday_temp      ! [-] Day factor mean temp
  real(dp), public, dimension(int(YearMonths,i4)) :: fnight_temp    ! [-] Night factor mean temp

  ! -------------------------------------------------------------------
  ! DEFINE OUTPUTS 
  ! -------------------------------------------------------------------
  !
  integer(i4)                      :: timeStep_model_outputs ! timestep for writing model outputs
  logical, dimension(nOutFlxState) :: outputFlxState         ! Define model outputs see "mhm_outputs.nml"
  !                                                            dim1 = number of output variables to be written 
  !
END MODULE mo_global_variables
