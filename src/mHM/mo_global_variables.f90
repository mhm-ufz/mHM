!> \file mo_global_variables.f90

!> \brief Global variables ONLY used in reading, writing and startup.

!> \details

!> \authors Luis Samaniego
!> \date Dec 2012

MODULE mo_global_variables

  ! This module provides

  !
  ! Written   Luis Samaniego,     Dec 2005
  ! Modified  Luis Samaniego,     Feb 2013 - new variable names, new modules, units
  !           Rohini Kumar,       Jul 2013 - fraction of perfectly sealed area within city added
  !           Rohini Kumar,       Aug 2013 - name changed "inputFormat" to "inputFormat_meteo_forcings"
  !           Rohini Kumar,       Aug 2013 - name changed from "L0_LAI" to "L0_LCover_LAI"
  !           Rohini Kumar,       Aug 2013 - added dirSoil_LUT and dirGeology_LUT
  !           Luis Samaniego,     Nov 2013 - documentation of dimensions
  !           Matthias Zink,      Nov 2013 - added "InflowGauge" and inflow gauge variabels in basin
  !           Rohini Kumar,       May 2014 - added options for the model run cordinate system
  !           Stephan Thober,     Jun 2014 - added timeStep_model_inputs and readPer
  !           Stephan Thober,     Jun 2014 - added perform_mpr, updated restart flags
  !           Cuntz M. & Mai J.,  Nov 2014 - LAI input from daily, monthly or yearly files
  !           Matthias Zink,      Dec 2014 - adopted inflow gauges to ignore headwater cells
  !           Matthias Zink,      Mar 2015 - added optional soil mositure readin: dirSoil_moisture, L1_sm
  !           Stephan Thober,     Aug 2015 - moved routing related variables to mRM
  !           Oldrich Rakovec,    Oct 2015 - added definition of basin averaged TWS data
  !           Rohini Kumar,       Mar 2016 - new variables for handling different soil databases
  !           Johann Brenner,     Feb 2017 - added optional evapotranspiration readin: dirEvapotranspiration, L1_et
  !           Zink M. Demirel C., Mar 2017 - added Jarvis soil water stress variable for SM process(3) 
  !           Demirel M.C.        May 2017 - added L1_petLAIcorFactor for PET correction
  !           O. Rakovec, R.Kumar Nov 2017 - added project description for the netcdf outputs
  
  USE mo_kind,             ONLY: i4, i8, dp
  use mo_common_variables, ONLY: period
  USE mo_mhm_constants,    ONLY: nOutFlxState, YearMonths, maxNoBasins, maxNLCovers

  IMPLICIT NONE

  ! Types have to be public
  PUBLIC :: period
  ! PUBLIC :: gaugingStation
  PUBLIC :: soilType
  PUBLIC :: gridGeoRef
  PUBLIC :: basinInfo

  ! -------------------------------------------------------------------
  ! PROJECT DESCRIPTION for the NETCDF output file
  ! -------------------------------------------------------------------
  character(1024),public                             :: project_details            ! project including funding instituion., PI, etc.
  character(1024),public                             :: setup_description          ! any specific description of simulation 
  character(1024),public                             :: simulation_type            ! e.g. seasonal forecast, climate projection, ...
  character(256),public                              :: Conventions                ! convention used for dataset
  character(1024),public                             :: contact                    ! contact details, incl. PI name
  character(1024),public                             :: mHM_details                ! developing institution, specific mHM revision
  character(1024),public                             :: history                    ! details on version/creation date

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mHM
  ! -------------------------------------------------------------------
  integer(i4),   public                              :: timeStep                   ! [h] simulation time step (= TS) in [h]
  integer(i4),   dimension(:), allocatable, public   :: timeStep_model_inputs      ! frequency for reading meteo input
  real(dp),      dimension(:), allocatable, public   :: resolutionHydrology        ! [m or degree] resolution of hydrology - Level 1
  real(dp),      dimension(:), allocatable, public   :: resolutionRouting          ! [m or degree] resolution of routing - Level 11
  integer(i4),   dimension(:), allocatable, public   :: L0_Basin
  logical,       public                              :: read_restart               ! flag
  logical,       public                              :: write_restart              ! flag
  logical,       public                              :: perform_mpr                ! switch for performing
  logical,       public                              :: read_meteo_weights         ! read weights for tavg and pet
  ! multiscale parameter regionalization
  character(256),public                              :: inputFormat_meteo_forcings ! format of meteo input data(bin or nc)
  ! LAI information
  character(256), public                             :: inputFormat_gridded_LAI    ! format of gridded LAI data(bin or nc)
  integer(i4),    public                             :: timeStep_LAI_input         ! time step of gridded LAI input
  ! Optional data
  integer(i4),    public                             :: timeStep_sm_input          ! time step of optional data: soil moisture sm
  integer(i4),    public                             :: timeStep_neutrons_input    ! time step of optional data: neutron counts
  integer(i4),    public                             :: timeStep_et_input          ! time step of optional data: evapotransp. et
  integer(i4),    public                             :: iFlag_cordinate_sys        ! options model for the run cordinate system
  integer(i4),    public                             :: iFlag_soilDB               ! options to handle different soil databases
  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nBasins
  character(256), dimension(:), allocatable, public :: dirMorpho          ! Directory where morphological files are located
  character(256), dimension(:), allocatable, public :: dirLCover          ! Directory where land cover files are located
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
  character(256), dimension(:), allocatable, public :: dirgridded_LAI     ! Directory where gridded LAI is located
                                                                          ! used when timeStep_LAI_input < 0
  character(256), dimension(:), allocatable, public :: fileLatLon          ! directory to lat lon files

  character(256), dimension(:), allocatable, public :: dirSoil_moisture        ! File of monthly soil moisture
  character(256), dimension(:), allocatable, public :: fileTWS                 ! File of tws data
  character(256), dimension(:), allocatable, public :: dirNeutrons             ! File of spatio-temporal neutron data
  character(256), dimension(:), allocatable, public :: dirEvapotranspiration   ! File of monthly soil moisture

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

  ! ------------------------------------------------------------------
  ! BASIN AVERAGED TOTAL WATER STORAGE DATA
  ! ------------------------------------------------------------------
  type TWSstructure
     integer(i4),    dimension(:),     allocatable :: basinId            ! Basin Id
     character(256), dimension(:),     allocatable :: fname              ! file name
     real(dp),       dimension(:,:),   allocatable :: TWS                ! [mm]
  end type TWSstructure
  type(TWSstructure), public                       :: basin_avg_TWS_obs   ! [mm] basin average TWS observational data

  real(dp), public, dimension(:,:), allocatable    :: basin_avg_TWS_sim  ! variable containing basin average TWS for each basin
  integer(i4), public                              :: nMeasPerDay_TWS    ! Number of WTS observations per day,
  !                                                                      ! e.g. 24 -> hourly, 1 -> daily

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
  integer(i4),    public                              :: nLCoverScene        ! Number of land cover scene
  character(256), public, dimension(:),   allocatable :: LCfilename          ! file names for the different land cover scenes
  integer(i4),    public, dimension(:,:), allocatable :: LCyearId            ! Mapping of landcover scenes (1, 2, ...) for each basin
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
  type(gridGeoRef), public                  :: level2       ! Reference of the metereological variables

  real(dp), dimension(:), allocatable, public :: L0_longitude  ! 1d longitude array
  real(dp), dimension(:), allocatable, public :: L0_latitude   ! 1d latitude  array
  real(dp), dimension(:), allocatable, public :: L1_longitude  ! 1d longitude array for active grid cells
  real(dp), dimension(:), allocatable, public :: L1_latitude   ! 1d latitude  array for active grid cells
  real(dp), dimension(:), allocatable, public :: L1_rect_longitude  ! 1d longitude array for whole basin rectangle
  real(dp), dimension(:), allocatable, public :: L1_rect_latitude   ! 1d latitude  array for whole basin rectangle

  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type(period), dimension(:), allocatable, public :: warmPer     ! time period for warming
  type(period), dimension(:), allocatable, public :: evalPer     ! time period for model evaluation
  type(period), dimension(:), allocatable, public :: simPer      ! warmPer + evalPer
  type(period),                            public :: readPer     ! start and end dates of read period
  integer(i4),  dimension(:), allocatable, public :: warmingDays ! number of days for warm up period

  ! -------------------------------------------------------------------
  ! BASIN general description
  ! -------------------------------------------------------------------
  integer(i4), public                           :: nBasins            ! Number of basins for multi-basin optimization
  type basinInfo
     ! dim1 = basin id
     ! for remapping
     integer(i4), dimension(:), allocatable     :: L0_iStart          ! Starting cell index of a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iEnd            ! Ending cell index of a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iStartMask      ! Starting cell index of mask a given basin at L0
     integer(i4), dimension(:), allocatable     :: L0_iEndMask        ! Ending cell index of mask a given basin at L0
     logical,     dimension(:), pointer         :: L0_mask            ! Mask of level0 based on DEM

     integer(i4), dimension(:), allocatable     :: L1_iStart          ! Starting cell index of a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iEnd            ! Ending cell index of a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iStartMask      ! Starting cell index of mask a given basin at L1
     integer(i4), dimension(:), allocatable     :: L1_iEndMask        ! Ending cell index of mask a given basin at L1
     logical,     dimension(:), allocatable     :: L1_mask            ! Mask of level1

     Integer(i4), dimension(:), allocatable     :: L2_iStart          ! Starting cell index of a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iEnd            ! Ending cell index of a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iStartMask      ! Starting cell index of mask a given basin at L2
     integer(i4), dimension(:), allocatable     :: L2_iEndMask        ! Ending cell index of mask a given basin at L2
     logical,     dimension(:), allocatable     :: L2_mask            ! Mask of level2

  end type basinInfo
  type(basinInfo), public                       :: basin              ! Basin structure
  logical, dimension(:), allocatable, target    :: L0_mask            ! target variable for mRM

  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain>
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  real(dp), public, dimension(:), allocatable, target :: L0_elev    ! [m]      Elevation (sinks removed)
  !                                                                 !           target variable for coupling to mRM
  real(dp), public, dimension(:), allocatable         :: L0_slope    ! [%]      Slope
  real(dp), public, dimension(:), allocatable         :: L0_asp     ! [degree]  Aspect degree
  integer(i4), public, dimension(:,:), allocatable    :: L0_soilId  !           soil id (iFlag_soilDB = 0)  
  !  [dim1=number grid cells, dim2=Number of soil horizons] note: for iFlag_soilDB=0, dim2=1
  integer(i4), public, dimension(:), allocatable   :: L0_geoUnit    !      Geologic formation (unit)

  ! input data - land cover
  integer(i4), public, dimension(:), allocatable           :: L0_LCover_LAI  ! Special landcover class for the LAI index
  integer(i4), public, dimension(:,:), allocatable, target :: L0_LCover      ! Classic mHM landcover class (upto 3 classes)
  !                                                                          ! dim1=number grid cells, dim2=Number of land cover scenes
  !                                                                          ! target variable for coupling to mRM

  ! mHM derived variables
  ! dim1 = number grid cells L0
  integer(i4), public                                :: L0_nCells            ! Number of valid cells
  integer(i4), public, dimension(:,:), allocatable   :: L0_cellCoor          ! Cell coordinates (row,col) for each grid cell, dim2=2
  integer(i4), public, dimension(:),   allocatable   :: L0_Id                ! Level-0 id
  real(dp),    public, dimension(:),   allocatable   :: L0_slope_emp         ! Empirical quantiles of slope
  !
  real(dp),    public, dimension(:,:), allocatable   :: L0_gridded_LAI       ! gridded AI data used when timeStep_LAI_input<0 or==1
  !                                                                          ! dim1=number of gridcells, dim2=number LAI timesteps

  real(dp),    public, dimension(:),   allocatable   :: L0_areaCell          ! [m2] Area of a cell at level-0

  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  ! dim2 = 2
  integer(i4), public                              :: L1_nCells        !       No. of cells at this level
  integer(i4), public, dimension(:,:), allocatable :: L1_cellCoor      !       Cell coordinates (row,col)
  !                                                                    !       -> <only for the domain> L1 modelling
  integer(i4), public, dimension(:), allocatable   :: L1_Id            !       Level-1 id
  real(dp),    public, dimension(:), allocatable   :: L1_areaCell      ! [km2] Effective area of cell at this level

  integer(i4), public, dimension(:), allocatable   :: L1_upBound_L0    ! Row start at finer level-0 scale
  integer(i4), public, dimension(:), allocatable   :: L1_downBound_L0  ! Row end at finer level-0 scale
  integer(i4), public, dimension(:), allocatable   :: L1_leftBound_L0  ! Col start at finer level-0 scale
  integer(i4), public, dimension(:), allocatable   :: L1_rightBound_L0 ! Col end at finer level-0 scale
  integer(i4), public, dimension(:), allocatable   :: L1_nTCells_L0    ! Total number of valid L0 cells
  !                                                                    ! in a given L1 cell
  ! Forcings
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  real(dp), public, dimension(:,:,:), allocatable  :: L1_temp_weights  ! hourly temperature weights for daily values
  real(dp), public, dimension(:,:,:), allocatable  :: L1_pet_weights   ! hourly pet weights for daily values
  real(dp), public, dimension(:,:,:), allocatable  :: L1_pre_weights   ! hourly pre weights for daily values
  real(dp), public, dimension(:,:), allocatable    :: L1_pre           ! [mm]    Precipitation
  real(dp), public, dimension(:,:), allocatable    :: L1_temp          ! [degC]  Air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_pet           ! [mm TST-1] Potential evapotranspiration
  real(dp), public, dimension(:,:), allocatable    :: L1_tmin          ! [degC]  minimum daily air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_tmax          ! [degC]  maximum daily air temperature
  real(dp), public, dimension(:,:), allocatable    :: L1_netrad        ! [W m2]  net radiation
  real(dp), public, dimension(:,:), allocatable    :: L1_absvappress   ! [Pa]    absolute vapour pressure
  real(dp), public, dimension(:,:), allocatable    :: L1_windspeed     ! [m s-1] windspeed

  ! optional data
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  ! soil moisture
  real(dp), public, dimension(:,:), allocatable    :: L1_sm                  ! [-] soil moisture input for optimization
  logical,  public, dimension(:,:), allocatable    :: L1_sm_mask             ! [-] mask for valid data in L1_sm
  integer(i4)                                      :: nTimeSteps_L1_sm       ! [-] number of time steps in L1_sm_mask
  integer(i4)                                      :: nSoilHorizons_sm_input ! No. of mhm soil horizons equivalent to sm input
  ! neutrons
  real(dp), public, dimension(:,:), allocatable    :: L1_neutronsdata            ! [cph] ground albedo neutrons input
  logical,  public, dimension(:,:), allocatable    :: L1_neutronsdata_mask       ! [cph] mask for valid data in L1_neutrons
  integer(i4)                                      :: nTimeSteps_L1_neutrons     ! [-] number of time steps in L1_neutrons_mask
  ! evapotranspiration
  real(dp), public, dimension(:,:), allocatable    :: L1_et                 ! [mm] Evapotranspiration input for optimization
  logical,  public, dimension(:,:), allocatable    :: L1_et_mask            ! [mm] mask for valid data in L1_neutrons
  integer(i4)                                      :: nTimeSteps_L1_et      ! [-] number of time steps in L1_sm_mask

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
  real(dp), public, dimension(:), allocatable   :: L1_neutrons     ! [mm]  Ground Albedo Neutrons

  ! Fluxes
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:),   allocatable :: L1_pet_calc     ! [mm TS-1] estimated/corrected potential evapotranspiration
  real(dp), public, dimension(:,:), allocatable :: L1_aETSoil      ! [mm TS-1] Actual ET from soil layers
  real(dp), public, dimension(:),   allocatable :: L1_aETCanopy    ! [mm TS-1] Real evaporation intensity from canopy
  real(dp), public, dimension(:),   allocatable :: L1_aETSealed    ! [mm TS-1] Real evap. from free water surfaces
  real(dp), public, dimension(:),   allocatable :: L1_baseflow     ! [mm TS-1] Baseflow
  real(dp), public, dimension(:,:), allocatable :: L1_infilSoil    ! [mm TS-1] Infiltration intensity each soil horizon
  real(dp), public, dimension(:),   allocatable :: L1_fastRunoff   ! [mm TS-1] Fast runoff component
  real(dp), public, dimension(:),   allocatable :: L1_melt         ! [mm TS-1] Melting snow depth.
  real(dp), public, dimension(:),   allocatable :: L1_percol       ! [mm TS-1] Percolation.
  real(dp), public, dimension(:),   allocatable :: L1_preEffect    ! [mm TS-1] Effective precip. depth (snow melt + rain)
  real(dp), public, dimension(:),   allocatable :: L1_rain         ! [mm TS-1] Rain precipitation depth
  real(dp), public, dimension(:),   allocatable :: L1_runoffSeal   ! [mm TS-1] Direct runoff from impervious areas
  real(dp), public, dimension(:),   allocatable :: L1_slowRunoff   ! [mm TS-1] Slow runoff component
  real(dp), public, dimension(:),   allocatable :: L1_snow         ! [mm TS-1] Snow precipitation depth
  real(dp), public, dimension(:),   allocatable :: L1_Throughfall  ! [mm TS-1] Throughfall.
  real(dp), public, dimension(:),   allocatable :: L1_total_runoff ! [m3 TS-1] Generated runoff

  ! Effective parameters
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable   :: L1_alpha                 ! [1]            Exponent for the upper reservoir
  real(dp), public, dimension(:), allocatable   :: L1_degDayInc             ! [d-1 degC-1]   Increase of the Degree-day factor
  !                                                                         !                per mm of increase in precipitation
  real(dp), public, dimension(:), allocatable   :: L1_degDayMax             ! [mm-1 degC-1]  Maximum Degree-day factor
  real(dp), public, dimension(:), allocatable   :: L1_degDayNoPre           ! [mm-1 degC-1]  Degree-day factor with no precipitation.
  real(dp), public, dimension(:), allocatable   :: L1_degDay                ! [mm d-1degC-1] Degree-day factor.
  real(dp), public, dimension(:), allocatable   :: L1_karstLoss             ! [1]    Karstic percolation loss
  real(dp), public, dimension(:), allocatable   :: L1_fAsp                  ! [1]    PET correction for aspect
  real(dp), public, dimension(:), allocatable   :: L1_petLAIcorFactor       ! [-]   PET correction based on LAI (KC by GEUS.dk)

  real(dp), public, dimension(:), allocatable   :: L1_HarSamCoeff           ! [1]    Hargreaves Samani coeffiecient
  real(dp), public, dimension(:,:), allocatable :: L1_PrieTayAlpha          ! [1]    Priestley Taylor coeffiecient
  real(dp), public, dimension(:,:), allocatable :: L1_aeroResist            ! [s m-1] aerodynamical resitance
  real(dp), public, dimension(:,:), allocatable :: L1_surfResist            ! [s m-1] bulk surface resitance
  !                                                                         ! dim1 = No cells for basin, dim2 = No of Months in year
  real(dp), public, dimension(:,:), allocatable :: L1_fRoots                ! [1]    Fraction of roots in soil horizons
  real(dp), public, dimension(:), allocatable   :: L1_maxInter              ! [mm]   Maximum interception
  

  
  real(dp), public, dimension(:), allocatable   :: L1_kfastFlow             ! [d-1]  Fast interflow recession coefficient
  real(dp), public, dimension(:), allocatable   :: L1_kSlowFlow             ! [d-1]  Slow interflow recession coefficient
  real(dp), public, dimension(:), allocatable   :: L1_kBaseFlow             ! [d-1]  Baseflow recession coefficient
  real(dp), public, dimension(:), allocatable   :: L1_kPerco                ! [d-1]  percolation coefficient
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistFC           ! [mm]   Soil moisture below which actual ET
  !                                                                         !        is reduced linearly till PWP
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistSat          ! [mm]   Saturation soil moisture for each horizon [mm]
  real(dp), public, dimension(:,:), allocatable :: L1_soilMoistExp          ! [1]    Exponential parameter to how non-linear
  !                                                                         !        is the soil water retention
  real(dp), public, dimension(:),   allocatable :: L1_jarvis_thresh_c1      ![1] jarvis critical value for normalized soil 
  !                                                                         !        water content 
  real(dp), public, dimension(:), allocatable   :: L1_tempThresh            ! [degC]   Threshold temperature for snow/rain
  real(dp), public, dimension(:), allocatable   :: L1_unsatThresh           ! [mm]  Threshold waterdepth controlling fast interflow
  real(dp), public, dimension(:), allocatable   :: L1_sealedThresh          ! [mm]  Threshold waterdepth for surface runoff
  !                                                                         !       in sealed surfaces
  real(dp), public, dimension(:,:), allocatable :: L1_wiltingPoint          ! [mm]  Permanent wilting point: below which neither
  !                                                                         !       plant can take water nor water can drain in
            
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

  ! -------------------------------------------------------------------
  ! AUXILIARY VARIABLES
  ! -------------------------------------------------------------------
  !

  real(dp), public, dimension(:), allocatable :: neutron_integral_AFast ! pre-calculated integrand for
                                                                        ! vertical projection of isotropic neutron flux

END MODULE mo_global_variables
