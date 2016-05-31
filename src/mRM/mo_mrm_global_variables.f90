!> \file mo_mrm_global_variables.f90

!> \brief Global variables for mRM only

!> \details 

!> \authors Luis Samaniego, Stephan Thober
!> \date Aug 2015

module mo_mrm_global_variables
  
  use mo_kind,             only: i4, i8, dp
  use mo_mrm_constants,    only: nOutFlxState
  use mo_common_variables, only: period
  
  implicit none
  
  PUBLIC :: gaugingStation

  ! -------------------------------------------------------------------
  ! General variables
  ! -------------------------------------------------------------------
  integer(i4) :: mrm_coupling_mode ! 0 = stand alone
  !                                ! 1 = general coupling to a model
  !                                ! 2 = specific coupling to mhm
  logical :: is_start              ! flag for first timestep for mpr

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mRM
  ! -------------------------------------------------------------------
  integer(i4), public :: timeStep ! [h] simulation time step (= TS) in [h]
  integer(i4), dimension(:), allocatable, public :: timeStep_model_inputs ! frequency for reading meteo input
  integer(i4), public :: iFlag_cordinate_sys ! options model for the run cordinate system
  integer(i4), dimension(:), allocatable, public :: L0_Basin
  real(dp), dimension(:), allocatable, public :: resolutionRouting ! [m or degree] resolution of routing - Level 11
  real(dp), dimension(:), allocatable, public :: resolutionHydrology ! [m or degree] resolution of routing - Level 11
  logical, public :: read_restart ! flag 
  logical, public :: write_restart ! flag 
  logical, public :: perform_mpr ! switch for performing multiscale parameter regionalization

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nBasins
  character(256), public :: dirConfigOut
  character(256), public :: dirCommonFiles ! directory where common input files should be located
  character(256), dimension(:), allocatable, public :: dirMorpho ! Directory where morphological files are located
  character(256), dimension(:), allocatable, public :: dirLCover ! Directory where land cover files are located
  character(256), dimension(:), allocatable, public :: dirGauges ! Directory where discharge files are located
  character(256), dimension(:), allocatable, public :: dirTotalRunoff ! Directory wherse simulated total runoff files are located
  character(256), dimension(:), allocatable, public :: dirOut ! Directory where output is written to
  character(256), dimension(:), allocatable, public :: dirRestartOut ! Directory where output of restart is written
  character(256), dimension(:), allocatable, public :: dirRestartIn! Directory where input of restart is read from
  character(256), dimension(:), allocatable, public :: fileLatLon ! Directory where the Lat Lon Files are located
  
  ! ------------------------------------------------------------------
  ! CONSTANT 
  ! ------------------------------------------------------------------
  integer(i4), public :: nTstepDay ! Number of time intervals per day
  !                                ! (was previously NAGG)

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type gridGeoRef
  !  dim1=basin id
     integer(i4), dimension(:), allocatable :: ncols ! Number of columns
     integer(i4), dimension(:), allocatable :: nrows ! Number of rows
     real(dp), dimension(:), allocatable :: xllcorner ! x coordinate of the lowerleft corner
     real(dp), dimension(:), allocatable :: yllcorner ! y coordinate of the lowerleft corner
     real(dp), dimension(:), allocatable :: cellsize ! Cellsize x = cellsize y
     real(dp), dimension(:), allocatable :: nodata_value ! Code to define the mask
  end type gridGeoRef
  type(gridGeoRef), public :: level0 ! grid information at morphological level (e.g., dem, fDir)
  type(gridGeoRef), public :: level1 ! grid information at runoff level
  type(gridGeoRef), public :: level11 ! Reference of the routing variables
  type(gridGeoRef), public :: level110 ! Reference of the routing variables at L0 scale (e.g., L0_floodPlain)
  !
  real(dp), dimension(:), allocatable, public :: L0_longitude  ! 1d longitude array
  real(dp), dimension(:), allocatable, public :: L0_latitude   ! 1d latitude  array
  real(dp), dimension(:), allocatable, public :: L11_rect_longitude  ! 1d longitude array for whole basin rectangle
  real(dp), dimension(:), allocatable, public :: L11_rect_latitude   ! 1d latitude  array for whole basin rectangle

  ! -----------------------------------------------------------------
  ! RUNOFF variable
  ! -----------------------------------------------------------------
  real(dp), dimension(:,:), allocatable, public :: mRM_runoff ! variable containing runoff for each basin and gauge
  
  ! -----------------------------------------------------------------
  ! GAUGED station data
  ! -----------------------------------------------------------------
  integer(i4), public :: nGaugesTotal ! Number of evaluation gauges for all basins 
  integer(i4), public :: nInflowGaugesTotal ! Number of evaluation gauges for all basins 
  integer(i4), public :: nMeasPerDay ! Number of observations per day,
  !                                  ! e.g. 24 -> hourly discharge, 1 -> daily discharge
  type gaugingStation
     integer(i4), dimension(:), allocatable :: basinId ! Basin Id
     integer(i4), dimension(:), allocatable :: gaugeId ! Gauge Id (e.g. 0000444)
     character(256), dimension(:), allocatable :: fname ! Name runoff file
     real(dp), dimension(:,:), allocatable :: Q ! [m3 s-1] observed daily mean discharge (simPer)
     !                                          ! dim1=number observations, dim2=number of gauges
  end type gaugingStation
  type(gaugingStation), public :: gauge ! Gauging station information
  type(gaugingStation), public :: InflowGauge ! inflow gauge information

  ! -----------------------------------------------------------------
  ! LAND COVER DATA
  ! -----------------------------------------------------------------
  ! Land cover information
  real(dp), public :: fracSealed_cityArea ! fraction of area within city assumed to be
  !                                       ! perfectly sealed [0-1] 
  integer(i4), public :: nLCoverScene ! Number of land cover scene
  character(256), public, dimension(:), allocatable :: LCfilename ! file names for the different land cover scenes
  integer(i4), public, dimension(:,:), allocatable :: LCyearId ! Mapping of landcover scenes (1, 2,..) for each basin

  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type(period), dimension(:), allocatable, public :: warmPer     ! time period for warming
  type(period), dimension(:), allocatable, public :: evalPer     ! time period for model evaluation
  type(period), dimension(:), allocatable, public :: simPer      ! warmPer + evalPer
  type(period), dimension(:), allocatable, public :: readPer     ! start and end dates of read period
  integer(i4), dimension(:), allocatable, public :: warmingDays_mrm ! number of days for warm up period

  ! -------------------------------------------------------------------
  ! BASIN general description
  ! -------------------------------------------------------------------
  integer(i4), public :: nBasins ! Number of basins for multi-basin optimization
  type basinInfo
     ! dim1 = basin id
     ! dim2 = maximum number of gauges in a given basin
     ! discharge measurement gauges
     integer(i4), dimension(:),   allocatable :: nGauges        ! Number of gauges within a basin
     integer(i4), dimension(:,:), allocatable :: gaugeIdList    ! Gauge Id list (e.g. 0000444 0000445)
     integer(i4), dimension(:,:), allocatable :: gaugeIndexList ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
     integer(i4), dimension(:,:), allocatable :: gaugeNodeList  ! Gauge node list at L11

     ! discharge inflow gauges (e.g if headwar bsins are missing)
     integer(i4), dimension(:),   allocatable :: nInflowGauges        ! Number of gauges within a basin
     integer(i4), dimension(:,:), allocatable :: InflowGaugeIdList    ! Gauge Id list (e.g. 0000444 0000445)
     integer(i4), dimension(:,:), allocatable :: InflowGaugeIndexList ! Gauge index list (e.g. 1 for 00444, 2 for 00445)
     integer(i4), dimension(:,:), allocatable :: InflowGaugeNodeList  ! Gauge node list at L11
     logical,     dimension(:,:), allocatable :: InflowGaugeHeadwater ! if headwater cells of inflow gauge will be considered

     ! basin outlet
     integer(i4), dimension(:),   allocatable :: L0_Noutlet
     integer(i4), dimension(:,:), allocatable :: L0_rowOutlet   ! Outlet locations in L0 
     integer(i4), dimension(:,:), allocatable :: L0_colOutlet   ! Outlet locations in L0 

     ! for remapping                                           
     integer(i4), dimension(:), allocatable :: L0_iStart      ! Starting cell index of a given basin at L0
     integer(i4), dimension(:), allocatable :: L0_iEnd        ! Ending cell index of a given basin at L0
     integer(i4), dimension(:), allocatable :: L0_iStartMask  ! Starting cell index of mask a given basin at L0
     integer(i4), dimension(:), allocatable :: L0_iEndMask    ! Ending cell index of mask a given basin at L0
     logical,     dimension(:), pointer     :: L0_mask        ! Mask of level0 based on DEM

     integer(i4), dimension(:), allocatable :: L1_iStart      ! Starting cell index of a given basin at L1
     integer(i4), dimension(:), allocatable :: L1_iEnd        ! Ending cell index of a given basin at L1
     integer(i4), dimension(:), allocatable :: L1_iStartMask  ! Starting cell index of mask a given basin at L1
     integer(i4), dimension(:), allocatable :: L1_iEndMask    ! Ending cell index of mask a given basin at L1
     logical,     dimension(:), allocatable :: L1_mask        ! Mask of level1

     integer(i4), dimension(:), allocatable :: L11_iStart ! Sarting cell index of a given basin at L11 = node
     integer(i4), dimension(:), allocatable :: L11_iEnd ! Ending cell index of a given basin at L11   = node
     integer(i4), dimension(:), allocatable :: L11_iStartMask ! Starting cell index of mask a given basin at L11
     integer(i4), dimension(:), allocatable :: L11_iEndMask   ! Ending cell index of mask a given basin at L11
     logical,     dimension(:), allocatable :: L11_mask       ! Mask of level11

     integer(i4), dimension(:), allocatable :: L110_iStart    ! Sarting cell index of L0_floodPlain 
     !                                                        ! at a given basin at L110 = node
     integer(i4), dimension(:), allocatable :: L110_iEnd      ! Ending cell index of L0_floodPlain 
     !                                                        ! at a given basin at L110   = node
  end type basinInfo
  type(basinInfo), public :: basin_mrm ! Basin structure
  logical, dimension(:), allocatable, target :: L0_mask_mRM ! global target variable in case mRM is not used to mHM

  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain> 
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  real(dp), public, dimension(:), allocatable, target  :: L0_elev_read ! only read if not coupled to mhm
  real(dp), public, dimension(:), pointer :: L0_elev_mRM ! [m] Elevation (sinks removed) variable used for routing
  integer(i4), public, dimension(:), allocatable :: L0_fAcc ! Flow accumulation
  integer(i4), public, dimension(:), allocatable :: L0_fDir ! Flow direction (standard ArcGIS)
  integer(i4), public, dimension(:,:), pointer :: L0_LCover_mRM ! Normal  landcover id (upto 3 classes) 
  !                                                             ! dim1=number grid cells, dim2=Number of land cover scenes
  integer(i4), public, dimension(:,:), allocatable, target :: L0_LCover_read ! read variable if not coupled to mhm

  integer(i4), public, dimension(:,:), allocatable :: L0_cellCoor ! Cell coordinates (row,col) for each grid cell, dim2=2
  integer(i4), public, dimension(:), allocatable :: L0_Id ! Level-0 id

  integer(i4), public, dimension(:), allocatable :: L0_gaugeLoc ! Location of gauges within the catchment
  integer(i4), public, dimension(:), allocatable :: L0_InflowGaugeLoc ! Location of inflow gauges within catchment
  integer(i4), public, dimension(:), allocatable :: L0_L11_Id ! Mapping of L11 Id on L0  
  !                                                           ! (sub-cat. id. == cell Id L11)
  real(dp), public, dimension(:), allocatable :: L0_areaCell ! [m2] Area of a cell at level-0 

  !
  ! mRM derived variables
  ! dim1 = number grid cells L0
  integer(i4), public, dimension(:), allocatable :: L0_draSC      ! Index of draining cell of each sub catchment 
  !                                                               ! i.e. a routing cell L11
  integer(i4), public, dimension(:), allocatable :: L0_draCell    ! Draining cell id at L11 of ith cell of L0
  integer(i4), public, dimension(:), allocatable :: L0_streamNet  ! Stream network
  integer(i4), public, dimension(:), allocatable :: L0_floodPlain ! Floodplains of stream i
  integer(i4) :: L0_nCells ! Number of cells at level 0
  
  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  integer(i4), public, dimension(:), allocatable :: L11_L1_ID  ! Mapping of L11 Id on L1
  real(dp),    public, dimension(:), allocatable :: L1_areaCell ! [km2] Effective area of cell at this level
  integer(i4), public, dimension(:), allocatable :: L1_Id ! Ids of grid at level-11           
  integer(i4) :: L1_nCells ! Number of cells at level 1
  ! -------------------------------------------------------------------
  ! L1 variables
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  ! dim2 = number of timesteps
  real(dp), public, dimension(:,:), allocatable :: L1_total_runoff_in

  ! -------------------------------------------------------------------
  ! L11 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L11
  ! dim2 = 2
  integer(i4), public :: L11_nCells ! No. of routing cells  (= nSC = nNodes)
  integer(i4), public, dimension(:,:), allocatable :: L11_cellCoor ! Cell coordinates (row,col)
  !                                                                ! -> <only domain> Routing
  integer(i4), public, dimension(:), allocatable :: L1_L11_ID  ! Mapping of L1 Id on L11
  real(dp),    public, dimension(:), allocatable :: L11_areaCell ! [km2] Effective area of cell at this level
  integer(i4), public, dimension(:), allocatable :: L11_Id ! Ids of grid at level-11           
  integer(i4), public, dimension(:), allocatable :: L11_fDir ! Flow direction (standard notation)
  integer(i4), public, dimension(:), allocatable :: L11_nOutlets

  ! Reference
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable :: L11_upBound_L0    ! Row start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable :: L11_downBound_L0  ! Row end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable :: L11_leftBound_L0  ! Col start at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable :: L11_rightBound_L0 ! Col end at finer level-0 scale 
  integer(i4), public, dimension(:), allocatable :: L11_upBound_L1    ! Row start at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable :: L11_downBound_L1  ! Row end at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable :: L11_leftBound_L1  ! Col start at finer level-1 scale 
  integer(i4), public, dimension(:), allocatable :: L11_rightBound_L1 ! Col end at finer level-1 scale 

  ! Constants
  ! dim1 = number grid cells L11
  integer(i4), public, dimension(:), allocatable :: L11_rowOut ! Grid vertical location of the Outlet
  integer(i4), public, dimension(:), allocatable :: L11_colOut ! Grid horizontal location  of the Outlet

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

  real(dp), public, dimension(:), allocatable     :: L11_FracFPimp   ! [1]     Fraction of the flood plain with
  integer(i4), public, dimension(:), allocatable  :: L11_fromN       !         From node (sinks are at the end)
  integer(i4), public, dimension(:), allocatable  :: L11_toN         !         To node (sinks are at the end)
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
  ! DEFINE OUTPUTS 
  ! -------------------------------------------------------------------
  !
  integer(i4)                      :: timeStep_model_outputs_mrm ! timestep for writing model outputs
  logical, dimension(nOutFlxState) :: outputFlxState_mrm         ! Define model outputs see "mhm_outputs.nml"
  !                                                              ! dim1 = number of output variables to be written 

end module mo_mrm_global_variables
