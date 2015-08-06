!> \file mo_global_variables_routing.f90
!> \brief Global variables for routing ONLY used in reading, writing and startup.
!> \details 
!> \authors Stephan Thober
!> \date Aug 2015
module mo_global_variables_routing
  use mo_kind, only: i4, i8, dp
  implicit none
  ! Types have to be public
  PUBLIC :: gaugingStation

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mRM
  ! -------------------------------------------------------------------
  real(dp), dimension(:), allocatable, public :: resolutionRouting ! [m or °] resolution of routing - Level 11
  
  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  character(256), dimension(:), allocatable, public :: dirGauges ! Directory where discharge files are located

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
  type(gridGeoRef), public :: level11 ! Reference of the routing variables
  type(gridGeoRef), public :: level110 ! Reference of the routing variables at L0 scale (e.g., L0_floodPlain)

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

  
  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain> 
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells
  ! input data - morphological variables
  integer(i4), public, dimension(:), allocatable :: L0_fAcc ! Flow accumulation
  integer(i4), public, dimension(:), allocatable :: L0_fDir ! Flow direction (standard ArcGIS)
  integer(i4), public, dimension(:), allocatable :: L0_gaugeLoc ! Location of gauges within the catchment
  integer(i4), public, dimension(:), allocatable :: L0_InflowGaugeLoc ! Location of inflow gauges within catchment
  integer(i4), public, dimension(:), allocatable :: L0_L11_Id ! Mapping of L11 Id on L0  
  !                                                           ! (sub-cat. id. == cell Id L11)
  !
  ! mRM derived variables
  ! dim1 = number grid cells L0
  integer(i4), public, dimension(:), allocatable   :: L0_draSC      !      Index of draining cell of each sub catchment 
  !                                                                 !      i.e. a routing cell L11
  integer(i4), public, dimension(:), allocatable   :: L0_draCell    !      Draining cell id at L11 of ith cell of L0
  integer(i4), public, dimension(:), allocatable   :: L0_streamNet  !      Stream network
  integer(i4), public, dimension(:), allocatable   :: L0_floodPlain !      Floodplains of stream i

  
  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! dim1 = number grid cells L1
  integer(i4), public, dimension(:), allocatable :: L1_L11_Id ! Mapping of L11 Id on L1

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

  real(dp), public, dimension(:), allocatable     :: L11_FracFPimp   ! [1]     Fraction of the flood plain with
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
  !                                                                  !         impervious cover
  real(dp), public, dimension(:), allocatable     :: L11_slope       ! [1]     Average slope of river link

  ! Parameters
  ! dim1 = number grid cells L11
  real(dp), public, dimension(:), allocatable     :: L11_K           ! [d]     kappa: Muskingum travel time parameter.
  real(dp), public, dimension(:), allocatable     :: L11_xi          ! [1]     xi:    Muskingum diffusion parameter
  !                                                                  !                (attenuation).
  real(dp), public, dimension(:), allocatable     :: L11_C1          ! [-]     Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
  real(dp), public, dimension(:), allocatable     :: L11_C2          ! [-]     Routing parameter C2 (")

end module mo_global_variables_routing
