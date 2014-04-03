Module mo_set_netcdf_restart
  !
  ! This module sets the netcdf structure V, which means it sets the attributes,
  ! the dimensions and the data arrays
  !
  ! number precision
  use mo_kind, only: i4, dp
  !
  ! use netcdfstruct V
  use mo_NCWrite, only: V, Dnc, Gatt, ndims, nVars
  !
  ! use constants of netcdf library
  use netcdf, only: NF90_CHAR, NF90_FLOAT, NF90_INT, NF90_UNLIMITED, NF90_DOUBLE
  !
  ! incorporate nodata value
  use mo_mhm_constants, only: nodata_dp, nodata_i4
  !
  implicit none
  !
  public :: set_state
  public :: set_L11_config
  public :: set_config
  !
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! provide global target L1 and L0 Variables <<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  integer(i4), dimension(:,:), allocatable, target :: L0_rowCoor_out    ! row cell Coordinates at Level 0
  integer(i4), dimension(:,:), allocatable, target :: L0_colCoor_out    ! column cell Coordinates at Level 0  
  integer(i4), dimension(:,:), allocatable, target :: L0_Id_out         ! Ids of grid at level-0 
  real(dp),    dimension(:,:), allocatable, target :: L0_areaCell_out   ! Ids of grid at level-0
  real(dp),    dimension(:,:), allocatable, target :: L0_slope_emp_out  ! Empirical quantiles of slope
  integer(i4), dimension(:,:), allocatable, target :: L1_basin_Mask_out ! Mask at Level 1
  integer(i4), dimension(:,:), allocatable, target :: L1_Id_out         ! Ids of grid at level-1
  integer(i4), dimension(:,:), allocatable, target :: L1_rowCoor_out    ! row cell Coordinates at Level 1
  integer(i4), dimension(:,:), allocatable, target :: L1_colCoor_out    ! column cell Coordinates at Level 1  
  integer(i4), dimension(:,:), allocatable, target :: L1_upBound_L0_out ! Row start at finer level-0 scale 
  integer(i4), dimension(:,:), allocatable, target :: L1_downBound_L0_out! Row end at finer level-0 scale
  integer(i4), dimension(:,:), allocatable, target :: L1_leftBound_L0_out! Col start at finer level-0 scale
  integer(i4), dimension(:,:), allocatable, target :: L1_rightBound_L0_out! Col end at finer level-0 scale
  real(dp),    dimension(:,:), allocatable, target :: L1_areaCell_out   ! [km2] Effective area of cell at this level
  integer(i4), dimension(:,:), allocatable, target :: L1_nTCells_L0_out ! Total number of valid L0 cells in a given L1 cell

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! provide global target L11 routing Variables <<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! dimension variables
  integer(i4), dimension(:),   allocatable, target :: nrows0             ! Number of rows at Level 0
  integer(i4), dimension(:),   allocatable, target :: ncols0             ! Number of colss at Level 0
  integer(i4), dimension(:),   allocatable, target :: nrows1             ! Number of rows at Level 1
  integer(i4), dimension(:),   allocatable, target :: ncols1             ! Number of cols at Level 1
  integer(i4), dimension(:),   allocatable, target :: nrows11            ! Number of rows at Level 11
  integer(i4), dimension(:),   allocatable, target :: ncols11            ! Number of cols at Level 11
  integer(i4), dimension(:),   allocatable, target :: NoutletCoord       ! Dimension of outlet coordiantes at Level 0
  integer(i4), dimension(:),   allocatable, target :: Ngauges            ! Number of evaluation gauges
  integer(i4), dimension(:),   allocatable, target :: nInflowGauges      ! Number of inflow gauges

  ! actual L11 Variables
  integer(i4), dimension(:,:), allocatable, target :: L11_basin_Mask_out ! Mask at Level 11
  integer(i4), dimension(:,:), allocatable, target :: L11_rowCoor_out    ! row cell Coordinates at Level 11
  integer(i4), dimension(:,:), allocatable, target :: L11_colCoor_out    ! column cell Coordinates at Level 11
  integer(i4), dimension(:,:), allocatable, target :: L11_Id_out         ! Ids of grid at level-11 
  integer(i4), dimension(2),                target :: L0_OutletCoord_out ! Coordinates of Outlet at Level 0
  integer(i4), dimension(:,:), allocatable, target :: L0_draSC_out       ! Index of draining cell of each sub catchment 
  integer(i4), dimension(:,:), allocatable, target :: L0_L11_Id_out      ! Mapping of L11 Id on L0 
  integer(i4), dimension(:,:), allocatable, target :: L1_L11_Id_out      ! Mapping of L11 Id on L1
  integer(i4), dimension(:,:), allocatable, target :: L11_fDir_out       ! Flow direction (standard notation)
  integer(i4), dimension(:,:), allocatable, target :: L11_rowOut_out     ! Grid vertical location of the Outlet
  integer(i4), dimension(:,:), allocatable, target :: L11_colOut_out     ! Grid horizontal location  of the Outlet
  integer(i4), dimension(:,:), allocatable, target :: L11_upBound_L0_out ! Row start at finer level-0 scale 
  integer(i4), dimension(:,:), allocatable, target :: L11_downBound_L0_out! Row end at finer level-0 scale
  integer(i4), dimension(:,:), allocatable, target :: L11_leftBound_L0_out! Col start at finer level-0 scale
  integer(i4), dimension(:,:), allocatable, target :: L11_rightBound_L0_out! Col end at finer level-0 scale
  integer(i4), dimension(:,:), allocatable, target :: L11_upBound_L1_out ! Row start at finer level-1 scale 
  integer(i4), dimension(:,:), allocatable, target :: L11_downBound_L1_out ! Row end at finer level-1 scale
  integer(i4), dimension(:,:), allocatable, target :: L11_leftBound_L1_out ! Col start at finer level-1 scale
  integer(i4), dimension(:,:), allocatable, target :: L11_rightBound_L1_out! Col end at finer level-1 scale 
  integer(i4), dimension(:,:), allocatable, target :: L11_fromN_out      ! From node
  integer(i4), dimension(:,:), allocatable, target :: L11_toN_out        ! To node
  integer(i4), dimension(:,:), allocatable, target :: L11_rOrder_out     ! Network routing order
  integer(i4), dimension(:,:), allocatable, target :: L11_label_out      ! Label Id [0='', 1=HeadWater, 2=Sink]
  integer(i4), dimension(:,:), allocatable, target :: L11_sink_out       ! .true. if sink node reached
  integer(i4), dimension(:,:), allocatable, target :: L11_netPerm_out    ! Routing sequence (permutation of L11_rOrder)
  integer(i4), dimension(:,:), allocatable, target :: L11_fRow_out       ! From row in L0 grid 
  integer(i4), dimension(:,:), allocatable, target :: L11_fCol_out       ! From col in L0 grid
  integer(i4), dimension(:,:), allocatable, target :: L11_tRow_out       ! To row in L0 grid
  integer(i4), dimension(:,:), allocatable, target :: L11_tCol_out       ! To col in L0 grid 
  integer(i4), dimension(:,:), allocatable, target :: L0_draCell_out     ! Draining cell id at L11 of ith cell of L0
  integer(i4), dimension(:),   allocatable, target :: gaugeNodeList_out  ! Id of evaluation gauges on L11
  integer(i4), dimension(:),   allocatable, target :: InflowGaugeNodeList_out ! Id of inflow gauges on L11
  integer(i4), dimension(:,:), allocatable, target :: L0_streamNet_out   ! Stream network
  integer(i4), dimension(:,:), allocatable, target :: L0_floodPlain_out  ! Floodplains of stream i
  real(dp),    dimension(:,:), allocatable, target :: L11_length_out     ! [m]     Total length of river link
  real(dp),    dimension(:,:), allocatable, target :: L11_aFloodPlain_out! [m2]    Area of the flood plain
  real(dp),    dimension(:,:), allocatable, target :: L11_slope_out      ! Average slope of river link

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! provide global target STATE and FLUXES Variables <<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! dimension Variables
  integer(i4), dimension(:), allocatable, target :: x1  ! number of rows at level 1
  integer(i4), dimension(:), allocatable, target :: y1  ! number of cols at level 1
  integer(i4), dimension(:), allocatable, target :: z1  ! number of soil layers at level 1
  integer(i4), dimension(:), allocatable, target :: zM  ! number of months per year at level 1
  integer(i4), dimension(:), allocatable, target :: x11 ! number of rows at level 11
  integer(i4), dimension(:), allocatable, target :: y11 ! number of cols at level 11
  integer(i4), dimension(:), allocatable, target :: z11 ! lag in routing at level 11

  ! provide global Land Cover variables at level 1
  real(dp), dimension(:,:), allocatable, target   :: L1_fSealed_out  ! fraction of Sealed area
  real(dp), dimension(:,:), allocatable, target   :: L1_fForest_out  ! fraction of Forest area
  real(dp), dimension(:,:), allocatable, target   :: L1_fPerm_out    ! fraction of permeable area

  ! provide global state Variables at level 1
  real(dp), dimension(:,:), allocatable, target   :: L1_Inter_out    ! Interception storage
  real(dp), dimension(:,:), allocatable, target   :: L1_snowPack_out ! Snowpack
  real(dp), dimension(:,:), allocatable, target   :: L1_sealSTW_out  ! Retention storage of
                                                                   ! impervious areas
  real(dp), dimension(:,:,:), allocatable, target :: L1_soilMoist_out ! soil moisture at level 1
  real(dp), dimension(:,:), allocatable, target   :: L1_unsatSTW_out ! upper soil storage
  real(dp), dimension(:,:), allocatable, target   :: L1_satSTW_out   ! groundwater storage

  ! provide global flux Variables at level 1
  real(dp), dimension(:,:,:), allocatable, target :: L1_aETSoil_out ! soil actual ET
  real(dp), dimension(:,:), allocatable, target   :: L1_aETCanopy_out ! canopy actual ET
  real(dp), dimension(:,:), allocatable, target   :: L1_aETSealed_out ! sealed actual ET
  real(dp), dimension(:,:), allocatable, target   :: L1_baseflow_out  ! baseflow
  real(dp), dimension(:,:,:), allocatable, target :: L1_infilSoil_out ! soil in-exfiltration
  real(dp), dimension(:,:), allocatable, target   :: L1_fastRunoff_out! fast runoff
  real(dp), dimension(:,:), allocatable, target   :: L1_melt_out      ! snow melt

  real(dp), dimension(:,:), allocatable, target   :: L1_percol_out ! percolation
  real(dp), dimension(:,:), allocatable, target   :: L1_preEffect_out ! effective precip. depth (snow melt + rain)
  real(dp), dimension(:,:), allocatable, target   :: L1_rain_out ! rain (liquid water)
  real(dp), dimension(:,:), allocatable, target   :: L1_runoffSeal_out ! runoff from impervious area
  real(dp), dimension(:,:), allocatable, target   :: L1_slowRunoff_out ! slow runoff 
  real(dp), dimension(:,:), allocatable, target   :: L1_snow_out ! snow (solid water) 
  real(dp), dimension(:,:), allocatable, target   :: L1_Throughfall_out ! throughfall 
  real(dp), dimension(:,:), allocatable, target   :: L1_total_runoff_out ! throughfall 

  !-------------------------------------------
  ! EFFECTIVE PARAMETERS
  !-------------------------------------------
  real(dp), dimension(:,:), allocatable, target   :: L1_alpha_out ! exponent for the upper reservoir
  real(dp), dimension(:,:), allocatable, target   :: L1_degDayInc_out ! increase of the Degree-day 
  !                                                                   ! factor per mm of increase in precipitation
  real(dp), dimension(:,:), allocatable, target   :: L1_degDayMax_out ! maximum degree-day factor 
  real(dp), dimension(:,:), allocatable, target   :: L1_degDayNoPre_out ! degree-day factor with no precipitation
  real(dp), dimension(:,:), allocatable, target   :: L1_degDay_out ! degree-day factor
  real(dp), dimension(:,:), allocatable, target   :: L1_karstLoss_out ! Karstic percolation loss
  real(dp), dimension(:,:), allocatable, target   :: L1_fAsp_out  ! PET correction factor due to terrain aspect
  real(dp), dimension(:,:,:), allocatable, target :: L1_fRoots_out ! Fraction of roots in soil horizons
  real(dp), dimension(:,:), allocatable, target   :: L1_maxInter_out ! Maximum interception 
  real(dp), dimension(:,:), allocatable, target   :: L1_kfastFlow_out ! fast interflow recession coefficient 
  real(dp), dimension(:,:), allocatable, target   :: L1_kSlowFlow_out ! slow interflow recession coefficient 
  real(dp), dimension(:,:), allocatable, target   :: L1_kBaseFlow_out ! baseflow recession coefficient 
  real(dp), dimension(:,:), allocatable, target   :: L1_kPerco_out ! percolation coefficient 
  real(dp), dimension(:,:,:), allocatable, target :: L1_soilMoistFC_out ! Soil moisture below which actual ET 
  !                                                                     ! is reduced linearly till PWP
  real(dp), dimension(:,:,:), allocatable, target :: L1_soilMoistSat_out ! Saturation soil moisture for each horizon [mm]
  real(dp), dimension(:,:,:), allocatable, target :: L1_soilMoistExp_out ! Exponential parameter to how 
  !                                                                      ! non-linear is the soil water retention
  real(dp), dimension(:,:), allocatable, target   :: L1_tempThresh_out ! Threshold temperature for snow/rain 
  real(dp), dimension(:,:), allocatable, target   :: L1_unsatThresh_out ! Threshhold water depth controlling fast interflow
  real(dp), dimension(:,:), allocatable, target   :: L1_sealedThresh_out ! Threshhold water depth for surface runoff 
  !                                                                      ! in sealed surfaces
  real(dp), dimension(:,:,:), allocatable, target :: L1_wiltingPoint_out ! Permanent wilting point
  !-------------------------------------------
  ! L11 ROUTING STATE VARIABLES, FLUXES AND
  !             PARAMETERS
  !-------------------------------------------
  real(dp), dimension(:,:), allocatable, target   :: L11_Qmod_out ! simulated discharge at each node
  real(dp), dimension(:,:), allocatable, target   :: L11_qOUT_out ! Total outflow from cells L11 at time tt
  real(dp), dimension(:,:,:), allocatable, target :: L11_qTIN_out ! Total discharge inputs at t-1 and t
  real(dp), dimension(:,:,:), allocatable, target :: L11_qTR_out ! Routed outflow leaving a node
  real(dp), dimension(:,:), allocatable, target   :: L11_K_out ! kappa: Muskingum travel time parameter.
  real(dp), dimension(:,:), allocatable, target   :: L11_xi_out ! xi: Muskingum diffusion parameter
  real(dp), dimension(:,:), allocatable, target   :: L11_C1_out ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
  real(dp), dimension(:,:), allocatable, target   :: L11_C2_out ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
  real(dp), dimension(:,:), allocatable, target   :: L11_FracFPimp_out ! Fraction of the flood plain with impervious cover

  ! 
  !
contains
  !
  !     NAME
  !        set_config

  !     PURPOSE
  !        sets the netcdf output variable for all configuration variables of L0 and L1. These
  !        variables are defined in L0_variable_init and L1_variable_init (mo_startup) and the corresponding
  !        target variables are set in config_variables_set in mo_restart

  !     HISTORY
  !        author Stephan Thober
  !        date Jul 2013
  subroutine set_config
    !
    use mo_string_utils, only: num2str
    !
    implicit none
    !
    integer(i4) :: ii
    !
    nDims = 4
    nVars = 15 + nDims
    !
    ! allocate arrays
    if (allocated(Dnc)   ) deallocate ( Dnc)
    if (allocated(V)     ) deallocate ( V  )
    allocate ( Dnc(nDims)  )
    allocate ( V(nVars) )
    !
    !
    allocate( nrows0( size(L0_rowCoor_out, 1) ) )
    do ii = 1, size( nrows0 )
       nrows0(ii) = ii
    end do
    !
    allocate( ncols0( size(L0_colCoor_out, 2) ) )
    do ii = 1, size( ncols0 )
       ncols0(ii) = ii
    end do
    !
    allocate( nrows1( size( L1_basin_Mask_out, 1 ) ) )
    do ii = 1, size( nrows1 )
       nrows1(ii) = ii
    end do
    !
    allocate( ncols1( size( L1_basin_Mask_out, 2 ) ) )
    do ii = 1, size( ncols1 )
       ncols1(ii) = ii
    end do
    !
    Dnc(1)%name      = "nrows0"
    Dnc(1)%len       = size(nrows0,1)
    Dnc(2)%name      = "ncols0"
    Dnc(2)%len       = size(ncols0,1)
    Dnc(3)%name      = "nrows1"
    Dnc(3)%len       = size(nrows1,1)
    Dnc(4)%name      = "ncols1"
    Dnc(4)%len       = size(ncols1,1)
    !
    ! DIMENSION VARIABLES ------------------------------------------------------
    !
    ! loop over dimension variables
    do ii = 1, size( DNC )
       V(ii)%name        =  Dnc(ii)%name
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  0
       V(ii)%nSubs       =  0
       V(ii)%nDims       =  1
       V(ii)%dimTypes    =  (/ii,0,0,0,0/)
       V(ii)%wFlag       = .true.
       !
       ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
       V(ii)%nAtt          = 2
       !
       V(ii)%att(1)%name   = "units"
       V(ii)%att(1)%xType  = NF90_CHAR
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values  = ""
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
    end do
    !
    ! nrows11
    V(1)%G1_i        => nrows0
    V(1)%att(2)%values = "Number of rows at level 0"   
    ! ncols11
    V(2)%G1_i        => ncols0
    V(2)%att(2)%values = "Number of cols at level 0"     
    ! nrows0
    V(3)%G1_i        => nrows1
    V(3)%att(2)%values = "Number of rows at level 1"   
    ! ncols0
    V(4)%G1_i        => ncols1
    V(4)%att(2)%values = "Number of cols at level 1"
    
    ! loop over level 0 variables
    do ii = size(DNC) + 1, size(DNC) + 3
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/1,2,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_INT
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    ! 
    ii = size(DNC) + 1
    V(ii)%name        =  "L0_rowCoor"
    V(ii)%G2_i        => L0_rowCoor_out
    V(ii)%att(2)%values = "row coordinate of cell at level 0"
    ! 
    ii = size(DNC) + 2
    V(ii)%name        =  "L0_colCoor"
    V(ii)%G2_i        => L0_colCoor_out
    V(ii)%att(2)%values = "col coordinate of cell at level 0"
    ! 
    ii = size(DNC) + 3
    V(ii)%name        =  "L0_Id"
    V(ii)%G2_i        => L0_Id_out
    V(ii)%att(2)%values = "cell IDs at level 0"
    !
    do ii = size(DNC) + 4, size(DNC) + 5
       V(ii)%xType       =  NF90_DOUBLE
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/1,2,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_DOUBLE
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    !
    ii = size(DNC) + 4
    V(ii)%name        =  "L0_areaCell"
    V(ii)%G2_d        => L0_areaCell_out
    V(ii)%att(2)%values = "Area of a cell at level-0 [m2]"
    !
    ii = size(DNC) + 5
    V(ii)%name        =  "L0_slope_emp"
    V(ii)%G2_d        => L0_slope_emp_out
    V(ii)%att(2)%values = "Empirical quantiles of slope"

    ! loop over level 1 variables
    do ii = size(DNC) + 6, size(DNC) + 15
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/3,4,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_INT
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    !
    ii = size(DNC) + 6
    V(ii)%name        =  "L1_basin_Mask"
    V(ii)%G2_i        => L1_basin_Mask_out
    V(ii)%att(2)%values = "Mask at Level 1"
    !
    ii = size(DNC) + 7
    V(ii)%name        =  "L1_Id"
    V(ii)%G2_i        => L1_Id_out
    V(ii)%att(2)%values = "Ids of grid at level-1"
    !
    ii = size(DNC) + 8
    V(ii)%name        =  "L1_rowCoor"
    V(ii)%G2_i        => L1_rowCoor_out
    V(ii)%att(2)%values = "row cell Coordinates at Level 1"
    !
    ii = size(DNC) + 9
    V(ii)%name        =  "L1_colCoor"
    V(ii)%G2_i        => L1_colCoor_out
    V(ii)%att(2)%values = "column cell Coordinates at Level 1"
    !
    ii = size(DNC) + 10
    V(ii)%name        =  "L1_upBound_L0"
    V(ii)%G2_i        => L1_upBound_L0_out
    V(ii)%att(2)%values = "Row start at finer level-0 scale"
    !
    ii = size(DNC) + 11
    V(ii)%name        =  "L1_downBound_L0"
    V(ii)%G2_i        => L1_downBound_L0_out
    V(ii)%att(2)%values = "Row end at finer level-0 scale"
    !
    ii = size(DNC) + 12
    V(ii)%name        =  "L1_leftBound_L0"
    V(ii)%G2_i        => L1_leftBound_L0_out
    V(ii)%att(2)%values = "Col start at finer level-0 scale"
    !
    ii = size(DNC) + 13
    V(ii)%name        =  "L1_rightBound_L0"
    V(ii)%G2_i        => L1_rightBound_L0_out
    V(ii)%att(2)%values = "Col end at finer level-0 scale"
    !
    ii = size(DNC) + 14
    V(ii)%name        =  "L1_nTCells_L0"
    V(ii)%G2_i        => L1_nTCells_L0_out
    V(ii)%att(2)%values = "Total number of valid L0 cells in a given L1 cell"
    !
    ii = size(DNC) + 15
    V(ii)%name        =  "L1_areaCell"
    V(ii)%xType       =  NF90_DOUBLE
    V(ii)%nLvls       =  1
    V(ii)%nSubs       =  1
    V(ii)%nDims       =  2
    V(ii)%dimTypes    =  (/3,4,0,0,0/)
    ! printing
    V(ii)%wFlag       =  .true.
    V(ii)%G2_d        => L1_areaCell_out
    !
    ! attributes 
    V(ii)%nAtt          = 2
    !
    V(ii)%att(2)%name   = "long_name"
    V(ii)%att(2)%xType  = NF90_CHAR
    V(ii)%att(2)%nValues= 1
    V(ii)%att(2)%values = "Effective area of cell at this level [km2]"
    !
    V(ii)%att(1)%name   = "_FillValue"
    V(ii)%att(1)%xType  = NF90_DOUBLE
    V(ii)%att(1)%nValues= 1
    V(ii)%att(1)%values = trim(num2str(nodata_i4))
    !
  end subroutine set_config

  !     NAME
  !        set_L11_config

  !     PURPOSE
  !        sets the netcdf output variable for all configuration variables at L11. These
  !        variables are defined in the subroutines in mo_net_startup and the corresponding
  !        target variables are set in L11_config_set in mo_restart

  !     HISTORY
  !        author    Stephan Thober
  !        date      Jul 2013

  !        Modified  Matthias Zink , Apr 2014 - added inflow gauge

  subroutine set_L11_config
    !
    use mo_string_utils, only: num2str
    !
    implicit none
    !
    integer(i4) :: ii
    !
    ! define parameters
    nDims = 9            ! Number of dimensions
    nVars = 37 + nDims   ! total number of Variables
    !
    ! allocate arrays
    if (allocated(Dnc)   ) deallocate ( Dnc)
    if (allocated(V)     ) deallocate ( V  )
    allocate ( Dnc(nDims)  )
    allocate ( V(nVars) )
    !
    allocate( nrows11( size( L11_basin_Mask_out, 1 ) ) )
    do ii = 1, size( nrows11 )
       nrows11(ii) = ii
    end do
    !
    allocate( ncols11( size( L11_basin_Mask_out, 2 ) ) )
    do ii = 1, size( ncols11 )
       ncols11(ii) = ii
    end do
    !
    allocate( nrows0( size(L0_draCell_out, 1) ) )
    do ii = 1, size( nrows0 )
       nrows0(ii) = ii
    end do
    !
    allocate( ncols0( size(L0_draCell_out, 2) ) )
    do ii = 1, size( ncols0 )
       ncols0(ii) = ii
    end do
    !
    allocate( nrows1( size(L1_L11_Id_out, 1) ) )
    do ii = 1, size( nrows1 )
       nrows1(ii) = ii
    end do
    !
    allocate( ncols1( size(L1_L11_Id_out, 2) ) )
    do ii = 1, size( ncols1 )
       ncols1(ii) = ii
    end do
    !
    allocate(NoutletCoord( size(L0_OutletCoord_out)))
    do ii = 1, size(NoutletCoord)
       NoutletCoord(ii) = ii
    end do
    !
    allocate(Ngauges( size(gaugeNodeList_out)))
    do ii = 1, size(Ngauges)
       Ngauges(ii) = ii
    end do
    !
    allocate(nInflowGauges( size(InflowGaugeNodeList_out)))
    do ii = 1, size(nInflowGauges)
       nInflowGauges(ii) = ii
    end do
    !
    ! define dimensions --------------------------------------------------------
    Dnc(1)%name      = "nrows11"
    Dnc(1)%len       = size(nrows11,1)
    Dnc(2)%name      = "ncols11"
    Dnc(2)%len       = size(ncols11,1)
    Dnc(3)%name      = "nrows0"
    Dnc(3)%len       = size(nrows0,1)
    Dnc(4)%name      = "ncols0"
    Dnc(4)%len       = size(ncols0,1)
    Dnc(5)%name      = "nrows1"
    Dnc(5)%len       = size(nrows1,1)
    Dnc(6)%name      = "ncols1"
    Dnc(6)%len       = size(ncols1,1)
    Dnc(7)%name      = "NoutletCoord"
    Dnc(7)%len       = size(NoutletCoord,1)
    Dnc(8)%name      = "Ngauges"
    Dnc(8)%len       = size(Ngauges)
    Dnc(9)%name      = "nInflowGauges"
    Dnc(9)%len       = size(nInflowGauges)
    !
    ! DIMENSION VARIABLES ------------------------------------------------------
    !
    ! loop over dimension variables
    do ii = 1, size( DNC )
       V(ii)%name        =  Dnc(ii)%name
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  0
       V(ii)%nSubs       =  0
       V(ii)%nDims       =  1
       V(ii)%dimTypes    =  (/ii,0,0,0,0/)
       V(ii)%wFlag       = .true.
       !
       ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
       V(ii)%nAtt          = 2
       !
       V(ii)%att(1)%name   = "units"
       V(ii)%att(1)%xType  = NF90_CHAR
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values  = ""
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
    end do
    !
    ! nrows11
    V(1)%G1_i        => nrows11
    V(1)%att(2)%values = "Number of rows at level 11"   
    ! ncols11
    V(2)%G1_i        => ncols11
    V(2)%att(2)%values = "Number of cols at level 11"     
    ! nrows0
    V(3)%G1_i        => nrows0
    V(3)%att(2)%values = "Number of rows at level 0"   
    ! ncols0
    V(4)%G1_i        => ncols0
    V(4)%att(2)%values = "Number of cols at level 0"
    ! nrows1
    V(5)%G1_i        => nrows1
    V(5)%att(2)%values = "Number of rows at level 1"
    ! ncols1
    V(6)%G1_i        => ncols1
    V(6)%att(2)%values = "Number of cols at level 1"
    ! Number of Outlet Coordinates
    V(7)%G1_i        => NoutletCoord
    V(7)%att(2)%values = "Number of Outlet Coordinates"
    ! Number of gauges
    V(8)%G1_i        => Ngauges
    V(8)%att(2)%values = "Number of evaluation gauges"
    ! Number of inflow gauges
    V(9)%G1_i        => nInflowGauges
    V(9)%att(2)%values = "Number of inflow gauges"
     
    !
    ! loop over Level 11 variables
    do ii = size(DNC) + 1, size(DNC) + 21
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/1,2,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_INT
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    !
    ii = size(DNC) + 1
    V(ii)%name        =  "L11_basin_Mask"
    V(ii)%G2_i        => L11_basin_Mask_out
    V(ii)%att(2)%values = "Mask at Level 11"
    !    
    ii = size(DNC) + 2
    V(ii)%name        =  "L11_rowCoor"
    V(ii)%G2_i        => L11_rowCoor_out
    V(ii)%att(2)%values = "row coordinates at Level 11"
    !    
    ii = size(DNC) + 3
    V(ii)%name        =  "L11_colCoor"
    V(ii)%G2_i        => L11_colCoor_out
    V(ii)%att(2)%values = "col coordinates at Level 11"
    !    
    ii = size(DNC) + 4
    V(ii)%name        =  "L11_Id"
    V(ii)%G2_i        => L11_Id_out
    V(ii)%att(2)%values = "cell Ids at Level 11"
    !    
    ii = size(DNC) + 5
    V(ii)%name        =  "L11_fDir"
    V(ii)%G2_i        => L11_fDir_out
    V(ii)%att(2)%values = "flow Direction at Level 11"
    !    
    ii = size(DNC) + 6
    V(ii)%name        =  "L11_rowOut"
    V(ii)%G2_i        => L11_rowOut_out
    V(ii)%att(2)%values = "Grid vertical location of the Outlet at Level 11"
    !    
    ii = size(DNC) + 7
    V(ii)%name        =  "L11_colOut"
    V(ii)%G2_i        => L11_colOut_out
    V(ii)%att(2)%values = "Grid horizontal location of the Outlet at Level 11"
    !    
    ii = size(DNC) + 8
    V(ii)%name        =  "L11_upBound_L0"
    V(ii)%G2_i        => L11_upBound_L0_out
    V(ii)%att(2)%values = "Row start at finer level-0 scale of Level 11 cell"
    !    
    ii = size(DNC) + 9
    V(ii)%name        =  "L11_downBound_L0"
    V(ii)%G2_i        => L11_downBound_L0_out
    V(ii)%att(2)%values = "Row end at finer level-0 scale of Level 11 cell"
    !    
    ii = size(DNC) + 10
    V(ii)%name        =  "L11_leftBound_L0"
    V(ii)%G2_i        => L11_leftBound_L0_out
    V(ii)%att(2)%values = "Col start at finer level-0 scale of Level 11 cell"
    !    
    ii = size(DNC) + 11
    V(ii)%name        =  "L11_rightBound_L0"
    V(ii)%G2_i        => L11_leftBound_L0_out
    V(ii)%att(2)%values = "Col end at finer level-0 scale of Level 11 cell"
    !    
    ii = size(DNC) + 12
    V(ii)%name        =  "L11_fromN"
    V(ii)%G2_i        => L11_fromN_out
    V(ii)%att(2)%values = "From Node"
    !    
    ii = size(DNC) + 13
    V(ii)%name        =  "L11_toN"
    V(ii)%G2_i        => L11_toN_out
    V(ii)%att(2)%values = "To Node"
    !    
    ii = size(DNC) + 14
    V(ii)%name        =  "L11_rOrder"
    V(ii)%G2_i        => L11_rOrder_out
    V(ii)%att(2)%values = "Network routing order at Level 11"
    !    
    ii = size(DNC) + 15
    V(ii)%name        =  "L11_label"
    V(ii)%G2_i        => L11_label_out
    V(ii)%att(2)%values = "Label Id [0='', 1=HeadWater, 2=Sink] at Level 11"
    !    
    ii = size(DNC) + 16
    V(ii)%name        =  "L11_sink"
    V(ii)%G2_i        => L11_sink_out
    V(ii)%att(2)%values = ".true. if sink node reached at Level 11"
    !    
    ii = size(DNC) + 17
    V(ii)%name        =  "L11_netPerm"
    V(ii)%G2_i        => L11_netPerm_out
    V(ii)%att(2)%values = "Routing sequence (permutation of L11_rOrder) at Level 11"
    !    
    ii = size(DNC) + 18
    V(ii)%name        =  "L11_fRow"
    V(ii)%G2_i        => L11_fRow_out
    V(ii)%att(2)%values = "From row in L0 grid at Level 11"
    !    
    ii = size(DNC) + 19
    V(ii)%name        =  "L11_fCol"
    V(ii)%G2_i        => L11_fCol_out
    V(ii)%att(2)%values = "From col in L0 grid at Level 11"
    !    
    ii = size(DNC) + 20
    V(ii)%name        =  "L11_tRow"
    V(ii)%G2_i        => L11_tRow_out
    V(ii)%att(2)%values = "To row in L0 grid at Level 11"
    !    
    ii = size(DNC) + 21
    V(ii)%name        =  "L11_tCol"
    V(ii)%G2_i        => L11_tCol_out
    V(ii)%att(2)%values = "To Col in L0 grid at Level 11"

    !
    ! loop over Level 11 variables
    do ii = size(DNC) + 22, size(DNC) + 24
       V(ii)%xType       =  NF90_DOUBLE
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/1,2,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_DOUBLE
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_dp))
       !
    end do
    !
    ii = size(DNC) + 22
    V(ii)%name        =  "L11_length"
    V(ii)%G2_d        => L11_length_out
    V(ii)%att(2)%values = "Total length of river link [m]"
    !
    ii = size(DNC) + 23
    V(ii)%name        =  "L11_aFloodPlain"
    V(ii)%G2_d        => L11_aFloodPlain_out
    V(ii)%att(2)%values = "Area of the flood plain [m2]"
    !
    ii = size(DNC) + 24
    V(ii)%name        =  "L11_slope"
    V(ii)%G2_d        => L11_slope_out
    V(ii)%att(2)%values = "Average slope of river link"

    !
    ! loop over Level 0 variables
    do ii = size(DNC) + 25, size(DNC) + 29
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/3,4,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_INT
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    !
    ii = size(DNC) + 25
    V(ii)%name        =  "L0_draCell"
    V(ii)%G2_i        => L0_draCell_out
    V(ii)%att(2)%values = "Draining cell id at L11 of ith cell of L0"
    !
    ii = size(DNC) + 26
    V(ii)%name        =  "L0_streamNet"
    V(ii)%G2_i        => L0_streamNet_out
    V(ii)%att(2)%values = "Stream network"
    !
    ii = size(DNC) + 27
    V(ii)%name        =  "L0_floodPlain"
    V(ii)%G2_i        => L0_floodPlain_out
    V(ii)%att(2)%values = "Floodplains of stream i"
    !
    ii = size(DNC) + 28
    V(ii)%name        =  "L0_draSC"
    V(ii)%G2_i        => L0_draSC_out
    V(ii)%att(2)%values = "Floodplains of stream i"
    !
    ii = size(DNC) + 29
    V(ii)%name        =  "L0_L11_Id"
    V(ii)%G2_i        => L0_L11_Id_out
    V(ii)%att(2)%values = "Mapping of L11 Id on L0"

    !
    ! loop over Level 1 variables
    do ii = size(DNC) + 30, size(DNC) + 34
       V(ii)%xType       =  NF90_INT
       V(ii)%nLvls       =  1
       V(ii)%nSubs       =  1
       V(ii)%nDims       =  2
       V(ii)%dimTypes    =  (/5,6,0,0,0/)
       ! printing
       V(ii)%wFlag       =  .true.
       !
       ! attributes 
       V(ii)%nAtt          = 2
       !
       V(ii)%att(2)%name   = "long_name"
       V(ii)%att(2)%xType  = NF90_CHAR
       V(ii)%att(2)%nValues= 1
       !
       V(ii)%att(1)%name   = "_FillValue"
       V(ii)%att(1)%xType  = NF90_INT
       V(ii)%att(1)%nValues= 1
       V(ii)%att(1)%values = trim(num2str(nodata_i4))
       !
    end do
    !
    ii = size(DNC) + 30
    V(ii)%name        =  "L1_L11_Id"
    V(ii)%G2_i        => L1_L11_Id_out
    V(ii)%att(2)%values = "Mapping of L11 Id on L1"
    !
    ii = size(DNC) + 31
    V(ii)%name        =  "L11_upBound_L1"
    V(ii)%G2_i        => L11_upBound_L1_out
    V(ii)%att(2)%values = "Row start at finer level-1 scale"
    !
    ii = size(DNC) + 32
    V(ii)%name        =  "L11_downBound_L1"
    V(ii)%G2_i        => L11_downBound_L1_out
    V(ii)%att(2)%values = "Row end at finer level-1 scale"
    !
    ii = size(DNC) + 33
    V(ii)%name        =  "L11_leftBound_L1"
    V(ii)%G2_i        => L11_leftBound_L1_out
    V(ii)%att(2)%values = "Col start at finer level-1 scale"
    !
    ii = size(DNC) + 34
    V(ii)%name        =  "L11_rightBound_L1"
    V(ii)%G2_i        => L11_rightBound_L1_out
    V(ii)%att(2)%values = "Col start at finer level-1 scale"

    ! outlet coordinates
    ii = size(DNC) + 35
    V(ii)%name        =  "L0_OutletCoord"
    V(ii)%xType       =  NF90_INT
    V(ii)%nLvls       =  1
    V(ii)%nSubs       =  1
    V(ii)%nDims       =  1
    V(ii)%dimTypes    =  (/7,0,0,0,0/)
    ! printing
    V(ii)%wFlag       =  .true.
    !
    ! pointer
    V(ii)%G1_i        => L0_OutletCoord_out
    !
    ! attributes 
    V(ii)%nAtt          = 2
    !
    V(ii)%att(2)%name   = "long_name"
    V(ii)%att(2)%xType  = NF90_CHAR
    V(ii)%att(2)%nValues= 1
    V(ii)%att(2)%values = "Outlet Coordinates at Level 0"
    !
    V(ii)%att(1)%name   = "_FillValue"
    V(ii)%att(1)%xType  = NF90_INT
    V(ii)%att(1)%nValues= 1
    V(ii)%att(1)%values = trim(num2str(nodata_i4))
    
    ! gauge node list
    ii = size(DNC) + 36
    V(ii)%name        =  "gaugeNodeList"
    V(ii)%xType       =  NF90_INT
    V(ii)%nLvls       =  1
    V(ii)%nSubs       =  1
    V(ii)%nDims       =  1
    V(ii)%dimTypes    =  (/8,0,0,0,0/)
    ! printing
    V(ii)%wFlag       =  .true.
    !
    ! pointer
    V(ii)%G1_i        => gaugeNodeList_out
    !
    ! attributes 
    V(ii)%nAtt          = 2
    !
    V(ii)%att(2)%name   = "long_name"
    V(ii)%att(2)%xType  = NF90_CHAR
    V(ii)%att(2)%nValues= 1
    V(ii)%att(2)%values = "cell ID of gauges"
    !
    V(ii)%att(1)%name   = "_FillValue"
    V(ii)%att(1)%xType  = NF90_INT
    V(ii)%att(1)%nValues= 1
    V(ii)%att(1)%values = trim(num2str(nodata_i4))
    !
    ! inflow gauge node list
    ii = size(DNC) + 37
    V(ii)%name        =  "InflowGaugeNodeList"
    V(ii)%xType       =  NF90_INT
    V(ii)%nLvls       =  1
    V(ii)%nSubs       =  1
    V(ii)%nDims       =  1
    V(ii)%dimTypes    =  (/9,0,0,0,0/)
    ! printing
    V(ii)%wFlag       =  .true.
    !
    ! pointer
    V(ii)%G1_i        => InflowGaugeNodeList_out
    !
    ! attributes 
    V(ii)%nAtt          = 2
    !
    V(ii)%att(2)%name   = "long_name"
    V(ii)%att(2)%xType  = NF90_CHAR
    V(ii)%att(2)%nValues= 1
    V(ii)%att(2)%values = "cell ID of inflow gauges"
    !
    V(ii)%att(1)%name   = "_FillValue"
    V(ii)%att(1)%xType  = NF90_INT
    V(ii)%att(1)%nValues= 1
    V(ii)%att(1)%values = trim(num2str(nodata_i4))
    !
  end subroutine set_L11_config
  !
  !     NAME
  !        set_config

  !     PURPOSE
  !        sets the netcdf output variable for all state variables and fluxes. These
  !        variables are defined in state_variables_init and the corresponding
  !        target variables are set in state_variables_set in mo_restart

  !     HISTORY
  !        author Stephan Thober
  !        date Jul 2013  
  !        modified Stephan Thober,     Nov 2013 - only set L11 variables if routing is switched on
  !
  subroutine set_state( L11_flag )
    !
    use mo_string_utils, only: num2str
    !
    implicit none
    !
    logical, intent(in)  :: L11_flag
    integer(i4)          :: i
    !
    !
    ! define parameters
    if ( L11_flag ) then
       ! nr. of dimensions and variables with routing
       nDims  = 6           
       nVars  = 53 + nDims  
    else
       ! nr. of dimensions and variables without routing
       NDims  = 3           
       nVars  = 44 + nDims
    end if
    !
    ! allocate arrays
    if ( allocated(Dnc)   ) deallocate ( Dnc)
    if ( allocated(V)     ) deallocate ( V  )
    allocate ( Dnc(nDims)  )
    allocate ( V(nVars)    )
    !
    allocate( x1( size( L1_fSealed_out, 1 ) ) )
    do i = 1, size(x1,1)
       x1(i) = i
    end do
    allocate( y1( size( L1_fSealed_out, 2 ) ) )
    do i = 1, size( y1,1)
       y1(i) = i
    end do
    !
    allocate( z1( size( L1_soilMoist_out, 3) ) )
    do i = 1, size( z1, 1)
       z1(i) = i
    end do
    !
    if ( L11_flag ) then
       allocate( x11( size( L11_Qmod_out, 1 ) ) )
       do i = 1, size(x11,1)
          x11(i) = i
       end do
       !
       allocate( y11( size( L11_Qmod_out, 2 ) ) )
       do i = 1, size( y11,1)
          y11(i) = i
       end do
       !
       allocate( z11( size( L11_qTR_out, 3) ) )
       do i = 1, size( z11,1)
          z11(i) = i
       end do
    end if
    !
    ! define dimensions --------------------------------------------------------
    Dnc(1)%name      = "nrows1"
    Dnc(1)%len       = size(x1,1)
    Dnc(2)%name      = "ncols1"
    Dnc(2)%len       = size(y1,1)
    Dnc(3)%name      = "L1_soilhorizons"
    Dnc(3)%len       = size(z1,1)
    if ( L11_flag ) then
       Dnc(4)%name      = "nrows11"
       Dnc(4)%len       = size(x11,1)
       Dnc(5)%name      = "ncols11"
       Dnc(5)%len       = size(y11,1)
       Dnc(6)%name      = "nIT"
       Dnc(6)%len       = size(z11,1)
    end if
    !
    ! DIMENSION VARIABLES ------------------------------------------------------
    !
    ! Number of Cells in one direction
    i                =  1
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/1,0,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer to actual data
    V(i)%G1_i        => x1
    !
    ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = ""
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "Number of rows at level 1"   
    !
    ! Number of Cells in one direction
    i                = i + 1_i4
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/2,0,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer to actual data
    V(i)%G1_i        => y1
    !
    ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = ""
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "Number of columns at level 1"   
    !
    ! Number of Cells in one direction
    i                = i + 1_i4
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/3,0,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer to actual data
    V(i)%G1_i        => z1
    !
    ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = ""
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "number of soil horizons at level 1"   
    !
    ! maximum Number of neighboring cells
    if ( L11_flag ) then
       i                = i + 1_i4
       V(i)%name        =  Dnc(i)%name
       V(i)%xType       =  NF90_INT
       V(i)%nLvls       =  0
       V(i)%nSubs       =  0
       V(i)%nDims       =  1
       V(i)%dimTypes    =  (/4,0,0,0,0/)
       V(i)%wFlag       = .true.
       !
       ! pointer to actual data
       V(i)%G1_i        => x11
       !
       ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "units"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values  = ""
       !
       V(i)%att(2)%name   = "long_name"
       V(i)%att(2)%xType  = NF90_CHAR
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = "Number of rows at level 11"
       !
       ! Number of Anchor cells
       i                = i + 1_i4
       V(i)%name        =  Dnc(i)%name
       V(i)%xType       =  NF90_INT
       V(i)%nLvls       =  0
       V(i)%nSubs       =  0
       V(i)%nDims       =  1
       V(i)%dimTypes    =  (/5,0,0,0,0/)
       V(i)%wFlag       = .true.
       !
       ! pointer to actual data
       V(i)%G1_i        => y11
       !
       ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "units"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values  = ""
       !
       V(i)%att(2)%name   = "long_name"
       V(i)%att(2)%xType  = NF90_CHAR
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = "Number of colums at level 11" 
       !
       ! Number of Anchor cells
       i                = i + 1_i4
       V(i)%name        =  Dnc(i)%name
       V(i)%xType       =  NF90_INT
       V(i)%nLvls       =  0
       V(i)%nSubs       =  0
       V(i)%nDims       =  1
       V(i)%dimTypes    =  (/6,0,0,0,0/)
       V(i)%wFlag       = .true.
       !
       ! pointer to actual data
       V(i)%G1_i        => z11
       !
       ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "units"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values  = ""
       !
       V(i)%att(2)%name   = "long_name"
       V(i)%att(2)%xType  = NF90_CHAR
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = "Number of timesteps routing is stored for" 
    end if
    !
    ! FIELD VARIABLES ----------------------------------------------------------

    ! provide global Land Cover variables at level 1
    i                = i + 1_i4
    V(i)%name        =  "L1_fSealed"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_fSealed_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "fraction of Sealed area at level 1"
    !
    V(i)%att(1)%name   = "_FillValue"
    V(i)%att(1)%xType  = NF90_DOUBLE
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = trim(num2str(nodata_dp))

    i                = i + 1_i4
    V(i)%name        =  "L1_fForest"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_fForest_out 
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "fraction of Forest area at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)

    i                = i + 1_i4
    V(i)%name        =  "L1_fPerm"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_fPerm_out   
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "fraction of permeable area at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)

    ! provide global state Variables at level 1
    i                = i + 1_i4
    V(i)%name        =  "L1_Inter"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_Inter_out   
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Interception storage at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_snowPack"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_snowPack_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Snowpack at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_sealSTW"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_sealSTW_out 
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Retention storage of impervious areas at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_soilMoist"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_soilMoist_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "soil moisture at level 1 at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_unsatSTW"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_unsatSTW_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "upper soil storage at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_satSTW"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_satSTW_out  
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "groundwater storage at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    

  ! provide global flux Variables at level 1
    i                = i + 1_i4
    V(i)%name        =  "L1_aETSoil"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_aETSoil_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "soil actual ET at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_aETCanopy"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_aETCanopy_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "canopy actual ET at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_aETSealed"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_aETSealed_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "sealed actual ET at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_baseflow"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_baseflow_out 
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "baseflow at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_infilSoil"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_infilSoil_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "soil in-exfiltration at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_fastRunoff"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_fastRunoff_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "fast runoff"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)

    i                = i + 1_i4
    V(i)%name        =  "L1_melt"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_melt_out     
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "snow melt at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    

    i                = i + 1_i4
    V(i)%name        =  "L1_percol"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_percol_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "percolation at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_preEffect"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_preEffect_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "effective precip. depth (snow melt + rain) at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_rain"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_rain_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "rain (liquid water) at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_runoffSeal"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_runoffSeal_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "runoff from impervious area at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_slowRunoff"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_slowRunoff_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "slow runoff  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_snow"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_snow_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "snow (solid water)  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_Throughfall"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_Throughfall_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "throughfall  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_total_runoff"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_total_runoff_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "throughfall  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    

  !-------------------------------------------
  ! EFFECTIVE PARAMETERS
  !-------------------------------------------
    i                = i + 1_i4
    V(i)%name        =  "L1_alpha"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_alpha_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "exponent for the upper reservoir at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_degDayInc"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_degDayInc_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "increase of the Degree-day factor per mm of increase in precipitation at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_degDayMax"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_degDayMax_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "maximum degree-day factor  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_degDayNoPre"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_degDayNoPre_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "degree-day factor with no precipitation at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_degDay"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_degDay_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "degree-day factor at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_karstLoss"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_karstLoss_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Karstic percolation loss at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_fAsp"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_fAsp_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "PET correction factor due to terrain aspect at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_fRoots"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_fRoots_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Fraction of roots in soil horizons at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_maxInter"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_maxInter_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Maximum interception  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_kfastFlow"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_kfastFlow_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "fast interflow recession coefficient  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_kSlowFlow"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_kSlowFlow_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "slow interflow recession coefficient  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_kBaseFlow"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_kBaseFlow_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "baseflow recession coefficient  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_kPerco"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_kPerco_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "percolation coefficient  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_soilMoistFC"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_soilMoistFC_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Soil moisture below which actual ET is reduced linearly till PWP at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_soilMoistSat"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_soilMoistSat_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Saturation soil moisture for each horizon [mm] at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_soilMoistExp"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_soilMoistExp_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Exponential parameter to how non-linear is the soil water retention at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_tempThresh"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_tempThresh_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Threshold temperature for snow/rain  at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_unsatThresh"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_unsatThresh_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Threshhold water depth controlling fast interflow at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_sealedThresh"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => L1_sealedThresh_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Threshhold water depth for surface runoff in sealed surfaces at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
    i                = i + 1_i4
    V(i)%name        =  "L1_wiltingPoint"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G3_d        => L1_wiltingPoint_out
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Permanent wilting point at level 1"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_DOUBLE
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = num2str(nodata_dp)
    
  !-------------------------------------------
  ! L11 ROUTING STATE VARIABLES, FLUXES AND
  !             PARAMETERS
  !-------------------------------------------
    if ( L11_flag ) then
       i                = i + 1_i4
       V(i)%name        =  "L11_Qmod"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_Qmod_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "simulated discharge at each node at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_qOUT"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_qOUT_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Total outflow from cells L11 at time tt at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_qTIN"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  3
       V(i)%dimTypes    =  (/1,2,3,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G3_d        => L11_qTIN_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Total discharge inputs at t-1 and t at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_qTR"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  3
       V(i)%dimTypes    =  (/1,2,3,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G3_d        => L11_qTR_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Routed outflow leaving a node at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_K"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_K_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "kappa: Muskingum travel time parameter. at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_xi"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_xi_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "xi: Muskingum diffusion parameter at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_C1"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_C1_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Routing parameter C1=f(K,xi, DT) (Chow, 25-41) at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_C2"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_C2_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Routing parameter C2 =f(K,xi, DT) (Chow, 25-41) at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)

       i                = i + 1_i4
       V(i)%name        =  "L11_FracFPimp"
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  2
       V(i)%dimTypes    =  (/4,5,0,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer  
       V(i)%G2_d        => L11_FracFPimp_out
       !
       ! attributes 
       V(i)%nAtt          = 2
       !
       V(i)%att(1)%name   = "long_name"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       V(i)%att(1)%values = "Fraction of the flood plain with impervious cover at level 11"
       !
       V(i)%att(2)%name   = "_FillValue"
       V(i)%att(2)%xType  = NF90_DOUBLE
       V(i)%att(2)%nValues= 1
       V(i)%att(2)%values = num2str(nodata_dp)
    end if
    
    ! File Attributes <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !
    Gatt(1)%name       = "Title"
    Gatt(1)%values     = "restart file for the state variables of mHM"
    !
    Gatt(2)%name       = "history"
    Gatt(2)%values     = ""
    !
  end subroutine  set_state
end Module mo_set_netcdf_restart
