!> \file mo_mrm_net_startup.f90

!> \brief Startup drainage network for mHM.

!> \details This module initializes the drainage network at L11 in mHM.\n
!>          - Delineation of drainage network at level 11.    \n
!>          - Setting network topology (i.e. nodes and link). \n
!>          - Determining routing order.                      \n
!>          - Determining cell locations for network links.   \n
!>          - Find drainage outlet.                           \n
!>          - Determine stream (links) features.              \n

!> \authors Luis Samaniego
!> \date Dec 2012
!         Modified
!         Rohini Kumar, May 2014   - cell area calulation based on a regular lat-lon grid or 
!                                     on a regular X-Y coordinate system

module mo_mrm_net_startup
  use mo_kind, only: i4, dp
  implicit none
  PUBLIC :: L11_variable_init
  PUBLIC :: L11_flow_direction
  PUBLIC :: L11_set_network_topology
  PUBLIC :: L11_routing_order
  PUBLIC :: L11_link_location
  PUBLIC :: L11_set_drain_outlet_gauges
  PUBLIC :: L11_stream_features
  PUBLIC :: L11_fraction_sealed_floodplain
  PUBLIC :: get_distance_two_lat_lon_points
contains
  ! --------------------------------------------------------------------------

  !     NAME
  !         L11_variable_init
  !     PURPOSE
  !>        \brief Cell numbering at ROUTING LEVEL-11

  !>        \details Cell numbering at ROUTING LEVEL-11  \n
  !>                 List of Level- 0 and 1 cells contained within a given Level-11 cell.\n

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    ::  iBasin"        Basin Id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None                 

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call L11_variable_init(1)

  !     LITERATURE
  !         None
  
  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Stephan Thober, Aug 2015 - ported to mRM
  !                  Stephan Thober, Sep 2015 - create L11 mask based on Level 0 not based on Level 1
  ! --------------------------------------------------------------------------

  subroutine L11_variable_init(iBasin)
    use mo_utils, only: ge
    use mo_mrm_constants, only: nodata_i4, nodata_dp
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm, calculate_grid_properties
    use mo_mrm_global_variables, only: &
         level1, &
         basin_mrm, &
         level11, resolutionRouting, &
         nBasins, &
         L0_areaCell,       &
         L1_Id,             &
         L1_L11_ID,         & ! INOUT: mapping of L1 Id on L11 in case L11 > L1
         L11_L1_ID,         & ! INOUT: mapping of L11 Id on L1 in case L1 > L11
         L11_areaCell,      & ! INOUT: cell area in [km2] at L11
         L11_cellCoor,      & ! cell coordinates (row,col)
         L11_nCells,        & ! Total No. of routing cells  (= nNodes)
         L11_Id,            & ! ids of grid at level-11    
         L11_upBound_L1,    & ! INOUT: row start at finer level-1 scale 
         L11_downBound_L1,  & ! INOUT: row end at finer level-1 scale 
         L11_leftBound_L1,  & ! INOUT: col start at finer level-1 scale 
         L11_rightBound_L1    ! INOUT: col end at finer level-1 scale 

    implicit none

    integer(i4), intent(in)                  :: iBasin      ! basin 

    ! local
    integer(i4)                              :: nrows0, ncols0
    real(dp)                                 :: xllcorner0, yllcorner0
    real(dp)                                 :: cellsize0
    !
    integer(i4)                              :: nrows1, ncols1
    integer(i4)                              :: nrows11, ncols11
    integer(i4)                              :: ncells
    integer(i4)                              :: iStart0, iEnd0
    integer(i4)                              :: iStart1, iEnd1
    integer(i4)                              :: iStartMask1, iEndMask1
    logical,     dimension(:,:), allocatable :: mask0, mask1, mask11
    integer(i4), dimension(:),   allocatable :: upBound1, downBound1, leftBound1, rightBound1 
    integer(i4), dimension(:,:), allocatable :: cellCoor
    integer(i4)                              :: kk
    integer(i4)                              :: ic, jc, icc, jcc
    integer(i4)                              :: iu, id, jl, jr
    integer(i4), dimension(:,:), allocatable :: L11Id_on_L1 ! mapping of L11 Id on L1
    integer(i4), dimension(:,:), allocatable :: L1Id_on_L11 ! mapping of L1 Id on L11
    integer(i4), dimension(:,:), allocatable :: Id11        ! ids of grid at level-11     
    integer(i4), dimension(:,:), allocatable :: Id1         ! ids of grid at level-1
    real(dp), dimension(:,:), allocatable    :: areaCell0_2D
    real(dp), dimension(:), allocatable      :: areaCell
    real(dp)                                 :: cellFactorR
    real(dp)                                 :: cellFactorRbyH
    integer(i4)                              :: cellFactorRbyH_inv
    !--------------------------------------------------------
    ! STEPS::
    ! 1) Estimate each variable locally for a given basin
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! level-0 information
    call get_basin_info_mrm( iBasin, 0, nrows0, ncols0, iStart=iStart0, iEnd=iEnd0, &
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0, mask=mask0) 

    if(iBasin == 1) then
       ! allocate
       allocate( level11%nrows        (nBasins) )
       allocate( level11%ncols        (nBasins) )
       allocate( level11%xllcorner    (nBasins) )
       allocate( level11%yllcorner    (nBasins) )
       allocate( level11%cellsize     (nBasins) )
       allocate( level11%nodata_value (nBasins) )

       ! initialize
       level11%nrows(:)        = nodata_i4   
       level11%ncols(:)        = nodata_i4
       level11%xllcorner(:)    = nodata_dp
       level11%yllcorner(:)    = nodata_dp
       level11%cellsize(:)     = nodata_dp
       level11%nodata_value(:) = nodata_dp
    end if

    ! grid information
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp, &
         resolutionRouting(iBasin) ,                                                              &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin),                 &
         level11%yllcorner(iBasin), level11%cellsize(iBasin), level11%nodata_value(iBasin) )
    ! level-1 information
    call get_basin_info_mrm (iBasin, 1, nrows1, ncols1, iStart=iStart1, iEnd=iEnd1,                   &
         iStartMask=iStartMask1, iEndMask=iEndMask1, mask=mask1 ) 

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11) 

    ! allocate and initialize
    allocate( mask11(nrows11, ncols11) )
    mask11(:,:) = .FALSE.

    cellFactorR    = level11%cellsize(iBasin) / cellsize0
    ! create a mask: Id with respect to Level 0
    ! create a mask at Level 11: Id with respect to Level 0
    do jc = 1, ncols0
       jcc = ceiling(real(jc, dp) / cellFactorR)
       do ic = 1, nrows0
          if (.NOT. mask0(ic, jc)) cycle
          icc = ceiling(real(ic,dp) / cellFactorR)
          mask11(icc, jcc) = .TRUE.
       end do
    end do

    ! set number of cells (equals nNodes)
    ncells = count( mask11 )

    ! initialize areacell
    allocate ( areacell(ncells) )
    areacell = nodata_i4
    ! allocate bounds
    allocate ( upBound1    (ncells) )
    allocate ( downBound1  (ncells) )
    allocate ( leftBound1  (ncells) )
    allocate ( rightBound1 (ncells) )
    upBound1(:)    = nodata_i4   
    downBound1(:)  = nodata_i4 
    leftBound1(:)  = nodata_i4 
    rightBound1(:) = nodata_i4 

    ! allocate variables for mapping L11 Ids and L1 Ids
    allocate ( L11Id_on_L1  (nrows1, ncols1 ) )
    allocate ( L1Id_on_L11  (nrows11, ncols11 ) )
    allocate ( Id11         (nrows11, ncols11 ) )
    allocate ( Id1          (nrows1, ncols1 ) )
    L11Id_on_L1(:,:) = nodata_i4
    L1Id_on_L11(:,:) = nodata_i4
    Id11(:,:)        = nodata_i4
    Id1(:,:)         = nodata_i4

    ! allocate
    allocate ( cellCoor(nCells,2) )

    ! initialize
    cellCoor(:,:) = nodata_i4

    ! counting valid cells at level 11
    kk = 0
    do jcc = 1, ncols11
       do icc = 1, nrows11
          if ( .not. mask11(icc,jcc) ) cycle
          kk = kk + 1
          Id11(icc,jcc)  = kk
          cellCoor(kk,1) = icc
          cellCoor(kk,2) = jcc
       end do
    end do

    !--------------------------------------------------------
    ! UPDATE BASIN_MRM VARIABLE
    !--------------------------------------------------------
    if(iBasin == 1) then

       ! allocate
       allocate(basin_mrm%L11_iStart     (nBasins))
       allocate(basin_mrm%L11_iEnd       (nBasins))
       allocate(basin_mrm%L11_iStartMask (nBasins))
       allocate(basin_mrm%L11_iEndMask   (nBasins))

       ! initialize   
       basin_mrm%L11_iStart(:)     = nodata_i4  
       basin_mrm%L11_iEnd(:)       = nodata_i4 
       basin_mrm%L11_iStartMask(:) = nodata_i4
       basin_mrm%L11_iEndMask(:)   = nodata_i4

       ! basin information
       basin_mrm%L11_iStart(iBasin) = 1
       basin_mrm%L11_iEnd  (iBasin) = basin_mrm%L11_iStart(iBasin) + nCells - 1

       basin_mrm%L11_iStartMask(iBasin) = 1
       basin_mrm%L11_iEndMask  (iBasin) = basin_mrm%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    else

       ! basin information
       basin_mrm%L11_iStart(iBasin) = basin_mrm%L11_iEnd(iBasin-1) + 1
       basin_mrm%L11_iEnd  (iBasin) = basin_mrm%L11_iStart(iBasin) + nCells - 1

       basin_mrm%L11_iStartMask(iBasin) = basin_mrm%L11_iEndMask(iBasin-1) + 1
       basin_mrm%L11_iEndMask  (iBasin) = basin_mrm%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    end if

    !--------------------------------------------------------
    ! CALCULATE L11_AREACELL AND CELL ID MAPPING WITH L1
    !--------------------------------------------------------
    ! level-0 cell area
    allocate( areaCell0_2D(nrows0,ncols0) )
    areaCell0_2D(:,:) = UNPACK( L0_areaCell(iStart0:iEnd0), mask0, nodata_dp )

    ! set cell factor for routing
    cellFactorRbyH = level11%cellsize(iBasin) / level1%cellsize(iBasin)
    cellFactorR    = level11%cellsize(iBasin) / cellsize0
    
    kk = 0
    do jcc = 1, ncols11
       do icc = 1, nrows11
          if( .not. mask11(icc,jcc)) cycle
          kk = kk + 1

          ! coord. of all corners L11 -> of finer scale level-0
          iu = (icc-1) * nint(cellFactorR,i4) + 1
          id =     icc * nint(cellFactorR,i4)
          jl = (jcc-1) * nint(cellFactorR,i4) + 1
          jr =     jcc * nint(cellFactorR,i4)
          ! constrain the range of up, down, left, and right boundaries
          if(iu < 1     ) iu = 1
          if(id > nrows0) id = nrows0
          if(jl < 1     ) jl = 1
          if(jr > ncols0) jr = ncols0

          ! effective area [km2] & total no. of L0 cells within a given L1 cell
          areaCell(kk) = sum( areacell0_2D(iu:id, jl:jr), mask0(iu:id, jl:jr) )*1.0E-6

          ! coord. of all corners L11 -> of finer scale level-1
          iu = (icc-1) * nint(cellFactorRbyH,i4) + 1
          id =     icc * nint(cellFactorRbyH,i4)
          jl = (jcc-1) * nint(cellFactorRbyH,i4) + 1
          jr =     jcc * nint(cellFactorRbyH,i4)

          ! constrain the range of up, down, left, and right boundaries
          if( iu < 1   ) iu =  1
          if( id > nrows1 ) id =  nrows1
          if( jl < 1   ) jl =  1
          if( jr > ncols1 ) jr =  ncols1
          
          upBound1   (kk) = iu
          downBound1 (kk) = id
          leftBound1 (kk) = jl
          rightBound1(kk) = jr

          ! set mapping
          if (ge(cellFactorRbyH, 1._dp)) then 
             ! Delimitation of level-11 cells on level-1 for L11 resolution lower than L1 resolution
             L11Id_on_L1(iu:id, jl:jr) = Id11(icc, jcc)
          end if
       end do
    end do

    ! create mapping between L11 and L1 for L11 resolution higher than L1 resolution
    if (cellFactorRbyH .lt. 1._dp) then
       cellFactorRbyH_inv = int(1. / cellFactorRbyH, i4)
       kk = 0
       do jcc = 1, ncols1
          do icc = 1, nrows1
             if( .not. mask1(icc,jcc)) cycle
             kk = kk + 1
             !
             iu = (icc - 1) * cellFactorRbyH_inv + 1
             id =       icc * cellFactorRbyH_inv
             jl = (jcc - 1) * cellFactorRbyH_inv + 1
             jr =       jcc * cellFactorRbyH_inv
             !
             Id1(icc,jcc) = kk
             L1Id_on_L11(iu:id, jl:jr) = merge(Id1(icc,jcc), nodata_i4, mask11(iu:id, jl:jr))
          end do
       end do
    end if


    ! L1 data sets
    call append( L1_id, pack ( Id1(:,:), mask1 ) )
    call append( L1_L11_Id, pack ( L11Id_on_L1(:,:), mask1 ) )

    ! L11 data sets
    call append( L11_L1_Id, PACK ( L1Id_on_L11(:,:), mask11)  )
    call append( basin_mrm%L11_Mask,  RESHAPE( mask11, (/nrows11*ncols11/)  )  )
    ! other L11 data sets
    call append( L11_cellCoor, cellCoor )
    call append( L11_Id, pack(Id11, mask11) )
    call append( L11_areaCell, areacell)
    call append( L11_upBound_L1, upBound1(:) )
    call append( L11_downBound_L1, downBound1(:) )
    call append( L11_leftBound_L1, leftBound1(:) )
    call append( L11_rightBound_L1, rightBound1(:) )

    L11_nCells = size( L11_Id, 1 )

    ! free space
    deallocate(Id11, mask1, mask11, cellCoor, areacell, &
         upBound1, downBound1, leftBound1, rightBound1, &
         L11Id_on_L1, L1Id_on_L11)
    
  end subroutine L11_variable_init

  ! --------------------------------------------------------------------------

  !     NAME
  !         L11_flow_direction

  !     PURPOSE

  !>       \brief Determine the flow direction of the upscaled river
  !>    network at level L11.

  !>       \details The hydrographs generated at each cell are routed
  !>                through the drainage network at level-11 towards their 
  !>                outlets. The drainage network at level-11 is conceptualized as a
  !>                graph whose nodes are hypothetically located at the center of
  !>                each grid cell connected by links that represent the river
  !>                reaches. The flow direction of a link correspond to the
  !>                direction towards a neighboring cell in which the net flow
  !>                accumulation (outflows minus inflows) attains its maximum
  !>                value. The net flow accumulation across a cell's boundary at
  !>                level-11 is estimated based on flow direction and flow
  !>                accumulation obtained at level-0 (\ref fig_routing "Routing
  !>                Network"). Note: level-1 denotes the modeling level, whereas
  !>                level-L11 is at least as coarse as level-1. Experience has
  !>                shown that routing can be done at a coarser resolution as
  !>                level-1, hence the level-11 was introduced.

  !>                 \image html  routing.png "Upscaling routing network from L0 to L1 (or L11)"
  !>                \anchor fig_routing \image latex routing.pdf "Upscaling routing network from L0 to L1 (or L11)" width=14cm

  !>                The left panel depicts a schematic derivation of a drainage
  !>                network at the level-11 based on level-0 flow direction and
  !>                flow accumulation. The dotted line circle denotes the point
  !>                with the highest flow accumulation within a grid cell. The
  !>                topology of a tipical drainage routing network at level-11 is
  !>                shown in the right panel. Gray color areas denote the flood
  !>                plains estimated in mo_init_mrm, where the network
  !>                upscaling is also carried out.

  !>                For the sake of simplicity, it is assumed that all runoff leaving
  !>                a given cell would exit through a major direction.

  !>                Note that multiple outlets can exist within the modelling domain.

  !>                If a variable is added or removed here, then it also has to 
  !>                be added or removed in the subroutine L11_config_set in
  !>                module mo_restart and in the subroutine set_L11_config in module
  !>                mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"             Basin Id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Rohini Kumar,   Apr 2014 - Case of L0 is same as L11 implemented
  !                  Stephan Thober, Aug 2015 - ported to mRM
  !                  Stephan Thober, Sep 2015 - create mapping between L11 and L1 if L11 resolution
  !                                             is higher than L1 resolution
  !                  Stephan Thober, May 2016 - introducing multiple outlets
  ! --------------------------------------------------------------------------
  subroutine L11_flow_direction(iBasin)
    use mo_append, only: append
    use mo_message, only: message
    use mo_mrm_constants, only: nodata_i4
    use mo_mrm_global_variables, only: &
         basin_mrm, &
         level0, &
         nBasins, &
         level11,   &
         L0_fAcc, L0_fDir,  &
         L0_draSC,          & ! INOUT: draining cell of each sub catchment (== cell L11)
         L0_cellCoor,       &
         L0_id,             &
         L0_L11_Id,         & ! INOUT: mapping of L11 Id on L0
         L11_Id,            &
         L11_cellCoor,      &
         L11_rowOut,        & ! INOUT: grid vertical location of the Outlet
         L11_colOut,        & ! INOUT: grid horizontal location  of the Outlet
         L11_fDir,          & ! INOUT: flow direction at L11 (standard notation)
         L11_upBound_L0,    & ! INOUT: row start at finer level-0 scale 
         L11_downBound_L0,  & ! INOUT: row end at finer level-0 scale 
         L11_leftBound_L0,  & ! INOUT: col start at finer level-0 scale 
         L11_rightBound_L0, & ! INOUT: col end at finer level-0 scale
         L11_nOutlets
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_string_utils, only: num2str

    implicit none

    integer(i4), intent(in)                  :: iBasin         ! basin 

    ! local
    integer(i4)                              :: nCells0
    integer(i4)                              :: nrows0, ncols0
    integer(i4)                              :: nrows11, ncols11
    integer(i4)                              :: nNodes      ! =  ncells11
    integer(i4)                              :: iStart0, iEnd0
    integer(i4)                              :: iStart11, iEnd11
    logical,     dimension(:,:), allocatable :: mask0, mask11
    integer(i4), dimension(:),   allocatable :: upBound0, downBound0, leftBound0, rightBound0 
    real(dp)                                 :: cellFactorR
    integer(i4)                              :: icc, jcc
    integer(i4)                              :: ii, jj, kk, ic, jc 
    integer(i4)                              :: iu, id
    integer(i4)                              :: jl, jr
    integer(i4)                              :: iRow, jCol
    integer(i4), dimension(:,:), allocatable :: L11Id_on_L0 ! mapping of L11 Id on L0
    integer(i4), dimension(:,:), allocatable :: iD0            
    integer(i4), dimension(:,:), allocatable :: fDir0          
    integer(i4), dimension(:,:), allocatable :: fAcc0          
    integer(i4), dimension(:,:), allocatable :: Id11        ! ids of grid at level-11     
    integer(i4), dimension(:,:), allocatable :: fDir11         
    integer(i4), dimension(:,:), allocatable :: cellCoor0
    integer(i4), dimension(:,:), allocatable :: cellCoor11
    integer(i4), dimension(:),   allocatable :: rowOut      ! northing cell loc. of the Outlet
    integer(i4), dimension(:),   allocatable :: colOut      ! easting cell loc. of the Outlet
    integer(i4), dimension(:,:), allocatable :: draSC0         
    integer(i4), dimension(:,:), allocatable :: oLoc        ! output location in L0
    integer(i4)                              :: side
    integer(i4)                              :: fAccMax, idMax
    integer(i4)                              :: Noutlet     ! Number of outlet found
    integer(i4)                              :: old_Noutlet ! Number of outlets before this basin
    integer(i4), dimension(:,:), allocatable :: dummy       ! helping variable for storing L0 outlet coordinates
    logical                                  :: is_outlet   ! flag whether outlet is found

    !--------------------------------------------------------
    ! STEPS:
    ! 1) Estimate each variable locally for a given basin
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------

    ! level-0 information
    call get_basin_info_mrm (iBasin, 0, nrows0, ncols0, ncells=nCells0,   &
         iStart=iStart0, iEnd=iEnd0, mask=mask0) 

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11, ncells=nNodes, &
         iStart=iStart11, iEnd=iEnd11, mask=mask11)

    ! set cell factors
    cellFactorR = level11%cellsize(iBasin) / level0%cellsize(iBasin)

    ! allocate
    allocate ( upBound0    (nNodes) )
    allocate ( downBound0  (nNodes) )
    allocate ( leftBound0  (nNodes) )
    allocate ( rightBound0 (nNodes) )

    allocate ( L11Id_on_L0  (nrows0, ncols0 ) )
    allocate ( Id11(nrows11, ncols11) )

    ! initialize
    Noutlet        = 0_i4
    upBound0(:)    = nodata_i4
    downBound0(:)  = nodata_i4
    leftBound0(:)  = nodata_i4
    rightBound0(:) = nodata_i4

    L11Id_on_L0(:,:) = nodata_i4
    Id11(:,:)        = nodata_i4
    
    ! get Ids of L11 
    Id11(:,:) = unpack( L11_Id(iStart11:iEnd11),  mask11, nodata_i4 )

    kk   = 0
    do jcc = 1, ncols11
       do icc = 1, nrows11
          if( .not. mask11(icc,jcc)) cycle
          kk = kk + 1

          ! coord. of all corners L11 -> on finer scale level-0
          iu = (icc-1) * nint(cellFactorR,i4) + 1
          id =    icc  * nint(cellFactorR,i4)
          jl = (jcc-1) * nint(cellFactorR,i4) + 1
          jr =    jcc  * nint(cellFactorR,i4)

          ! constrain the range of up, down, left, and right boundaries
          if( iu < 1   ) iu =  1
          if( id > nrows0 ) id =  nrows0
          if( jl < 1   ) jl =  1
          if( jr > ncols0 ) jr =  ncols0

          upBound0   (kk) = iu
          downBound0 (kk) = id
          leftBound0 (kk) = jl
          rightBound0(kk) = jr

          ! Delimitation of level-11 cells on level-0
          L11Id_on_L0(iu:id, jl:jr) = Id11(icc, jcc)

       end do
    end do

    !------------------------------------------------------------------
    !                Set Flow Direction at Level 11
    !                       Searching order
    !                             jl    jr
    !                          iu  +  4  +
    !                              3     1
    !                          id  +  2  +
    !------------------------------------------------------------------

    ! flow direction at level-11

    ! allocate
    allocate ( iD0         ( nrows0, ncols0 ) )
    allocate ( fAcc0       ( nrows0, ncols0 ) )
    allocate ( fDir0       ( nrows0, ncols0 ) )
    allocate ( draSC0      ( nrows0, ncols0 ) )
    allocate ( cellCoor0   ( nCells0, 2 ) )  
    allocate ( cellCoor11  ( nNodes,  2 ) )  
    allocate ( fDir11      ( nrows11, ncols11 ) )
    allocate ( rowOut      ( nNodes ) )
    allocate ( colOut      ( nNodes ) )
    allocate ( oLoc        ( 1, 2 ) )

    ! initialize
    iD0(:,:)        = nodata_i4   
    fAcc0(:,:)      = nodata_i4
    fDir0(:,:)      = nodata_i4
    draSC0(:,:)     = nodata_i4    
    cellCoor0(:,:)  = nodata_i4
    cellCoor11(:,:) = nodata_i4
    fDir11(:,:)     = nodata_i4 
    rowOut(:)       = nodata_i4 
    colOut(:)       = nodata_i4
    oLoc(:,:)       = nodata_i4

    ! get iD, fAcc, fDir at L0
    iD0(:,:)   = UNPACK( L0_Id   (iStart0:iEnd0),  mask0, nodata_i4 )
    fAcc0(:,:) = UNPACK( L0_fAcc (iStart0:iEnd0),  mask0, nodata_i4 )
    fDir0(:,:) = UNPACK( L0_fDir (iStart0:iEnd0),  mask0, nodata_i4 )

    cellCoor0(:,:)  = L0_cellCoor  (iStart0 : iEnd0,  :)
    cellCoor11(:,:) = L11_cellCoor (iStart11: iEnd11, :)

    ! case where routing and input data scale is similar
    IF(nCells0 .EQ. nNodes) THEN
      oLoc(1, :) = maxloc( fAcc0, mask0 )
      kk   = L11Id_on_L0( oLoc(1, 1), oLoc(1, 2) )
      ! for a single node model run
      if(nCells0 .EQ. 1) then
       fDir11(1,1) = fDir0(oLoc(1, 1), oLoc(1, 2)) 
      else
        fDir11(:,:) = fDir0(:,:)
      end if
      fDir11 ( cellCoor11(kk,1), cellCoor11(kk,2) ) = 0
      ! set location of main outlet in L11
      do kk = 1, nNodes
         ii = cellCoor11( kk, 1 )
         jj = cellCoor11( kk, 2 )
         rowOut(kk) = ii
         colOut(kk) = jj
      end do
      do kk = 1, ncells0 
         ii = cellCoor0( kk, 1 )
         jj = cellCoor0( kk, 2 )
         draSC0(ii, jj) = kk
      end do

      ! case where routing and input data scale differs 
   ELSE
      ! =======================================================================
      ! ST: find all cells whose downstream cells are outside the domain
      ! =======================================================================
      do ii = 1, nCells0
         iRow = cellCoor0(ii, 1)
         jCol = cellCoor0(ii, 2)
         call moveDownOneCell(fDir0(iRow, jCol), iRow, jCol)
         ! check whether new location is inside bound
         is_outlet = .False.
         if ((iRow .le. 0_i4) .or. (iRow .gt. nrows0) .or. &
             (jCol .le. 0_i4) .or. (jCol .gt. ncols0)) then
            is_outlet = .True.
         else
            if (fdir0(iRow, jCol) .lt. 0) is_outlet = .True.
         end if
         !
         if (is_outlet) then
            Noutlet = Noutlet + 1_i4
            ! cell is an outlet
            if (Noutlet .eq. 1) then
               oLoc(1, :) = cellCoor0(ii, :)
            else
               call append(oLoc, cellCoor0(ii:ii, :))
            end if
            ! drain this cell into corresponding L11 cell
            kk   = L11Id_on_L0(oLoc(Noutlet, 1), oLoc(Noutlet, 2))
            draSC0(oLoc(Noutlet, 1), oLoc(Noutlet, 2)) = kk
            ! check whether cell has maximum flow accumulation
            ! coord. of all corners
            iu = upBound0   (kk)
            id = downBound0 (kk)
            jl = leftBound0 (kk)
            jr = rightBound0(kk)
            if (maxval(facc0(iu: id, jl: jr)) .eq. facc0(oLoc(Noutlet, 1), oLoc(Noutlet, 2)))  then
               ! set location of outlet at L11
               rowOut(kk) = oLoc(Noutlet, 1)
               colOut(kk) = oLoc(Noutlet, 2)
               fdir11(cellCoor11(kk,1), cellCoor11(kk,2)) = 0
            end if
         end if
      end do
      
      ! finding cell L11 outlets -  using L0_fAcc

      do kk = 1, nNodes

         ! exclude outlet L11
         if ( rowOut(kk) > 0 ) cycle

         ic = cellCoor11(kk,1)
         jc = cellCoor11(kk,2)

         ! coord. of all corners
         iu = upBound0   (kk)
         id = downBound0 (kk)
         jl = leftBound0 (kk)
         jr = rightBound0(kk)

         fAccMax = -9
         idMax   =  0
         side    = -1
         ! searching on side 4
         do jj = jl,jr
            if ( ( fAcc0(iu,jj) > fAccMax       )  .and. &
                 ( fDir0(iu,jj) ==  32 .or.  &
                 fDir0(iu,jj) ==  64 .or.  &
                 fDir0(iu,jj) == 128       )        ) then
               fAccMax = fAcc0(iu,jj)
               idMax   =   id0(iu,jj)
               side    = 4
            end if
         end do

         ! searching on side 1
         do ii = iu,id
            if ( ( fAcc0(ii,jr) > fAccMax       )  .and. &
                 ( fDir0(ii,jr) ==   1 .or.  &
                 fDir0(ii,jr) ==   2 .or.  &
                 fDir0(ii,jr) == 128       )        ) then
               fAccMax = fAcc0(ii,jr)
               idMax   =   id0(ii,jr)
               side    = 1
            end if
         end do

         ! searching on side 2
         do jj = jl,jr
            if ( ( fAcc0(id,jj) > fAccMax       )  .and. &
                 ( fDir0(id,jj) ==   2 .or.  &
                 fDir0(id,jj) ==   4 .or.  &
                 fDir0(id,jj) ==   8       )        ) then
               fAccMax = fAcc0(id,jj)
               idMax   =   id0(id,jj)
               side    = 2
            end if
         end do

         ! searching on side 3
         do ii = iu,id  
            if ( ( fAcc0(ii,jl) > fAccMax       )  .and. &
                 ( fDir0(ii,jl) ==   8 .or.  &
                 fDir0(ii,jl) ==  16 .or.  &
                 fDir0(ii,jl) ==  32       )        ) then
               fAccMax = fAcc0(ii,jl)
               idMax   =   id0(ii,jl)
               side    = 3
            end if
         end do

         ! set location of the cell-outlet (row, col) in L0
         ii = cellCoor0( idMax, 1 )
         jj = cellCoor0( idMax, 2 )
         rowOut(kk) = ii
         colOut(kk) = jj
         draSC0(ii,jj) = kk

         ! set fDir at L11
         if     ( ii == iu .and.  jj == jl ) then
            select case ( fDir0(ii,jj) )
            case (8,16)
               fDir11(ic,jc) = 16
            case (32)
               fDir11(ic,jc) = 32
            case (64,128)
               fDir11(ic,jc) = 64
            end select
         elseif ( ii == iu .and.  jj == jr ) then
            select case ( fDir0(ii,jj) )
            case (32,64)
               fDir11(ic,jc) = 64
            case (128)
               fDir11(ic,jc) = 128
            case (1,2)
               fDir11(ic,jc) = 1
            end select
         elseif ( ii == id .and.  jj == jl ) then
            select case ( fDir0(ii,jj) )
            case (2,4)
               fDir11(ic,jc) = 4
            case (8)
               fDir11(ic,jc) = 8
            case (16,32)
               fDir11(ic,jc) = 16
            end select
         elseif ( ii == id .and.  jj == jr ) then
            select case ( fDir0(ii,jj) )
            case (128,1)
               fDir11(ic,jc) = 1
            case (2)
               fDir11(ic,jc) = 2
            case (4,8)
               fDir11(ic,jc) = 4
            end select
         else
            ! cell on one side
            select case (side)
            case (1)
               fDir11(ic,jc) = 1
            case (2)
               fDir11(ic,jc) = 4
            case (3)
               fDir11(ic,jc) = 16
            case (4)
               fDir11(ic,jc) = 64
            case default
               stop 'Error L11_flow_direction: side = -1'
            end select
         endif

      end do
      
   END IF
   !--------------------------------------------------------
   ! Start padding up local variables to global variables
   !--------------------------------------------------------
   
   ! allocate space for row and col Outlet
   if(iBasin .eq. 1) then
      allocate( basin_mrm%L0_Noutlet(nBasins) )
      allocate( basin_mrm%L0_rowOutlet(1, nBasins) ) 
      allocate( basin_mrm%L0_colOutlet(1, nBasins) )
      basin_mrm%L0_Noutlet = nodata_i4
      basin_mrm%L0_rowOutlet = nodata_i4
      basin_mrm%L0_colOutlet = nodata_i4
   end if

   ! L0 data sets
   call append( L0_draSC, PACK ( draSC0(:,:),  mask0)  ) 
   call append( L0_L11_Id, PACK ( L11Id_on_L0(:,:), mask0)  )
   basin_mrm%L0_Noutlet(iBasin) = Noutlet
   ! set L0 outlet coordinates
   old_Noutlet = size(basin_mrm%L0_rowOutlet, dim=1)
   if (Noutlet .le. old_Noutlet) then
      basin_mrm%L0_rowOutlet(:Noutlet, iBasin) = oLoc(:, 1)
      basin_mrm%L0_colOutlet(:Noutlet, iBasin) = oLoc(:, 2)
   else
      ! store up to size of old_Noutlet
      basin_mrm%L0_rowOutlet(:old_Noutlet, iBasin) = oLoc(:old_Noutlet, 1)
      basin_mrm%L0_colOutlet(:old_Noutlet, iBasin) = oLoc(:old_Noutlet, 2)
      ! enlarge rowOutlet and colOutlet in basin_mrm structure
      allocate(dummy(Noutlet - old_Noutlet, nBasins))
      dummy = nodata_i4
      dummy(:, iBasin) = oLoc(old_Noutlet + 1:, 1)
      call append(basin_mrm%L0_rowOutlet, dummy)
      dummy(:, iBasin) = oLoc(old_Noutlet + 1:, 2)
      call append(basin_mrm%L0_colOutlet, dummy)
      deallocate(dummy)
   end if
   
   ! L11 data sets
   call append( L11_nOutlets,          count(fdir11 .eq. 0_i4) )
   call append( L11_fDir,     PACK ( fDir11(:,:),      mask11) )
   call append( L11_rowOut          ,  rowOut(:)               )
   call append( L11_colOut          ,  colOut(:)               )
   call append( L11_upBound_L0      ,  upBound0(:)             )
   call append( L11_downBound_L0    ,  downBound0(:)           )
   call append( L11_leftBound_L0    ,  leftBound0(:)           )
   call append( L11_rightBound_L0   ,  rightBound0(:)          )

   ! communicate
   call message('      Number of outlets found at Level 0:.. '//num2str(Noutlet, '(i7)'))
   call message('      Number of outlets found at Level 11:. '//num2str(count(fdir11 .eq. 0_i4), '(i7)'))

   ! free space
   deallocate(mask0, mask11, &
        upBound0, downBound0, leftBound0, rightBound0, & 
        L11Id_on_L0, Id11,                             &
        iD0, fDir0, fAcc0, fDir11, cellCoor0,          &
        cellCoor11, rowOut, colOut, draSC0         )   
   
 end subroutine L11_flow_direction

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_set_network_topology

  !     PURPOSE
  !>        \brief Set network topology

  !>        \details Set network topology from and to node for all links
  !>                 at level-11 (\ref fig_routing "Routing Network"). \n

  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"             Basin Id

  !     INTENT(INOUT)
  !         None           

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Stephan Thober, Aug 2015 - ported to mRM
  !                  Stephan Thober, May 2016 - moved calculation of sink here
  ! ------------------------------------------------------------------

  subroutine L11_set_network_topology(iBasin)
    use mo_mrm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         L11_Id, L11_cellCoor,     &
         L11_fDir,                 &
         L11_fromN,                & ! INOUT: from node 
         L11_toN                     ! INOUT: to node

    implicit none

    integer(i4), intent(in)                   :: iBasin         ! basin 

    ! local
    integer(i4)                               :: nNodes
    integer(i4)                               :: nrows11, ncols11
    integer(i4)                               :: iStart11, iEnd11
    logical,     dimension(:,:), allocatable  :: mask11
    integer(i4), dimension(:,:), allocatable  :: cellCoor11
    integer(i4), dimension(:,:), allocatable  :: Id11           ! ids of grid at level-11     
    integer(i4), dimension(:,:), allocatable  :: fDir11         
    integer(i4)                               :: jj, kk, ic, jc 
    integer(i4)                               :: fn, tn

    integer(i4), dimension(:), allocatable    :: nLinkFromN, nLinkToN 

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11, ncells=nNodes, &
         iStart=iStart11, iEnd=iEnd11, mask=mask11)

    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a basin. 

    ! allocate
    allocate ( nLinkFromN ( nNodes    ) )  ! valid from (1 : nLinks)
    allocate ( nLinkToN   ( nNodes    ) )  ! "
    allocate ( cellCoor11 ( nNodes, 2 ) )  
    allocate ( Id11       ( nrows11, ncols11 ) )
    allocate ( fDir11     ( nrows11, ncols11 ) )

    ! initialize
    nLinkFromN(:)   = nodata_i4
    nLinkToN(:)     = nodata_i4
    cellCoor11(:,:) = nodata_i4
    Id11(:,:)       = nodata_i4
    fDir11(:,:)     = nodata_i4

    ! get grids of L11 
    Id11(:,:) =    UNPACK( L11_Id   ( iStart11 : iEnd11),  mask11, nodata_i4 )
    fDir11(:,:) =  UNPACK( L11_fDir ( iStart11 : iEnd11),  mask11, nodata_i4 )
    cellCoor11(:,:) = L11_cellCoor ( iStart11 : iEnd11, : )

    ! ------------------------------------------------------------------
    !  network topology
    ! ------------------------------------------------------------------

    jj = 0
    do kk = 1, nNodes
       ic = cellCoor11(kk,1)
       jc = cellCoor11(kk,2)
       fn = kk
       call moveDownOneCell(fDir11(ic, jc), ic, jc)
       tn = Id11(ic,jc)
       if (fn == tn) cycle 
       jj = jj + 1
       nLinkFromN(jj) = fn
       nLinkToN(jj)   = tn
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L11 data sets
    call append( L11_fromN, nLinkFromN(:) ) ! sinks are at the end 
    call append( L11_toN,   nLinkToN(:)   )

    ! free space
    deallocate (mask11, cellCoor11, Id11, fDir11, nLinkFromN, nLinkToN)   

  end subroutine L11_set_network_topology

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_routing_order

  !     PURPOSE
  !>        \brief Find routing order, headwater cells and sink

  !>        \details Find routing order, headwater cells and sink. \n
  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"             Basin Id         

  !     INTENT(INOUT)
  !         None            

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None
  
  !     EXAMPLE
  !         None
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Sa. Ku.         Jan 2015 - corrected initialization of nLinkSink
  !                  Stephan Thober, Aug 2015 - ported to mRM
  ! ------------------------------------------------------------------

  subroutine L11_routing_order(iBasin)
    use mo_mrm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         L11_fromN,                & ! IN:    from node 
         L11_toN,                  & ! IN:    to node
         L11_fDir,                 & ! IN:    flow direction to identify sink
         L11_nOutlets,             & ! IN:    number of sinks/outlets
         L11_sink,                 & ! IN: == .true. if sink node reached
         L11_rOrder,               & ! INOUT: network routing order
         L11_label,                & ! INOUT: label Id [0='', 1=HeadWater, 2=Sink]
         L11_netPerm                 ! INOUT: routing order (permutation)

    implicit none

    integer(i4), intent(in)                   :: iBasin         ! basin 

    ! local
    integer(i4)                               :: nNodes
    integer(i4)                               :: nLinks
    integer(i4)                               :: nrows11, ncols11
    integer(i4)                               :: iStart11, iEnd11
    integer(i4), dimension(:), allocatable    :: nLinkFromN      ! from node
    integer(i4), dimension(:), allocatable    :: nLinkToN        ! to node
    integer(i4), dimension(:), allocatable    :: nLinkROrder     ! network routing order  
    integer(i4), dimension(:), allocatable    :: nLinkLabel      ! label Id [0='', 1=HeadWater, 2=Sink]
    logical,     dimension(:), allocatable    :: nLinkSink       ! == .true. if sink node reached
    integer(i4), dimension(:), allocatable    :: netPerm         ! routing order (permutation)
    integer(i4)                               :: ii, jj, kk
    logical                                   :: flag

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)
    
    nLinks  = nNodes - L11_nOutlets(iBasin)
    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a basin. 

    ! allocate
    allocate ( nLinkFromN  ( nNodes ) )  ! all vectors valid from (1 : nLinks)
    allocate ( nLinkToN    ( nNodes ) )
    allocate ( nLinkROrder ( nNodes ) )
    allocate ( nLinkLabel  ( nNodes ) )
    allocate ( nLinkSink   ( nNodes ) )
    allocate ( netPerm     ( nNodes ) )
    ! initialize
    nLinkFromN(:)         = nodata_i4
    nLinkToN(:)           = nodata_i4
    nLinkROrder(1:nLinks) = 1
    nLinkROrder(nNodes)   = nodata_i4
    nLinkLabel(1:nLinks)  =  0           
    nLinkLabel(nNodes)    = nodata_i4
    nLinkSink(:)          = .FALSE.
    netPerm(:)            = nodata_i4

    ! for a single node model run
    if(nNodes .GT. 1) then
      ! get network vectors of L11 
      nLinkFromN(:) = L11_fromN ( iStart11 : iEnd11 )
      nLinkToN(:)   = L11_toN   ( iStart11 : iEnd11 )

      loop1: do ii = 1, nLinks
         loop2: do jj = 1, nLinks
            if ( jj == ii ) cycle loop2
            if ( nLinkFromN(ii) == nLinkToN(jj) ) then
               nLinkROrder(ii) = -9
            end if
            if ( nLinkROrder(ii) == -9 ) cycle loop1
         end do loop2
      end do loop1
      ! counting headwaters
      kk = 0
      do ii = 1, nLinks
         if ( nLinkROrder(ii) == 1) then
            kk = kk + 1
            nLinkROrder(ii) = kk
            nLinkLabel(ii)  = 1  ! 'Head Water'
         end if
      end do
      ! counting downstream
      do while ( minval( nLinkROrder( 1 : nLinks ) ) < 0 )
       !!  print *, count(nLinkROrder .lt. 0), minval(nLinkROrder)
         loop3: do ii = 1, nLinks
            if ( .NOT. nLinkROrder(ii) == -9 ) cycle loop3
            flag = .TRUE.
            loop4: do jj = 1, nLinks
               if ( jj == ii .OR. nLinkFromN(ii)  /=  nLinkToN(jj) ) then
                  cycle loop4
               else if (.NOT. (  nLinkFromN(ii)  == nLinkToN(jj)  .AND. nLinkROrder(jj) > 0 )) then
                  flag = .FALSE.
                  exit loop4
               else
               end if
            end do loop4

            if (flag) then
               kk = kk + 1
               nLinkROrder(ii) = kk
            end if
         end do loop3
      end do
     
      ! identify sink cells
      do ii = 1, nLinks
         if (L11_fdir(iStart11 + nLinkToN(ii) - 1_i4) .eq. 0_i4) nlinksink(ii) = .True. 
      end do
      where(nlinksink) nLinkLabel = 2 !  'Sink'

      ! keep routing order
      do ii = 1, nLinks
         netPerm(nLinkROrder(ii)) = ii
      end do
     
      ! end of multi-node network design loop
    end if
   
    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L11 network data sets 
    call append( L11_rOrder,  nLinkROrder(:) )
    call append( L11_label,   nLinkLabel(:)  )
    call append( L11_sink,    nLinkSink(:)   )
    call append( L11_netPerm, netPerm(:)     )

    ! free space
    deallocate (nLinkFromN, nLinkToN, nLinkROrder, nLinkLabel, nLinkSink, netPerm )   

  end subroutine L11_routing_order

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_link_location

  !     PURPOSE
  !>        \brief Estimate the LO (row,col) location for each routing link at level L11

  !>        \details If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"        Basin Id 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note Cell location  can ONLY be called after routing order is done.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Stephan Thober, Aug 2015 - ported to mRM
  ! ------------------------------------------------------------------

  subroutine L11_link_location(iBasin)
    use mo_mrm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         basin_mrm,   & ! IN
         L0_fDir,     & ! IN:    flow direction (standard notation) L0
         L11_nOutlets,& ! IN:    Number of Outlets/Sinks
         L0_draSC,    & ! IN:    Index of draining cell of each sub catchment (== cell L11)
         L11_fromN,   & ! IN:    from node 
         L11_rowOut,  & ! IN:    grid vertical location of the Outlet
         L11_colOut,  & ! IN:    grid horizontal location  of the Outlet
         L11_netPerm, & ! IN:    routing order (permutation)
         L11_fRow,    & ! INOUT: from row in L0 grid 
         L11_fCol,    & ! INOUT: from col in L0 grid
         L11_tRow,    & ! INOUT: to row in L0 grid
         L11_tCol       ! INOUT: to col in L0 grid 

    implicit none

    integer(i4), intent(in)                   :: iBasin         ! basin 

    ! local
    integer(i4)                               :: nNodes
    integer(i4)                               :: nLinks
    integer(i4)                               :: nrows0, ncols0
    integer(i4)                               :: nrows110, ncols110
    integer(i4)                               :: iStart0, iEnd0
    integer(i4)                               :: iStart110, iEnd110
    integer(i4)                               :: nrows11, ncols11
    integer(i4)                               :: iStart11, iEnd11
    integer(i4), dimension(:), allocatable    :: rowOut         ! northing cell loc. of the Outlet
    integer(i4), dimension(:), allocatable    :: colOut         ! easting cell loc. of the Outlet
    integer(i4), dimension(:), allocatable    :: nLinkFromN     
    integer(i4), dimension(:), allocatable    :: netPerm
    integer(i4), dimension(:), allocatable    :: nLinkFromRow   
    integer(i4), dimension(:), allocatable    :: nLinkFromCol   
    integer(i4), dimension(:), allocatable    :: nLinkToRow     
    integer(i4), dimension(:), allocatable    :: nLinkToCol     
    logical,   dimension(:,:), allocatable    :: mask0
    integer(i4), dimension(:,:), allocatable  :: fDir0          
    integer(i4), dimension(:,:), allocatable  :: draSC0
    integer(i4)                               :: ii, rr, kk
    integer(i4)                               :: iNode, iRow, jCol
    integer(i4), dimension(:,:), allocatable  :: oLoc           ! output location in L0
    integer(i4)                               :: nOutlets       ! number of outlets in basin
    logical                                   :: is_outlet      ! flag for finding outlet

    ! level-0 information
    call get_basin_info_mrm (iBasin, 0, nrows0, ncols0, iStart=iStart0, iEnd=iEnd0, mask=mask0) 

    ! level-110 information
    call get_basin_info_mrm (iBasin, 110, nrows110, ncols110, iStart=iStart110, iEnd=iEnd110) 

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)
    nOutlets = L11_nOutlets(iBasin)

    nLinks  = nNodes - nOutlets

    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a basin. 
    ! allocate
    allocate ( rowOut        ( nNodes ) )
    allocate ( colOut        ( nNodes ) )
    allocate ( nLinkFromN    ( nNodes ) )  ! all network vectors valid from (1 : nLinks)
    allocate ( netPerm       ( nNodes ) )  
    allocate ( nLinkFromRow  ( nNodes ) )
    allocate ( nLinkFromCol  ( nNodes ) )
    allocate ( nLinkToRow    ( nNodes ) )
    allocate ( nLinkToCol    ( nNodes ) )
    allocate ( fDir0         ( nrows0, ncols0 ) )
    allocate ( draSC0        ( nrows0, ncols0 ) )

    ! initialize
    rowOut       = nodata_i4    
    colOut       = nodata_i4    
    nLinkFromN   = nodata_i4    
    netPerm      = nodata_i4    
    nLinkFromRow = nodata_i4    
    nLinkFromCol = nodata_i4    
    nLinkToRow   = nodata_i4    
    nLinkToCol   = nodata_i4    
    fDir0        = nodata_i4    
    draSC0       = nodata_i4    

    ! for a single node model run
    if(nNodes .GT. 1) then
      ! get fDir at L0
      fDir0(:,:) =   UNPACK( L0_fDir  (iStart0:iEnd0),  mask0, nodata_i4 )
      draSC0(:,:) =  UNPACK( L0_draSC (iStart110:iEnd110),  mask0, nodata_i4 )

      ! get network vectors of L11 
      nLinkFromN(:) = L11_fromN   ( iStart11 : iEnd11 )
      netPerm(:)    = L11_netPerm ( iStart11 : iEnd11 )
      rowOut(:)     = L11_rowOut  ( iStart11 : iEnd11 )
      colOut(:)     = L11_colOut  ( iStart11 : iEnd11 )  

      ! finding main outlet (row, col) in L0
      allocate(oLoc(Noutlets, 2))
      oLoc(:, 1) = basin_mrm%L0_rowOutlet(:Noutlets, iBasin)
      oLoc(:, 2) = basin_mrm%L0_colOutlet(:Noutlets, iBasin) 

      ! Location of the stream-joint cells  (row, col)
      do rr = 1, nLinks

         ii = netPerm(rr)
         iNode = nLinkFromN(ii)
         iRow = rowOut(iNode)
         jCol = colOut(iNode) 
         call moveDownOneCell( fDir0(iRow,jcol), iRow, jcol ) 
         ! set "from" cell
         nLinkFromRow(ii) = iRow
         nLinkFromCol(ii) = jCol

         ! check whether this location is an outlet
         is_outlet = .False.
         do kk = 1, Noutlets
            if (iRow .eq. oLoc(kk, 1) .and. jCol .eq. oLoc(kk, 2)) is_outlet = .True.
         end do

         if (is_outlet) then

            nLinkToRow(ii) = iRow
            nLinkToCol(ii) = jCol

         else

            do while ( .not. ( draSC0(iRow,jCol) > 0 ) )
               call moveDownOneCell( fDir0(iRow,jcol), iRow, jCol )
               ! check whether this location is an outlet and exit
               do kk = 1, Noutlets
                  if (iRow .eq. oLoc(kk, 1) .and. jCol .eq. oLoc(kk, 2)) exit
               end do
               ! if ( iRow == oLoc(1) .and. jCol == oLoc(2)) exit
            end do
            ! set "to" cell (when an outlet is reached)
            nLinkToRow(ii) = iRow
            nLinkToCol(ii) = jCol

         end if
      end do

      ! end of multi-node network design loop
    end if
    
    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L11 network data sets 
    call append( L11_fRow,   nLinkFromRow(:) )
    call append( L11_fCol,   nLinkFromCol(:) )
    call append( L11_tRow,   nLinkToRow(:)   )
    call append( L11_tCol,   nLinkToCol(:)   )

    ! free space
    deallocate ( rowOut, colOut, nLinkFromN, netPerm, nLinkFromRow, &
         nLinkFromCol,  nLinkToRow, nLinkToCol, fDir0, draSC0       )  

  end subroutine L11_link_location

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_set_drain_outlet_gauges
  
  !     PURPOSE
  !>        \brief Draining cell identification and Set gauging node

  !>        \details Perform the following tasks: \n
  !>                 - Draining cell identification (cell at L0 to draining cell outlet at L11). \n
  !>                 - Set gauging nodes \n
  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"        Basin Id 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Matthias Zink , Mar 2014 - bugfix, added inflow gauge
  !                  Rohini Kumar  , Apr 2014 - variable index is changed to index_gauge 
  !                  Stephan Thober, Aug 2015 - ported to mRM
  ! ------------------------------------------------------------------
  subroutine L11_set_drain_outlet_gauges(iBasin)
    use mo_mrm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         basin_mrm,   & 
         L0_fDir,     &       ! IN: flow direction (standard notation) L0
         L0_draSC,    &       ! IN: Index of draining cell of each sub catchment (== cell L11)
         L0_gaugeLoc, &       ! IN: location of gauges (read with gauge Id then 
         !                    !     transformed into gauge running ID => [1,nGaugesTotal]
         L0_cellCoor, &       ! IN: cell coordinates (row,col) -> <only domain> input data
         L0_InflowgaugeLoc, & ! IN: location of gauges (read with gauge Id then  
         !                    !     transformed into gauge running ID => [1,nGaugesTotal]
         L0_L11_Id,   &       ! IN: mapping of L11 Id on L0
         L0_draCell           ! INOUT: draining cell id at L11 of ith cell of L0

    implicit none

    integer(i4), intent(in)                   :: iBasin         ! basin 

    ! local
    integer(i4)                               :: nCells0
    integer(i4)                               :: nrows0, ncols0
    integer(i4)                               :: nrows110, ncols110
    integer(i4)                               :: iStart0, iEnd0
    integer(i4)                               :: iStart110, iEnd110
    logical,     dimension(:,:), allocatable  :: mask0
    integer(i4), dimension(:,:), allocatable  :: cellCoor0
    integer(i4), dimension(:,:), allocatable  :: draSC0         
    integer(i4), dimension(:,:), allocatable  :: fDir0
    integer(i4), dimension(:,:), allocatable  :: gaugeLoc0   
    integer(i4), dimension(:,:), allocatable  :: InflowGaugeLoc0   
    integer(i4), dimension(:,:), allocatable  :: draCell0
    integer(i4), dimension(:,:), allocatable  :: L11Id_on_L0
    integer(i4)                               :: ii, jj, kk, ll
    integer(i4)                               :: iSc
    integer(i4)                               :: iRow, jCol

    ! level-0 information
    call get_basin_info_mrm ( iBasin, 0, nrows0, ncols0, ncells=nCells0, &
         iStart=iStart0, iEnd=iEnd0, mask=mask0     ) 

    ! level-110 information (nrows110,ncols110) always equal to (nrows0,ncols0)
    call get_basin_info_mrm ( iBasin, 110, nrows110, ncols110, &
         iStart=iStart110, iEnd=iEnd110 ) 

    ! allocate
    allocate ( cellCoor0       ( nCells0, 2 ) )  
    allocate ( draSC0          ( nrows0, ncols0 ) )
    allocate ( fDir0           ( nrows0, ncols0 ) )
    allocate ( gaugeLoc0       ( nrows0, ncols0 ) )
    allocate ( InflowGaugeLoc0 ( nrows0, ncols0 ) )
    allocate ( draCell0        ( nrows0, ncols0 ) )
    allocate ( L11Id_on_L0     ( nrows0, ncols0 ) )

    ! initialize
    cellCoor0(:,:)         = nodata_i4  
    draSC0(:,:)            = nodata_i4
    fDir0(:,:)             = nodata_i4
    gaugeLoc0(:,:)         = nodata_i4
    InflowGaugeLoc0(:,:)   = nodata_i4
    draCell0(:,:)          = nodata_i4
    L11Id_on_L0(:,:)       = nodata_i4

    ! get L0 fields
    cellCoor0(:,:)       = L0_cellCoor(iStart0 : iEnd0, :)

    draSC0(:,:)          = UNPACK( L0_draSC          (iStart110:iEnd110), mask0, nodata_i4 )
    fDir0(:,:)           = UNPACK( L0_fDir           (iStart0:iEnd0),     mask0, nodata_i4 )
    gaugeLoc0(:,:)       = UNPACK( L0_gaugeLoc       (iStart0:iEnd0),     mask0, nodata_i4 )
    InflowGaugeLoc0(:,:) = UNPACK( L0_InflowgaugeLoc (iStart0:iEnd0),     mask0, nodata_i4 )
    L11Id_on_L0(:,:)     = UNPACK( L0_L11_Id         (iStart110:iEnd110), mask0, nodata_i4 ) 


    do kk = 1, nCells0
       ii   = cellCoor0(kk,1)
       jj   = cellCoor0(kk,2)
       iSc  = draSC0(ii,jj)
       ! find drainage path
       iRow = ii
       jCol = jj
       do while ( .NOT. iSC > 0 )
          ! move downstream
          call moveDownOneCell( fDir0(iRow,jCol), iRow, jCol )
          iSC = draSC0(iRow,jCol)
       end do
       draCell0(ii,jj) = iSC

       ! find cell at L11 corresponding to gauges in basin at L0 !>> L11Id_on_L0 is Id of
       ! the routing cell at level-11
        if ( gaugeLoc0(ii,jj) .NE. nodata_i4 ) then 
          ! evaluation gauges
          do ll = 1, basin_mrm%nGauges(iBasin)
             ! search for gaugeID in L0 grid and save ID on L11
             if (basin_mrm%gaugeIdList(iBasin, ll) .EQ. gaugeLoc0(ii,jj)) then
                basin_mrm%gaugeNodeList(iBasin, ll) = L11Id_on_L0(ii, jj)
             end if
          end do
       end if

       if ( InflowGaugeLoc0(ii,jj) .NE. nodata_i4 ) then 
          ! inflow gauges
          do ll = 1, basin_mrm%nInflowGauges(iBasin)
             ! search for gaugeID in L0 grid and save ID on L11
             if ( basin_mrm%InflowGaugeIdList(iBasin, ll) .EQ. InflowGaugeLoc0(ii,jj)) &
                  basin_mrm%InflowGaugeNodeList( iBasin, ll ) = L11Id_on_L0(ii,jj)
          end do
       end if
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L0 data sets 
    call append(L0_draCell, PACK(draCell0(:,:),  mask0)) 

    ! free space
    deallocate ( mask0, cellCoor0, draSC0, fDir0, gaugeLoc0, draCell0, L11Id_on_L0) 

  end subroutine  L11_set_drain_outlet_gauges

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_stream_features
  
  !     PURPOSE
  !>        \brief Stream features (stream network and floodplain)

  !>        \details Stream features (stream network and floodplain)\n
  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"        Basin Id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  R. Kumar      , Oct 2013 - stack size increased from nNodes to 100 
  !                  Stephan Thober, Aug 2015 - ported to mRM
  ! ------------------------------------------------------------------

  subroutine L11_stream_features(iBasin)
    use mo_mrm_constants, only: nodata_i4, nodata_dp
    use mo_append, only: append
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         L0_elev_mRM,         & ! IN:    elevation (sinks removed)  [m]
         iFlag_cordinate_sys, & ! IN:    coordinate system
         L0_Id,           & ! IN:    level-0 id
         L0_fDir,         & ! IN:    flow direction (standard notation) L0
         L0_areaCell,     & ! IN:    area of a cell at level-0, -> is same for all basin [m2]
         L11_fRow,        & ! IN:    from row in L0 grid 
         L11_fCol,        & ! IN:    from col in L0 grid
         L11_tRow,        & ! IN:    to row in L0 grid 
         L11_tCol,        & ! IN:    to col in L0 grid 
         L11_netPerm,     & ! IN:    routing order (permutation)
         L0_streamNet,    & ! IN:    stream network
         L0_floodPlain,   & ! IN:    floodplains of stream i
         L11_length,      & ! IN:    total length [m] 
         L11_aFloodPlain, & ! IN:    area of the flood plain [m2]
         L11_nOutlets,    & ! IN:    Number of Outlets/Sinks
         L11_slope          ! INOUT: normalized average slope

    implicit none

    integer(i4), intent(in)                  :: iBasin         ! basin 

    ! local
    integer(i4)                              :: nCells0
    integer(i4)                              :: nNodes
    integer(i4)                              :: nLinks
    integer(i4)                              :: nrows0, ncols0
    integer(i4)                              :: iStart0, iEnd0
    integer(i4)                              :: nrows11, ncols11
    integer(i4)                              :: iStart11, iEnd11
    logical,     dimension(:,:), allocatable :: mask0
    integer(i4), dimension(:,:), allocatable :: iD0
    integer(i4), dimension(:,:), allocatable :: fDir0
    real(dp),    dimension(:,:), allocatable :: elev0
    real(dp),    dimension(:,:), allocatable :: areaCell0
    integer(i4), dimension(:,:), allocatable :: streamNet0
    integer(i4), dimension(:,:), allocatable :: floodPlain0
    integer(i4), dimension(:),   allocatable :: netPerm         ! routing order (permutation)
    integer(i4), dimension(:),   allocatable :: nLinkFromRow   
    integer(i4), dimension(:),   allocatable :: nLinkFromCol
    integer(i4), dimension(:),   allocatable :: nLinkToRow     
    integer(i4), dimension(:),   allocatable :: nLinkToCol  
    real(dp),    dimension(:),   allocatable :: nLinkLength
    real(dp),    dimension(:),   allocatable :: nLinkAFloodPlain
    real(dp),    dimension(:),   allocatable :: nLinkSlope
    integer(i4)                              :: ii, rr, ns
    integer(i4)                              :: frow, fcol
    integer(i4)                              :: fId,  tId
    integer(i4), dimension(:,:), allocatable :: stack, append_chunk
    integer(i4), dimension(:),   allocatable :: dummy_1d
    real(dp)                                 :: length
    integer(i4), dimension(:,:), allocatable :: nodata_i4_tmp
    real(dp),    dimension(:,:), allocatable :: nodata_dp_tmp

    ! level-0 information
    call get_basin_info_mrm ( iBasin, 0, nrows0, ncols0, ncells=nCells0, &
         iStart=iStart0, iEnd=iEnd0, mask=mask0     ) 

    ! level-11 information
    call get_basin_info_mrm (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)

    nLinks  = nNodes - L11_nOutlets(iBasin)

    ! allocate
    allocate ( iD0           ( nrows0, ncols0 ) )
    allocate ( elev0         ( nrows0, ncols0 ) )
    allocate ( fDir0         ( nrows0, ncols0 ) )
    allocate ( areaCell0     ( nrows0, ncols0 ) )
    allocate ( streamNet0    ( nrows0, ncols0 ) )
    allocate ( floodPlain0   ( nrows0, ncols0 ) )

    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a basin.
    allocate ( stack             ( nNodes, 2 ) ) !>> stack(nNodes, 2)
    allocate ( dummy_1d          ( 2 ))
    allocate ( append_chunk      ( 8,      2 ) )
    allocate ( netPerm           ( nNodes ) )  
    allocate ( nLinkFromRow      ( nNodes ) )
    allocate ( nLinkFromCol      ( nNodes ) )
    allocate ( nLinkToRow        ( nNodes ) )  
    allocate ( nLinkToCol        ( nNodes ) ) 
    allocate ( nLinkLength       ( nNodes ) )
    allocate ( nLinkAFloodPlain  ( nNodes ) )
    allocate ( nLinkSlope        ( nNodes ) )

    allocate (nodata_i4_tmp      ( nrows0, ncols0 ) )
    allocate (nodata_dp_tmp      ( nrows0, ncols0 ) )

    ! initialize
    iD0(:,:)             = nodata_i4
    elev0(:,:)           = nodata_dp
    fDir0(:,:)           = nodata_i4
    areaCell0(:,:)       = nodata_dp
    streamNet0(:,:)      = nodata_i4
    floodPlain0(:,:)     = nodata_i4

    stack(:,:)           = nodata_i4
    append_chunk(:,:)    = nodata_i4
    netPerm(:)           = nodata_i4
    nLinkFromRow(:)      = nodata_i4
    nLinkFromCol(:)      = nodata_i4
    nLinkToRow(:)        = nodata_i4
    nLinkToCol(:)        = nodata_i4
    nLinkLength(:)       = nodata_dp
    nLinkAFloodPlain(:)  = nodata_dp
    nLinkSlope(:)        = nodata_dp

    nodata_i4_tmp(:,:)   = nodata_i4
    nodata_dp_tmp(:,:)   = nodata_dp

    ! for a single node model run
    if(nNodes .GT. 1) then
      ! get L0 fields
      iD0(:,:) =         UNPACK( L0_Id   (iStart0:iEnd0),  mask0, nodata_i4_tmp )
      elev0(:,:) =       UNPACK( L0_elev_mRM (iStart0:iEnd0),  mask0, nodata_dp_tmp )
      fDir0(:,:) =       UNPACK( L0_fDir (iStart0:iEnd0),  mask0, nodata_i4_tmp )
      areaCell0(:,:) =   UNPACK( L0_areaCell (iStart0:iEnd0),  mask0, nodata_dp_tmp )

      ! get network vectors of L11 
      netPerm(:)      = L11_netPerm ( iStart11 : iEnd11 )
      nLinkFromRow(:) = L11_fRow    ( iStart11 : iEnd11 )
      nLinkFromCol(:) = L11_fCol    ( iStart11 : iEnd11 )
      nLinkToRow(:)   = L11_tRow    ( iStart11 : iEnd11 )
      nLinkToCol(:)   = L11_tCol    ( iStart11 : iEnd11 )

      ! Flood plains:  stream network delineation
      streamNet0(:,:)  = nodata_i4
      floodPlain0(:,:) = nodata_i4

      do rr = 1, nLinks

         ii    = netPerm(rr)
         frow = nLinkFromRow(ii)
         fcol = nLinkFromCol(ii)

         ! Init
         streamNet0( frow, fcol) = ii
         floodPlain0(frow, fcol) = ii
         stack = 0
         append_chunk = 0
         ns    = 1
         stack(ns,1) = frow
         stack(ns,2) = fcol

         call cellLength(iBasin, fDir0(frow,fcol), fRow, fCol, iFlag_cordinate_sys, nLinkLength(ii) )
         nLinkSlope(ii) = elev0(frow, fcol)

         fId = iD0( frow, fcol )
         tId = iD0( nLinkToRow(ii) , nLinkToCol(ii) )

         do while ( .NOT. (fId == tId))
            ! Search flood plain from point(frow,fcol) upwards, keep co-ordinates in STACK
            do while (ns > 0)
               if (ns + 8 .gt. size(stack,1)) then 
                  call append(stack,append_chunk)
               end if
               call moveUp( elev0, fDir0, frow, fcol, stack, ns )
               stack(1,1) = 0
               stack(1,2) = 0
               ! stack = cshift(stack, SHIFT = 1, DIM = 1)
               ! substitute cshift <<<
               dummy_1d = stack(1, :)
               stack(:size(stack, dim=1) - 1, :) = stack(2:, :)
               stack(size(stack, dim=1), :) = dummy_1d
               ! substitute cshift >>>
               if (stack(1,1) > 0 .and. stack(1,2) > 0 ) floodPlain0( stack(1,1), stack(1,2) ) = ii
               ns = count( stack > 0 ) / 2
            end do

            ! move downstream
            call moveDownOneCell( fDir0(frow,fcol), frow, fcol )
            streamNet0(frow, fcol)  = ii
            floodPlain0(frow, fcol) = ii
            fId = iD0(frow, fcol)
            stack = 0
            ns = 1
            stack(ns,1) = frow
            stack(ns,2) = fcol
            call cellLength(iBasin, fDir0(fRow,fCol), fRow, fCol, iFlag_cordinate_sys, length )
            nLinkLength(ii) = nLinkLength(ii) + length

         end do

         ! stream bed slope
         nLinkSlope(ii) = ( nLinkSlope(ii) - elev0(frow, fcol) ) / nLinkLength(ii)

         if ( nLinkSlope(ii) < 0.0001_dp) nLinkSlope(ii) = 0.0001_dp

         ! calculate area of floodplains (avoid overwriting)
         nLinkAFloodPlain(ii) = sum ( areaCell0(:,:),  mask = ( floodPlain0(:,:) == ii ) )
         !  old > real( count( floodPlain0(:,:,) == i), dp ) * areaCell0

      end do

      ! end of multi-node network design loop
    end if

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L0 data sets 
    call append( L0_streamNet,    PACK ( streamNet0(:,:),  mask0)   ) 
    call append( L0_floodPlain,   PACK ( floodPlain0(:,:),  mask0)  ) 

    ! L11 network data sets 
    call append( L11_length,       nLinkLength(:)      )
    call append( L11_aFloodPlain,  nLinkAFloodPlain(:) )
    call append( L11_slope,        nLinkSlope(:)       )

    ! free space
    deallocate (&
         mask0, iD0, elev0, fDir0, areaCell0, streamNet0, floodPlain0,      &
         stack, netPerm, nLinkFromRow, nLinkFromCol, nLinkToRow, nLinkToCol, &       
         nLinkLength, nLinkAFloodPlain, nLinkSlope, dummy_1d) 
    deallocate(nodata_i4_tmp,nodata_dp_tmp)    

  end subroutine L11_stream_features

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_fraction_sealed_floodplain

  !     PURPOSE
  !         \brief Fraction of the flood plain with impervious cover

  !>        \details Fraction of the flood plain with impervious cover (\ref fig_routing "Routing
  !>                 Network"). This proportion is used to regionalize the Muskingum parameters.
  !>                 Samaniego et al. \cite SB05 found out that this fraction is one of the statistically
  !>                 significant predictor variables of peak discharge in mesoscale basins.\n

  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: nLinks"           number of links for a given basin
  !>        \param[in] "integer(i4)        :: LCover0"          land cover id field (basin)
  !>        \param[in] "integer(i4)        :: floodPlain0"      floodplains of stream i (basin)
  !>        \param[in] "real(dp)           :: areaCell0"        area of a cell at level-0 [m2]
  !>        \param[in] "real(dp)           :: nLinkAFloodPlain" area of the flood plain at level-11 [m2] 
  !>        \param[in] "integer(i4)        :: LCClassImp"       Impervious land cover class Id, e.g. = 2 (old code)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)         :: nLinkFracFPimp"   Fraction of the flood plain with impervious cover

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note - Call only after L0 and L11 initialization routines\n
  !>             - All spatial input variables are 2D for a given basin, unlike the L11_, L0_ vectors

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  !                  Stephan Thober, Aug 2015 - ported to mRM
  ! ------------------------------------------------------------------
  subroutine L11_fraction_sealed_floodplain( &
       nLinks, LCover0, floodPlain0, & ! INTENT IN
       areaCell0, nLinkAFloodPlain,  & ! INTENT IN
       LCClassImp,                   & ! INTENT IN
       nLinkFracFPimp    )             ! INTENT OUT
    use mo_mrm_constants, only: nodata_dp
    implicit none

    integer(i4),                intent(in)  :: nLinks
    integer(i4), dimension(:),  intent(in)  :: LCover0
    integer(i4), dimension(:),  intent(in)  :: floodPlain0
    real(dp),    dimension(:),  intent(in)  :: areaCell0
    real(dp),    dimension(:),  intent(in)  :: nLinkAFloodPlain
    integer(i4),                intent(in)  :: LCClassImp         ! e.g. = 2 (old code)
    real(dp), dimension(nLinks),intent(out) :: nLinkFracFPimp  

    ! local
    integer(i4) :: ii

    ! initalization
    nLinkFracFPimp(1:nLinks) = nodata_dp

    ! for a single node model run
    if(nLinks .GT. 0) then
      do ii = 1, nLinks
         nLinkFracFPimp(ii) =  sum ( areaCell0(:),  & 
                              mask = ( floodPlain0(:) == ii .and. LCover0(:) == LCClassImp ) ) &
                              /  nLinkAFloodPlain(ii)
       end do
    end if
     
  end subroutine L11_fraction_sealed_floodplain

  ! ------------------------------------------------------------------
  !  MOVE UPSTREAM FROM-TO
  ! ------------------------------------------------------------------
  subroutine moveUp(elev0, fDir0, fi, fj, ss, nn)

    use mo_mrm_constants,    only: deltaH
    use mo_utils,            only: le, ge

    implicit none

    real(dp),    dimension(:,:), allocatable, intent(IN)    :: elev0
    integer(i4), dimension(:,:), allocatable, intent(IN)    :: fDir0
    integer(i4),                              intent(IN)    :: fi, fj  ! co-ordinate of the stream bed
    integer(i4), dimension(:,:),              intent(INOUT) :: ss
    integer(i4),                              intent(INOUT) :: nn

    ! local
    integer(i4) :: ii ,jj, ip, im, jp, jm
    integer(i4) :: nrows, ncols

    ii = ss(1,1)
    jj = ss(1,2)
    ip = ii+1
    im = ii-1
    jp = jj+1
    jm = jj-1

    nrows = size(fDir0, 1)
    ncols = size(fDir0, 2)

    !E
    if   (jp                <= ncols ) then
       if ( (   fdir0(ii,jp) == 16        )                   .and. &
            ( le(( elev0(ii,jp)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = ii
          ss(nn,2) = jp
          !print *, i,jp
       end if
    end if
    
    !SE
    if ( ( ip                <= nrows ) .and. &
         ( jp                <= ncols )      ) then
       if ( (   fdir0(ip,jp) == 32        )                   .and. &
            ( le(( elev0(ip,jp)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = ip
          ss(nn,2) = jp
          !print *, ip,jp
       end if
    end if

    !S
    if ( ( ip               <= nrows )  .and. &
         ( jp               <= ncols )     ) then
       if ( (   fdir0(ip,jj) == 64        )                 .and. &
            ( le(( elev0(ip,jj)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = ip
          ss(nn,2) = jj
          !print *, ip,j
       end if
    end if

    !SW
    if ( ( ip                <= nrows ) .and. &
         ( jp                <= ncols ) .and. &
         ( jm                >= 1         )     ) then
       if ( (   fdir0(ip,jm) == 128       )                 .and. &
            ( le(( elev0(ip,jm)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = ip
          ss(nn,2) = jm
          !print *, ip,jm
       end if
    end if

    !W
    if ( ( jm                 >= 1         ) .and. &
         (jp                 <= ncols      ) ) then
       if ( (   fdir0(ii,jm)  == 1         )                 .and. &
            ( le(( elev0(ii,jm)   - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)   - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = ii
          ss(nn,2) = jm
          !print *, i,jm
       end if
    end if

    !NW
    if ( ( im                >= 1         ) .and. &
         ( jp                <= ncols     ) .and. &
         ( jm                >= 1         )      )  then
       if ( (   fdir0(im,jm) == 2         )                 .and. &
            ( le(( elev0(im,jm)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = im
          ss(nn,2) = jm
          !print *, im,jm
       end if
    end if

    !N
    if ( (  im                >= 1         ) .and. &
         ( jp                 <= ncols     ) ) then
       if ( (   fdir0(im,jj)  == 4         )                 .and. &
            ( le(( elev0(im,jj)   - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)   - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = im
          ss(nn,2) = jj
          !print *, im,j
       end if
    end if

    !NE
    if ( ( im                >= 1       ) .and. &
         ( jp                <= ncols   )         )  then
       if ( (   fdir0(im,jp) == 8           )               .and. &
            ( le(( elev0(im,jp)  - elev0(fi,fj) ), deltaH) )  .and. &
            ( ge(( elev0(ii,jp)  - elev0(fi,fj) ), 0.0_dp) )        &
            ) then
          nn = nn + 1
          ss(nn,1) = im
          ss(nn,2) = jp
          !print *, im,jp
       end if
    end if

  end subroutine moveUp

  ! ------------------------------------------------------------------
  !  MOVE DOWNSTREAM
  ! ------------------------------------------------------------------
  subroutine moveDownOneCell(fDir, iRow, jCol)

    implicit none

    integer(i4), intent(IN)     :: fDir
    integer(i4), intent(INOUT)  :: iRow, jCol

    select case (fDir)
    case(1)   !E
       jCol = jCol + 1
    case(2)   !SE
       iRow = iRow + 1
       jCol = jCol + 1
    case(4)   !S
       iRow = iRow + 1
    case(8)   !SW
       iRow = iRow + 1
       jCol = jCol - 1
    case(16)  !W
       jCol = jCol - 1
    case(32)  !NW
       iRow = iRow - 1
       jCol = jCol - 1
    case(64)  !N
       iRow = iRow - 1
    case(128) !NE
       iRow = iRow - 1
       jCol = jCol + 1
    case default !sink
       ! do nothing
    end select

  end subroutine moveDownOneCell

  ! ------------------------------------------------------------------
  !  CELL LENGTH
  ! ------------------------------------------------------------------
  subroutine cellLength(iBasin, fDir, iRow, jCol, iCoorSystem, length)

    use mo_constants, only: SQRT2_dp
    use mo_mrm_global_variables, only: level0

    implicit none

    integer(i4), intent(IN)  :: iBasin
    integer(i4), intent(IN)  :: fDir
    integer(i4), intent(IN)  :: iRow
    integer(i4), intent(IN)  :: jCol
    integer(i4), intent(IN)  :: iCoorSystem
    real(dp),    intent(OUT) :: length
    
    ! local variables
    integer(i4)              :: iRow_to, jCol_to
    real(dp)                 :: lat_1, long_1, lat_2, long_2
    

    ! regular X-Y cordinate system
    IF(iCoorSystem .EQ. 0) THEN
    
       select case (fDir)
           case(1, 4, 16, 64)       ! E, S, W, N
              length = 1.0_dp
           case(2, 8, 32, 128)      ! SE, SW, NW, NE
              length = SQRT2_dp
        end select
        length = length * level0%cellsize(iBasin)
        
    ! regular lat-lon cordinate system
    ELSE IF(iCoorSystem .EQ. 1) THEN
        iRow_to = iRow
        jCol_to = jCol
        
        ! move in the direction of flow
        call moveDownOneCell(fDir, iRow_to, jCol_to)
        
        ! estimate lat-lon points
        lat_1  = level0%yllcorner(iBasin) + real( (level0%ncols(iBasin)-jCol),dp)*level0%cellsize(iBasin) + &
                                            0.5_dp*level0%cellsize(iBasin)
        long_1 = level0%xllcorner(iBasin) + real( (iRow-1)                   ,dp)*level0%cellsize(iBasin) + &
                                            0.5_dp*level0%cellsize(iBasin)
        
        lat_2  = level0%yllcorner(iBasin) + real( (level0%ncols(iBasin)-jCol_to),dp)*level0%cellsize(iBasin) + &
                                            0.5_dp*level0%cellsize(iBasin)
        long_2 = level0%xllcorner(iBasin) + real( (iRow_to-1)                   ,dp)*level0%cellsize(iBasin) + &
                                            0.5_dp*level0%cellsize(iBasin)
        ! get distance between two points
        call get_distance_two_lat_lon_points(lat_1, long_1, lat_2, long_2, length)
        
    END IF
    !
  end subroutine cellLength


  ! --------------------------------------------------------------------------

  !     NAME
  !         get_distance_two_lat_lon_points
  
  !     PURPOSE
  !>        \brief estimate distance in [m] between two points in a lat-lon
  
  !>        \details estimate distance in [m] between two points in a lat-lon
  
  !     INTENT(IN)
  !>        \param[in] "real(dp)    :: lat1"    latitude  of point-1
  !>        \param[in] "real(dp)    :: long1"   longitude of point-1
  !>        \param[in] "real(dp)    :: lat2"    latitude  of point-2
  !>        \param[in] "real(dp)    :: long2"   longitude of point-2

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)    :: distance_out"    distance between two points [m]

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call L11_variable_init(1)

  !     LITERATURE
  !         Code is based on one that is implemented in the VIC-3L model 

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date   May 2014
  !         Modified,
  !                  Stephan Thober, Aug 2015 - ported to mRM

  ! --------------------------------------------------------------------------
  subroutine get_distance_two_lat_lon_points(lat1, long1, lat2, long2, distance_out)

    use mo_constants,     only: TWOPI_dp, RadiusEarth_dp
    implicit none

    real(dp), intent(in)             :: lat1, long1, lat2, long2
    real(dp), intent(out)            :: distance_out
   
    ! local variables
    real(dp)                         :: theta1
    real(dp)                         :: phi1
    real(dp)                         :: theta2
    real(dp)                         :: phi2
    real(dp)                         :: dtor
    real(dp)                         :: term1
    real(dp)                         :: term2
    real(dp)                         :: term3
    real(dp)                         :: temp

    dtor   = TWOPI_dp/360.0_dp
    theta1 = dtor*long1
    phi1   = dtor*lat1
    theta2 = dtor*long2
    phi2   = dtor*lat2
      
    term1  = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2)
    term2  = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2)
    term3  = sin(phi1)*sin(phi2)
    temp   = term1+term2+term3
    if(temp .GT. 1.0_dp) temp = 1.0_dp

    distance_out = RadiusEarth_dp*acos(temp);

  end subroutine get_distance_two_lat_lon_points


end module mo_mrm_net_startup
