!> \file mo_net_startup.f90

!> \brief Startup drainage network for mHM.

!> \details This module initializes the drainage network at L11 in mHM.\n
!>  - Delineation of drainage network at level 11.    \n
!>  - Setting network topology (i.e. nodes and link). \n
!>  - Determining routing order.                      \n
!>  - Determining cell locations for network links.   \n
!>  - Find drainage outlet.                           \n
!>  - Determine stream (links) features.              \n

!> \authors Luis Samaniego
!> \date Dec 2012

MODULE mo_net_startup

  ! This module sets the river network characteristics and routing order.

  ! Written  Luis Samaniego, Mar 2005

  USE mo_kind,          ONLY: i4, dp
  USE mo_init_states,   ONLY: get_basin_info
  USE mo_mhm_constants, ONLY: nodata_i4, nodata_dp
  USE mo_append,        ONLY: append

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L11_variable_init
  PUBLIC :: L11_flow_direction
  PUBLIC :: L11_set_network_topology
  PUBLIC :: L11_routing_order
  PUBLIC :: L11_link_location
  PUBLIC :: L11_set_drain_outlet_gauges
  PUBLIC :: L11_stream_features
  PUBLIC :: L11_fraction_sealed_floodplain
  PUBLIC :: routing_dummy_alloc

CONTAINS

  ! --------------------------------------------------------------------------

  !     NAME
  !         L11_variable_init
  !     PURPOSE
  !>        \brief Cell numbering at ROUTING LEVEL-11

  !>        \details Cell numbering at ROUTING LEVEL-11  \n
  !>        List of Level- 0 and 1 cells contained within a given Level-11 cell.\n
  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  ! --------------------------------------------------------------------------

  subroutine L11_variable_init(iBasin)

    use mo_global_variables, only :          &
         nBasins, basin,                     &
         level1, level11, resolutionRouting, &
         L11_cellCoor,                       & ! cell coordinates (row,col)
         L11_nCells,                         & ! Total No. of routing cells  (= nNodes)
         L11_Id                                ! ids of grid at level-11    
    use  mo_init_states,  only : calculate_grid_properties

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
    integer(i4)                              :: iStart1, iEnd1
    integer(i4)                              :: iStartMask1, iEndMask1
    logical,     dimension(:,:), allocatable :: mask1, mask11
    integer(i4), dimension(:,:), allocatable :: cellCoor
    integer(i4)                              :: kk
    integer(i4)                              :: ic, jc, icc, jcc
    real(dp)                                 :: cellFactorRbyH
    integer(i4), dimension(:),   allocatable :: Id                         ! old Id11  ids of grid at level-11               
    !--------------------------------------------------------
    ! STEPS::
    ! 1) Estimate each variable locally for a given basin
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0, xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0) 

    if(iBasin == 1) then
       allocate( level11%nrows     (nBasins) )
       allocate( level11%ncols     (nBasins) )
       allocate( level11%xllcorner (nBasins) )
       allocate( level11%yllcorner (nBasins) )
    end if

    ! grid information
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp, &
         resolutionRouting ,                                                                      &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin),                 &
         level11%yllcorner(iBasin), level11%cellsize, level11%nodata_value        ) 
    ! level-1 information
    call get_basin_info (iBasin, 1, nrows1, ncols1, iStart=iStart1, iEnd=iEnd1,                   &
         iStartMask=iStartMask1, iEndMask=iEndMask1, mask=mask1 ) 

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11) 

    allocate( mask11(nrows11, ncols11) )
    mask11(:,:) = .FALSE.

    cellFactorRbyH = level11%cellsize / level1%cellsize

    ! create a mask: Id
    do jc = 1, ncols1
       jcc = ceiling ( real(jc, dp)/cellFactorRbyH )
       do ic = 1, nrows1
          if ( .not. mask1(ic,jc) ) cycle
          ! Identify grids (of level-1) which will take part at the routing level-11
          icc = ceiling ( real(ic, dp)/cellFactorRbyH )
          mask11(icc, jcc) = .TRUE.
       end do
    end do

    ncells = count( mask11 )
    allocate ( cellCoor(nCells,2) )
    allocate ( Id( ncells) )

    ! counting valid cells at level 11
    kk = 0
    do jcc = 1, ncols11
       do icc = 1, nrows11
          if ( .not. mask11(icc,jcc) ) cycle
          kk = kk + 1
          Id(kk)         = kk
          cellCoor(kk,1) = icc
          cellCoor(kk,2) = jcc
       end do
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    if(iBasin == 1) then

       !
       allocate(basin%L11_iStart     (nBasins))
       allocate(basin%L11_iEnd       (nBasins))
       allocate(basin%L11_iStartMask (nBasins))
       allocate(basin%L11_iEndMask   (nBasins))    

       ! basin information
       basin%L11_iStart(iBasin) = 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells - 1

       basin%L11_iStartMask(iBasin) = 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    else

       ! basin information
       basin%L11_iStart(iBasin) = basin%L11_iEnd(iBasin-1) + 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells - 1

       basin%L11_iStartMask(iBasin) = basin%L11_iEndMask(iBasin-1) + 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    end if

    call append( basin%L11_Mask,  RESHAPE( mask11, (/nrows11*ncols11/)  )  )
    ! other L11 data sets
    call append( L11_cellCoor, cellCoor )
    call append( L11_Id, Id )

    L11_nCells = size( L11_Id, 1 )

    ! free space
    deallocate(Id, mask1, mask11, cellCoor)

  end subroutine L11_variable_init

  ! --------------------------------------------------------------------------

  !     NAME
  !         L11_flow_direction

  !     PURPOSE

  !>       \brief Determine the flow direction of the upscaled river
  !>    network at level L11.

  !>       \details The hydrographs generated at each cell are routed
  !>    through the drainage network at level-11 towards the basin's
  !>    outlet. The drainage network at level-11 is conceptualized as a
  !>    graph whose nodes are hypothetically located at the center of
  !>    each grid cell connected by links that represent the river
  !>    reaches. The flow direction of a link correspond to the
  !>    direction towards a neighboring cell in which the net flow
  !>    accumulation (outflows minus inflows) attains its maximum
  !>    value. The net flow accumulation across a cell's boundary at
  !>    level-11 is estimated based on flow direction and flow
  !>    accumulation obtained at level-0 (\ref fig_routing "Routing
  !>    Network"). Note: level-1 denotes the modeling level, whereas
  !>    level-L11 is at least as coarse as level-1. Experience has
  !>    shown that routing can be done at a coarser resolution as
  !>    level-1, hence the level-11 was introduced.

  !>     \image html  routing.png "Upscaling routing network from L0 to L1 (or L11)"
  !>    \anchor fig_routing \image latex routing.pdf "Upscaling routing network from L0 to L1 (or L11)" width=14cm

  !>    The left panel depicts a schematic derivation of a drainage
  !>    network at the level-11 based on level-0 flow direction and
  !>    flow accumulation. The dotted line circle denotes the point
  !>    with the highest flow accumulation within a grid cell. The
  !>    topology of a tipical drainage routing network at level-11 is
  !>    shown in the right panel. Gray color areas denote the flood
  !>    plains estimated in mo_net_startup, where the network
  !>    upscaling is also carried out.

  !>    For the sake of simplicity, it is assumed that all runoff leaving
  !>    a given cell would exit through a major direction.

  !>    If a variable is added or removed here, then it also has to 
  !>    be added or removed in the subroutine L11_config_set in
  !>    module mo_restart and in the subroutine set_L11_config in module
  !>    mo_set_netcdf_restart

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

  !     LITERATURE

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  ! --------------------------------------------------------------------------
  subroutine L11_flow_direction(iBasin)

    use mo_global_variables, only : &
         nBasins, basin, level0, level1, level11,   &
         L0_fAcc, L0_fDir, L0_id, &
         L0_draSC,          & ! INOUT: draining cell of each sub catchment (== cell L11)
         L0_cellCoor,       &
         L0_L11_Id,         & ! INOUT: mapping of L11 Id on L0
         L1_L11_Id,         & ! INOUT: mapping of L11 Id on L1
         L11_Id,            &
         L11_cellCoor,      &
         L11_rowOut,        & ! INOUT: grid vertical location of the Outlet
         L11_colOut,        & ! INOUT: grid horizontal location  of the Outlet
         L11_fDir,          & ! INOUT: flow direction at L11 (standard notation)
         L11_upBound_L0,    & ! INOUT: row start at finer level-0 scale 
         L11_downBound_L0,  & ! INOUT: row end at finer level-0 scale 
         L11_leftBound_L0,  & ! INOUT: col start at finer level-0 scale 
         L11_rightBound_L0, & ! INOUT: col end at finer level-0 scale 
         L11_upBound_L1,    & ! INOUT: row start at finer level-1 scale 
         L11_downBound_L1,  & ! INOUT: row end at finer level-1 scale 
         L11_leftBound_L1,  & ! INOUT: col start at finer level-1 scale 
         L11_rightBound_L1    ! INOUT: col end at finer level-1 scale 

    implicit none

    integer(i4), intent(in)                  :: iBasin         ! basin 

    ! local
    integer(i4)                              :: nCells0
    integer(i4)                              :: nCells1
    integer(i4)                              :: nNodes      ! =  ncells11
    integer(i4)                              :: nrows0, ncols0
    integer(i4)                              :: nrows1, ncols1
    integer(i4)                              :: nrows11, ncols11
    integer(i4)                              :: iStart0, iEnd0
    integer(i4)                              :: iStart1, iEnd1
    integer(i4)                              :: iStart11, iEnd11
    logical,     dimension(:,:), allocatable :: mask0, mask1, mask11
    integer(i4), dimension(:),   allocatable :: upBound1, downBound1, leftBound1, rightBound1 
    integer(i4), dimension(:),   allocatable :: upBound0, downBound0, leftBound0, rightBound0 
    real(dp)                                 :: cellFactorR, cellFactorRbyH
    integer(i4)                              :: icc, jcc
    integer(i4)                              :: ii, jj, kk, ic, jc 
    integer(i4)                              :: iu, id
    integer(i4)                              :: jl, jr
    integer(i4), dimension(:,:), allocatable :: L11Id_on_L0 ! mapping of L11 Id on L0
    integer(i4), dimension(:,:), allocatable :: L11Id_on_L1 ! mapping of L11 Id on L1
    integer(i4), dimension(:,:), allocatable :: Id11        ! ids of grid at level-11     
    integer(i4), dimension(:,:), allocatable :: iD0            
    integer(i4), dimension(:,:), allocatable :: fDir0          
    integer(i4), dimension(:,:), allocatable :: fAcc0          
    integer(i4), dimension(:,:), allocatable :: fDir11         
    integer(i4), dimension(:,:), allocatable :: cellCoor0
    integer(i4), dimension(:,:), allocatable :: cellCoor11
    integer(i4), dimension(:),   allocatable :: rowOut      ! northing cell loc. of the Outlet
    integer(i4), dimension(:),   allocatable :: colOut      ! easting cell loc. of the Outlet
    integer(i4), dimension(:,:), allocatable :: draSC0         
    integer(i4), dimension(2)                :: oLoc        ! output location in L0
    integer(i4)                              :: side
    integer(i4)                              :: fAccMax, idMax

    !--------------------------------------------------------
    ! STEPS:
    ! 1) Estimate each variable locally for a given basin
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------

    ! level-0 information
    call get_basin_info (iBasin, 0, nrows0, ncols0, ncells=nCells0,   &
         iStart=iStart0, iEnd=iEnd0, mask=mask0) 

    ! level-1 information
    call get_basin_info (iBasin, 1, nrows1, ncols1, ncells=nCells1,   &
         iStart=iStart1, iEnd=iEnd1, mask=mask1) 

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nNodes, &
         iStart=iStart11, iEnd=iEnd11, mask=mask11)

    allocate ( upBound1    (nCells1) )
    allocate ( downBound1  (nCells1) )
    allocate ( leftBound1  (nCells1) )
    allocate ( rightBound1 (nCells1) )

    allocate ( upBound0    (nNodes) )
    allocate ( downBound0  (nNodes) )
    allocate ( leftBound0  (nNodes) )
    allocate ( rightBound0 (nNodes) )

    allocate ( L11Id_on_L0  (nrows0, ncols0 ) )
    allocate ( L11Id_on_L1  (nrows1, ncols1 ) )
    allocate ( Id11     (nrows11, ncols11 ) )

    cellFactorR    = level11%cellsize / level0%cellsize
    cellFactorRbyH = level11%cellsize / level1%cellsize

    ! get Ids of L11 
    Id11(:,:) =  UNPACK( L11_Id(iStart11:iEnd11),  mask11, nodata_i4 )

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

          ! Delimitation of level-11 cells on level-1
          L11Id_on_L1(iu:id, jl:jr) = Id11(icc, jcc)

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

    allocate ( iD0         ( nrows0, ncols0 ) )
    allocate ( fAcc0       ( nrows0, ncols0 ) )
    allocate ( fDir0       ( nrows0, ncols0 ) )
    allocate ( draSC0      ( nrows0, ncols0 ) )
    allocate ( cellCoor0   ( nCells0, 2 ) )  
    allocate ( cellCoor11  ( nNodes,  2 ) )  
    allocate ( fDir11      ( nrows11, ncols11 ) )
    allocate ( rowOut      ( nNodes ) )
    allocate ( colOut      ( nNodes ) )

    fDir11 = nodata_i4
    draSC0 = nodata_i4
    rowOut = 0
    colOut = 0

    ! get iD, fAcc, fDir at L0
    iD0(:,:)   =  UNPACK( L0_Id   (iStart0:iEnd0),  mask0, nodata_i4 )
    fAcc0(:,:) =  UNPACK( L0_fAcc (iStart0:iEnd0),  mask0, nodata_i4 )
    fDir0(:,:) =  UNPACK( L0_fDir (iStart0:iEnd0),  mask0, nodata_i4 )

    cellCoor0(:,:)  = L0_cellCoor  (iStart0 : iEnd0,  :)
    cellCoor11(:,:) = L11_cellCoor (iStart11: iEnd11, :)

    ! finding main outlet (row, col) in L11
    oLoc = maxloc ( fAcc0, mask0 )
    kk    = L11Id_on_L0( oLoc(1), oLoc(2) )
    fDir11 ( cellCoor11(kk,1), cellCoor11(kk,2) ) = 0

    ! set location of main outlet in L11
    rowOut(kk) = oLoc(1)
    colOut(kk) = oLoc(2)
    draSC0 ( oLoc(1), oLoc(2) ) = kk

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

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! allocate space for row and col Outlet
    if(iBasin .eq. 1) then
       allocate( basin%L0_rowOutlet(nBasins) ) 
       allocate( basin%L0_colOutlet(nBasins) )
    end if


    ! L0 data sets

    basin%L0_rowOutlet(iBasin) = oLoc(1)
    basin%L0_colOutlet(iBasin) = oLoc(2)

    call append( L0_draSC,     PACK ( draSC0(:,:),  mask0)  ) 
    call append( L0_L11_Id,    PACK ( L11Id_on_L0(:,:), mask0)  )

    ! L1 data sets
    call append( L1_L11_Id,    PACK ( L11Id_on_L1(:,:), mask1)  )

    ! L11 data sets
    call append( L11_fDir,     PACK ( fDir11(:,:),      mask11) )
    call append( L11_rowOut          ,  rowOut(:)               )
    call append( L11_colOut          ,  colOut(:)               )
    call append( L11_upBound_L0      ,  upBound0(:)             )
    call append( L11_downBound_L0    ,  downBound0(:)           )
    call append( L11_leftBound_L0    ,  leftBound0(:)           )
    call append( L11_rightBound_L0   ,  rightBound0(:)          )
    call append( L11_upBound_L1      ,  upBound1(:)             )
    call append( L11_downBound_L1    ,  downBound1(:)           )
    call append( L11_leftBound_L1    ,  leftBound1(:)           )
    call append( L11_rightBound_L1   ,  rightBound1(:)          )

    ! free space
    deallocate(mask0, mask1, mask11, &
         upBound1, downBound1, leftBound1, rightBound1, & 
         upBound0, downBound0, leftBound0, rightBound0, & 
         L11Id_on_L0, L11Id_on_L1, Id11,                &      
         iD0, fDir0, fAcc0, fDir11, cellCoor0,          &
         cellCoor11, rowOut, colOut, draSC0         )   

  end subroutine L11_flow_direction

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_set_network_topology

  !     PURPOSE
  !>        \brief Set network topology

  !>        \details Set network topology from and to node for all links
  !>        at level-11 (\ref fig_routing "Routing Network") \n

  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart.

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

  !     LITERATURE

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  ! ------------------------------------------------------------------

  subroutine L11_set_network_topology(iBasin)

    use mo_global_variables, only: &
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
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nNodes, &
         iStart=iStart11, iEnd=iEnd11, mask=mask11)

    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a basin. 

    allocate ( nLinkFromN ( nNodes    ) )  ! valid from (1 : nLinks)
    allocate ( nLinkToN   ( nNodes    ) )  ! "
    allocate ( cellCoor11 ( nNodes, 2 ) )  
    allocate ( Id11       ( nrows11, ncols11 ) )
    allocate ( fDir11     ( nrows11, ncols11 ) )

    ! get grids of L11 
    Id11(:,:) =    UNPACK( L11_Id   ( iStart11 : iEnd11),  mask11, nodata_i4 )
    fDir11(:,:) =  UNPACK( L11_fDir ( iStart11 : iEnd11),  mask11, nodata_i4 )
    cellCoor11(:,:) = L11_cellCoor ( iStart11 : iEnd11, : )

    ! initialize
    nLinkFromN(:) = nodata_i4
    nLinkToN(:)   = nodata_i4

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
    call append( L11_fromN, nLinkFromN(:) )
    call append( L11_toN, nLinkToN(:)   )

    ! free space
    deallocate (mask11, cellCoor11, Id11, fDir11, nLinkFromN, nLinkToN )   

  end subroutine L11_set_network_topology

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_routing_order

  !     PURPOSE
  !>        \brief Find routing order, headwater cells and sink

  !>        \details Find routing order, headwater cells and sink. \n

  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego, Jan 2013 - modular version
  ! ------------------------------------------------------------------

  subroutine L11_routing_order(iBasin)

    use mo_global_variables, only: &
         L11_fromN,                & ! IN:    from node 
         L11_toN,                  & ! IN:    to node
         L11_rOrder,               & ! INOUT: network routing order
         L11_label,                & ! INOUT: label Id [0='', 1=HeadWater, 2=Sink]
         L11_sink,                 & ! INOUT: == .true. if sink node reached
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
    integer(i4), dimension(1)                 :: iSink
    logical                                   :: flag

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)

    nLinks  = nNodes - 1

    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a basin. 

    allocate ( nLinkFromN  ( nNodes ) )  ! all vectors valid from (1 : nLinks)
    allocate ( nLinkToN    ( nNodes ) )
    allocate ( nLinkROrder ( nNodes ) )
    allocate ( nLinkLabel  ( nNodes ) )
    allocate ( nLinkSink   ( nNodes ) )
    allocate ( netPerm     ( nNodes ) )

    ! get network vectors of L11 
    nLinkFromN(:) = L11_fromN ( iStart11 : iEnd11 )
    nLinkToN(:)   = L11_toN   ( iStart11 : iEnd11 )

    ! initialize
    nLinkROrder(1:nLinks) = 1
    nLinkROrder(nNodes)   = nodata_i4
    netPerm(:)            = nodata_i4
    nLinkSink(:)          = .FALSE.

    loop1: do ii = 1, nLinks
       loop2: do jj = 1, nLinks
          if ( jj == ii ) cycle loop2
          if ( nLinkFromN(ii) == nLinkToN(jj) ) then
             nLinkROrder(ii) = -9
          end if
          if ( nLinkROrder(ii) == -9 ) cycle loop1
       end do loop2
    end do loop1

    nLinkLabel(:) = 0  ! ''

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

    ! identify sink cell
    iSink = maxloc ( nLinkROrder( 1 : nLinks ) )
    nLinkLabel( iSink ) = 2    !  'Sink'
    nLinkSink(  iSink ) = .TRUE.

    ! keep routing order
    do ii = 1, nLinks
       netPerm( nLinkROrder(ii) ) = ii
    end do

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

  !>        \details
  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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
  ! ------------------------------------------------------------------

  subroutine L11_link_location(iBasin)

    use mo_global_variables, only: &
         basin,       &
         L0_fDir,     & ! IN:    flow direction (standard notation) L0
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
    integer(i4)                               :: iStart0, iEnd0
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
    integer(i4)                               :: ii, rr
    integer(i4)                               :: iNode, iRow, jCol
    integer(i4), dimension(2)                 :: oLoc           ! output location in L0

    ! level-0 information
    call get_basin_info (iBasin, 0, nrows0, ncols0, iStart=iStart0, iEnd=iEnd0, mask=mask0) 

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)

    nLinks  = nNodes - 1

    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a basin. 

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

    ! get fDir at L0
    fDir0(:,:) =   UNPACK( L0_fDir  (iStart0:iEnd0),  mask0, nodata_i4 )
    draSC0(:,:) =  UNPACK( L0_draSC (iStart0:iEnd0),  mask0, nodata_i4 )

    ! get network vectors of L11 
    nLinkFromN(:) = L11_fromN   ( iStart11 : iEnd11 )
    netPerm(:)    = L11_netPerm ( iStart11 : iEnd11 )
    rowOut(:)     = L11_rowOut  ( iStart11 : iEnd11 )
    colOut(:)     = L11_colOut  ( iStart11 : iEnd11 )  

    ! finding main outlet (row, col) in L0
    oLoc(1) = basin%L0_rowOutlet(iBasin)
    oLoc(2) = basin%L0_colOutlet(iBasin) 

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

       if(iRow == oLoc(1) .and. jCol == oLoc(2)) then

          nLinkToRow(ii) = iRow
          nLinkToCol(ii) = jCol

       else

          do while ( .not. ( draSC0(iRow,jCol) > 0 ) )
             call moveDownOneCell( fDir0(iRow,jcol), iRow, jCol )
             if ( iRow == oLoc(1) .and. jCol == oLoc(2)) exit
          end do
          ! set "to" cell (when an outlet is reached)
          nLinkToRow(ii) = iRow
          nLinkToCol(ii) = jCol

       end if
    end do

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
  !>        \brief

  !>        \details Perform the following tasks: \n
  !>        - Draining cell identification (cell at L0 to draining cell outlet at L11). 
  !>        - Set gauging nodes

  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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
  ! ------------------------------------------------------------------

  subroutine L11_set_drain_outlet_gauges(iBasin)

    use mo_global_variables, only: &
         basin,       & 
         L0_fDir,     & ! IN: flow direction (standard notation) L0
         L0_draSC,    & ! IN: Index of draining cell of each sub catchment (== cell L11)
         L0_cellCoor, & ! IN: cell coordinates (row,col) -> <only domain> input data
         L0_gaugeLoc, & ! IN: location of gauges (read with gauge Id then 
         !              !     transformed into gauge running ID => [1,nGaugesTotal]
         L0_L11_Id,   & ! IN: mapping of L11 Id on L0
         L0_draCell     ! INOUT: draining cell id at L11 of ith cell of L0

    implicit none

    integer(i4), intent(in)                   :: iBasin         ! basin 

    ! local
    integer(i4)                               :: nCells0
    integer(i4)                               :: nrows0, ncols0
    integer(i4)                               :: iStart0, iEnd0
    logical,     dimension(:,:), allocatable  :: mask0
    integer(i4), dimension(:,:), allocatable  :: cellCoor0
    integer(i4), dimension(:,:), allocatable  :: draSC0         
    integer(i4), dimension(:,:), allocatable  :: fDir0
    integer(i4), dimension(:,:), allocatable  :: gaugeLoc0      
    integer(i4), dimension(:,:), allocatable  :: draCell0
    integer(i4), dimension(:,:), allocatable  :: L11Id_on_L0
    integer(i4)                               :: ii, jj, kk
    integer(i4)                               :: iSc
    integer(i4)                               :: iRow, jCol
    integer(i4)                               :: gaugeCounter

    ! level-0 information
    call get_basin_info ( iBasin, 0, nrows0, ncols0, ncells=nCells0, &
         iStart=iStart0, iEnd=iEnd0, mask=mask0     ) 

    allocate ( cellCoor0   ( nCells0, 2 ) )  
    allocate ( draSC0      ( nrows0, ncols0 ) )
    allocate ( fDir0       ( nrows0, ncols0 ) )
    allocate ( gaugeLoc0   ( nrows0, ncols0 ) )
    allocate ( draCell0    ( nrows0, ncols0 ) )
    allocate ( L11Id_on_L0 ( nrows0, ncols0 ) )

    ! get L0 fields
    cellCoor0(:,:)  = L0_cellCoor(iStart0 : iEnd0, :)

    draSC0(:,:) =      UNPACK( L0_draSC    (iStart0:iEnd0),  mask0, nodata_i4 )
    fDir0(:,:) =       UNPACK( L0_fDir     (iStart0:iEnd0),  mask0, nodata_i4 )
    gaugeLoc0(:,:) =   UNPACK( L0_gaugeLoc (iStart0:iEnd0),  mask0, nodata_i4 )
    L11Id_on_L0(:,:) = UNPACK( L0_L11_Id   (iStart0:iEnd0),  mask0, nodata_i4 ) 

    draCell0(:,:) = nodata_i4

    gaugeCounter = 0

    do kk = 1, nCells0

       ii   = cellCoor0(kk,1)
       jj   = cellCoor0(kk,2)
       iSc = draSC0(ii,jj)
       ! find drainage path
       iRow = ii
       jCol = jj
       do while ( .NOT. iSC > 0 )
          ! move downstream
          call moveDownOneCell( fDir0(iRow,jCol), iRow, jCol )
          iSC = draSC0(iRow,jCol)
       end do
       draCell0(ii,jj) = iSC

       ! set gauging nodes !>> G0%ScId is Id of the routing cell at level-11
       if ( gaugeLoc0(ii,jj) /= nodata_i4 ) then 
          gaugeCounter = gaugeCounter + 1
          basin%gaugeNodeList( iBasin, gaugeCounter ) = L11Id_on_L0(ii,jj)
       end if
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L0 data sets 
    call append( L0_draCell,     PACK ( draCell0(:,:),  mask0)  ) 

    ! free space
    deallocate ( mask0, cellCoor0, draSC0, fDir0, gaugeLoc0, draCell0, L11Id_on_L0) 

  end subroutine  L11_set_drain_outlet_gauges

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_stream_features
  !     PURPOSE
  !>        \brief

  !>        \details Stream features (stream network and floodplain)\n

  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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
  ! ------------------------------------------------------------------

  subroutine L11_stream_features(iBasin)

    use mo_global_variables, only: &
         L0_Id,           & ! IN:    level-0 id
         L0_elev,         & ! IN:    elevation (sinks removed)  [m]
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
         L11_slope          ! INOUT: normalized average slope
    use mo_mhm_constants, only: nodata_i4, nodata_dp

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
    integer(i4), dimension(:,:), allocatable :: stack,append_chunk
    real(dp)                                 :: length
    integer(i4), dimension(:,:), allocatable :: nodata_i4_tmp
    real(dp),    dimension(:,:), allocatable :: nodata_dp_tmp

    ! level-0 information
    call get_basin_info ( iBasin, 0, nrows0, ncols0, ncells=nCells0, &
         iStart=iStart0, iEnd=iEnd0, mask=mask0     ) 

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nNodes, iStart=iStart11, iEnd=iEnd11)

    nLinks  = nNodes - 1

    allocate ( iD0           ( nrows0, ncols0 ) )
    allocate ( elev0         ( nrows0, ncols0 ) )
    allocate ( fDir0         ( nrows0, ncols0 ) )
    allocate ( areaCell0     ( nrows0, ncols0 ) )
    allocate ( streamNet0    ( nrows0, ncols0 ) )
    allocate ( floodPlain0   ( nrows0, ncols0 ) )

    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a basin.
    allocate ( stack             ( nNodes, 2 ) ) !>> stack(nNodes, 2)
    allocate ( append_chunk      ( 8,      2 ) )
    allocate ( netPerm           ( nNodes ) )  
    allocate ( nLinkFromRow      ( nNodes ) )
    allocate ( nLinkFromCol      ( nNodes ) )
    allocate ( nLinkToRow        ( nNodes ) )  
    allocate ( nLinkToCol        ( nNodes ) ) 
    allocate ( nLinkLength       ( nNodes ) )
    allocate ( nLinkAFloodPlain  ( nNodes ) )
    allocate ( nLinkSlope        ( nNodes ) )

    allocate(nodata_i4_tmp(nrows0,ncols0))
    allocate(nodata_dp_tmp(nrows0,ncols0))
    nodata_i4_tmp = nodata_i4
    nodata_dp_tmp = nodata_dp

    ! get L0 fields
    iD0(:,:) =         UNPACK( L0_Id   (iStart0:iEnd0),  mask0, nodata_i4_tmp )
    elev0(:,:) =       UNPACK( L0_elev (iStart0:iEnd0),  mask0, nodata_dp_tmp )
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

       call cellLength(fDir0(frow,fcol),  nLinkLength(ii) )
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
             stack = cshift(stack, SHIFT = 1, DIM = 1)
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
          call cellLength( fDir0(fRow,fCol), length )
          nLinkLength(ii) = nLinkLength(ii) + length

       end do

       ! stream bed slope
       nLinkSlope(ii) = ( nLinkSlope(ii) - elev0(frow, fcol) ) / nLinkLength(ii)

       if ( nLinkSlope(ii) <= 0.0001_dp) nLinkSlope(ii) = 0.0001_dp

       ! calculate area of floodplains (avoid overwriting)
       nLinkAFloodPlain(ii) = sum ( areaCell0(:,:),  mask = ( floodPlain0(:,:) == ii ) )
       !  old > real( count( floodPlain0(:,:,) == i), dp ) * areaCell0

    end do

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
         nLinkLength, nLinkAFloodPlain, nLinkSlope) 
    deallocate(nodata_i4_tmp,nodata_dp_tmp)    

  end subroutine L11_stream_features

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_fraction_sealed_floodplain

  !     PURPOSE
  !         \brief

  !>        \details Fraction of the flood plain with impervious cover (\ref fig_routing "Routing
  !>        Network"). This proportion is used to regionalize the Muskingum parameters.
  !>        Samaniego et al. \cite SB05 found out that this fraction is one of the statistically
  !>        significant predictor variables of peak discharge in mesoscale basins.\n

  !>        If a variable is added or removed here, then it also has to 
  !>        be added or removed in the subroutine L11_config_set in
  !>        module mo_restart and in the subroutine set_L11_config in module
  !>        mo_set_netcdf_restart

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
  ! ------------------------------------------------------------------
  subroutine L11_fraction_sealed_floodplain( &
       nLinks, LCover0, floodPlain0, & ! INTENT IN
       areaCell0, nLinkAFloodPlain,  & ! INTENT IN
       LCClassImp,                   & ! INTENT IN
       nLinkFracFPimp    )             ! INTENT OUT
    implicit none

    integer(i4),                intent(in)  :: nLinks
    integer(i4), dimension(:),  intent(in)  :: LCover0
    integer(i4), dimension(:),  intent(in)  :: floodPlain0
    real(dp),    dimension(:),  intent(in)  :: areaCell0
    real(dp),    dimension(:),  intent(in)  :: nLinkAFloodPlain
    integer(i4),                intent(in)  :: LCClassImp         ! e.g. = 2 (old code)
    real(dp),    dimension(:),  intent(out) :: nLinkFracFPimp  

    ! local
    integer(i4) :: ii

    do ii = 1, nLinks
       nLinkFracFPimp(ii) =  sum ( areaCell0(:),  & 
            mask = ( floodPlain0(:) == ii .and. LCover0(:) == LCClassImp ) ) &
            /  nLinkAFloodPlain(ii)
    end do

  end subroutine L11_fraction_sealed_floodplain

  ! ------------------------------------------------------------------

  !     NAME
  !         routing_dummy_alloc

  !     PURPOSE
  !>        \brief routing_dummy_alloc related to routing

  !>        \details Allocate L0 variable that are initialized
  !>        for routing, when routing is switched off. This is a dummy
  !>        allocation required for the mhm call. No initialization is
  !>        performed. These variables are all the variables that would
  !>        would have been initialized by the net startup or restart
  !>        L11 config, if routing would be switched on.

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
  !>        \author  Stephan Thober
  !>        \date    Sep 2013

  ! ------------------------------------------------------------------
  subroutine routing_dummy_alloc( iBasin )

    use mo_kind,             only: i4
    use mo_init_states,      only: get_basin_info
    use mo_global_variables, only: &
         L11_cellCoor,      & ! cell Coordinates at Level 11
         L11_Id,            & ! cell Ids at Level 11
         L11_nCells,        & ! Number of Cells at Level 11
         L0_draSC,          &
         L0_L11_Id,         &
         L1_L11_Id,         &
         L11_fDir,          &
         L11_rowOut,        &
         L11_colOut,        &
         L11_upBound_L0,    &
         L11_downBound_L0,  &
         L11_leftBound_L0,  &
         L11_rightBound_L0, &
         L11_upBound_L1,    &
         L11_downBound_L1,  &
         L11_leftBound_L1,  &
         L11_rightBound_L1, &
         L11_fromN,         &
         L11_toN,           &
         L11_rOrder,        &
         L11_label,         &
         L11_sink,          &
         L11_netPerm,       &
         L11_fRow,          &
         L11_fCol,          &
         L11_tRow,          &
         L11_tCol,          &
         L0_draCell,        &
         L0_streamNet,      &
         L0_floodPlain,     &
         L11_length,        &
         L11_aFloodPlain,   &
         L11_slope

    implicit none

    integer(i4), intent(in) :: iBasin

    ! local
    integer(i4) :: ncols0
    integer(i4) :: nrows0
    integer(i4) :: ncells0
    integer(i4) :: ncols1
    integer(i4) :: nrows1
    integer(i4) :: ncells1

    ! get L0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0, ncells=ncells0 )

    ! get L1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1 )

    ! L0 variables ---------------------------------------------------
    call extend( L0_draSC, ncells0 )
    call extend( L0_L11_Id, ncells0 )
    call extend( L0_draCell, ncells0 )
    call extend( L0_streamnet, ncells0 )
    call extend( L0_floodplain, ncells0 )
    ! L1 variables ---------------------------------------------------
    call extend( L1_L11_Id, ncells1 )
    call extend( L11_upBound_L1, ncells1 )
    call extend( L11_downBound_L1, ncells1 )
    call extend( L11_leftBound_L1, ncells1 )
    call extend( L11_rightBound_L1, ncells1 )

    ! L11 variables --------------------------------------------------
    if ( .not. allocated( L11_cellCoor) )        allocate( L11_cellCoor( 1, 1) )
    if ( .not. allocated( L11_Id ) )             allocate( L11_Id( 1 ) )
    L11_nCells = 1_i4
    if ( .not. allocated( L11_fDir ) )           allocate( L11_fDir( 1 ) )
    if ( .not. allocated( L11_rowOut ) )         allocate( L11_rowOut( 1 ) )
    if ( .not. allocated( L11_colOut ) )         allocate( L11_colOut( 1 ) )
    if ( .not. allocated( L11_upBound_L0 ) )     allocate( L11_upBound_L0( 1 ) )
    if ( .not. allocated( L11_downBound_L0 ) )   allocate( L11_downBound_L0( 1 ) )
    if ( .not. allocated( L11_leftBound_L0 ) )   allocate( L11_leftBound_L0( 1 ) )
    if ( .not. allocated( L11_rightBound_L0 ) )  allocate( L11_rightBound_L0( 1 ) )
    if ( .not. allocated( L11_fromN ) )          allocate( L11_fromN( 1 ) )
    if ( .not. allocated( L11_toN ) )            allocate( L11_toN( 1 ) )
    if ( .not. allocated( L11_rOrder) )          allocate( L11_rOrder( 1 ) )
    if ( .not. allocated( L11_label) )           allocate( L11_label( 1 ) )
    if ( .not. allocated( L11_sink) )            allocate( L11_sink( 1 ) )
    if ( .not. allocated( L11_netPerm) )         allocate( L11_netPerm( 1 ) )
    if ( .not. allocated( L11_fRow) )            allocate( L11_fRow( 1 ) )
    if ( .not. allocated( L11_fCol) )            allocate( L11_fCol( 1 ) )
    if ( .not. allocated( L11_tRow) )            allocate( L11_tRow( 1 ) )
    if ( .not. allocated( L11_tCol) )            allocate( L11_tCol( 1 ) )
    if ( .not. allocated( L11_length) )          allocate( L11_length( 1) )
    if ( .not. allocated( L11_aFloodPlain) )     allocate( L11_aFloodPlain( 1) )
    if ( .not. allocated( L11_slope) )           allocate( L11_slope( 1) )

  end subroutine routing_dummy_alloc

  ! -------------------------------------------------------------------
  ! extend allocated arrays
  ! -------------------------------------------------------------------
  subroutine extend( arr, ext_size)

    use mo_kind, only: i4

    implicit none

    integer(i4), dimension(:), allocatable, intent(inout) :: arr
    integer(i4),                            intent(in)    :: ext_size
    integer(i4)                                           :: old_size

    old_size = 0_i4
    if ( allocated( arr ) ) then
       old_size = size( arr )
       deallocate( arr )
    end if
    allocate( arr( old_size + ext_size ) )

  end subroutine extend

  ! ------------------------------------------------------------------
  !  MOVE UPSTREAM FROM-TO
  ! ------------------------------------------------------------------
  subroutine moveUp(elev0, fDir0, fi, fj, ss, nn)

    use mo_mhm_constants,    only: deltaH

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
       if ( (   fdir0(ii,jp) == 16        )                 .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
       if ( (   fdir0(ip,jp) == 32        )                 .and. &
            ( ( elev0(ip,jp)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(ip,jj)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(ip,jm)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(ii,jm)   - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)   - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(im,jm)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(im,jj)   - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)   - elev0(fi,fj) ) >= 0.0_dp )        &
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
            ( ( elev0(im,jp)  - elev0(fi,fj) ) <= deltaH )  .and. &
            ( ( elev0(ii,jp)  - elev0(fi,fj) ) >= 0.0_dp )        &
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
  subroutine cellLength(fDir, length)

    use mo_constants,        only: SQRT2_dp
    use mo_global_variables, only: level0

    implicit none

    integer(i4), intent(IN)  :: fDir
    real(dp),    intent(OUT) :: length

    select case (fDir)
    case(1, 4, 16, 64)       ! E, S, W, N
       length = 1.0_dp
    case(2, 8, 32, 128)      ! SE, SW, NW, NE
       length = SQRT2_dp
    end select

    length = length * level0%cellsize

  end subroutine cellLength

END MODULE mo_net_startup
