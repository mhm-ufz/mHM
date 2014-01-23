!> \file mo_restart.f90

!> \brief reading and writing states, fluxes and configuration for restart of mHM.

!> \details routines are seperated for reading and writing variables for:\n
!>          - states and fluxes, and \n
!>          - configuration.\n
!>          Reading of L11 configuration is also seperated from the rest, 
!>          since it is only required when routing is activated.

!> \authors Stephan Thober
!> \date Jul 2013

MODULE mo_restart

  ! This module is a restart for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Stephan Thober, Apr 2011

  IMPLICIT NONE

  PUBLIC :: read_restart_states     ! read restart files for state variables from a given path
  PUBLIC :: write_restart_states    ! write restart files for state variables to a given path

  PUBLIC :: read_restart_config     ! read restart files for configuration from a given path
  PUBLIC :: read_restart_L11_config ! read L11 configuration
  PUBLIC :: write_restart_config    ! write restart files for configuration to a given path

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------

  !      NAME
  !         read_restart_L11_config

  !     PURPOSE
  !>        \brief reads Level 11 configuration from a restart directory

  !>        \details read Level 11 configuration variables from a given restart
  !>        directory and initializes all Level 11 configuration variables,
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

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

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_config 

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013

  subroutine read_restart_L11_config( iBasin, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_kind,             only: i4, dp
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nodata_dp
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: basin, & ! basin database
         level11,           &
         resolutionRouting, &
         nBasins,           & ! Number of Basins
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

    integer(i4),    intent(in) :: iBasin
    character(256), intent(in) :: InPath ! list of Output paths per Basin

    ! local variables
    character(256)             :: Fname

    ! local variables
    integer(i4)                                          :: nrows0   ! Number of rows at level 0
    integer(i4)                                          :: ncols0   ! Number of cols at level 0
    logical, dimension(:,:), allocatable                 :: mask0    ! Mask at Level 0
    integer(i4)                                          :: nrows1   ! Number of rows at level 1
    integer(i4)                                          :: ncols1   ! Number of cols at level 1
    logical, dimension(:,:), allocatable                 :: mask1    ! Mask at Level 1
    integer(i4)                                          :: nrows11  ! Number of rows at level 11
    integer(i4)                                          :: ncols11  ! Number of cols at level 11
    logical, dimension(:,:), allocatable                 :: mask11   ! Mask at Level 11
    integer(i4)                                          :: nCells11 ! Number of cells at lev. 11
    real(dp)                                             :: xllcorner0, yllcorner0
    real(dp)                                             :: cellsize0

    ! DUMMY variables
    integer(i4), dimension(:),     allocatable           :: dummyI1  ! dummy, 1 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI2  ! dummy, 2 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp),    dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension DP

    ! set file name
    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_L11_config.nc'
    call message('    Reading Restart-file: ', trim(adjustl(Fname)),' ...')

    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0,&
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0, mask=mask0 )

    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1, mask=mask1 )

    ! calculate l11 grid resolutionRouting
    if(iBasin .eq. 1) then
       allocate( level11%nrows     (nBasins) )
       allocate( level11%ncols     (nBasins) )
       allocate( level11%xllcorner (nBasins) )
       allocate( level11%yllcorner (nBasins) )
    end if
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,            &
         resolutionRouting , &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin), &
         level11%yllcorner(iBasin), level11%cellsize, level11%nodata_value        )

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11)

    ! read L11 mask
    allocate( dummyI2( nrows11, ncols11 ), mask11( nrows11, ncols11) )
    call Get_NcVar( Fname,  'L11_basin_Mask', dummyI2 )
    mask11 = (dummyI2 .eq. 1_i4)
    call append( basin%L11_Mask, reshape( mask11, (/nrows11*ncols11/)))

    ! get Number of cells
    nCells11 = count( dummyI2 .eq. 1_i4 )
    deallocate( dummyI2 )

    ! update basin database
    if (iBasin .eq. 1) then

       !
       allocate(basin%L11_iStart     (nBasins))
       allocate(basin%L11_iEnd       (nBasins))
       allocate(basin%L11_iStartMask (nBasins))
       allocate(basin%L11_iEndMask   (nBasins))    

       ! basin information
       basin%L11_iStart(iBasin) = 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells11 - 1

       basin%L11_iStartMask(iBasin) = 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    else

       ! basin information
       basin%L11_iStart(iBasin) = basin%L11_iEnd(iBasin-1) + 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells11 - 1

       basin%L11_iStartMask(iBasin) = basin%L11_iEndMask(iBasin-1) + 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    end if

    ! read L11 cellCoor
    allocate( dummyI2( nrows11, ncols11 ) )
    allocate( dummyI22( count(mask11), 2 ))
    call Get_NcVar( Fname, 'L11_rowCoor', dummyI2 )
    dummyI22(:,1) = pack( dummyI2, mask11 )
    call Get_NcVar( Fname, 'L11_colCoor', dummyI2 )
    dummyI22(:,2) = pack( dummyI2, mask11 )
    call append( L11_cellCoor, dummyI22)
    deallocate( dummyI22 )

    ! read L11 IDs
    call Get_NcVar( Fname, 'L11_Id', dummyI2 )
    call append( L11_Id, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! update Number of cells at Level 11
    L11_nCells = size( L11_Id, 1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L0 outlet Coordinates
    allocate(dummyI1(2))
    call Get_NcVar( Fname, 'L0_OutletCoord', dummyI1)

    ! allocate space for row and col Outlet
    if (iBasin .eq. 1) then
       allocate( basin%L0_rowOutlet(nBasins) ) 
       allocate( basin%L0_colOutlet(nBasins) )
    end if

    ! L0 data sets
    basin%L0_rowOutlet(iBasin) = dummyI1(1)
    basin%L0_colOutlet(iBasin) = dummyI1(2)
    deallocate(dummyI1)

    ! read L0 draining cell index
    allocate(dummyI2(nrows0,ncols0))
    call Get_NcVar( Fname, 'L0_draSC', dummyI2)
    call append( L0_draSC,     PACK ( dummyI2,  mask0)  )

    ! Mapping of L11 Id on L0
    call Get_NcVar( Fname, 'L0_L11_Id', dummyI2 )
    call append( L0_L11_Id,    PACK ( dummyI2, mask0)  )
    deallocate( dummyI2 )

    ! L1 data sets
    ! Mapping of L11 Id on L1
    allocate( dummyI2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L1_L11_Id', dummyI2 )
    call append( L1_L11_Id,    PACK ( dummyI2, mask1)  )
    deallocate( dummyI2 )

    ! L11 data sets
    ! Flow direction (standard notation)
    allocate( dummyI2( nrows11, ncols11) )
    call Get_NcVar( Fname, 'L11_fDir', dummyI2 )
    call append( L11_fDir,     PACK ( dummyI2,      mask11) )
    deallocate( dummyI2 )

    ! Grid vertical location of the Outlet
    allocate( dummyI2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_rowOut', dummyI2 )
    call append( L11_rowOut, pack( dummyI2, mask11 ) )

    ! Grid horizontal location  of the Outlet
    call Get_NcVar( Fname, 'L11_colOut', dummyI2 )
    call append( L11_colOut, pack( dummyI2, mask11 ) )

    ! Row start at finer level-0 scale
    call Get_NcVar( Fname, 'L11_upBound_L0', dummyI2 )
    call append( L11_upBound_L0, pack( dummyI2, mask11 ) )

    ! Row end at finer level-0 scale
    call Get_NcVar( Fname, 'L11_downBound_L0', dummyI2 )
    call append( L11_downBound_L0, pack( dummyI2, mask11) )

    ! Col start at finer level-0 scale
    call Get_NcVar( Fname, 'L11_leftBound_L0', dummyI2 )
    call append( L11_leftBound_L0, pack( dummyI2, mask11) )

    ! Col end at finer level-0 scale 
    call Get_NcVar( Fname, 'L11_rightBound_L0', dummyI2 )
    call append( L11_rightBound_L0, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! Row start at finer level-1 scale
    allocate( dummyI2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L11_upBound_L1', dummyI2 )
    call append( L11_upBound_L1, pack( dummyI2, mask1) )

    ! Row end at finer level-1 scale
    call Get_NcVar( Fname, 'L11_downBound_L1', dummyI2 )
    call append( L11_downBound_L1, pack( dummyI2, mask1) )

    ! Col start at finer level-1 scale
    call Get_NcVar( Fname, 'L11_leftBound_L1', dummyI2 )
    call append( L11_leftBound_L1, pack( dummyI2, mask1) )

    ! Col end at finer level-1 scale 
    call Get_NcVar( Fname, 'L11_rightBound_L1', dummyI2 )
    call append( L11_rightBound_L1, pack(dummyI2, mask1) ) 
    deallocate( dummyI2 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From Node
    allocate( dummyI2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_fromN', dummyI2)
    call append( L11_fromN, pack( dummyI2, mask11) )

    ! To Node
    call Get_NcVar( Fname, 'L11_toN', dummyI2 )
    call append( L11_toN, pack( dummyI2, mask11) )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Network routing order
    call Get_NcVar( Fname, 'L11_rOrder', dummyI2 )
    call append( L11_rOrder, pack(dummyI2, mask11) )

    ! Label Id [0='', 1=HeadWater, 2=Sink]
    call Get_NcVar( Fname, 'L11_label', dummyI2 )
    call append( L11_label, pack( dummyI2, mask11) )

    ! .true. if sink node reached
    call Get_NcVar( Fname, 'L11_sink', dummyI2 )
    call append( L11_sink, (pack(dummyI2,mask11) .eq. 1_i4) )

    ! Routing sequence (permutation of L11_rOrder)
    call Get_NcVar( Fname, 'L11_netPerm', dummyI2 )
    call append( L11_netPerm, pack( dummyI2, mask11) )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From row in L0 grid
    call Get_NcVar( Fname, 'L11_fRow', dummyI2 )
    call append( L11_fRow, pack( dummyI2, mask11) )

    ! From col in L0 grid
    call Get_NcVar( Fname, 'L11_fCol', dummyI2 )
    call append( L11_fCol, pack( dummyI2, mask11) )

    ! To row in L0 grid
    call Get_NcVar( Fname, 'L11_tRow', dummyI2 )
    call append( L11_tRow, pack( dummyI2, mask11) )

    ! To col in L0 grid
    call Get_NcVar( Fname, 'L11_tCol', dummyI2 )
    call append( L11_tCol, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate( dummyI2( nrows0, ncols0 ) )
    call Get_NcVar( Fname, 'L0_draCell', dummyI2 )
    call append( L0_draCell,     PACK ( dummyI2,  mask0)  ) 

    ! read gaugenodelist
    allocate(dummyI1( size(basin%gaugeNodeList(iBasin,:))))
    call Get_NcVar( Fname, 'gaugeNodeList', dummyI1)
    basin%gaugeNodeList( iBasin, : ) = dummyI1
    deallocate(dummyI1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! L0 data sets 
    ! Stream network
    call Get_NcVar( Fname, 'L0_streamNet', dummyI2 )
    call append( L0_streamNet,    PACK ( dummyI2,  mask0)   ) 

    ! Floodplains of stream i
    call Get_NcVar( Fname, 'L0_floodPlain', dummyI2 )
    call append( L0_floodPlain,   PACK ( dummyI2,  mask0)  ) 
    deallocate( dummyI2 )

    ! L11 network data sets
    ! [m]     Total length of river link
    allocate( dummyD2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_length', dummyD2 )
    call append( L11_length, pack( dummyD2, mask11 ) )

    ! [m2]    Area of the flood plain
    call Get_NcVar( Fname, 'L11_aFloodPlain', dummyD2 )
    call append( L11_aFloodPlain, pack( dummyD2, mask11) )

    ! Average slope of river link
    call Get_NcVar( Fname, 'L11_slope', dummyD2 )
    call append( L11_slope, pack( dummyD2, mask11 ) )
    deallocate( dummyD2 )

  end subroutine read_restart_L11_config

  ! ------------------------------------------------------------------

  !      NAME
  !         read_restart_config

  !     PURPOSE
  !>        \brief reads configuration apart from Level 11 configuration
  !>        from a restart directory

  !>        \details read configuration variables from a given restart
  !>        directory and initializes all configuration variables,
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

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

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_config 

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013

  subroutine read_restart_config( iBasin, soilId_isPresent, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nodata_dp
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: L0_cellCoor    , & 
         L0_Id         , & ! Ids of grid at level-0 
         L0_areaCell   , & ! Ids of grid at level-0
         L0_slope_emp  , & ! Empirical quantiles of slope
         basin, & 
         nBasins, &
         level1, &
         L0_soilId, &
         L0_nCells, &
         nSoilTypes, &
         resolutionHydrology, &
         L1_nCells,      &
         L1_Id         , & ! Ids of grid at level-1
         L1_cellCoor   , &
         L1_upBound_L0 , & ! Row start at finer level-0 scale 
         L1_downBound_L0, & ! Row end at finer level-0 scale
         L1_leftBound_L0, & ! Col start at finer level-0 scale
         L1_rightBound_L0, & ! Col end at finer level-0 scale
         L1_areaCell   , & ! [km2] Effective area of cell at this level
         L1_nTCells_L0     ! Total number of valid L0 cells in a given L1 cell

    implicit none

    !
    integer(i4), intent(in) :: iBasin
    integer(i4), dimension(:), allocatable, intent(inout) :: soilId_isPresent
    character(256), intent(in)     :: InPath             ! list of Output paths per Basin

    !
    integer(i4)                                          :: nrows0   ! Number of rows at level 0
    integer(i4)                                          :: ncols0   ! Number of cols at level 
    integer(i4)                                          :: iStart0, iEnd0
    logical, dimension(:,:), allocatable                 :: mask0    ! Mask at Level 0
    integer(i4)                                          :: nrows1   ! Number of rows at level 1
    integer(i4)                                          :: ncols1   ! Number of cols at level 1
    logical, dimension(:,:), allocatable                 :: mask1    ! Mask at Level 1
    real(dp)                                             :: xllcorner0, yllcorner0
    real(dp)                                             :: cellsize0
    !
    ! Dummy Variables
    integer(i4)                                          :: k, j
    integer(i4), dimension(:,:),   allocatable           :: dummyI2  ! dummy, 2 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp),    dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension DP 

    ! local variables
    character(256) :: Fname

    ! read config
    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_config.nc'
    call message('    Reading Restart-file: ', trim(adjustl(Fname)),' ...')
 
    !
    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0, iStart= iStart0, iEnd=iEnd0, mask=mask0, &
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0  )
    !
    if (iBasin .eq. 1) then
       allocate( level1%nrows     (nBasins) )
       allocate( level1%ncols     (nBasins) )
       allocate( level1%xllcorner (nBasins) )
       allocate( level1%yllcorner (nBasins) )
    end if
    ! grid properties
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,         &
         resolutionHydrology , &
         level1%nrows(iBasin), level1%ncols(iBasin), level1%xllcorner(iBasin), &
         level1%yllcorner(iBasin), level1%cellsize, level1%nodata_value        )
    !
    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L0 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! read L0 cellCoor
    allocate( dummyI2( nrows0, ncols0 ) )
    allocate( dummyI22( count(mask0), 2 ))
    call Get_NcVar( Fname, 'L0_rowCoor', dummyI2 )
    dummyI22(:,1) = pack( dummyI2, mask0 )
    call Get_NcVar( Fname, 'L0_colCoor', dummyI2 )
    dummyI22(:,2) = pack( dummyI2, mask0 )
    call append( L0_cellCoor, dummyI22)
    deallocate( dummyI22 )
    !
    call Get_NcVar( Fname, 'L0_Id', dummyI2)
    call append( L0_Id, pack(dummyI2, mask0) )
    deallocate( dummyI2 )
    !
    allocate( dummyD2( nrows0, ncols0 ) )
    call Get_NcVar( Fname, 'L0_areaCell', dummyD2 )
    call append( L0_areaCell, pack(dummyD2, mask0) )
    !
    call Get_NcVar( Fname, 'L0_slope_emp', dummyD2 )
    call append( L0_slope_emp, pack(dummyD2, mask0) )
    deallocate( dummyD2 )

    L0_nCells = size(L0_Id,1)

    !------------------------------------------------------
    ! Assign whether a given soil type is present or not
    !------------------------------------------------------
    if ( iBasin .eq. 1 ) then
       allocate( soilId_isPresent(nSoilTypes) )
       soilId_isPresent(:) = 0
    end if

    do k = iStart0, iEnd0
       j = L0_soilId(k)
       soilId_isPresent(j) = 1
    end do
    !
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L1 mask
    allocate( dummyI2( nrows1, ncols1 ), mask1( nrows1, ncols1) )
    call Get_NcVar( Fname,  'L1_basin_Mask', dummyI2 )
    mask1 = (dummyI2 .eq. 1_i4)
    call append( basin%L1_Mask, reshape( mask1, (/nrows1*ncols1/)))
    !
    ! update basin database
    if( iBasin .eq. 1 ) then

       allocate(basin%L1_iStart(nBasins))
       allocate(basin%L1_iEnd  (nBasins))
       allocate(basin%L1_iStartMask(nBasins))
       allocate(basin%L1_iEndMask   (nBasins))    

       ! basin information
       basin%L1_iStart(iBasin) = 1
       basin%L1_iEnd  (iBasin) = basin%L1_iStart(iBasin) + count(mask1) - 1

       basin%L1_iStartMask(iBasin) = 1
       basin%L1_iEndMask  (iBasin) = basin%L1_iStartMask(iBasin) + nrows1*ncols1 - 1

    else

       ! basin information
       basin%L1_iStart(iBasin) = basin%L1_iEnd(iBasin-1) + 1
       basin%L1_iEnd  (iBasin) = basin%L1_iStart(iBasin) + count(mask1) - 1

       basin%L1_iStartMask(iBasin) = basin%L1_iEndMask(iBasin-1) + 1
       basin%L1_iEndMask  (iBasin) = basin%L1_iStartMask(iBasin) + nrows1*ncols1 - 1

    end if
    !
    ! L1 cell Ids
    call Get_NcVar( Fname, 'L1_Id', dummyI2)
    call append( L1_Id, pack(dummyI2, mask1) )
    ! L1 cell coordinates
    allocate( dummyI22( count(mask1), 2 ))
    call Get_NcVar( Fname, 'L1_rowCoor', dummyI2 )
    dummyI22(:,1) = pack( dummyI2, mask1 )
    call Get_NcVar( Fname, 'L1_colCoor', dummyI2 )
    dummyI22(:,2) = pack( dummyI2, mask1 )
    call append( L1_cellCoor, dummyI22)
    deallocate( dummyI22 )
    ! 
    call Get_NcVar( Fname, 'L1_upBound_L0', dummyI2)
    call append( L1_upBound_L0   , pack( dummyI2, mask1) )
    !
    call Get_NcVar( Fname, 'L1_downBound_L0', dummyI2 )
    call append( L1_downBound_L0 , pack( dummyI2, mask1)  )
    !
    call Get_NcVar( Fname, 'L1_leftBound_L0', dummyI2 )
    call append( L1_leftBound_L0 , pack( dummyI2, mask1)  )
    !
    call Get_NcVar( Fname, 'L1_rightBound_L0', dummyI2 )
    call append( L1_rightBound_L0, pack( dummyI2, mask1) )
    !
    call Get_NcVar( Fname, 'L1_nTCells_L0', dummyI2 )
    call append( L1_nTCells_L0   , pack( dummyI2, mask1)    )
    deallocate( dummyI2)
    !
    allocate( dummyD2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L1_areaCell', dummyD2 )
    call append( L1_areaCell     , pack( dummyD2, mask1)   )

    L1_nCells = size( L1_Id, 1 )

  end subroutine read_restart_config

  ! ------------------------------------------------------------------

  !      NAME
  !         read_restart_states

  !     PURPOSE
  !>        \brief reads fluxes and state variables from file

  !>        \details read fluxes and state variables from given 
  !>        restart directory and initialises all state variables
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

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

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_states 

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013
  !         Modified   R. Kumar, J. Mai    Sep. 2013  - Splitting allocation and initialization of arrays

  subroutine read_restart_states( iBasin, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nRoutingStates
    use mo_global_variables, only: processMatrix, &
         L1_fSealed, &
         L1_fForest, &
         L1_fPerm, &
         L1_Inter, &
         L1_snowPack, &
         L1_sealSTW, &
         L1_soilMoist, &
         L1_unsatSTW, &
         L1_satSTW, &
         L1_aETSoil, &
         L1_aETCanopy, &
         L1_aETSealed, &
         L1_baseflow, &
         L1_infilSoil, &
         L1_fastRunoff, &
         L1_melt, &
         L1_percol, &
         L1_preEffect, &
         L1_rain, &
         L1_runoffSeal, &
         L1_slowRunoff, &
         L1_snow, &
         L1_Throughfall, &
         L1_total_runoff, &
         L1_alpha, &
         L1_degDayInc, &
         L1_degDayMax, &
         L1_degDayNoPre, &
         L1_degDay, &
         L1_karstLoss, &
         L1_fAsp, &
         L1_fRoots, &
         L1_maxInter, &
         L1_kfastFlow, &
         L1_kSlowFlow, &
         L1_kBaseFlow, &
         L1_kPerco, &
         L1_soilMoistFC, &
         L1_soilMoistSat, &
         L1_soilMoistExp, &
         L1_tempThresh, &
         L1_unsatThresh, &
         L1_sealedThresh, &
         L1_wiltingPoint, &
         L11_Qmod, &
         L11_qOUT, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp, &
         nSoilHorizons_mHM

    implicit none

    integer(i4),    intent(in) :: iBasin
    character(256), intent(in) :: InPath ! list of Output paths per Basin

    character(256)                                    :: Fname
    integer(i4)                                       :: ii       ! loop index
    integer(i4)                                       :: s1       ! start index at level 1
    integer(i4)                                       :: e1       ! end index at level 1
    integer(i4)                                       :: ncols1   ! number of colums at level 1
    integer(i4)                                       :: nrows1   ! number of rows at level 1
    integer(i4)                                       :: ncells1  ! number of cells at level 1
    logical, dimension(:,:), allocatable              :: mask1    ! mask at level 1
    integer(i4)                                       :: s11      ! start index at level 11
    integer(i4)                                       :: e11      ! end index at level 11
    integer(i4)                                       :: ncols11  ! number of colums at level 11
    integer(i4)                                       :: nrows11  ! number of rows at level 11
    integer(i4)                                       :: ncells11 ! number of cells at level 11
    logical, dimension(:,:), allocatable              :: mask11   ! mask at level 11

    real(dp), dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable           :: dummyD3  ! dummy, 3 dimension

    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_states.nc'
    call message('    Reading Restart-file: ', trim(adjustl(Fname)),' ...')

    ! get basin information at level 1
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 )

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------
    allocate( dummyD2( nrows1, ncols1 ) )

    call Get_NcVar( Fname,  'L1_fSealed', dummyD2 )
    L1_fSealed(s1:e1) = pack( dummyD2, mask1 )

    call Get_NcVar( Fname,  'L1_fForest', dummyD2 )
    L1_fForest(s1:e1) = pack( dummyD2, mask1 )

    call Get_NcVar( Fname,  'L1_fPerm', dummyD2 )
    L1_fPerm(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------

    ! Interception
    call Get_NcVar( Fname,  'L1_Inter', dummyD2 )
    L1_inter(s1:e1) = pack( dummyD2, mask1 )

    ! Snowpack
    call Get_NcVar( Fname,  'L1_snowPack', dummyD2 )
    L1_snowPack(s1:e1) = pack( dummyD2, mask1 )

    ! Retention storage of impervious areas
    call Get_NcVar( Fname,  'L1_sealSTW', dummyD2 )
    L1_sealSTW(s1:e1) = pack( dummyD2, mask1 )

    ! upper soil storage
    call Get_NcVar( Fname,  'L1_unsatSTW', dummyD2 )
    L1_unsatSTW(s1:e1) = pack( dummyD2, mask1 )

    ! groundwater storage
    call Get_NcVar( Fname,  'L1_satSTW', dummyD2 )
    L1_satSTW(s1:e1) = pack( dummyD2, mask1 )

    ! Soil moisture of each horizon
    deallocate( dummyD2 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1, nSoilHorizons_mHM ) )

    call Get_NcVar( Fname,  'L1_soilMoist', dummyD3 )

    do ii = 1, nSoilHorizons_mHM
       L1_soilMoist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------   

    !  soil actual ET
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_aETSoil', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_aETSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do


    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! canopy actual ET
    call Get_NcVar( Fname,  'L1_aETCanopy', dummyD2 )
    L1_aETCanopy(s1:e1) = pack( dummyD2, mask1 ) 

    ! sealed area actual ET
    call Get_NcVar( Fname,  'L1_aETSealed', dummyD2 )
    L1_aETSealed(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow
    call Get_NcVar( Fname,  'L1_baseflow', dummyD2 )
    L1_baseflow(s1:e1) = pack( dummyD2, mask1 ) 

    ! soil in-exfiltration
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_infilSoil', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_infilSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! fast runoff
    call Get_NcVar( Fname,  'L1_fastRunoff', dummyD2 )
    L1_fastRunoff(s1:e1) = pack( dummyD2, mask1 )

    ! snow melt
    call Get_NcVar( Fname,  'L1_melt', dummyD2 )
    L1_melt(s1:e1) = pack( dummyD2, mask1 )

    ! percolation
    call Get_NcVar( Fname,  'L1_percol', dummyD2 )
    L1_percol(s1:e1) = pack( dummyD2, mask1 ) 

    ! effective precip. depth (snow melt + rain)
    call Get_NcVar( Fname,  'L1_preEffect', dummyD2 )
    L1_preEffect(s1:e1) = pack( dummyD2, mask1 )

    ! rain (liquid water)
    call Get_NcVar( Fname,  'L1_rain', dummyD2 )
    L1_rain(s1:e1) = pack( dummyD2, mask1 ) 

    ! runoff from impervious area
    call Get_NcVar( Fname,  'L1_runoffSeal', dummyD2 )
    L1_runoffSeal(s1:e1) = pack( dummyD2, mask1 )

    ! slow runoff
    call Get_NcVar( Fname,  'L1_slowRunoff', dummyD2 )
    L1_slowRunoff(s1:e1) = pack( dummyD2, mask1 ) 

    ! snow (solid water)
    call Get_NcVar( Fname,  'L1_snow', dummyD2 )
    L1_snow(s1:e1) = pack( dummyD2, mask1 )

    ! throughfall 
    call Get_NcVar( Fname,  'L1_Throughfall', dummyD2 )
    L1_Throughfall(s1:e1) = pack( dummyD2, mask1 )

    ! total runoff
    call Get_NcVar( Fname,  'L1_total_runoff', dummyD2 )
    L1_total_runoff(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    call Get_NcVar( Fname,  'L1_alpha', dummyD2 )
    L1_alpha(s1:e1) = pack( dummyD2, mask1 ) 

    ! increase of the Degree-day factor per mm of increase in precipitation
    call Get_NcVar( Fname,  'L1_degDayInc', dummyD2 )
    L1_degDayInc(s1:e1) = pack( dummyD2, mask1 ) 

    ! maximum degree-day factor 
    call Get_NcVar( Fname,  'L1_degDayMax', dummyD2 )
    L1_degDayMax(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor with no precipitation
    call Get_NcVar( Fname,  'L1_degDayNoPre', dummyD2 )
    L1_degDayNoPre(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor
    call Get_NcVar( Fname,  'L1_degDay', dummyD2 )
    L1_degDay(s1:e1) = pack( dummyD2, mask1 ) 

    ! Karstic percolation loss
    call Get_NcVar( Fname,  'L1_karstLoss', dummyD2 )
    L1_karstLoss(s1:e1) = pack( dummyD2, mask1 ) 

    ! PET correction factor due to terrain aspect
    call Get_NcVar( Fname,  'L1_fAsp', dummyD2 )
    L1_fAsp(s1:e1) = pack( dummyD2, mask1 ) 

    ! Fraction of roots in soil horizons    
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_fRoots', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_fRoots(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Maximum interception 
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD2( nrows1, ncols1 ) )

    call Get_NcVar( Fname, 'L1_maxInter', dummyD2 )
    L1_maxInter(s1:e1) = pack(dummyD2, mask1) 

    ! fast interflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kfastFlow', dummyD2 )
    L1_kfastFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! slow interflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kSlowFlow', dummyD2 )
    L1_kSlowFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kBaseFlow', dummyD2 )
    L1_kBaseFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! percolation coefficient
    call Get_NcVar( Fname,  'L1_kPerco', dummyD2 )
    L1_kPerco(s1:e1) = pack( dummyD2, mask1 ) 

    ! Soil moisture below which actual ET is reduced linearly till PWP
    deallocate( dummyD2 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistFC', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistFC(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Saturation soil moisture for each horizon [mm]
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistSat', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistSat(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Exponential parameter to how non-linear is the soil water retention
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistExp', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistExp(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do


    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! Threshold temperature for snow/rain 
    call Get_NcVar( Fname,  'L1_tempThresh', dummyD2 )
    L1_tempThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Threshhold water depth controlling fast interflow
    call Get_NcVar( Fname,  'L1_unsatThresh', dummyD2 )
    L1_unsatThresh(s1:e1) = pack( dummyD2, mask1 )

    ! Threshhold water depth for surface runoff in sealed surfaces
    call Get_NcVar( Fname,  'L1_sealedThresh', dummyD2 )
    L1_sealedThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Permanent wilting point
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_wiltingPoint', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_wiltingPoint(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do


    deallocate( dummyD2, dummyD3 )

    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------
    if (  processMatrix(8, 1) /= 0 ) then

       ! level-11 information
       call get_basin_info( iBasin, 11, nrows11, ncols11, ncells=ncells11, &
            iStart=s11, iEnd=e11, mask=mask11 )
       allocate( dummyD2( nrows11, ncols11 ) )

       ! simulated discharge at each node
       call Get_NcVar( Fname,  'L11_Qmod', dummyD2 )
       L11_Qmod(s11:e11) = pack( dummyD2, mask11 ) 

       ! Total outflow from cells L11 at time tt
       call Get_NcVar( Fname,  'L11_qOUT', dummyD2 )
       L11_qOUT(s11:e11) = pack( dummyD2, mask11 )

       ! Total discharge inputs at t-1 and t
       deallocate( dummyD2 )
       allocate( dummyD3( nrows1, ncols1, nRoutingStates ) )
       allocate( dummyD2( ncells11, nRoutingStates ) )

       call Get_NcVar( Fname, 'L11_qTIN', dummyD3 )
       do ii = 1, nRoutingStates
          L11_qTIN(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
       end do


       !  Routed outflow leaving a node
       deallocate( dummyD2, dummyD3 )
       allocate( dummyD3( nrows1, ncols1, nRoutingStates ) )
       allocate( dummyD2( ncells11,  nRoutingStates ) )

       call Get_NcVar( Fname, 'L11_qTR', dummyD3 )
       do ii = 1, nRoutingStates
          L11_qTR(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
       end do

       deallocate(dummyD2)
       allocate( dummyD2( nrows11, ncols11 ) )

       ! kappa: Muskingum travel time parameter.
       call Get_NcVar( Fname,  'L11_K', dummyD2 )
       L11_K(s11:e11) = pack( dummyD2, mask11 )

       !  xi:    Muskingum diffusion parameter
       call Get_NcVar( Fname,  'L11_xi', dummyD2 )
       L11_xi(s11:e11) = pack( dummyD2, mask11 )

       ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
       call Get_NcVar( Fname,  'L11_C1', dummyD2 )
       L11_C1(s11:e11) = pack( dummyD2, mask11 )

       ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
       call Get_NcVar( Fname,  'L11_C2', dummyD2 )
       L11_C2(s11:e11) = pack( dummyD2, mask11 )

       ! Fraction of the flood plain with impervious cover
       call Get_NcVar( Fname,  'L11_FracFPimp', dummyD2 )
       L11_FracFPimp(s11:e11) = pack( dummyD2, mask11 )

       ! free memory
       deallocate( dummyD2, dummyD3 )

    end if

  end subroutine read_restart_states

  ! ------------------------------------------------------------------

  !      NAME
  !         write_restart_states

  !     PURPOSE
  !>        \brief writes states to a restart file

  !>        \details writes all state variables and fluxes to file in 
  !>        a given restart directory

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath" Output Path including trailing slash
  !>                                                             for each basin (size(OutPath)==nBasins)

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

  !     RESTRICTIONS
  !>        \note file name cannot be specified, it is <number of basin>_states.nc, with
  !>        number of basins as three digits

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013
  !         modified Stephan Thober,    Nov 2013 - only writting L11 variables if routing is switched on

  subroutine write_restart_states( OutPath )

    use mo_kind,             only: i4
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_ncwrite,          only: create_netcdf, write_static_netcdf, close_netcdf
    use mo_global_variables, only: processMatrix

    implicit none

    character(256)                           :: Fname
    character(256), dimension(:), intent(in) :: OutPath ! list of Output paths per Basin
    integer(i4)                              :: hnc     ! handle for netcdf file
    integer(i4)                              :: i

    do i = 1, size(OutPath) 

       Fname = trim(OutPath(i)) // trim(num2str(i, '(i3.3)')) // '_states.nc'
       call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
       call state_variables_set( i )
       ! call create_netcdf( trim(Fname), hnc, lfs = .true., netcdf4 = .false.)
       call create_netcdf( trim(Fname), hnc, netcdf4 = .true.)
       call write_static_netcdf( hnc)
       call close_netcdf(hnc)

       ! free memory
       call free_memory_states( (processMatrix(8,1) .ne. 0))

    end do

  end subroutine write_restart_states

  ! ------------------------------------------------------------------

  !      NAME
  !         write_restart_config

  !     PURPOSE
  !>        \brief writes configuration variables to a restart file

  !>        \details writes configuration variables EXCLUDING state
  !>        variables to a file in a given restart directory. Two files
  !>        are written: 1) for L11_configuration containing configuration
  !>        variables for routing and 2) other configuration variables

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath" Output Path including trailing slash
  !>                                                             for each basin (size(OutPath)==nBasins)

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

  !     RESTRICTIONS
  !>        \note file names cannot be specified, it is <number of basin>_L11_config.nc 
  !>        and <number of basin>_config.nc, respectively. L11_config is only written, 
  !>        when routing is switched on in the process matrix

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Jul 2013
  subroutine write_restart_config( OutPath )

    use mo_kind,         only: i4
    use mo_message,      only: message
    use mo_string_utils, only: num2str
    use mo_ncwrite,      only: create_netcdf, write_static_netcdf, close_netcdf
    use mo_global_variables, only: processMatrix, nBasins

    implicit none

    character(256)                           :: Fname
    character(256), dimension(:), intent(in) :: OutPath ! list of Output paths per Basin
    integer(i4)                              :: hnc     ! handle for netcdf file
    integer(i4)                              :: i

    if ( size( OutPath ) .ne. nBasins ) then
       stop 'write_restart_config: size of OutPath does not match Number of Basins'
    End if

    do i = 1, nBasins

       Fname = trim(OutPath(i)) // trim(num2str(i, '(i3.3)')) // '_config.nc'
       call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
       call config_variables_set( i )
       call create_netcdf( trim(Fname), hnc, netcdf4 = .true.)
       call write_static_netcdf( hnc)
       call close_netcdf(hnc)

       ! free memory
       call free_memory_config

       ! if L11 routing is activated, then write it to an OutFile
       ! L11: network initialization
       if ( processMatrix(8, 1) .ne. 0 ) then

          Fname = trim(OutPath(i)) // trim(num2str(i, '(i3.3)')) // '_L11_config.nc'
          call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
          call L11_config_set( i )
          call create_netcdf( trim(Fname), hnc, netcdf4 = .true.)
          call write_static_netcdf( hnc)
          call close_netcdf(hnc)

          ! free memory
          call free_memory_L11_config

       end if

    end do

  end subroutine write_restart_config
  ! ------------------------------------------------------------------

  !      NAME
  !         config_variables_set

  !     PURPOSE
  !        sets netcdf target variables for all mHM configure variables. These target variables
  !        are defined in the module mo_set_netcdf_restart. The actual setting of the netcdf
  !        output variable V, that points to the target variables here, takes place in the 
  !        subroutine set_config contained in mo_set_netcdf_restart. In general each variable
  !        defined in the subroutine L0_variable_init and L1_variable_init (mo_startup) 
  !        requires a netcdf target variable.

  !     INTENT(IN)
  !         None

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

  !     RESTRICTIONS
  !         None

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         author Stephan Thober
  !         date Jul 2013
  subroutine config_variables_set( iBasin )
    !
    use mo_kind,             only: i4
    use mo_mhm_constants,    only: nodata_i4, nodata_dp
    USE mo_init_states,      ONLY: get_basin_info
    use mo_set_netcdf_restart, only: set_config
    use mo_global_variables, only: &
         L0_cellCoor    ,          & 
         L0_Id         ,           & ! Ids of grid at level-0 
         L0_areaCell   ,           & ! Ids of grid at level-0
         L0_slope_emp  ,           & ! Empirical quantiles of slope
         L1_Id         ,           & ! Ids of grid at level-1
         L1_cellCoor    ,          &
         L1_upBound_L0 ,           & ! Row start at finer level-0 scale 
         L1_downBound_L0,          & ! Row end at finer level-0 scale
         L1_leftBound_L0,          & ! Col start at finer level-0 scale
         L1_rightBound_L0,         & ! Col end at finer level-0 scale
         L1_areaCell   ,           & ! [km2] Effective area of cell at this level
         L1_nTCells_L0               ! Total number of valid L0 cells in a given L1 cell
    use mo_set_netcdf_restart,     &
         only: L0_rowCoor_out    , & ! row cell Coordinates at Level 0
         L0_colCoor_out    ,       & ! column cell Coordinates at Level 0  
         L0_Id_out         ,       & ! Ids of grid at level-0 
         L0_areaCell_out   ,       & ! Ids of grid at level-0
         L0_slope_emp_out  ,       & ! Empirical quantiles of slope
         L1_basin_Mask_out ,       & ! Mask at Level 1
         L1_Id_out         ,       & ! Ids of grid at level-1
         L1_rowCoor_out    ,       & ! row cell Coordinates at Level 1
         L1_colCoor_out    ,       & ! column cell Coordinates at Level 1  
         L1_upBound_L0_out ,       & ! Row start at finer level-0 scale 
         L1_downBound_L0_out,      & ! Row end at finer level-0 scale
         L1_leftBound_L0_out,      & ! Col start at finer level-0 scale
         L1_rightBound_L0_out,     & ! Col end at finer level-0 scale
         L1_areaCell_out   ,       & ! [km2] Effective area of cell at this level
         L1_nTCells_L0_out           ! Total number of valid L0 cells in a given L1 cell
    !
    implicit none 
    !
    ! Input Variables
    integer(i4), intent(in) :: iBasin
    !
    ! Local variables
    integer(i4)                                :: nrows0    ! Number of rows at Level 0
    integer(i4)                                :: ncols0    ! Number of cols at Level 0
    integer(i4)                                :: iStart0   ! Start of Basin at Level 0
    integer(i4)                                :: iEnd0     ! End of Basin at Level 0
    logical, dimension(:,:), allocatable       :: mask0     ! Mask at Level 0
    integer(i4)                                :: nrows1    ! Number of rows at Level 1
    integer(i4)                                :: ncols1    ! Number of cols at Level 1
    integer(i4)                                :: iStart1   ! Start of Basin at Level 1
    integer(i4)                                :: iEnd1     ! End of Basin at Level 1
    logical, dimension(:,:), allocatable       :: mask1     ! Mask at Level 1
    !
    ! get level-0 information
    call get_basin_info(iBasin, 0, nrows0, ncols0, iStart=iStart0, &
         iEnd = iEnd0, mask = mask0 )

    ! get level-1 information
    call get_basin_info(iBasin, 1, nrows1, ncols1, iStart=iStart1, &
         iEnd = iEnd1, mask = mask1 )  
    !
    ! Variables set by L0_variable_init <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate( L0_rowCoor_out( nrows0, ncols0 ) )
    L0_rowCoor_out = unpack( L0_cellCoor(iStart0:iEnd0,1), mask0, nodata_i4 )
    !
    allocate( L0_colCoor_out( nrows0, ncols0 ) )
    L0_colCoor_out = unpack( L0_cellCoor(iStart0:iEnd0,2), mask0, nodata_i4 )
    !
    allocate( L0_Id_out( nrows0, ncols0 ) )
    L0_Id_out = unpack( L0_Id(iStart0:iEnd0), mask0, nodata_i4 )
    !
    allocate( L0_areaCell_out( nrows0, ncols0 ) )
    L0_areaCell_out = unpack( L0_areaCell(iStart0:iEnd0), mask0, nodata_dp )
    !
    allocate( L0_slope_emp_out( nrows0, ncols0 ) )
    L0_slope_emp_out = unpack( L0_slope_emp(iStart0:iEnd0), mask0, nodata_dp )

    ! variables set by L1_variable_init <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate( L1_basin_Mask_out( nrows1, ncols1 ) )
    L1_basin_Mask_out = merge( 1_i4, 0_i4,  mask1 )
    !
    allocate( L1_Id_out( nrows1, ncols1 ) )
    L1_Id_out = unpack( L1_Id(iStart1:iEnd1), mask1, nodata_i4 )
    !
    allocate( L1_rowCoor_out( nrows1, ncols1 ) )
    L1_rowCoor_out = unpack( L1_cellCoor(iStart1:iEnd1,1), mask1, nodata_i4 )
    !
    allocate( L1_colCoor_out( nrows1, ncols1 ) )
    L1_colCoor_out = unpack( L1_cellCoor(iStart1:iEnd1,2), mask1, nodata_i4 )
    !
    allocate( L1_upBound_L0_out( nrows1, ncols1 ) )
    L1_upBound_L0_out = unpack( L1_upBound_L0(iStart1:iEnd1), mask1, nodata_i4 )
    !
    allocate( L1_downBound_L0_out( nrows1, ncols1 ) )
    L1_downBound_L0_out = unpack( L1_downBound_L0(iStart1:iEnd1), mask1, nodata_i4 ) 
    !
    allocate( L1_leftBound_L0_out( nrows1, ncols1 ) )
    L1_leftBound_L0_out = unpack( L1_leftBound_L0(iStart1:iEnd1), mask1, nodata_i4 )
    !
    allocate( L1_rightBound_L0_out( nrows1, ncols1 ) )
    L1_rightBound_L0_out = unpack( L1_rightBound_L0(iStart1:iEnd1), mask1, nodata_i4 )
    !
    allocate( L1_areaCell_out( nrows1, ncols1 ) )
    L1_areaCell_out = unpack( L1_areaCell(iStart1:iEnd1), mask1, nodata_dp )
    !
    allocate( L1_nTCells_L0_out( nrows1, ncols1 ) )
    L1_nTCells_L0_out = unpack( L1_nTCells_L0(iStart1:iEnd1), mask1, nodata_i4) 

    ! allocate V
    call set_config

  end subroutine config_variables_set
  ! ------------------------------------------------------------------

  !      NAME
  !         L11_config_set

  !     PURPOSE
  !        sets netcdf target variables for all mHM L11 configuration variables. These target variables
  !        are defined in the module mo_set_netcdf_restart. The actual setting of the netcdf
  !        output variable V, that points to the target variables here, takes place in the 
  !        subroutine set_L11_config contained in mo_set_netcdf_restart. In general each variable
  !        defined in the subroutines contained in the module mo_net_startup 
  !        requires a netcdf target variable.

  !     INTENT(IN)
  !         None

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

  !     RESTRICTIONS
  !         None

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         author Stephan Thober
  !         date Jul 2013
  subroutine L11_config_set( iBasin )

    use mo_kind,               only: i4
    use mo_mhm_constants,      only: nodata_i4, nodata_dp
    USE mo_init_states,        ONLY: get_basin_info
    use mo_set_netcdf_restart, only: set_L11_config
    use mo_global_variables,   only: basin, &
         L11_cellCoor, &
         L11_Id, &
         L0_L11_Id, &
         L1_L11_Id, &
         L11_rowOut, &
         L11_colOut, &
         L11_fDir, &
         L11_upBound_L0, &
         L11_downBound_L0, &
         L11_leftBound_L0, &
         L11_rightBound_L0, &
         L11_upBound_L1, &
         L11_downBound_L1, &
         L11_leftBound_L1, &
         L11_rightBound_L1, &
         L11_fDir, &
         L11_fromN, &
         L11_toN, &
         L11_rOrder, &
         L11_label, &
         L11_sink, &
         L11_netPerm, &
         L11_rowOut, &
         L11_colOut, &
         L11_fRow, &
         L11_fCol, &
         L11_tRow, &
         L11_tCol, &
         L0_draSC, &
         L0_draCell, &
         L0_streamNet, &
         L0_floodPlain, &
         L11_length, &
         L11_aFloodPlain, &
         L11_slope
    use mo_set_netcdf_restart, only: L11_basin_Mask_out, &
         L11_rowCoor_out, &
         L11_colCoor_out, &
         L11_Id_out, &
         L0_OutletCoord_out, &
         L0_draSC_out, &
         L0_L11_Id_out, &
         L1_L11_Id_out, &
         L11_rowOut_out, &
         L11_colOut_out, &
         L11_fDir_out, &
         L11_upBound_L0_out, &
         L11_downBound_L0_out, &
         L11_leftBound_L0_out, &
         L11_rightBound_L0_out, &
         L11_upBound_L1_out, &
         L11_downBound_L1_out, &
         L11_leftBound_L1_out, &
         L11_rightBound_L1_out, &
         L11_fromN_out, &
         L11_toN_out, &
         L11_rOrder_out, &
         L11_label_out, &
         L11_sink_out, &
         L11_netPerm_out, &
         L11_fRow_out, &
         L11_fCol_out, &
         L11_tRow_out, &
         L11_tCol_out, &
         L0_draCell_out, &
         gaugeNodeList_out, &
         L0_streamNet_out, &
         L0_floodPlain_out, &
         L11_length_out, &
         L11_aFloodPlain_out, &
         L11_slope_out

    implicit none

    ! Input variables
    integer(i4), intent(in)                    :: iBasin    ! Basin number

    ! Local variables
    integer(i4)                                :: nrows0    ! Number of rows at Level 0
    integer(i4)                                :: ncols0    ! Number of cols at Level 0
    integer(i4)                                :: iStart0   ! Start of Basin at Level 0
    integer(i4)                                :: iEnd0     ! End of Basin at Level 0
    logical, dimension(:,:), allocatable       :: mask0     ! Mask at Level 0
    integer(i4)                                :: nrows1    ! Number of rows at Level 1
    integer(i4)                                :: ncols1    ! Number of cols at Level 1
    integer(i4)                                :: iStart1   ! Start of Basin at Level 1
    integer(i4)                                :: iEnd1     ! End of Basin at Level 1
    logical, dimension(:,:), allocatable       :: mask1     ! Mask at Level 1
    integer(i4)                                :: nCells11  ! Number of cells at Level 11
    integer(i4)                                :: nrows11   ! Number of rows at Level 11
    integer(i4)                                :: ncols11   ! Number of cols at Level 11
    integer(i4)                                :: iStart11  ! Start of Basin at Level 11
    integer(i4)                                :: iEnd11    ! End of Basin at level 11
    logical, dimension(:,:), allocatable       :: mask11    ! Mask at Level 11
    integer(i4)                                :: startMask ! starting position of mask at L11
    integer(i4)                                :: endMask   ! ending position of mask at L11

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_variable_init(iBasin) sets
    ! level11,resolutionRouting, L11_cellCoor, L11_nCells, L11_Id

    ! get level-0 information
    call get_basin_info(iBasin, 0, nrows0, ncols0, iStart=iStart0, &
         iEnd = iEnd0, mask = mask0 )

    ! get level-11 information
    call get_basin_info(iBasin, 1, nrows1, ncols1, iStart=iStart1, &
         iEnd = iEnd1, mask = mask1 )

    ! get level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11, ncells=nCells11, &
         iStart=iStart11, iEnd=iEnd11, mask=mask11)
    startMask = basin%L11_iStartMask(iBasin)
    endMask   = basin%L11_iEndMask(iBasin)
    nCells11  = count( basin%L11_Mask(startMask:endMask) )

    ! Level 11 mask
    allocate( L11_basin_Mask_out( nrows11, ncols11))
    L11_basin_Mask_out = merge( 1_i4, 0_i4,  &
         reshape(basin%L11_Mask(startMask:endMask),(/nrows11,ncols11/)) )

    ! Level 11 cell coordinates
    allocate(L11_rowCoor_out(nrows11,ncols11))
    L11_rowCoor_out = unpack( L11_cellCoor(iStart11:iEnd11,1), mask11, nodata_i4 )
    allocate(L11_colCoor_out(nrows11,ncols11))
    L11_colCoor_out = unpack( L11_cellCoor(iStart11:iEnd11,2), mask11, nodata_i4 )


    ! Level 11 id
    allocate(L11_Id_out(nrows11,ncols11))
    L11_Id_out = unpack( L11_Id(iStart11:iEnd11), mask11, nodata_i4 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_flow_direction(iBasin) sets
    ! basin%L0_rowOutlet, basin%L0_colOutlet
    ! L0_drasC, L0_L11_Id, L1_L11_Id, L11_Id, L11_rowOut, 
    ! L11_colOut, L11_fDir, L11_upBound_L0, L11_downBound_L0, L11_leftBound_L0,
    ! L11_rightBound_L0, L11_upBound_L1, L11_downBound_L1, L11_leftBound_L1, 
    ! L11_rightBound_L1

    ! L0 data sets <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! L0 Outlet coordinates
    L0_OutletCoord_out = (/ basin%L0_rowOutlet(iBasin), basin%L0_colOutlet(iBasin) /)

    ! L0 draining cell index
    allocate( L0_draSC_out( nrows0,ncols0 ) )
    L0_draSC_out = unpack( L0_draSC(iStart0:iEnd0), mask0, nodata_i4 )

    ! Mapping of L11 Id on L0  
    allocate( L0_L11_Id_out( nrows0,ncols0 ) )
    L0_L11_Id_out = unpack( L0_L11_Id(iStart0:iEnd0), mask0, nodata_i4)

    ! L1 data sets <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Mapping of L11 Id on L1
    allocate( L1_L11_Id_out( nrows1, ncols1) )
    L1_L11_Id_out = unpack( L1_L11_Id(iStart1:iEnd1), mask1, nodata_i4)

    ! L11 data sets <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Flow direction (standard notation)
    allocate(L11_fDir_out(nrows11, ncols11))
    L11_fDir_out = unpack( L11_fDir(iStart11:iEnd11), mask11, nodata_i4)

    ! Grid vertical location of the Outlet
    allocate(L11_rowOut_out(nrows11, ncols11))
    L11_rowOut_out = unpack( L11_rowOut(iStart11:iEnd11), mask11, nodata_i4 )

    ! Grid horizontal location  of the Outlet
    allocate(L11_colOut_out(nrows11, ncols11))
    L11_colOut_out = unpack( L11_colOut(iStart11:iEnd11), mask11, nodata_i4 )

    ! Row start at finer level-0 scale 
    allocate(L11_upBound_L0_out(nrows11, ncols11))
    L11_upBound_L0_out = unpack( L11_upBound_L0(iStart11:iEnd11), mask11, nodata_i4 )

    ! Row end at finer level-0 scale
    allocate(L11_downBound_L0_out(nrows11, ncols11))
    L11_downBound_L0_out = unpack( L11_downBound_L0(iStart11:iEnd11), mask11, nodata_i4 )

    ! Col start at finer level-0 scale 
    allocate(L11_leftBound_L0_out(nrows11, ncols11))
    L11_leftBound_L0_out = unpack( L11_leftBound_L0(iStart11:iEnd11), mask11, nodata_i4 )

    ! Col end at finer level-0 scale 
    allocate(L11_rightBound_L0_out(nrows11, ncols11))
    L11_rightBound_L0_out = unpack( L11_rightBound_L0(iStart11:iEnd11), mask11, nodata_i4 )

    ! Row start at finer level-1 scale 
    allocate(L11_upBound_L1_out(nrows1, ncols1))
    L11_upBound_L1_out = unpack( L11_upBound_L1(iStart1:iEnd1), mask1, nodata_i4 )

    ! Row end at finer level-1 scale
    allocate(L11_downBound_L1_out(nrows1, ncols1))
    L11_downBound_L1_out = unpack( L11_downBound_L1(iStart1:iEnd1), mask1, nodata_i4 )

    ! Col start at finer level-1 scale 
    allocate(L11_leftBound_L1_out(nrows1, ncols1))
    L11_leftBound_L1_out = unpack( L11_leftBound_L1(iStart1:iEnd1), mask1, nodata_i4 )

    ! Col end at finer level-1 scale 
    allocate(L11_rightBound_L1_out(nrows1, ncols1))
    L11_rightBound_L1_out = unpack( L11_rightBound_L1(iStart1:iEnd1), mask1, nodata_i4 )


    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_set_network_topology(iBasin) sets
    ! L11_fromN, L11_toN

    ! From node
    allocate( L11_fromN_out( nrows11, ncols11) )
    L11_fromN_out = unpack( L11_fromN( iStart11: iEnd11), mask11, nodata_i4 )

    ! To node
    allocate( L11_toN_out( nrows11, ncols11 ) )
    L11_toN_out = unpack( L11_toN( iStart11: iEnd11 ), mask11, nodata_i4 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_routing_order(iBasin) sets
    ! L11_rOrder, L11_label, L11_sink, L11_netPerm

    ! Network routing order
    allocate( L11_rOrder_out( nrows11, ncols11 ) )
    L11_rOrder_out = unpack( L11_rOrder( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! Label Id [0='', 1=HeadWater, 2=Sink]
    allocate( L11_label_out( nrows11, ncols11 ) )
    L11_label_out = unpack( L11_label( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! .true. if sink node reached
    allocate( L11_sink_out( nrows11, ncols11 ) )
    L11_sink_out = unpack(merge( 1_i4, 0_i4, L11_sink( iStart11 : iEnd11 ) ), mask11, nodata_i4)

    ! Routing sequence (permutation of L11_rOrder)
    allocate( L11_netPerm_out( nrows11, ncols11 ) )
    L11_netPerm_out = unpack( L11_netPerm( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_link_location(iBasin) sets
    ! L11_fRow, L11_fCol, L11_tRow, L11_tCol

    ! From row in L0 grid
    allocate( L11_fRow_out( nrows11, ncols11 ) )
    L11_fRow_out = unpack( L11_fRow( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! From col in L0 grid
    allocate( L11_fCol_out( nrows11, ncols11 ) )
    L11_fCol_out = unpack( L11_fCol( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! To row in L0 grid 
    allocate( L11_tRow_out( nrows11, ncols11 ) )
    L11_tRow_out = unpack( L11_tRow( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! To col in L0 grid 
    allocate( L11_tCol_out( nrows11, ncols11 ) )
    L11_tCol_out = unpack( L11_tCol( iStart11 : iEnd11 ), mask11, nodata_i4 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! call L11_set_drain_outlet_gauges(iBasin) sets
    ! L0_draCell

    ! Draining cell id at L11 of ith cell of L0
    allocate( L0_draCell_out( nrows0, ncols0 ) )
    L0_draCell_out = unpack( L0_draCell( iStart0 : iEnd0), mask0, nodata_i4 )

    allocate( gaugeNodeList_out( size( basin%gaugeNodeList(iBasin,:), dim=1)))
    gaugeNodeList_out = basin%gaugeNodeList(iBasin,:)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! stream characteristics
    ! call L11_stream_features(iBasin)
    ! L0_streamNet, L0_floodPlain, L11_length, L11_aFloodPlain, L11_slope

    ! L0 data sets
    ! Stream network
    allocate( L0_streamNet_out( nrows0, ncols0 ) )
    L0_streamNet_out = unpack( L0_streamNet( iStart0 : iEnd0 ), mask0, nodata_i4 )

    ! Floodplains of stream i
    allocate( L0_floodPlain_out( nrows0, ncols0 ) )
    L0_floodPlain_out = unpack( L0_floodPlain( iStart0 : iEnd0 ), mask0, nodata_i4 )

    ! L11 network data sets
    ! [m]     Total length of river link
    allocate( L11_length_out( nrows11, ncols11 ) )
    L11_length_out = unpack( L11_length( iStart11 : iEnd11 ), mask11, nodata_dp )

    ! [m2]    Area of the flood plain
    allocate( L11_aFloodPlain_out( nrows11, ncols11 ) )
    L11_aFloodPlain_out = unpack( L11_aFloodPlain( iStart11 : iEnd11 ), mask11, nodata_dp )

    ! Average slope of river link
    allocate( L11_slope_out( nrows11, ncols11 ) )
    L11_slope_out = unpack( L11_slope( iStart11 : iEnd11 ), mask11, nodata_dp )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! initialize netcdf output structure V
    call set_L11_config

  end subroutine L11_config_set
  ! ---------------------------------------------------------------------------

  !     NAME
  !        state_variables_set

  !     PURPOSE

  !        sets netcdf target variables for all mHM state variables. These target variables
  !        are defined in the module mo_set_netcdf_restart. The actual setting of the netcdf
  !        output variable V, that points to the target variables here, takes place in the 
  !        subroutine set_state contained in mo_set_netcdf_restart. In general each variable
  !        defined in the subroutine variables_default_init (mo_init_states) requires a netcdf 
  !        target variable.
  !        

  !     INTENT(IN)
  !        integer(i4)    :: iBasin"   id of the Basin to write state variables for

  !     INTENT(INOUT)
  !        None

  !     INTENT(OUT)
  !        None

  !     INTENT(IN), OPTIONAL
  !        None

  !     INTENT(INOUT), OPTIONAL
  !        None

  !     INTENT(OUT), OPTIONAL
  !        None

  !     RETURN
  !        None

  !     RESTRICTIONS
  !        None

  !     HISTORY
  !        author  Stephan Thober
  !        date    Apr 2013
  ! --------------------------------------------------------------------------
  subroutine state_variables_set( iBasin )

    use mo_kind,               only: i4
    use mo_init_states,        only: get_basin_info
    use mo_mhm_constants,      only: nodata_dp
    use mo_set_netcdf_restart, only: set_state, & 
         L1_fSealed_out, & ! netcdf output variables
         L1_fForest_out, &
         L1_fPerm_out, &
         L1_Inter_out, &
         L1_snowPack_out, &
         L1_sealSTW_out, &
         L1_soilMoist_out, &
         L1_unsatSTW_out, &
         L1_satSTW_out, &
         L1_aETSoil_out, &
         L1_aETCanopy_out, &
         L1_aETSealed_out, &
         L1_baseflow_out, &
         L1_infilSoil_out, &
         L1_fastRunoff_out, &
         L1_melt_out, &
         L1_percol_out, &
         L1_preEffect_out, &
         L1_rain_out, &
         L1_runoffSeal_out, &
         L1_slowRunoff_out, &
         L1_snow_out, &
         L1_Throughfall_out, &
         L1_total_runoff_out, &
         L1_alpha_out, &
         L1_degDayInc_out, &
         L1_degDayMax_out, &
         L1_degDayNoPre_out, &
         L1_degDay_out, &
         L1_karstLoss_out, &
         L1_fAsp_out, &
         L1_fRoots_out, &
         L1_maxInter_out, &
         L1_kfastFlow_out, &
         L1_kSlowFlow_out, &
         L1_kBaseFlow_out, &
         L1_kPerco_out, &
         L1_soilMoistFC_out, &
         L1_soilMoistSat_out, &
         L1_soilMoistExp_out, &
         L1_tempThresh_out, &
         L1_unsatThresh_out, &
         L1_sealedThresh_out, &
         L1_wiltingPoint_out, &
         L11_Qmod_out, &
         L11_qOUT_out, &
         L11_qTIN_out, &
         L11_qTR_out, &
         L11_K_out, &
         L11_xi_out, &
         L11_C1_out, &
         L11_C2_out, &
         L11_FracFPimp_out
    use mo_global_variables, only: processMatrix, & ! mhm_state_variables
         L1_fSealed, &
         L1_fForest, &
         L1_fPerm, &
         L1_Inter, &
         L1_snowPack, &
         L1_sealSTW, &
         L1_soilMoist, &
         L1_unsatSTW, &
         L1_satSTW, &
         L1_aETSoil, &
         L1_aETCanopy, &
         L1_aETSealed, &
         L1_baseflow, &
         L1_infilSoil, &
         L1_fastRunoff, &
         L1_melt, &
         L1_percol, &
         L1_preEffect, &
         L1_rain, &
         L1_runoffSeal, &
         L1_slowRunoff, &
         L1_snow, &
         L1_Throughfall, &
         L1_total_runoff, &
         L1_alpha, &
         L1_degDayInc, &
         L1_degDayMax, &
         L1_degDayNoPre, &
         L1_degDay, &
         L1_karstLoss, &
         L1_fAsp, &
         L1_fRoots, &
         L1_maxInter, &
         L1_kfastFlow, &
         L1_kSlowFlow, &
         L1_kBaseFlow, &
         L1_kPerco, &
         L1_soilMoistFC, &
         L1_soilMoistSat, &
         L1_soilMoistExp, &
         L1_tempThresh, &
         L1_unsatThresh, &
         L1_sealedThresh, &
         L1_wiltingPoint, &
         L11_Qmod, &
         L11_qOUT, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp

    implicit none

    ! Input and Output variables
    integer(i4),    intent(in) :: iBasin

    ! local variables
    integer(i4)                          :: i
    integer(i4)                          :: s1       ! start index at level 1
    integer(i4)                          :: e1       ! end index at level 1
    integer(i4)                          :: ncols1   ! number of colums at level 1
    integer(i4)                          :: nrows1   ! number of rows at level 1
    integer(i4)                          :: ncells1  ! number of cells at level 1
    logical, dimension(:,:), allocatable :: mask1    ! mask at level 1
    integer(i4)                          :: s11      ! start index at level 11
    integer(i4)                          :: e11      ! end index at level 11
    integer(i4)                          :: ncols11  ! number of colums at level 11
    integer(i4)                          :: nrows11  ! number of rows at level 11
    integer(i4)                          :: ncells11 ! number of cells at level 11
    logical, dimension(:,:), allocatable :: mask11   ! mask at level 11

    ! get Level1 information about the basin
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 ) 

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------
    allocate( L1_fSealed_out( nrows1, ncols1 ) )
    allocate( L1_fForest_out( nrows1, ncols1 ) )
    allocate( L1_fPerm_out(   nrows1, ncols1 ) )

    L1_fSealed_out = unpack( L1_fSealed(s1:e1), mask1, nodata_dp )
    L1_fForest_out = unpack( L1_fForest(s1:e1), mask1, nodata_dp )
    L1_fPerm_out   = unpack( L1_fPerm(s1:e1), mask1, nodata_dp )

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    allocate( L1_inter_out( nrows1, ncols1 ) )
    L1_inter_out = unpack( L1_inter(s1:e1), mask1, nodata_dp )

    !Snowpack
    allocate( L1_snowPack_out( nrows1, ncols1 ) )
    L1_snowPack_out = unpack( L1_snowPack(s1:e1), mask1, nodata_dp )

    !Retention storage of impervious areas
    allocate( L1_sealSTW_out( nrows1, ncols1 ) )
    L1_sealSTW_out = unpack( L1_sealSTW(s1:e1), mask1, nodata_dp )

    ! Soil moisture of each horizon    
    allocate( L1_soilMoist_out( nrows1, ncols1, size(L1_soilMoist,2) ) )
    do i = 1, size(L1_soilMoist_out, 3)
       L1_soilMoist_out(:,:,i) = unpack( L1_soilMoist(s1:e1,i), mask1, nodata_dp )
    end do

    ! upper soil storage
    allocate( L1_unsatSTW_out( nrows1, ncols1 ) )
    L1_unsatSTW_out = unpack( L1_unsatSTW(s1:e1), mask1, nodata_dp )

    ! groundwater storage
    allocate( L1_satSTW_out( nrows1, ncols1 ) )
    L1_satSTW_out = unpack( L1_satSTW(s1:e1), mask1, nodata_dp )

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------

    !  soil actual ET
    allocate( L1_aETSoil_out( nrows1, ncols1, size(L1_aETSoil, 2) ) )
    do i = 1, size( L1_aETSoil_out, 3 )
       L1_aETSoil_out(:,:,i) = unpack( L1_aETSoil(s1:e1,i), mask1, nodata_dp )
    end do

    ! canopy actual ET
    allocate( L1_aETCanopy_out( nrows1, ncols1 ) )
    L1_aETCanopy_out = unpack( L1_aETCanopy(s1:e1), mask1, nodata_dp )

    ! sealed area actual ET
    allocate( L1_aETSealed_out( nrows1, ncols1 ) )
    L1_aETSealed_out = unpack( L1_aETSealed(s1:e1), mask1, nodata_dp )

    ! baseflow
    allocate( L1_baseflow_out( nrows1, ncols1 ) )
    L1_baseflow_out = unpack( L1_baseflow(s1:e1), mask1, nodata_dp )

    !  soil in-exfiltration
    allocate( L1_infilSoil_out( nrows1, ncols1, size(L1_infilSoil,2) ) )
    do i = 1, size( L1_infilSoil_out, 3 )
       L1_infilSoil_out(:,:,i) = unpack( L1_infilSoil(s1:e1,i), mask1, nodata_dp)
    end do

    ! fast runoff
    allocate( L1_fastRunoff_out( nrows1, ncols1 ) )
    L1_fastRunoff_out = unpack( L1_fastRunoff(s1:e1), mask1, nodata_dp )

    ! snow melt
    allocate( L1_melt_out( nrows1, ncols1 ) )
    L1_melt_out = unpack( L1_melt(s1:e1), mask1, nodata_dp )

    ! percolation
    allocate( L1_percol_out( nrows1, ncols1 ) )
    L1_percol_out = unpack( L1_percol(s1:e1), mask1, nodata_dp )

    ! effective precip. depth (snow melt + rain)
    allocate( L1_preEffect_out( nrows1, ncols1 ) )
    L1_preEffect_out = unpack( L1_preEffect(s1:e1), mask1, nodata_dp )

    ! rain (liquid water)
    allocate( L1_rain_out( nrows1, ncols1 ) )
    L1_rain_out = unpack( L1_rain(s1:e1), mask1, nodata_dp )

    ! runoff from impervious area
    allocate( L1_runoffSeal_out( nrows1, ncols1 ) )
    L1_runoffSeal_out = unpack( L1_runoffSeal(s1:e1), mask1, nodata_dp )

    ! slow runoff 
    allocate( L1_slowRunoff_out( nrows1, ncols1 ) )
    L1_slowRunoff_out = unpack( L1_slowRunoff(s1:e1), mask1, nodata_dp )

    ! snow (solid water) 
    allocate( L1_snow_out( nrows1, ncols1 ) )
    L1_snow_out = unpack( L1_snow(s1:e1), mask1, nodata_dp )

    ! throughfall 
    allocate( L1_Throughfall_out( nrows1, ncols1 ) )
    L1_Throughfall_out = unpack( L1_Throughfall(s1:e1), mask1, nodata_dp )

    ! total runoff
    allocate( L1_total_runoff_out( nrows1, ncols1 ) )
    L1_total_runoff_out = unpack( L1_total_runoff(s1:e1), mask1, nodata_dp )

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    allocate( L1_alpha_out( nrows1, ncols1 ) )
    L1_alpha_out = unpack( L1_alpha(s1:e1), mask1, nodata_dp )

    ! increase of the Degree-day factor per mm of increase in precipitation
    allocate( L1_degDayInc_out( nrows1, ncols1 ) )
    L1_degDayInc_out = unpack( L1_degDayInc(s1:e1), mask1, nodata_dp )

    ! maximum degree-day factor 
    allocate( L1_degDayMax_out( nrows1, ncols1 ) )
    L1_degDayMax_out = unpack( L1_degDayMax(s1:e1), mask1, nodata_dp )

    ! degree-day factor with no precipitation
    allocate( L1_degDayNoPre_out( nrows1, ncols1 ) )
    L1_degDayNoPre_out = unpack( L1_degDayNoPre(s1:e1), mask1, nodata_dp )

    ! degree-day factor
    allocate( L1_degDay_out( nrows1, ncols1 ) )
    L1_degDay_out = unpack( L1_degDay(s1:e1), mask1, nodata_dp )

    ! Karstic percolation loss
    allocate( L1_karstLoss_out( nrows1, ncols1 ) )
    L1_karstLoss_out = unpack( L1_karstLoss(s1:e1), mask1, nodata_dp )

    ! PET correction factor due to terrain aspect
    allocate( L1_fAsp_out( nrows1, ncols1 ) )
    L1_fAsp_out = unpack( L1_fAsp(s1:e1), mask1, nodata_dp )

    ! Fraction of roots in soil horizons
    allocate( L1_fRoots_out( nrows1, ncols1, size(L1_fRoots, 2 ) ) )
    do i = 1, size( L1_fRoots_out, 3)
       L1_fRoots_out(:,:,i) = unpack( L1_fRoots(s1:e1,i), mask1, nodata_dp )
    end do

    ! Maximum interception 
    allocate( L1_maxInter_out( nrows1, ncols1 ) )
    L1_maxInter_out = unpack( L1_maxInter(s1:e1), mask1, nodata_dp )

    ! fast interflow recession coefficient 
    allocate( L1_kfastFlow_out( nrows1, ncols1 ) )
    L1_kfastFlow_out = unpack( L1_kfastFlow(s1:e1), mask1, nodata_dp )

    ! slow interflow recession coefficient 
    allocate( L1_kSlowFlow_out( nrows1, ncols1 ) )
    L1_kSlowFlow_out = unpack( L1_kSlowFlow(s1:e1), mask1, nodata_dp )

    ! baseflow recession coefficient 
    allocate( L1_kBaseFlow_out( nrows1, ncols1 ) )
    L1_kBaseFlow_out = unpack( L1_kBaseFlow(s1:e1), mask1, nodata_dp )

    ! percolation coefficient 
    allocate( L1_kPerco_out( nrows1, ncols1 ) )
    L1_kPerco_out = unpack( L1_kPerco(s1:e1), mask1, nodata_dp )

    ! Soil moisture below which actual ET is reduced linearly till PWP
    allocate( L1_soilMoistFC_out( nrows1, ncols1, size( L1_soilMoistFC, 2 ) ) )
    do i = 1, size( L1_soilMoistFC_out, 3 )
       L1_soilMoistFC_out(:,:,i) = unpack( L1_soilMoistFC(s1:e1, i), mask1, nodata_dp )
    end do

    ! Saturation soil moisture for each horizon [mm]
    allocate( L1_soilMoistSat_out( nrows1, ncols1, size( L1_soilMoistSat, 2 ) ) )
    do i = 1, size( L1_soilMoistSat_out, 3 )
       L1_soilMoistSat_out(:,:,i) = unpack( L1_soilMoistSat(s1:e1, i), mask1, nodata_dp )
    end do

    ! Exponential parameter to how non-linear is the soil water retention
    allocate( L1_soilMoistExp_out( nrows1, ncols1, size( L1_soilMoistExp, 2 ) ) )
    do i = 1, size( L1_soilMoistExp_out, 3 )
       L1_soilMoistExp_out(:,:,i) = unpack( L1_soilMoistExp(s1:e1, i), mask1, nodata_dp )
    end do

    ! Threshold temperature for snow/rain 
    allocate( L1_tempThresh_out( nrows1, ncols1 ) )
    L1_tempThresh_out = unpack( L1_tempThresh(s1:e1), mask1, nodata_dp )

    ! Threshhold water depth controlling fast interflow
    allocate( L1_unsatThresh_out( nrows1, ncols1 ) )
    L1_unsatThresh_out = unpack( L1_unsatThresh(s1:e1), mask1, nodata_dp )

    ! Threshhold water depth for surface runoff in sealed surfaces
    allocate( L1_sealedThresh_out( nrows1, ncols1 ) )
    L1_sealedThresh_out = unpack( L1_sealedThresh(s1:e1), mask1, nodata_dp )

    ! Permanent wilting point
    allocate( L1_wiltingPoint_out( nrows1, ncols1, size( L1_wiltingPoint, 2 ) ) )
    do i = 1, size( L1_wiltingPoint_out, 3 )
       L1_wiltingPoint_out(:,:,i) = unpack( L1_wiltingPoint(s1:e1, i), mask1, nodata_dp )
    end do

    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------
    if (  processMatrix(8, 1) /= 0 ) then

       ! level-11 information
       call get_basin_info( iBasin, 11, nrows11, ncols11, ncells=ncells11, &
            iStart=s11, iEnd=e11, mask=mask11 ) 

       ! simulated discharge at each node
       allocate( L11_Qmod_out( nrows11, ncols11 ) )
       L11_Qmod_out = unpack( L11_Qmod(s11:e11), mask11, nodata_dp )

       ! Total outflow from cells L11 at time tt
       allocate( L11_qOUT_out( nrows11, ncols11 ) )
       L11_qOUT_out = unpack( L11_qOUT(s11:e11), mask11, nodata_dp )

       ! Total discharge inputs at t-1 and t
       allocate( L11_qTIN_out( nrows11, ncols11, size( L11_qTIN, 2 ) ) )
       do i = 1, size( L11_qTIN_out, 3 )
          L11_qTIN_out(:,:,i) = unpack( L11_qTIN(s11:e11, i), mask11, nodata_dp )
       end do

       !  Routed outflow leaving a node
       allocate( L11_qTR_out( nrows11, ncols11, size( L11_qTR, 2 ) ) )
       do i = 1, size( L11_qTR_out, 3 )
          L11_qTR_out(:,:,i) = unpack( L11_qTR(s11:e11, i), mask11, nodata_dp )
       end do

       ! kappa: Muskingum travel time parameter.
       allocate( L11_K_out( nrows11, ncols11 ) )
       L11_K_out = unpack( L11_K(s11:e11), mask11, nodata_dp )

       ! xi:    Muskingum diffusion parameter
       allocate( L11_xi_out( nrows11, ncols11 ) )
       L11_xi_out = unpack( L11_xi(s11:e11), mask11, nodata_dp )

       ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
       allocate( L11_C1_out( nrows11, ncols11 ) )
       L11_C1_out = unpack( L11_C1(s11:e11), mask11, nodata_dp )

       ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
       allocate( L11_C2_out( nrows11, ncols11 ) )
       L11_C2_out = unpack( L11_C2(s11:e11), mask11, nodata_dp )

       ! Fraction of the flood plain with impervious cover
       allocate( L11_FracFPimp_out( nrows11, ncols11 ) )
       L11_FracFPimp_out = unpack( L11_FracFPimp(s11:e11), mask11, nodata_dp )

    end if

    ! initialize netcdf output structure V
    call set_state( (processMatrix( 8, 1) /= 0_i4 ) )

  end subroutine state_variables_set
  !
  ! ----------------------------------------------------------------------------
  !
  ! subroutine free memory config
  !
  ! frees all variables that are required by create_netcdf and that were set by
  ! config_set
  !
  ! author: Stephan Thober
  ! created: 11.07.2013
  !
  ! ----------------------------------------------------------------------------
  subroutine free_memory_config
    !
    use mo_set_netcdf_restart , only: nrows0, & ! Number of rows at Level 0
         ncols0               , & ! Number of colss at Level 0
         nrows1               , & ! Number of rows at Level 1
         ncols1               , & ! Number of cols at Level 1
         L0_rowCoor_out       , & ! row cell Coordinates at Level 0
         L0_colCoor_out       , & ! column cell Coordinates at Level 0  
         L0_Id_out            , & ! Ids of grid at level-0 
         L0_areaCell_out      , & ! Ids of grid at level-0
         L0_slope_emp_out     , & ! Empirical quantiles of slope
         L1_basin_Mask_out    , & ! Mask at Level 1
         L1_Id_out            , & ! Ids of grid at level-1
         L1_rowCoor_out       , & ! row cell Coordinates at Level 1
         L1_colCoor_out       , & ! column cell Coordinates at Level 1  
         L1_upBound_L0_out    , & ! Row start at finer level-0 scale 
         L1_downBound_L0_out  , & ! Row end at finer level-0 scale
         L1_leftBound_L0_out  , & ! Col start at finer level-0 scale
         L1_rightBound_L0_out , & ! Col end at finer level-0 scale
         L1_areaCell_out      , & ! [km2] Effective area of cell at this level
         L1_nTCells_L0_out        ! Total number of valid L0 cells in a given L1 cell
    !
    implicit none
    !
    deallocate(  nrows0,       & ! Number of rows at Level 0
         ncols0             ,  & ! Number of colss at Level 0
         nrows1             ,  & ! Number of rows at Level 1
         ncols1             ,  & ! Number of cols at Level 1
         L0_rowCoor_out    ,   & ! row cell Coordinates at Level 0
         L0_colCoor_out    ,   & ! column cell Coordinates at Level 0  
         L0_Id_out         ,   & ! Ids of grid at level-0 
         L0_areaCell_out   ,   & ! Ids of grid at level-0
         L0_slope_emp_out  ,   & ! Empirical quantiles of slope
         L1_basin_Mask_out ,   & ! Mask at Level 1
         L1_Id_out         ,   & ! Ids of grid at level-1
         L1_rowCoor_out    ,   & ! row cell Coordinates at Level 1
         L1_colCoor_out    ,   & ! column cell Coordinates at Level 1  
         L1_upBound_L0_out ,   & ! Row start at finer level-0 scale 
         L1_downBound_L0_out,  & ! Row end at finer level-0 scale
         L1_leftBound_L0_out,  & ! Col start at finer level-0 scale
         L1_rightBound_L0_out, & ! Col end at finer level-0 scale
         L1_areaCell_out   ,   & ! [km2] Effective area of cell at this level
         L1_nTCells_L0_out )     ! Total number of valid L0 cells in a given L1 cell
    !
  end subroutine free_memory_config
  !
  ! ----------------------------------------------------------------------------
  !
  ! subroutine free memory L11 config
  !
  ! frees all variables that are required by create_netcdf and that were set by
  ! L11_config_set
  !
  ! author: Stephan Thober
  ! created: 10.07.2013
  !
  ! ----------------------------------------------------------------------------
  subroutine free_memory_L11_config
    use mo_set_netcdf_restart, only: nrows0, & ! Number of rows at Level 0
         ncols0            ,                 & ! Number of colss at Level 0
         nrows1            ,                 & ! Number of rows at Level 1
         ncols1            ,                 & ! Number of cols at Level 1
         nrows11           ,                 & ! Number of rows at Level 11
         ncols11           ,                 & ! Number of cols at Level 11
         NoutletCoord      ,                 & ! Dimension of outlet coordiantes at Level 0
         Ngauges, &
         L11_basin_Mask_out,                 & ! Mask at Level 11
         L11_rowCoor_out   ,                 & ! row cell Coordinates at Level 11
         L11_colCoor_out   ,                 & ! column cell Coordinates at Level 11
         L11_Id_out        ,                 & ! Ids of grid at level-11 
         L0_draSC_out      ,                 & ! Index of draining cell of each sub catchment 
         L0_L11_Id_out     ,                 & ! Mapping of L11 Id on L0 
         L1_L11_Id_out     ,                 & ! Mapping of L11 Id on L1
         L11_fDir_out      ,                 & ! Flow direction (standard notation)
         L11_rowOut_out    ,                 & ! Grid vertical location of the Outlet
         L11_colOut_out    ,                 & ! Grid horizontal location  of the Outlet
         L11_upBound_L0_out,                 & ! Row start at finer level-0 scale 
         L11_downBound_L0_out,               & ! Row end at finer level-0 scale
         L11_leftBound_L0_out,               & ! Col start at finer level-0 scale
         L11_rightBound_L0_out,              & ! Col end at finer level-0 scale
         L11_upBound_L1_out,                 & ! Row start at finer level-1 scale 
         L11_downBound_L1_out,               & ! Row end at finer level-1 scale
         L11_leftBound_L1_out,               & ! Col start at finer level-1 scale
         L11_rightBound_L1_out,              & ! Col end at finer level-1 scale 
         L11_fromN_out     ,                 & ! From node
         L11_toN_out       ,                 & ! To node
         L11_rOrder_out    ,                 & ! Network routing order
         L11_label_out     ,                 & ! Label Id [0='', 1=HeadWater, 2=Sink]
         L11_sink_out      ,                 & ! .true. if sink node reached
         L11_netPerm_out   ,                 & ! Routing sequence (permutation of L11_rOrder)
         L11_fRow_out      ,                 & ! From row in L0 grid 
         L11_fCol_out      ,                 & ! From col in L0 grid
         L11_tRow_out      ,                 & ! To row in L0 grid
         L11_tCol_out      ,                 & ! To col in L0 grid 
         L0_draCell_out    ,                 & ! Draining cell id at L11 of ith cell of L0
         gaugeNodeList_out, &
         L0_streamNet_out  ,                 & ! Stream network
         L0_floodPlain_out ,                 & ! Floodplains of stream i
         L11_length_out    ,                 & ! [m]     Total length of river link
         L11_aFloodPlain_out,                & ! [m2]    Area of the flood plain
         L11_slope_out                         ! Average slope of river link

    implicit none

    ! free memory
    deallocate( nrows0, & ! Number of rows at Level 0
         ncols0            ,                 & ! Number of colss at Level 0
         nrows1            ,                 & ! Number of rows at Level 1
         ncols1            ,                 & ! Number of cols at Level 1
         nrows11           ,                 & ! Number of rows at Level 11
         ncols11           ,                 & ! Number of cols at Level 11
         NoutletCoord      ,                 & ! Dimension of outlet coordiantes at Level 0
         Ngauges, &
         L11_basin_Mask_out,                 & ! Mask at Level 11
         L11_rowCoor_out   ,                 & ! row cell Coordinates at Level 11
         L11_colCoor_out   ,                 & ! column cell Coordinates at Level 11
         L11_Id_out        ,                 & ! Ids of grid at level-11 
         L0_draSC_out      ,                 & ! Index of draining cell of each sub catchment 
         L0_L11_Id_out     ,                 & ! Mapping of L11 Id on L0 
         L1_L11_Id_out     ,                 & ! Mapping of L11 Id on L1
         L11_fDir_out      ,                 & ! Flow direction (standard notation)
         L11_rowOut_out    ,                 & ! Grid vertical location of the Outlet
         L11_colOut_out    ,                 & ! Grid horizontal location  of the Outlet
         L11_upBound_L0_out,                 & ! Row start at finer level-0 scale 
         L11_downBound_L0_out,               & ! Row end at finer level-0 scale
         L11_leftBound_L0_out,               & ! Col start at finer level-0 scale
         L11_rightBound_L0_out,              & ! Col end at finer level-0 scale
         L11_upBound_L1_out,                 & ! Row start at finer level-1 scale 
         L11_downBound_L1_out,               & ! Row end at finer level-1 scale
         L11_leftBound_L1_out,               & ! Col start at finer level-1 scale
         L11_rightBound_L1_out,              & ! Col end at finer level-1 scale 
         L11_fromN_out     ,                 & ! From node
         L11_toN_out       ,                 & ! To node
         L11_rOrder_out    ,                 & ! Network routing order
         L11_label_out     ,                 & ! Label Id [0='', 1=HeadWater, 2=Sink]
         L11_sink_out      ,                 & ! .true. if sink node reached
         L11_netPerm_out   ,                 & ! Routing sequence (permutation of L11_rOrder)
         L11_fRow_out      ,                 & ! From row in L0 grid 
         L11_fCol_out      ,                 & ! From col in L0 grid
         L11_tRow_out      ,                 & ! To row in L0 grid
         L11_tCol_out      ,                 & ! To col in L0 grid 
         L0_draCell_out    ,                 & ! Draining cell id at L11 of ith cell of L0
         gaugeNodeList_out, &
         L0_streamNet_out  ,                 & ! Stream network
         L0_floodPlain_out ,                 & ! Floodplains of stream i
         L11_length_out    ,                 & ! [m]     Total length of river link
         L11_aFloodPlain_out,                & ! [m2]    Area of the flood plain
         L11_slope_out )                       ! Average slope of river link

  end subroutine free_memory_L11_config
  !
  ! ----------------------------------------------------------------------------
  !
  ! subroutine free memory
  !
  ! frees all variables that are required by create_netcdf and that were set by 
  ! set_state
  !
  ! author: Stephan Thober
  ! created: 24.4.2013
  !
  ! ----------------------------------------------------------------------------
  subroutine free_memory_states( L11_flag )

    use mo_set_netcdf_restart, only: x1, y1, z1, x11, y11, z11, DNC, &
         V, &
         L1_fSealed_out, &
         L1_fForest_out, &
         L1_fPerm_out, &
         L1_Inter_out, &
         L1_snowPack_out, &
         L1_sealSTW_out, &
         L1_soilMoist_out, &
         L1_unsatSTW_out, &
         L1_satSTW_out, &
         L1_aETSoil_out, &
         L1_aETCanopy_out, &
         L1_aETSealed_out, &
         L1_baseflow_out, &
         L1_infilSoil_out, &
         L1_fastRunoff_out, &
         L1_melt_out, &
         L1_percol_out, &
         L1_preEffect_out, &
         L1_rain_out, &
         L1_runoffSeal_out, &
         L1_slowRunoff_out, &
         L1_snow_out, &
         L1_Throughfall_out, &
         L1_total_runoff_out, &
         L1_alpha_out, &
         L1_degDayInc_out, &
         L1_degDayMax_out, &
         L1_degDayNoPre_out, &
         L1_degDay_out, &
         L1_karstLoss_out, &
         L1_fAsp_out, &
         L1_fRoots_out, &
         L1_maxInter_out, &
         L1_kfastFlow_out, &
         L1_kSlowFlow_out, &
         L1_kBaseFlow_out, &
         L1_kPerco_out, &
         L1_soilMoistFC_out, &
         L1_soilMoistSat_out, &
         L1_soilMoistExp_out, &
         L1_tempThresh_out, &
         L1_unsatThresh_out, &
         L1_sealedThresh_out, &
         L1_wiltingPoint_out, &
         L11_Qmod_out, &
         L11_qOUT_out, &
         L11_qTIN_out, &
         L11_qTR_out, &
         L11_K_out, &
         L11_xi_out, &
         L11_C1_out, &
         L11_C2_out, &
         L11_FracFPimp_out

    implicit none

    logical, intent(in) :: L11_flag ! flag for L11_variables

    deallocate( x1, y1, z1, DNC, &
         V, &
         L1_fSealed_out, &
         L1_fForest_out, &
         L1_fPerm_out, &
         L1_Inter_out, &
         L1_snowPack_out, &
         L1_sealSTW_out, &
         L1_soilMoist_out, &
         L1_unsatSTW_out, &
         L1_satSTW_out, &
         L1_aETSoil_out, &
         L1_aETCanopy_out, &
         L1_aETSealed_out, &
         L1_baseflow_out, &
         L1_infilSoil_out, &
         L1_fastRunoff_out, &
         L1_melt_out, &
         L1_percol_out, &
         L1_preEffect_out, &
         L1_rain_out, &
         L1_runoffSeal_out, &
         L1_slowRunoff_out, &
         L1_snow_out, &
         L1_Throughfall_out, &
         L1_total_runoff_out, &
         L1_alpha_out, &
         L1_degDayInc_out, &
         L1_degDayMax_out, &
         L1_degDayNoPre_out, &
         L1_degDay_out, &
         L1_karstLoss_out, &
         L1_fAsp_out, &
         L1_fRoots_out, &
         L1_maxInter_out, &
         L1_kfastFlow_out, &
         L1_kSlowFlow_out, &
         L1_kBaseFlow_out, &
         L1_kPerco_out, &
         L1_soilMoistFC_out, &
         L1_soilMoistSat_out, &
         L1_soilMoistExp_out, &
         L1_tempThresh_out, &
         L1_unsatThresh_out, &
         L1_sealedThresh_out, &
         L1_wiltingPoint_out )

    if ( L11_flag ) then
       deallocate( x11, y11, z11, &
         L11_Qmod_out, &
         L11_qOUT_out, &
         L11_qTIN_out, &
         L11_qTR_out, &
         L11_K_out, &
         L11_xi_out, &
         L11_C1_out, &
         L11_C2_out, &
         L11_FracFPimp_out )
    end if           

  end subroutine free_memory_states

END MODULE mo_restart
