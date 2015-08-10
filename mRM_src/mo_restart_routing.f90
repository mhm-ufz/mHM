module mo_restart_routing
  use mo_kind, only: i4, dp
  implicit none
  public :: read_restart_routing
  public :: write_restart_routing
  public :: read_restart_L11_config
contains
  subroutine write_restart_routing(iBasin, OutPath)
    use mo_message, only: message
    use mo_ncwrite, only: var2nc
    use mo_string_utils, only: num2str
    use mo_global_variables_routing, only: &
         basin_mrm, &
         L11_Qmod, &
         L11_qOUT, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp, &
         L11_cellCoor, &
         L11_Id, &
         L0_draSC, &
         L0_draCell, &
         L0_streamNet, &
         L0_floodPlain, &
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
         L11_length, &
         L11_aFloodPlain, &
         L11_slope
    use mo_mrm_constants, only: nRoutingStates, nodata_dp, nodata_i4
    !ST: the following dependency has to be removed
    use mo_init_states, only: get_basin_info
    use mo_global_variables, only: basin

    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    character(256), dimension(:), intent(in) :: OutPath ! list of Output paths per Basin
    ! local variables
    character(256) :: Fname
    integer(i4) :: ii
    integer(i4) :: s0 ! start index at level 0
    integer(i4) :: e0 ! end index at level 0
    integer(i4) :: ncols0 ! number of colums at level 0
    integer(i4) :: nrows0 ! number of rows at level 0
    logical, dimension(:,:), allocatable :: mask0 ! mask at level 0
    integer(i4) :: s1 ! start index at level 1
    integer(i4) :: e1 ! end index at level 1
    integer(i4) :: ncols1 ! number of colums at level 1
    integer(i4) :: nrows1 ! number of rows at level 1
    logical, dimension(:,:), allocatable :: mask1 ! mask at level 1
    integer(i4) :: s11 ! start index at level 11
    integer(i4) :: e11 ! end index at level 11
    integer(i4) :: ncols11 ! number of colums at level 11
    integer(i4) :: nrows11 ! number of rows at level 11
    logical, dimension(:,:), allocatable :: mask11 ! mask at level 11
    integer(i4) :: s110 ! start index at pseudo level 110 
    integer(i4) :: e110 ! end index at pseudo level 110
    integer(i4) :: ncols110 ! number of colums at pseudo level 110
    integer(i4) :: nrows110 ! number of rows at pseudo level 110
    real(dp), dimension(:,:,:), allocatable  :: dummy_d3 ! dummy variable
    character(256), dimension(2) :: dims_L0 ! dimension names for L0 states
    character(256), dimension(4) :: dims_L1 ! dimension names for L1 states
    character(256), dimension(3) :: dims_L11 ! dimension names for L11 states
    character(256), dimension(1) :: dims_Links ! dimension names for L11 links
    character(256), dimension(1) :: dims_outlet ! dimension name  for outlet Coordinates
    character(256), dimension(1) :: dims_gauges ! dimension name  for number of gauges
    character(256), dimension(1) :: dims_inflow ! dimension name  for inflow gauge

    ! inititialize dimension names
    dims_L0(1)     = 'nrows0'
    dims_L0(2)     = 'ncols0'
    dims_L1(1)     = 'nrows1'
    dims_L1(2)     = 'ncols1'
    dims_L1(3)     = 'L1_soilhorizons'
    dims_L1(4)     = 'MonthsPerYear'
    dims_L11(1)    = 'nrows11'
    dims_L11(2)    = 'ncols11'
    dims_L11(3)    = 'nIT'
    dims_outlet(1) = 'NoutletCoord'
    dims_gauges(1) = 'Ngauges'
    dims_inflow(1) = 'nInflowGauges'
    dims_Links     = 'nLinks'
    
    ! get Level0 information about the basin
    call get_basin_info( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )
    ! get Level1 information about the basin
    call get_basin_info( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )
    ! get Level11 information about the basin
    call get_basin_info( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )
    ! get Level110 information about the basin
    call get_basin_info( iBasin, 110, nrows110, ncols110, iStart=s110, iEnd=e110)

    !ST set variable name for compitability with mHM, should be changed
    Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_states.nc'

    call var2nc( Fname, unpack( L11_Qmod(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_Qmod', &
         long_name = 'simulated discharge at each node at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_qOUT(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_qOUT', &
         long_name = 'Total outflow from cells L11 at time tt at level 11', missing_value = nodata_dp)

    allocate( dummy_d3( nrows11, ncols11, nRoutingStates ) )
    do ii = 1, size( dummy_d3, 3 )
       dummy_d3(:,:,ii) = unpack( L11_qTIN(s11:e11,ii), mask11, nodata_dp )
    end do
    call var2nc( Fname, dummy_d3, &
         dims_L11, 'L11_qTIN', &
         long_name = 'Total discharge inputs at t-1 and t at level 11', missing_value = nodata_dp)

    do ii = 1, size( dummy_d3, 3 )
       dummy_d3(:,:,ii) = unpack( L11_qTR(s11:e11,ii), mask11, nodata_dp )
    end do
    call var2nc( Fname, dummy_d3, &
         dims_L11, 'L11_qTR', &
         long_name = 'Routed outflow leaving a node at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_K(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_K', &
         long_name = 'kappa: Muskingum travel time parameter at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_xi(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_xi', &
         long_name = 'xi: Muskingum diffusion parameter at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_C1(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_C1', &
         long_name = 'Routing parameter C1=f(K,xi, DT) (Chow, 25-41) at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_C2(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_C2', &
         long_name = 'Routing parameter C2=f(K,xi, DT) (Chow, 25-41) at level 11', missing_value = nodata_dp)

    call var2nc( Fname, unpack( L11_FracFPimp(s11:e11), mask11, nodata_dp ), &
         dims_L11(1:2), 'L11_FracFPimp', &
         long_name = 'Fraction of the flood plain with impervious cover at level 11', missing_value = nodata_dp)

    ! ----------------------------------------------------------
    ! L11 config set - create new file
    ! ----------------------------------------------------------
    Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_L11_config.nc'
    call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
    call var2nc( Fname, &
         merge( 1_i4, 0_i4,  &
         reshape(basin%L11_Mask(basin%L11_iStartMask(iBasin):basin%L11_iEndMask(iBasin)),&
         (/nrows11,ncols11/)) ),&
         dims_L11(1:2), 'L11_basin_Mask', &
         long_name = 'Mask at Level 11', missing_value = nodata_i4, create = .true. )

    call var2nc( Fname, unpack( L11_cellCoor(s11:e11,1), mask11, nodata_i4 ), &
         dims_L11(1:2),'L11_rowCoor', &
         long_name = 'row coordinates at Level 11', missing_value = nodata_i4 )

    call var2nc( Fname, unpack( L11_cellCoor(s11:e11,2), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_colCoor', &
         long_name = 'col coordinates at Level 11', missing_value = nodata_i4 )

    call var2nc( Fname, unpack( L11_Id(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_Id', &
         long_name = 'cell Ids at Level 11', missing_value = nodata_i4 )

    call var2nc( Fname, unpack( L11_fDir(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_fDir', &
         long_name = 'flow Direction at Level 11', missing_value = nodata_i4 )     

    call var2nc( Fname, unpack( L11_rowOut(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_rowOut', &
         long_name = 'Grid vertical location of the Outlet at Level 11', missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_colOut(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_colOut', &
         long_name = 'Grid horizontal location of the Outlet at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_upBound_L0(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_upBound_L0', &
         long_name = 'Row start at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_downBound_L0(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_downBound_L0', &
         long_name = 'Row end at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_leftBound_L0(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_leftBound_L0', &
         long_name = 'Col start at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_rightBound_L0(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_rightBound_L0', &
         long_name = 'Col end at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

    call var2nc( Fname, L11_fromN(s11:e11), &
         dims_Links(1:1), 'L11_fromN', &
         long_name = 'From Node',missing_value=nodata_i4)

    call var2nc( Fname, L11_toN(s11:e11), &
         dims_Links(1:1), 'L11_toN', &
         long_name = 'To Node',missing_value=nodata_i4)

    call var2nc( Fname, L11_rOrder(s11:e11), &
         dims_Links(1:1), 'L11_rOrder', &
         long_name = 'Network routing order at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_label(s11:e11), &
         dims_Links(1:1), 'L11_label', &
         long_name = 'Label Id [0='', 1=HeadWater, 2=Sink] at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, merge( 1_i4, 0_i4, L11_sink(s11:e11)), &
         dims_Links(1:1), 'L11_sink', &
         long_name = '.true. if sink node reached at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_netPerm(s11:e11), &
         dims_Links(1:1), 'L11_netPerm', &
         long_name = 'Routing sequence (permutation of L11_rOrder) at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_fRow(s11:e11), &
         dims_Links(1:1), 'L11_fRow', &
         long_name = 'From row in L0 grid at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_fCol(s11:e11), &
         dims_Links(1:1), 'L11_fCol', &
         long_name = 'From col in L0 grid at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_tRow(s11:e11), &
         dims_Links(1:1), 'L11_tRow', &
         long_name = 'To row in L0 grid at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_tCol(s11:e11), &
         dims_Links(1:1), 'L11_tCol', &
         long_name = 'To Col in L0 grid at Level 11',missing_value=nodata_i4)

    call var2nc( Fname, L11_length(s11:e11), &
         dims_Links(1:1), 'L11_length', &
         long_name = 'Total length of river link [m]',missing_value=nodata_dp)

    call var2nc( Fname, L11_aFloodPlain(s11:e11), &
         dims_Links(1:1), 'L11_aFloodPlain', &
         long_name = 'Area of the flood plain [m2]',missing_value=nodata_dp)

    call var2nc( Fname, L11_slope(s11:e11), &
         dims_Links(1:1), 'L11_slope', &
         long_name = 'Average slope of river link',missing_value=nodata_dp)

    call var2nc( Fname, unpack( L0_draCell(s110:e110), mask0, nodata_i4 ), &
         dims_L0, 'L0_draCell', &
         long_name = 'Draining cell id at L11 of ith cell of L0',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L0_streamNet(s110:e110), mask0, nodata_i4 ), &
         dims_L0, 'L0_streamNet', &
         long_name = 'Stream network',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L0_floodPlain(s110:e110), mask0, nodata_i4 ), &
         dims_L0, 'L0_floodPlain', &
         long_name = 'Floodplains of stream i',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L0_draSC(s110:e110), mask0, nodata_i4 ), &
         dims_L0, 'L0_draSC', &
         long_name = 'Floodplains of stream i',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L0_L11_Id(s110:e110), mask0, nodata_i4 ), &
         dims_L0, 'L0_L11_Id', &
         long_name = 'Mapping of L11 Id on L0',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L1_L11_Id(s1:e1), mask1, nodata_i4 ), &
         dims_L1(1:2), 'L1_L11_Id', &
         long_name = 'Mapping of L11 Id on L1',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_upBound_L1(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_upBound_L1', &
         long_name = 'Row start at finer level-1 scale',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_downBound_L1(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_downBound_L1', &
         long_name = 'Row end at finer level-1 scale',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_leftBound_L1(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_leftBound_L1', &
         long_name = 'Col start at finer level-1 scale',missing_value=nodata_i4)

    call var2nc( Fname, unpack( L11_rightBound_L1(s11:e11), mask11, nodata_i4 ), &
         dims_L11(1:2), 'L11_rightBound_L1', &
         long_name = 'Col start at finer level-1 scale',missing_value=nodata_i4)

    call var2nc( Fname, (/ basin%L0_rowOutlet(iBasin), basin%L0_colOutlet(iBasin) /), &
         dims_outlet, 'L0_OutletCoord', &
         long_name = 'Outlet Coordinates at Level 0',missing_value=nodata_i4)

    call var2nc( Fname, basin_mrm%gaugeNodeList(iBasin,:), &
         dims_gauges, 'gaugeNodeList', &
         long_name = 'cell ID of gauges',missing_value=nodata_i4)

    call var2nc( Fname, basin_mrm%InflowGaugeNodeList(iBasin,:), &
         dims_inflow, 'InflowGaugeNodeList', &
         long_name = 'cell ID of gauges',missing_value=nodata_i4)

    ! free dummy variables
    deallocate( dummy_d3 )

  end subroutine write_restart_routing
  !
  subroutine read_restart_routing(iBasin, dirRestart)
    use mo_ncread, only: Get_NcVar
    use mo_string_utils, only: num2str
    use mo_mrm_constants, only: nRoutingStates
    use mo_global_variables_routing, only: &
         L11_Qmod, &
         L11_Qout, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp
    !ST: The following dependency has to be removed
    use mo_init_states, only: get_basin_info

    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    character(256), intent(in) :: dirRestart
    ! local variables
    integer(i4) :: ii
    integer(i4) :: s11 ! start index at level 11
    integer(i4) :: e11 ! end index at level 11
    integer(i4) :: ncols11 ! number of colums at level 11
    integer(i4) :: nrows11 ! number of rows at level 11
    integer(i4) :: ncells11 ! number of cells at level 11
    logical, dimension(:,:), allocatable :: mask11 ! mask at level 11
    real(dp), dimension(:,:), allocatable :: dummyD2 ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable :: dummyD3 ! dummy, 3 dimension
    character(256) :: Fname

    ! set file name
    Fname = trim(dirRestart) // trim(num2str(iBasin, '(i3.3)')) // '_states.nc'! '_restart.nc'
    
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
    allocate( dummyD3( nrows11, ncols11, nRoutingStates ) )
    call Get_NcVar( Fname, 'L11_qTIN', dummyD3 )
    do ii = 1, nRoutingStates
       L11_qTIN(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
    end do

    !  Routed outflow leaving a node
    deallocate( dummyD3 )
    allocate( dummyD3( nrows11, ncols11, nRoutingStates ) )

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

  end subroutine read_restart_routing
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
  !>        it would have been written by subroutine write_restart_files 

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Apr 2013

  !         Modified  Matthias Zink , Apr 2014 - added inflow gauge

  subroutine read_restart_L11_config( iBasin, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_kind,             only: i4, dp
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_ncread,           only: Get_NcVar
    use mo_mrm_constants,    only: nodata_dp
    use mo_global_variables_routing, only: &
         nBasins,           & ! Number of Basins
         basin_mrm,         &
         resolutionRouting, &
         level11,           &
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
    !ST The following dependencies have to be removed
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: basin ! basin database

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
    real(dp),    dimension(:),     allocatable           :: dummyD1  ! dummy, 1 dimension DP

    ! set file name
    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_L11_config.nc' ! '_restart.nc'
    call message('    Reading L11_config from ', trim(adjustl(Fname)),' ...')

    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0,&
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0, mask=mask0 )

    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1, mask=mask1 )

    ! calculate l11 grid resolutionRouting
    if(iBasin .eq. 1) then
       allocate( level11%nrows        (nBasins) )
       allocate( level11%ncols        (nBasins) )
       allocate( level11%xllcorner    (nBasins) )
       allocate( level11%yllcorner    (nBasins) )
       allocate( level11%cellsize     (nBasins) )
       allocate( level11%nodata_value (nBasins) )
    end if
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,            &
         resolutionRouting(iBasin) , &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin), &
         level11%yllcorner(iBasin), level11%cellsize(iBasin), level11%nodata_value(iBasin)        )

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
!    deallocate( dummyI2 )

    ! Row start at finer level-11 scale
!    allocate( dummyI2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L11_upBound_L1', dummyI2 )
    call append( L11_upBound_L1, pack( dummyI2, mask11) )

    ! Row end at finer level-11 scale
    call Get_NcVar( Fname, 'L11_downBound_L1', dummyI2 )
    call append( L11_downBound_L1, pack( dummyI2, mask11) )

    ! Col start at finer level-11 scale
    call Get_NcVar( Fname, 'L11_leftBound_L1', dummyI2 )
    call append( L11_leftBound_L1, pack( dummyI2, mask11) )

    ! Col end at finer level-11 scale 
    call Get_NcVar( Fname, 'L11_rightBound_L1', dummyI2 )
    call append( L11_rightBound_L1, pack(dummyI2, mask11) ) 
    deallocate( dummyI2 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From Node
    allocate( dummyI1( ncells11 ) )
    call Get_NcVar( Fname, 'L11_fromN', dummyI1)
    call append( L11_fromN, dummyI1 )

    ! To Node
    call Get_NcVar( Fname, 'L11_toN', dummyI1 )
    call append( L11_toN, dummyI1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Network routing order
    call Get_NcVar( Fname, 'L11_rOrder', dummyI1 )
    call append( L11_rOrder, dummyI1 )

    ! Label Id [0='', 1=HeadWater, 2=Sink]
    call Get_NcVar( Fname, 'L11_label', dummyI1 )
    call append( L11_label, dummyI1 )

    ! .true. if sink node reached
    call Get_NcVar( Fname, 'L11_sink', dummyI1 )
    call append( L11_sink, (dummyI1 .eq. 1_i4) )

    ! Routing sequence (permutation of L11_rOrder)
    call Get_NcVar( Fname, 'L11_netPerm', dummyI1 )
    call append( L11_netPerm, dummyI1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From row in L0 grid
    call Get_NcVar( Fname, 'L11_fRow', dummyI1 )
    call append( L11_fRow, dummyI1 )

    ! From col in L0 grid
    call Get_NcVar( Fname, 'L11_fCol', dummyI1 )
    call append( L11_fCol, dummyI1 )

    ! To row in L0 grid
    call Get_NcVar( Fname, 'L11_tRow', dummyI1 )
    call append( L11_tRow, dummyI1 )

    ! To col in L0 grid
    call Get_NcVar( Fname, 'L11_tCol', dummyI1 )
    call append( L11_tCol, dummyI1 )
    deallocate( dummyI1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate( dummyI2( nrows0, ncols0 ) )
    call Get_NcVar( Fname, 'L0_draCell', dummyI2 )
    call append( L0_draCell,     PACK ( dummyI2,  mask0)  ) 

    ! read gaugenodelist
    allocate(dummyI1( size(basin_mrm%gaugeNodeList(iBasin,:))))
    call Get_NcVar( Fname, 'gaugeNodeList', dummyI1)
    basin_mrm%gaugeNodeList( iBasin, : ) = dummyI1
    deallocate(dummyI1)

    ! read InflowGaugeNodelist
    if (basin_mrm%nInflowGauges(iBasin) > 0) then 
       allocate(dummyI1( size(basin_mrm%InflowGaugeNodeList(iBasin,:))))
       call Get_NcVar( Fname, 'InflowGaugeNodeList', dummyI1)
       basin_mrm%InflowgaugeNodeList( iBasin, : ) = dummyI1
       deallocate(dummyI1)
    end if

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
    allocate ( dummyD1( nCells11 ) )
    call Get_NcVar( Fname, 'L11_length', dummyD1 )
    call append( L11_length, dummyD1 )

    ! [m2]    Area of the flood plain
    call Get_NcVar( Fname, 'L11_aFloodPlain', dummyD1 )
    call append( L11_aFloodPlain, dummyD1  )

    ! Average slope of river link
    call Get_NcVar( Fname, 'L11_slope', dummyD1 )
    call append( L11_slope, dummyD1 )
    deallocate( dummyD1 )

  end subroutine read_restart_L11_config
end module mo_restart_routing
