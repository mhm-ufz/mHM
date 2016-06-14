!> \file mo_mrm_restart.f90

!> \brief Restart routines

!> \details This module contains the subroutines for reading and writing
!> routing related variables to file.

!> \authors Stephan Thober
!> \date Aug 2015
module mo_mrm_restart
  use mo_kind, only: i4, dp
  implicit none
  public :: mrm_write_restart
  public :: mrm_read_restart_states
  public :: mrm_read_restart_config
contains
  
  ! ------------------------------------------------------------------

  !      NAME
  !         mrm_write_restart

  !     PURPOSE
  !>        \brief write routing states and configuration

  !>        \details write configuration and state variables to a given restart
  !>        directory.

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
  !>        \note can only be used after mHM write_restart has been called because
  !>        state variables are added to the file containing the state variables of mHM.
  !>        This file must exist.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         modified, Sep 2015, Stephan Thober - added all write restart commands in this subroutine
  !                   Sep 2015, Stephan Thober - added L11_areaCell L1_ID and L1_L11_Id for routing
  !                                              resolution higher than hydrology resolution
  !                   Nov 2015, David Schaefer - mo_netcdf
  !                   May 2016, Stephan Thober - split L0_OutletCoord into L0_rowOutlet & L0_colOutlet
  !                                              because multiple outlets could exist
  subroutine mrm_write_restart(iBasin, OutPath)
    use mo_message, only: message
    use mo_netcdf, only: NcDataset, NcDimension, NcVariable
    use mo_string_utils, only: num2str
    use mo_mrm_constants, only: nRoutingStates, nodata_dp, nodata_i4
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         basin_mrm, &
         L0_cellCoor, &
         L0_Id, &
         L0_areaCell, &
         L1_Id, &
         L1_areaCell, &
         L1_L11_Id, &
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
         L11_areaCell, &
         L0_draSC, &
         L0_draCell, &
         L0_streamNet, &
         L0_floodPlain, &
         L0_L11_Id, &
         L11_L1_Id, &
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
    integer(i4) :: Noutlet ! number of outlets at level 0

    type(NcDataset)   :: nc
    type(NcDimension) :: rows0, cols0, rows1, cols1, rows11, cols11, it11, links, nout
    type(NcVariable)  :: var
    
    ! get Level0 information about the basin
    call get_basin_info_mrm( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )
    ! get Level1 information about the basin
    call get_basin_info_mrm( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )
    ! get Level11 information about the basin
    call get_basin_info_mrm( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )
    ! get Level110 information about the basin
    call get_basin_info_mrm( iBasin, 110, nrows110, ncols110, iStart=s110, iEnd=e110)

    ! set restart file name
    Fname = trim(OutPath(iBasin)) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc' 
    ! set number of outlets
    Noutlet = basin_mrm%L0_Noutlet(iBasin)

    call message('    Writing mRM restart file to ' // trim(Fname) // ' ...')

    nc     = NcDataset(fname,"w")
    rows0  = nc%setDimension("nrows0", nrows0)
    cols0  = nc%setDimension("ncols0", ncols0)
    rows1  = nc%setDimension("nrows1", nrows1)
    cols1  = nc%setDimension("ncols1", ncols1)
    rows11 = nc%setDimension("nrows11", nrows11)
    cols11 = nc%setDimension("ncols11", ncols11)
    it11   = nc%setDimension("nIT", nRoutingStates)
    links  = nc%setDimension("nLinks", size(L11_length(s11:e11)))
    nout   = nc%setDimension("Noutlet", Noutlet)

    
    var = nc%setVariable("L0_rowCoor","i32",(/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_cellCoor(s0:e0,1), mask0, nodata_i4))
    call var%setAttribute("long_name","row coordinates at Level 0")

    var = nc%setVariable("L0_colCoor", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_cellCoor(s0:e0,2), mask0, nodata_i4))
    call var%setAttribute("long_name", "col coordinates at Level 0")

    var = nc%setVariable("L0_Id", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_Id(s0:e0), mask0, nodata_i4))
    call var%setAttribute("long_name", "cell IDs at level 0")

    var = nc%setVariable("L0_areaCell", "f64", (/rows0, cols0/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L0_areaCell(s0:e0), mask0, nodata_dp))
    call var%setAttribute("long_name", "Area of a cell at level-0 [m2]")

    var = nc%setVariable("L1_basin_Mask", "i32", (/rows1, cols1/))
    call var%setFillValue(nodata_i4)
    call var%setData(merge( 1_i4, 0_i4,  &
         reshape(basin_mrm%L1_Mask(basin_mrm%L1_iStartMask(iBasin):basin_mrm%L1_iEndMask(iBasin)),&
         (/nrows1,ncols1/))))
    call var%setAttribute("long_name", "Mask at Level 1")

    var = nc%setVariable("L1_Id", "i32", (/rows1, cols1/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L1_Id(s1:e1), mask1, nodata_i4))
    call var%setAttribute("long_name", "cell IDs at level 1")

    var = nc%setVariable("L1_areaCell", "f64", (/rows1, cols1/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L1_areaCell(s1:e1), mask1, nodata_dp))
    call var%setAttribute("long_name", "Effective area of cell at this level [km2]")
    
    var = nc%setVariable("L1_L11_Id", "i32", (/rows1, cols1/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L1_L11_Id(s1:e1), mask1, nodata_i4))
    call var%setAttribute("long_name", "Mapping of L1 Id on L11")

    var = nc%setVariable("L11_Qmod", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_Qmod(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "simulated discharge at each node at level 11")

    var = nc%setVariable("L11_qOUT", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_qOUT(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Total outflow from cells L11 at time tt at level 11")

    allocate( dummy_d3( nrows11, ncols11, nRoutingStates ) )
    do ii = 1, size( dummy_d3, 3 )
       dummy_d3(:,:,ii) = unpack( L11_qTIN(s11:e11,ii), mask11, nodata_dp )
    end do
    var = nc%setVariable("L11_qTIN", "f64", (/rows11, cols11, it11/))
    call var%setFillValue(nodata_dp)
    call var%setData(dummy_d3)
    call var%setAttribute("long_name", "Total discharge inputs at t-1 and t at level 11")

    do ii = 1, size( dummy_d3, 3 )
       dummy_d3(:,:,ii) = unpack( L11_qTR(s11:e11,ii), mask11, nodata_dp )
    end do
    var = nc%setVariable("L11_qTR", "f64", (/rows11, cols11, it11/))
    call var%setFillValue(nodata_dp)
    call var%setData(dummy_d3)
    call var%setAttribute("long_name", "Routed outflow leaving a node at level 11")

    var = nc%setVariable("L11_K", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_K(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "kappa: Muskingum travel time parameter at level 11")

    var = nc%setVariable("L11_xi", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_xi(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "xi: Muskingum diffusion parameter at level 11")

    var = nc%setVariable("L11_C1", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_C1(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Routing parameter C1=f(K,xi, DT) (Chow, 25-41) at level 11")

    var = nc%setVariable("L11_C2", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_C2(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Routing parameter C2=f(K,xi, DT) (Chow, 25-41) at level 11")

    var = nc%setVariable("L11_FracFPimp", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_FracFPimp(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Fraction of the flood plain with impervious cover at level 11")

    ! ----------------------------------------------------------
    ! L11 config set
    ! ----------------------------------------------------------
    var = nc%setVariable("L11_basin_Mask", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(merge(1_i4, 0_i4,  &
         reshape(basin_mrm%L11_Mask(basin_mrm%L11_iStartMask(iBasin):basin_mrm%L11_iEndMask(iBasin)),&
         (/nrows11,ncols11/))))
    call var%setAttribute("long_name", "Mask at Level 11")

    var = nc%setVariable("L11_rowCoor", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_cellCoor(s11:e11,1), mask11, nodata_i4))
    call var%setAttribute("long_name", "row coordinates at Level 11")

    var = nc%setVariable("L11_colCoor", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_cellCoor(s11:e11,2), mask11, nodata_i4))
    call var%setAttribute("long_name", "col coordinates at Level 11")

    var = nc%setVariable("L11_Id", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_Id(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "cell Ids at Level 11")

    var = nc%setVariable("L11_areaCell", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_areaCell(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "cell area at Level 11")

    var = nc%setVariable("L11_fDir", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_fDir(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "flow Direction at Level 11")

    var = nc%setVariable("L11_rowOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_rowOut(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Grid vertical location of the Outlet at Level 11")

    var = nc%setVariable("L11_colOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_colOut(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Grid horizontal location of the Outlet at Level 11")

    var = nc%setVariable("L11_upBound_L0", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_upBound_L0(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Row start at finer level-0 scale of Level 11 cell")

    var = nc%setVariable("L11_downBound_L0", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_downBound_L0(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Row end at finer level-0 scale of Level 11 cell")

    var = nc%setVariable("L11_leftBound_L0", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_leftBound_L0(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Col start at finer level-0 scale of Level 11 cell")

    var = nc%setVariable("L11_rightBound_L0", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_rightBound_L0(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Col end at finer level-0 scale of Level 11 cell")

    var = nc%setVariable("L11_fromN", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_fromN(s11:e11))
    call var%setAttribute("long_name", "From Node")

    var = nc%setVariable("L11_toN", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_toN(s11:e11))
    call var%setAttribute("long_name", "To Node")

    var = nc%setVariable("L11_rOrder", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_rOrder(s11:e11))
    call var%setAttribute("long_name", "Network routing order at Level 11")

    var = nc%setVariable("L11_label", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_label(s11:e11))
    call var%setAttribute("long_name", "Label Id [0='', 1=HeadWater, 2=Sink] at Level 11")

    var = nc%setVariable("L11_sink", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(merge( 1_i4, 0_i4, L11_sink(s11:e11)))
    call var%setAttribute("long_name", ".true. if sink node reached at Level 11")

    var = nc%setVariable("L11_netPerm", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_netPerm(s11:e11))
    call var%setAttribute("long_name", "Routing sequence (permutation of L11_rOrder) at Level 11")

    var = nc%setVariable("L11_fRow", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_fRow(s11:e11))
    call var%setAttribute("long_name", "From row in L0 grid at Level 11")

    var = nc%setVariable("L11_fCol", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_fCol(s11:e11))
    call var%setAttribute("long_name", "From col in L0 grid at Level 11")

    var = nc%setVariable("L11_tRow", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_tRow(s11:e11))
    call var%setAttribute("long_name", "To row in L0 grid at Level 11")

    var = nc%setVariable("L11_tCol", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_tCol(s11:e11))
    call var%setAttribute("long_name", "To Col in L0 grid at Level 11")

    var = nc%setVariable("L11_length", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_length(s11:e11))
    call var%setAttribute("long_name", "Total length of river link [m]")

    var = nc%setVariable("L11_aFloodPlain", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_aFloodPlain(s11:e11))
    call var%setAttribute("long_name", "Area of the flood plain [m2]")

    var = nc%setVariable("L11_slope", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_slope(s11:e11))
    call var%setAttribute("long_name", "Average slope of river link")

    var = nc%setVariable("L0_draCell", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_draCell(s110:e110), mask0, nodata_i4))
    call var%setAttribute("long_name", "Draining cell id at L11 of ith cell of L0")

    var = nc%setVariable("L0_streamNet", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_streamNet(s110:e110), mask0, nodata_i4))
    call var%setAttribute("long_name", "Stream network")

    var = nc%setVariable("L0_floodPlain", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_floodPlain(s110:e110), mask0, nodata_i4))
    call var%setAttribute("long_name", "Floodplains of stream i")

    var = nc%setVariable("L0_draSC", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_draSC(s110:e110), mask0, nodata_i4))
    call var%setAttribute("long_name", "Floodplains of stream i")

    var = nc%setVariable("L0_L11_Id", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_L11_Id(s110:e110), mask0, nodata_i4))
    call var%setAttribute("long_name", "Mapping of L11 Id on L0")

    var = nc%setVariable("L11_L1_Id", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_L1_Id(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Mapping of L1 Id on L11")

    var = nc%setVariable("L11_upBound_L1", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_upBound_L1(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Row start at finer level-1 scale")

    var = nc%setVariable("L11_downBound_L1", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack( L11_downBound_L1(s11:e11), mask11, nodata_i4 ))
    call var%setAttribute("long_name", "Row end at finer level-1 scale")

    var = nc%setVariable("L11_leftBound_L1", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_leftBound_L1(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Col start at finer level-1 scale")

    var = nc%setVariable("L11_rightBound_L1", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_rightBound_L1(s11:e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Col start at finer level-1 scale")

    var = nc%setVariable("L0_RowOutletCoord", "i32", (/nout/))
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm%L0_rowOutlet(:Noutlet, iBasin))
    call var%setAttribute("long_name", "Row outlet coordinates at level 0")

    var = nc%setVariable("L0_ColOutletCoord", "i32", (/nout/))
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm%L0_colOutlet(:Noutlet, iBasin))
    call var%setAttribute("long_name", "Column outlet coordinates at level 0")

    var = nc%setVariable("gaugeNodeList", "i32", &
         (/nc%setDimension("Ngauges", size(basin_mrm%gaugeNodeList(iBasin,:)))/) &
         )
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm%gaugeNodeList(iBasin,:))
    call var%setAttribute("long_name", "cell ID of gauges")

    var = nc%setVariable("InflowGaugeNodeList", "i32", &
         (/nc%setDimension("nInflowGauges", size(basin_mrm%InflowGaugeNodeList(iBasin,:)))/) &
         )
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm%InflowGaugeNodeList(iBasin,:))
    call var%setAttribute("long_name", "cell ID of gauges")

    ! free dummy variables
    deallocate( dummy_d3 )
    call nc%close()
    
  end subroutine mrm_write_restart

  ! ------------------------------------------------------------------

  !      NAME
  !         mrm_read_restart_states

  !     PURPOSE
  !>        \brief read routing states

  !>        \details This subroutine reads the routing states from
  !>        mRM_states_<basin_id>.nc that has to be in the given
  !>        path directory. This subroutine has to be called directly
  !>        each time the mHM_eval or mRM_eval is called such that the
  !>        the states are always the same at the first simulation time
  !>        step, crucial for optimization.

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
  !>        \note This subroutine has to be called directly
  !>        each time the mHM_eval or mRM_eval is called such that the
  !>        the states are always the same at the first simulation time
  !>        step, crucial for optimization.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Sep 2015
  !         Modified  Mar 2016, David Schaefer - mo_netcdf
  !                   May 2016, Stephan Thober - split L0_OutletCoord into L0_rowOutlet & L0_colOutlet
  !                                              because multiple outlets could exist
  !
  subroutine mrm_read_restart_states(iBasin, InPath)
    use mo_netcdf, only: NcDataset, NcVariable
    use mo_string_utils, only: num2str
    use mo_mrm_constants, only: nRoutingStates
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: &
         L11_Qmod, &
         L11_Qout, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp
    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    character(256), intent(in) :: InPath
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
    character(256)   :: fname
    type(NcDataset)   :: nc
    type(NcVariable) :: var
    
    ! set file name
    fname = trim(InPath) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc' 
    
    ! get basin info at L11 including ncells, start, end, and mask
    call get_basin_info_mrm(iBasin, 11, nrows11, ncols11, ncells=ncells11, iStart=s11, iEnd=e11, mask=mask11)

    nc = NcDataset(fname, "r")

    ! simulated discharge at each node
    var = nc%getVariable("L11_Qmod")
    call var%getData(dummyD2)
    L11_Qmod(s11:e11) = pack( dummyD2, mask11 ) 

    ! Total outflow from cells L11 at time tt
    var = nc%getVariable("L11_qOUT")
    call var%getData(dummyD2)
    L11_qOUT(s11:e11) = pack( dummyD2, mask11 )

    ! Total discharge inputs at t-1 and t
    var = nc%getVariable("L11_qTIN")
    call var%getData(dummyD3)
    do ii = 1, nRoutingStates
       L11_qTIN(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
    end do

    !  Routed outflow leaving a node
    var = nc%getVariable("L11_qTR")
    call var%getData(dummyD3)
    do ii = 1, nRoutingStates
       L11_qTR(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
    end do

    ! kappa: Muskingum travel time parameter.
    var = nc%getVariable("L11_K")
    call var%getData(dummyD2)
    L11_K(s11:e11) = pack( dummyD2, mask11 )

    !  xi:    Muskingum diffusion parameter
    var = nc%getVariable("L11_xi")
    call var%getData(dummyD2)
    L11_xi(s11:e11) = pack( dummyD2, mask11 )

    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    var = nc%getVariable("L11_C1")
    call var%getData(dummyD2)
    L11_C1(s11:e11) = pack( dummyD2, mask11 )

    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    var = nc%getVariable("L11_C2")
    call var%getData(dummyD2)
    L11_C2(s11:e11) = pack( dummyD2, mask11 )

    ! Fraction of the flood plain with impervious cover
    var = nc%getVariable("L11_FracFPimp")
    call var%getData(dummyD2)
    L11_FracFPimp(s11:e11) = pack( dummyD2, mask11 )

    ! free memory
    deallocate( dummyD2, dummyD3 )

  end subroutine mrm_read_restart_states

  ! ------------------------------------------------------------------

  !      NAME
  !         mrm_read_restart

  !     PURPOSE
  !>        \brief reads Level 11 configuration from a restart directory

  !>        \details read Level 11 configuration variables from a given restart
  !>        directory and initializes all Level 11 configuration variables,
  !>        that are initialized in L11_variable_init,
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

  !         Modified Apr 2014, Matthias Zink - added inflow gauge
  !                  Aug 2015, Stephan Thober - adapted for mRM
  !                  Sep 2015, Stephan Thober - added all read restart commands in this subroutine
  !                  Sep 2015, Stephan Thober - added L11_areaCell, L1_ID and L1_L11_Id for routing
  !                                             resolution higher than hydrology resolution
  !                  Mar 2016, David Schaefer - mo_netcdf

  subroutine mrm_read_restart_config( iBasin, InPath )
    use mo_kind, only: i4, dp
    use mo_message, only: message
    use mo_string_utils, only: num2str
    use mo_kind, only: i4, dp
    use mo_append, only: append
    use mo_netcdf, only: NcDataset, NcVariable
    use mo_mrm_constants, only: nodata_dp, nodata_i4
    use mo_mrm_tools, only: get_basin_info_mrm, calculate_grid_properties
    use mo_mrm_global_variables, only: &
         L0_Basin, & ! check whether L0_Basin should be read
         nBasins,           & ! Number of Basins
         basin_mrm,         &
         resolutionHydrology, &
         resolutionRouting, &
         level1, &
         level11,           &
         L0_cellCoor   , &
         L0_Id         , & ! Ids of grid at level-0
         L0_areaCell, & ! Ids of grid at level-0
         L1_Id, &
         L1_areaCell, & ! [km2] Effective area of cell at this level
         L1_L11_Id,         &
         L11_cellCoor,      & ! cell Coordinates at Level 11
         L11_Id,            & ! cell Ids at Level 11
         L11_areaCell,      &
         L11_nCells,        & ! Number of Cells at Level 11
         L1_nCells,         & ! Number of Cells at Level 1
         L0_nCells,         & ! Number of Cells at Level 0
         L0_draSC,          &
         L0_L11_Id,         &
         L11_L1_Id,         &
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
         L11_nOutlets,      &
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

    integer(i4), intent(in):: iBasin
    character(256), intent(in) :: InPath ! list of Output paths per Basin

    ! local variables
    character(256) :: fname

    ! local variables
    integer(i4) :: nrows0 ! Number of rows at level 0
    integer(i4) :: ncols0 ! Number of cols at level 0
    logical, dimension(:,:), allocatable :: mask0 ! Mask at Level 0
    integer(i4) :: nrows1 ! Number of rows at level 1
    integer(i4) :: ncols1 ! Number of cols at level 1
    logical, dimension(:,:), allocatable :: mask1    ! Mask at Level 1
    integer(i4) :: nrows11 ! Number of rows at level 11
    integer(i4) :: ncols11 ! Number of cols at level 11
    logical, dimension(:,:), allocatable :: mask11 ! Mask at Level 11
    integer(i4) :: nCells11 ! Number of cells at lev. 11
    real(dp) :: xllcorner0, yllcorner0
    real(dp) :: cellsize0
    integer(i4) :: iStart0, iEnd0

    ! DUMMY variables
    integer(i4) :: Noutlet
    integer(i4) :: old_Noutlet
    integer(i4), dimension(:), allocatable :: dummyI1 ! dummy, 1 dimension I4
    integer(i4), dimension(:,:), allocatable :: dummy
    integer(i4), dimension(:,:), allocatable :: dummyI2 ! dummy, 2 dimension I4
    integer(i4), dimension(:,:), allocatable :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp), dimension(:), allocatable :: dummyD1 ! dummy, 1 dimension DP
    real(dp), dimension(:,:), allocatable :: dummyD2 ! dummy, 2 dimension
    type(NcDataset) :: nc
    type(NcVariable) :: var

    ! set file name
    fname = trim(InPath) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'
    call message('        Reading mRM restart file:  ', trim(adjustl(Fname)),' ...')

    if (iBasin .eq. 1) then
       allocate(level1%nrows(nBasins))
       allocate(level1%ncols(nBasins))
       allocate(level1%xllcorner(nBasins))
       allocate(level1%yllcorner(nBasins))
       allocate(level1%cellsize(nBasins))
       allocate(level1%nodata_value(nBasins))

       allocate(level11%nrows(nBasins))
       allocate(level11%ncols(nBasins))
       allocate(level11%xllcorner(nBasins))
       allocate(level11%yllcorner(nBasins))
       allocate(level11%cellsize(nBasins))
       allocate(level11%nodata_value(nBasins))
 
       allocate(basin_mrm%L1_iStart(nBasins))
       allocate(basin_mrm%L1_iEnd  (nBasins))
       allocate(basin_mrm%L1_iStartMask(nBasins))
       allocate(basin_mrm%L1_iEndMask   (nBasins))

       allocate(basin_mrm%L11_iStart(nBasins))
       allocate(basin_mrm%L11_iEnd(nBasins))
       allocate(basin_mrm%L11_iStartMask(nBasins))
       allocate(basin_mrm%L11_iEndMask(nBasins))

       allocate(basin_mrm%L0_Noutlet(nBasins))
       allocate(basin_mrm%L0_rowOutlet(1, nBasins))
       allocate(basin_mrm%L0_colOutlet(1, nBasins))
       basin_mrm%L0_Noutlet = nodata_i4
       basin_mrm%L0_rowOutlet = nodata_i4
       basin_mrm%L0_colOutlet = nodata_i4
    end if

    ! grid properties
    call get_basin_info_mrm(iBasin, 0, nrows0, ncols0, iStart= iStart0, iEnd=iEnd0, mask=mask0, &
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0)
    
    !call get_basin_info_mrm( iBasin, 1, nrows1, ncols1)
    call calculate_grid_properties(nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp, &
         resolutionHydrology(iBasin) ,                                                           &
         level1%nrows(iBasin), level1%ncols(iBasin), level1%xllcorner(iBasin),                   &
         level1%yllcorner(iBasin), level1%cellsize(iBasin), level1%nodata_value(iBasin))
    !
    ! level-1 information
    call calculate_grid_properties(nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp, &
         resolutionRouting(iBasin),                                                              &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin),                &
         level11%yllcorner(iBasin), level11%cellsize(iBasin), level11%nodata_value(iBasin))

    call get_basin_info_mrm(iBasin, 1, nrows1, ncols1)

    call get_basin_info_mrm(iBasin, 11, nrows11, ncols11)

    nc = NcDataset(fname, "r")

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L0 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! check whether L0 data should be read
    if (iBasin .eq. 1) then
       ! read L0 cellCoor
       allocate(dummyI22(count(mask0),2))

       var = nc%getVariable("L0_rowCoor")
       call var%getData(dummyI2)
       dummyI22(:,1) = pack(dummyI2, mask0)

       var = nc%getVariable("L0_colCoor")
       call var%getData(dummyI2)
       dummyI22(:,2) = pack(dummyI2, mask0)

       call append(L0_cellCoor, dummyI22)
       deallocate(dummyI22)

       var = nc%getVariable("L0_Id")
       call var%getData(dummyI2)
       call append(L0_Id, pack(dummyI2, mask0))
       
       var = nc%getVariable("L0_areaCell")
       call var%getData(dummyD2)
       call append( L0_areaCell, pack(dummyD2, mask0) )
       
    else
       if ( L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1) ) then
          ! read L0 cellCoor
          allocate(dummyI22(count(mask0), 2))

          var = nc%getVariable("L0_rowCoor")
          call var%getData(dummyI2)
          dummyI22(:,1) = pack( dummyI2, mask0 )

          var = nc%getVariable("L0_colCoor")
          call var%getData(dummyI2)
          dummyI22(:,2) = pack(dummyI2, mask0)

          call append(L0_cellCoor, dummyI22)
          deallocate(dummyI22)
          
          var = nc%getVariable("L0_Id")
          call var%getData(dummyI2)
          call append(L0_Id, pack(dummyI2, mask0))
          
          var = nc%getVariable("L0_areaCell")
          call var%getData(dummyD2)
          call append(L0_areaCell, pack(dummyD2, mask0))
          
       end if
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L1 mask
    var = nc%getVariable("L1_basin_Mask")
    call var%getData(dummyI2)
    mask1 = (dummyI2 .eq. 1_i4)
    call append(basin_mrm%L1_Mask, reshape((dummyI2 .eq. 1_i4), (/nrows1*ncols1/)))

    ! update basin database
    if( iBasin .eq. 1 ) then
       basin_mrm%L1_iStart(iBasin) = 1
       basin_mrm%L1_iStartMask(iBasin) = 1
    else
       basin_mrm%L1_iStart(iBasin) = basin_mrm%L1_iEnd(iBasin-1) + 1
       basin_mrm%L1_iStartMask(iBasin) = basin_mrm%L1_iEndMask(iBasin-1) + 1
    end if
    basin_mrm%L1_iEnd(iBasin) = basin_mrm%L1_iStart(iBasin) + count(mask1) - 1
    basin_mrm%L1_iEndMask(iBasin) = basin_mrm%L1_iStartMask(iBasin) + nrows1*ncols1 - 1

    ! get level1 information including mask
    call get_basin_info_mrm(iBasin, 1, nrows1, ncols1, mask=mask1)
    !
    var = nc%getVariable("L1_Id")
    call var%getData(dummyI2)
    call append(L1_Id, pack(dummyI2, mask1))

    var = nc%getVariable("L1_L11_Id")
    call var%getData(dummyI2)
    call append(L1_L11_Id, pack(dummyI2, mask1))

    ! allocate( dummyD2( nrows1, ncols1 ) )
    var = nc%getVariable("L1_areaCell")
    call var%getData(dummyD2)
    call append(L1_areaCell, pack( dummyD2, mask1))

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L11 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L11 mask
    var = nc%getVariable("L1_basin_Mask")
    call var%getData(dummyI2)
    mask11 = (dummyI2 .eq. 1_i4)
    call append(basin_mrm%L11_Mask, reshape(mask11, (/nrows11*ncols11/)))

    ! get Number of cells
    nCells11 = count(mask11)

    ! update basin database
    if (iBasin .eq. 1) then
       basin_mrm%L11_iStart(iBasin) = 1
       basin_mrm%L11_iStartMask(iBasin) = 1
    else
       basin_mrm%L11_iStart(iBasin) = basin_mrm%L11_iEnd(iBasin-1) + 1
       basin_mrm%L11_iStartMask(iBasin) = basin_mrm%L11_iEndMask(iBasin-1) + 1
    end if
    basin_mrm%L11_iEnd(iBasin) = basin_mrm%L11_iStart(iBasin) + nCells11 - 1
    basin_mrm%L11_iEndMask(iBasin) = basin_mrm%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    ! read L11 cellCoor
    allocate(dummyI22(count(mask11), 2))
    var = nc%getVariable("L11_rowCoor")
    call var%getData(dummyI2)
    dummyI22(:,1) = pack(dummyI2, mask11)

    var = nc%getVariable("L11_colCoor")
    call var%getData(dummyI2)
    dummyI22(:,2) = pack(dummyI2, mask11)
    call append(L11_cellCoor, dummyI22)
    deallocate(dummyI22)

    ! read L11 IDs
    var = nc%getVariable("L11_Id")
    call var%getData(dummyI2)
    call append(L11_Id, pack(dummyI2, mask11))

    var = nc%getVariable("L11_areaCell")
    call var%getData(dummyD2)
    call append(L11_areaCell, pack(dummyD2, mask11))

    ! update Number of cells at Level 11
    L11_nCells = size(L11_Id, 1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L0 outlet Coordinates
    var = nc%getVariable("L0_RowOutletCoord")
    call var%getData(dummyI1)

    Noutlet = size(dummyI1, dim=1)
    basin_mrm%L0_Noutlet(iBasin) = Noutlet
    deallocate(dummyI2)
    allocate(dummyI2(Noutlet, 2))
    dummyI2(:, 1) = dummyI1

    var = nc%getVariable("L0_ColOutletCoord")
    call var%getData(dummyI1)
    dummyI2(:, 2) = dummyI1

    old_Noutlet = size(basin_mrm%L0_rowOutlet, dim=1)
    ! store outlet coordinates to basin structure
    if (Noutlet .le. old_Noutlet) then
       basin_mrm%L0_rowOutlet(:Noutlet, iBasin) = dummyI2(:, 1)
       basin_mrm%L0_colOutlet(:Noutlet, iBasin) = dummyI2(:, 2)
    else
       ! store up to size of old_Noutlet
       basin_mrm%L0_rowOutlet(:old_Noutlet, iBasin) = dummyI2(:old_Noutlet, 1)
       basin_mrm%L0_colOutlet(:old_Noutlet, iBasin) = dummyI2(:old_Noutlet, 2)
       ! enlarge rowOutlet and colOutlet in basin_mrm structure
       allocate(dummy(Noutlet - old_Noutlet, nBasins))
       dummy = nodata_i4
       dummy(:, iBasin) = dummyI2(old_Noutlet + 1:, 1)
       call append(basin_mrm%L0_rowOutlet, dummy)
       dummy(:, iBasin) = dummyI2(old_Noutlet + 1:, 2)
       call append(basin_mrm%L0_colOutlet, dummy)
       deallocate(dummy)
    end if
    deallocate(dummyI1, dummyI2)

    ! read L0 draining cell index
    var = nc%getVariable("L0_draSC")
    call var%getData(dummyI2)
    call append(L0_draSC, pack(dummyI2, mask0))

    ! Mapping of L11 Id on L0
    var = nc%getVariable("L0_L11_Id")
    call var%getData(dummyI2)
    call append(L0_L11_Id, pack(dummyI2, mask0))

    ! L11 data sets
    ! Mapping of L1 Ids on L1
    var = nc%getVariable("L11_L1_Id")
    call var%getData(dummyI2)
    call append(L11_L1_Id, pack(dummyI2, mask11))

    ! Flow direction (standard notation)
    var = nc%getVariable("L11_fDir")
    call var%getData(dummyI2)
    call append(L11_fDir, pack(dummyI2, mask11))

    ! Grid vertical location of the Outlet
    var = nc%getVariable("L11_rowOut")
    call var%getData(dummyI2)
    call append(L11_rowOut, pack(dummyI2, mask11))

    ! Grid horizontal location  of the Outlet
    var = nc%getVariable("L11_colOut")
    call var%getData(dummyI2)
    call append(L11_colOut, pack(dummyI2, mask11))

    ! Row start at finer level-0 scale
    var = nc%getVariable("L11_upBound_L0")
    call var%getData(dummyI2)
    call append(L11_upBound_L0, pack(dummyI2, mask11))

    ! Row end at finer level-0 scale
    var = nc%getVariable("L11_downBound_L0")
    call var%getData(dummyI2)
    call append(L11_downBound_L0, pack(dummyI2, mask11))

    ! Col start at finer level-0 scale
    var = nc%getVariable("L11_leftBound_L0")
    call var%getData(dummyI2)
    call append(L11_leftBound_L0, pack(dummyI2, mask11))

    ! Col end at finer level-0 scale
    var = nc%getVariable("L11_rightBound_L0")
    call var%getData(dummyI2)
    call append(L11_rightBound_L0, pack(dummyI2, mask11))

    ! Row start at finer level-11 scale
    var = nc%getVariable("L11_upBound_L1")
    call var%getData(dummyI2)
    call append(L11_upBound_L1, pack(dummyI2, mask11))

    ! Row end at finer level-11 scale
    var = nc%getVariable("L11_downBound_L1")
    call var%getData(dummyI2)
    call append(L11_downBound_L1, pack(dummyI2, mask11))

    ! Col start at finer level-11 scale
    var = nc%getVariable("L11_leftBound_L1")
    call var%getData(dummyI2)
    call append(L11_leftBound_L1, pack(dummyI2, mask11))

    ! Col end at finer level-11 scale
    var = nc%getVariable("L11_rightBound_L1")
    call var%getData(dummyI2)
    call append(L11_rightBound_L1, pack(dummyI2, mask11))
    
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From Node
    var = nc%getVariable("L11_fromN")
    call var%getData(dummyI1)
    call append(L11_fromN, dummyI1)

    ! To Node
    var = nc%getVariable("L11_toN")
    call var%getData(dummyI1)
    call append(L11_toN, dummyI1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Network routing order
    var = nc%getVariable("L11_rOrder")
    call var%getData(dummyI1)
    call append(L11_rOrder, dummyI1)

    ! Label Id [0='', 1=HeadWater, 2=Sink]
    var = nc%getVariable("L11_label")
    call var%getData(dummyI1)
    call append(L11_label, dummyI1)

    ! .true. if sink node reached
    var = nc%getVariable("L11_sink")
    call var%getData(dummyI1)
    call append(L11_sink, (dummyI1 .eq. 1_i4))
    ! append Number of Outlets at Level 11
    call append(L11_nOutlets, count((dummyI1 .eq. 1_i4)))

    ! Routing sequence (permutation of L11_rOrder)
    var = nc%getVariable("L11_netPerm")
    call var%getData(dummyI1)
    call append(L11_netPerm, dummyI1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From row in L0 grid
    var = nc%getVariable("L11_fRow")
    call var%getData(dummyI1)
    call append(L11_fRow, dummyI1)

    ! From col in L0 grid
    var = nc%getVariable("L11_fCol")
    call var%getData(dummyI1)
    call append(L11_fCol, dummyI1)

    ! To row in L0 grid
    var = nc%getVariable("L11_tRow")
    call var%getData(dummyI1)
    call append(L11_tRow, dummyI1)

    ! To col in L0 grid
    var = nc%getVariable("L11_tCol")
    call var%getData(dummyI1)
    call append(L11_tCol, dummyI1)


    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    var = nc%getVariable("L0_draCell")
    call var%getData(dummyI2)
    call append(L0_draCell, pack(dummyI2, mask0))

    ! read gaugenodelist
    var = nc%getVariable("gaugeNodeList")
    call var%getData(dummyI1)
    basin_mrm%gaugeNodeList(iBasin,:) = dummyI1

    ! read InflowGaugeNodelist
    if (basin_mrm%nInflowGauges(iBasin) > 0) then
       var = nc%getVariable("InflowGaugeNodeList")
       call var%getData(dummyI1)
       basin_mrm%InflowgaugeNodeList(iBasin,:) = dummyI1
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! L0 data sets
    ! Stream network
    var = nc%getVariable("L0_streamNet")
    call var%getData(dummyI2)
    call append(L0_streamNet, pack(dummyI2, mask0))

    ! Floodplains of stream i
    var = nc%getVariable("L0_floodPlain")
    call var%getData(dummyI2)
    call append(L0_floodPlain, pack(dummyI2, mask0))

    ! L11 network data sets
    ! [m]     Total length of river link
    var = nc%getVariable("L11_length")
    call var%getData(dummyD1)
    call append(L11_length, dummyD1)

    ! [m2]    Area of the flood plain
    var = nc%getVariable("L11_aFloodPlain")
    call var%getData(dummyD1)
    call append(L11_aFloodPlain, dummyD1)

    ! Average slope of river link
    var = nc%getVariable("L11_slope")
    call var%getData(dummyD1)
    call append(L11_slope, dummyD1)

    ! update Number of cells at Level 1 and Level 0
    L1_nCells = size(L0_areaCell, 1)
    L0_nCells = size(L1_areaCell, 1)
    
  end subroutine mrm_read_restart_config

end module mo_mrm_restart
