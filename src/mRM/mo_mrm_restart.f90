!>       \file mo_mrm_restart.f90

!>       \brief Restart routines

!>       \details This module contains the subroutines for reading and writing
!>       routing related variables to file.

!>       \authors Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_restart
  use mo_kind, only : i4, dp
  implicit none
  public :: mrm_write_restart
  public :: mrm_read_restart_states
  public :: mrm_read_restart_config
contains

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_write_restart

  !    PURPOSE
  !>       \brief write routing states and configuration

  !>       \details write configuration and state variables to a given restart
  !>       directory.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin"                   number of basin
  !>       \param[in] "character(256), dimension(:) :: OutPath" list of Output paths per Basin

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober    Sep 2015 - added all write restart commands in this subroutine
  ! Stephan Thober    Sep 2015 - added L11_areaCell L1_ID and L1_L11_Id for routing
  !                              resolution higher than hydrology resolution
  ! David Schaefer    Nov 2015 - mo_netcdf
  ! Stephan Thober    May 2016 - split L0_OutletCoord into L0_rowOutlet & L0_colOutlet
  !                              because multiple outlets could exist
  ! Stephan Thober    Nov 2016 - added L11_TSrout, ProcessMatrix
  ! Matthias Kelbling Aug 2017 - added L11_fAcc, L0_slope, L0_celerity, L11_celerity, L11_meandering
  ! Robert Schweppe   Jun 2018 - refactoring and reformatting

  subroutine mrm_write_restart(iBasin, OutPath)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_restart, only : write_grid_info
    use mo_common_variables, only : level1, nLCoverScene, processMatrix
    use mo_message, only : message
    use mo_mrm_constants, only : nRoutingStates
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, L11_L1_Id, L11_Qmod, &
                                        L11_TSrout, L11_aFloodPlain, L11_colOut, L11_colOut, L11_fCol, L11_fDir, &
                                        L11_fDir, L11_fRow, L11_fromN, L11_label, L11_length, L11_nLinkFracFPimp, &
                                        L11_netPerm, L11_qOUT, L11_qTIN, L11_qTR, L11_rOrder, L11_rowOut, L11_rowOut, &
                                        L11_sink, L11_slope, L11_tCol, L11_tRow, L11_toN, L11_xi, L1_L11_Id, basin_mrm, &
                                        level11
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of basin
    integer(i4), intent(in) :: iBasin

    ! list of Output paths per Basin
    character(256), dimension(:), intent(in) :: OutPath

    character(256) :: Fname

    integer(i4) :: ii

    ! start index at level 1
    integer(i4) :: s1

    ! end index at level 1
    integer(i4) :: e1

    ! mask at level 1
    logical, dimension(:, :), allocatable :: mask1

    ! start index at level 11
    integer(i4) :: s11

    ! end index at level 11
    integer(i4) :: e11

    ! number of colums at level 11
    integer(i4) :: ncols11

    ! number of rows at level 11
    integer(i4) :: nrows11

    ! mask at level 11
    logical, dimension(:, :), allocatable :: mask11

    ! dummy variable
    real(dp), dimension(:, :, :), allocatable :: dummy_d3

    type(NcDataset) :: nc

    type(NcDimension) :: rows1, cols1, rows11, cols11, it11, lcscenes

    type(NcDimension) :: links, nts, nproc

    type(NcVariable) :: var


    ! get Level1 and Level11 information about the basin
    s1 = level1(iBasin)%iStart
    e1 = level1(iBasin)%iEnd
    mask1 = level1(iBasin)%mask
    s11 = level11(iBasin)%iStart
    e11 = level11(iBasin)%iEnd
    mask11 = level11(iBasin)%mask
    ncols11 = level11(iBasin)%ncols
    nrows11 = level11(iBasin)%nrows
    allocate(dummy_d3(nrows11, ncols11, nRoutingStates))

    ! set restart file name
    Fname = trim(OutPath(iBasin)) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'

    call message('    Writing mRM restart file to ' // trim(Fname) // ' ...')

    nc = NcDataset(fname, "w")

    call write_grid_info(level1(iBasin), "1", nc)
    call write_grid_info(level11(iBasin), "11", nc)

    rows1 = nc%getDimension("nrows1")
    cols1 = nc%getDimension("ncols1")
    rows11 = nc%getDimension("nrows11")
    cols11 = nc%getDimension("ncols11")

    it11 = nc%setDimension("nIT", nRoutingStates)
    links = nc%setDimension("nLinks", size(L11_length(s11 : e11)))
    nts = nc%setDimension("TS", 1)
    nproc = nc%setDimension("Nprocesses", size(processMatrix, dim = 1))
    lcscenes = nc%setDimension("LCoverScenes", nLCoverScene)


    var = nc%setVariable("L0_slope_mRM", "f64", (/rows0, cols0/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L0_slope_mRM(s0:e0), mask0, nodata_dp))
    call var%setAttribute("long_name", "slope at Level 0 [%]")

    var = nc%setVariable("L0_celerity", "f64", (/rows0, cols0/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L0_celerity(s0:e0), mask0, nodata_dp))
    call var%setAttribute("long_name", "celerity at Level 0 [m/s]")

    var = nc%setVariable("L0_fAcc", "i32", (/rows0, cols0/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L0_fAcc(s0:e0), mask0, nodata_i4))
    call var%setAttribute("long_name", "fAcc at Level 0 (n cells)")

    ! add processMatrix
    var = nc%setVariable("ProcessMatrix", "i32", (/nproc/))
    call var%setFillValue(nodata_i4)
    call var%setData(processMatrix(:, 1))
    call var%setAttribute("long_name", "Process Matrix")

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
    call var%setData(unpack(L1_L11_Id(s1 : e1), mask1, nodata_i4))
    call var%setAttribute("long_name", "Mapping of L1 Id on L11")

    var = nc%setVariable("L11_Qmod", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_Qmod(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "simulated discharge at each node at level 11")

    var = nc%setVariable("L11_qOUT", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_qOUT(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Total outflow from cells L11 at time tt at level 11")

    do ii = 1, size(dummy_d3, 3)
      dummy_d3(:, :, ii) = unpack(L11_qTIN(s11 : e11, ii), mask11, nodata_dp)
    end do
    var = nc%setVariable("L11_qTIN", "f64", (/rows11, cols11, it11/))
    call var%setFillValue(nodata_dp)
    call var%setData(dummy_d3)
    call var%setAttribute("long_name", "Total discharge inputs at t-1 and t at level 11")

    do ii = 1, size(dummy_d3, 3)
      dummy_d3(:, :, ii) = unpack(L11_qTR(s11 : e11, ii), mask11, nodata_dp)
    end do
    var = nc%setVariable("L11_qTR", "f64", (/rows11, cols11, it11/))
    call var%setFillValue(nodata_dp)
    call var%setData(dummy_d3)
    call var%setAttribute("long_name", "Routed outflow leaving a node at level 11")

    var = nc%setVariable("L11_K", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_K(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "kappa: Muskingum travel time parameter at level 11")

    var = nc%setVariable("L11_xi", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_xi(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "xi: Muskingum diffusion parameter at level 11")

    var = nc%setVariable("L11_C1", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_C1(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Routing parameter C1=f(K,xi, DT) (Chow, 25-41) at level 11")

    var = nc%setVariable("L11_C2", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_C2(s11 : e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "Routing parameter C2=f(K,xi, DT) (Chow, 25-41) at level 11")

    deallocate(dummy_d3)
    allocate(dummy_d3(nrows11, ncols11, nLCoverScene))
    do ii = 1, size(dummy_d3, 3)
      dummy_d3(:, :, ii) = unpack(L11_nLinkFracFPimp(s11 : e11, ii), mask11, nodata_dp)
    end do
    var = nc%setVariable("L11_nLinkFracFPimp", "f64", (/rows11, cols11, lcscenes/))
    call var%setFillValue(nodata_dp)
    call var%setData(dummy_d3)
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

    var = nc%setVariable("L11_TSrout", "i32", (/nts/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_TSrout(iBasin))
    call var%setAttribute("long_name", "routing resolution at Level 11")
    call var%setAttribute("units", "s")
    
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

    var = nc%setVariable("L11_fAcc", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_fAcc(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "flow accumulation at Level 11")

    var = nc%setVariable("L11_fDir", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_fDir(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "flow Direction at Level 11")

    var = nc%setVariable("L11_rowOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_rowOut(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Grid vertical location of the Outlet at Level 11")

    var = nc%setVariable("L11_colOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_colOut(s11 : e11), mask11, nodata_i4))
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
    call var%setData(L11_fromN(s11 : e11))
    call var%setAttribute("long_name", "From Node")

    var = nc%setVariable("L11_toN", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_toN(s11 : e11))
    call var%setAttribute("long_name", "To Node")

    var = nc%setVariable("L11_rOrder", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_rOrder(s11 : e11))
    call var%setAttribute("long_name", "Network routing order at Level 11")

    var = nc%setVariable("L11_label", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_label(s11 : e11))
    call var%setAttribute("long_name", "Label Id [0='', 1=HeadWater, 2=Sink] at Level 11")

    var = nc%setVariable("L11_sink", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(merge(1_i4, 0_i4, L11_sink(s11 : e11)))
    call var%setAttribute("long_name", ".true. if sink node reached at Level 11")

    var = nc%setVariable("L11_netPerm", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_netPerm(s11 : e11))
    call var%setAttribute("long_name", "Routing sequence (permutation of L11_rOrder) at Level 11")

    var = nc%setVariable("L11_fRow", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_fRow(s11 : e11))
    call var%setAttribute("long_name", "From row in L0 grid at Level 11")

    var = nc%setVariable("L11_fCol", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_fCol(s11 : e11))
    call var%setAttribute("long_name", "From col in L0 grid at Level 11")

    var = nc%setVariable("L11_tRow", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_tRow(s11 : e11))
    call var%setAttribute("long_name", "To row in L0 grid at Level 11")

    var = nc%setVariable("L11_tCol", "i32", (/links/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_tCol(s11 : e11))
    call var%setAttribute("long_name", "To Col in L0 grid at Level 11")

    var = nc%setVariable("L11_length", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_length(s11 : e11))
    call var%setAttribute("long_name", "Total length of river link [m]")

    var = nc%setVariable("L11_aFloodPlain", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_aFloodPlain(s11 : e11))
    call var%setAttribute("long_name", "Area of the flood plain [m2]")

    var = nc%setVariable("L11_slope", "f64", (/links/))
    call var%setFillValue(nodata_dp)
    call var%setData(L11_slope(s11 : e11))
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
    call var%setData(unpack(L11_L1_Id(s11 : e11), mask11, nodata_i4))
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
            (/nc%setDimension("Ngauges", size(basin_mrm(iBasin)%gaugeNodeList(:)))/) &
            )
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm(iBasin)%gaugeNodeList(:))
    call var%setAttribute("long_name", "cell ID of gauges")

    var = nc%setVariable("InflowGaugeNodeList", "i32", &
            (/nc%setDimension("nInflowGauges", size(basin_mrm(iBasin)%InflowGaugeNodeList(:)))/) &
            )
    call var%setFillValue(nodata_i4)
    call var%setData(basin_mrm(iBasin)%InflowGaugeNodeList(:))
    call var%setAttribute("long_name", "cell ID of gauges")

    if (processMatrix(8, 1) .eq. 3) then
      var = nc%setVariable("L11_celerity", "f64", (/links/))   ! celerity
      call var%setFillValue(nodata_dp)
      call var%setData(L11_celerity(s11:e11))
      call var%setAttribute("long_name", "celerity at Level 11")
    end if

    ! free dummy variables
    deallocate(dummy_d3)
    call nc%close()

  end subroutine mrm_write_restart

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_read_restart_states

  !    PURPOSE
  !>       \brief read routing states

  !>       \details This subroutine reads the routing states from
  !>       mRM_states_<basin_id>.nc that has to be in the given
  !>       path directory. This subroutine has to be called directly
  !>       each time the mHM_eval or mRM_eval is called such that the
  !>       the states are always the same at the first simulation time
  !>       step, crucial for optimization.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin"    number of basin
  !>       \param[in] "character(256) :: InPath" Input Path including trailing slash

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Sep 2015

  ! Modifications:
  ! David Schaefer Mar 2016 - mo_netcdf
  ! Stephan Thober May 2016 - split L0_OutletCoord into L0_rowOutlet & L0_colOutlet because multiple outlets could exist 
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_restart_states(iBasin, InPath)

    use mo_common_variables, only : nLCoverScene
    use mo_mrm_constants, only : nRoutingStates
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, L11_Qmod, L11_Qout, L11_nLinkFracFPimp, L11_qTIN, L11_qTR, &
                                        L11_xi, level11
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of basin
    integer(i4), intent(in) :: iBasin

    ! Input Path including trailing slash
    character(256), intent(in) :: InPath

    integer(i4) :: ii

    ! start index at level 11
    integer(i4) :: s11

    ! end index at level 11
    integer(i4) :: e11

    ! mask at level 11
    logical, dimension(:, :), allocatable :: mask11

    ! dummy, 2 dimension
    real(dp), dimension(:, :), allocatable :: dummyD2

    ! dummy, 3 dimension
    real(dp), dimension(:, :, :), allocatable :: dummyD3

    character(256) :: fname

    type(NcDataset) :: nc

    type(NcVariable) :: var


    ! set file name
    fname = trim(InPath) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'

    ! get basin info at L11 including ncells, start, end, and mask
    s11 = level11(iBasin)%iStart
    e11 = level11(iBasin)%iEnd
    mask11 = level11(iBasin)%mask

    nc = NcDataset(fname, "r")

    ! simulated discharge at each node
    var = nc%getVariable("L11_Qmod")
    call var%getData(dummyD2)
    L11_Qmod(s11 : e11) = pack(dummyD2, mask11)

    ! Total outflow from cells L11 at time tt
    var = nc%getVariable("L11_qOUT")
    call var%getData(dummyD2)
    L11_qOUT(s11 : e11) = pack(dummyD2, mask11)

    ! Total discharge inputs at t-1 and t
    var = nc%getVariable("L11_qTIN")
    call var%getData(dummyD3)
    do ii = 1, nRoutingStates
      L11_qTIN(s11 : e11, ii) = pack(dummyD3(:, :, ii), mask11)
    end do

    !  Routed outflow leaving a node
    var = nc%getVariable("L11_qTR")
    call var%getData(dummyD3)
    do ii = 1, nRoutingStates
      L11_qTR(s11 : e11, ii) = pack(dummyD3(:, :, ii), mask11)
    end do

    ! kappa: Muskingum travel time parameter.
    var = nc%getVariable("L11_K")
    call var%getData(dummyD2)
    L11_K(s11 : e11) = pack(dummyD2, mask11)

    !  xi:    Muskingum diffusion parameter
    var = nc%getVariable("L11_xi")
    call var%getData(dummyD2)
    L11_xi(s11 : e11) = pack(dummyD2, mask11)

    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    var = nc%getVariable("L11_C1")
    call var%getData(dummyD2)
    L11_C1(s11 : e11) = pack(dummyD2, mask11)

    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    var = nc%getVariable("L11_C2")
    call var%getData(dummyD2)
    L11_C2(s11 : e11) = pack(dummyD2, mask11)

    ! Fraction of the flood plain with impervious cover
    var = nc%getVariable("L11_nLinkFracFPimp")
    deallocate(dummyD3)
    call var%getData(dummyD3)
    do ii = 1, nLCoverScene
      L11_nLinkFracFPimp(s11 : e11, ii) = pack(dummyD3(:, :, ii), mask11)
    end do

    ! free memory
    deallocate(dummyD2, dummyD3)

  end subroutine mrm_read_restart_states

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_read_restart_config

  !    PURPOSE
  !>       \brief reads Level 11 configuration from a restart directory

  !>       \details read Level 11 configuration variables from a given restart
  !>       directory and initializes all Level 11 configuration variables,
  !>       that are initialized in L11_variable_init,
  !>       contained in module mo_startup.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin"    number of basin
  !>       \param[in] "character(256) :: InPath" Input Path including trailing slash

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink  Apr 2014 - added inflow gauge
  ! Stephan Thober Aug 2015 - adapted for mRM
  ! Stephan Thober Sep 2015 - added all read restart commands in this subroutine
  ! Stephan Thober Sep 2015 - added L11_areaCell, L1_ID and L1_L11_Id for routing resolution higher than hydrology resolution
  ! David Schaefer Mar 2016 - mo_netcdf
  ! Stephan Thober Nov 2016 - added L11_TSrout, ProcessMatrix
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_restart_config(iBasin, InPath)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, nBasins, processMatrix
    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_mrm_global_variables, only : L11_L1_Id, L11_TSrout, L11_aFloodPlain, L11_colOut, L11_fCol, &
                                        L11_fDir, L11_fRow, L11_fromN, L11_label, L11_length, L11_nOutlets, L11_netPerm, &
                                        L11_rOrder, L11_rowOut, L11_sink, L11_slope, L11_tCol, L11_tRow, L11_toN, &
                                        L1_L11_Id, basin_mrm, level11
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of basin
    integer(i4), intent(in) :: iBasin

    ! Input Path including trailing slash
    character(256), intent(in) :: InPath

    character(256) :: fname

    ! Mask at Level 1
    logical, allocatable, dimension(:, :) :: mask1

    ! Mask at Level 11
    logical, allocatable, dimension(:, :) :: mask11

    ! dummy, 1 dimension I4
    integer(i4), allocatable, dimension(:) :: dummyI1

    ! dummy, 2 dimension I4
    integer(i4), allocatable, dimension(:, :) :: dummyI2

    ! dummy, 1 dimension DP
    real(dp), allocatable, dimension(:) :: dummyD1

    type(NcDataset) :: nc

    type(NcVariable) :: var


    ! set file name
    fname = trim(InPath) // 'mRM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'
    call message('        Reading mRM restart file:  ', trim(adjustl(Fname)), ' ...')

    ! get basin info at L11 mask
    mask1 = level1(iBasin)%mask
    mask11 = level11(iBasin)%mask

    nc = NcDataset(fname, "r")

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read Process Matrix for check <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    var = nc%getVariable("ProcessMatrix")
    allocate(dummyI1(size(processMatrix, dim = 1)))
    call var%getData(dummyI1)
    if (dummyI1(8) .ne. processMatrix(8, 1)) then
      call message('***ERROR: process description for routing')
      call message('***ERROR: given in restart file does not match')
      call message('***ERROR: that in namelist')
      call message('***ERROR: restart file value:. ' // num2str(dummyI1(8), '(i2)'))
      call message('***ERROR: namelist value:..... ' // num2str(processMatrix(8, 1), '(i2)'))
      stop 'ERROR: mrm_read_restart_config'
    end if
    deallocate(dummyI1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    var = nc%getVariable("L1_L11_Id")
    call var%getData(dummyI2)
    call append(L1_L11_Id, pack(dummyI2, mask1))

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L11 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L11_TSrout
    if (iBasin .eq. 1) then
      allocate(L11_TSrout(nBasins))
      L11_tsRout = nodata_dp
    end if
    var = nc%getVariable("L11_TSrout")
    call var%getData(L11_TSrout(iBasin))

    ! L11 data sets
    ! Mapping of L1 Ids on L1
    var = nc%getVariable("L11_L1_Id")
    call var%getData(dummyI2)
    call append(L11_L1_Id, pack(dummyI2, mask11))

    ! Flow direction (standard notation)
    var = nc%getVariable("L11_fDir")
    call var%getData(dummyI2)
    call append(L11_fDir, pack(dummyI2, mask11))
    ! append Number of Outlets at Level 11 (where facc == 0 )
    call append(L11_nOutlets, count((dummyI2 .eq. 0_i4)))

    ! Grid vertical location of the Outlet
    var = nc%getVariable("L11_rowOut")
    call var%getData(dummyI2)
    call append(L11_rowOut, pack(dummyI2, mask11))

    ! Grid horizontal location  of the Outlet
    var = nc%getVariable("L11_colOut")
    call var%getData(dummyI2)
    call append(L11_colOut, pack(dummyI2, mask11))

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
    ! read gaugenodelist
    var = nc%getVariable("gaugeNodeList")
    call var%getData(dummyI1)
    basin_mrm(iBasin)%gaugeNodeList(:) = dummyI1

    ! read InflowGaugeNodelist
    if (basin_mrm(iBasin)%nInflowGauges > 0) then
      var = nc%getVariable("InflowGaugeNodeList")
      call var%getData(dummyI1)
      basin_mrm(iBasin)%InflowgaugeNodeList(:) = dummyI1
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

  end subroutine mrm_read_restart_config

end module mo_mrm_restart
