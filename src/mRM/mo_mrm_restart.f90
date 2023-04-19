!> \file mo_mrm_restart.f90
!> \brief \copybrief mo_mrm_restart
!> \details \copydetails mo_mrm_restart

!> \brief Restart routines
!> \details This module contains the subroutines for reading and writing routing related variables to file.
!> \authors Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_restart
  use mo_kind, only : i4, dp
  use mo_message, only : message, error_message

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
  !>       \param[in] "integer(i4) :: iDomain"                   number of domain
  !>       \param[in] "character(256), dimension(:) :: OutFile" list of Output paths per Domain

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
  ! Stephan Thober    Jun 2018 - including varying celerity functionality
  ! Stephan Thober    May 2019 - added L0 info required for Process 3

  subroutine mrm_write_restart(iDomain, domainID, OutFile)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_restart, only : write_grid_info
    use mo_common_variables, only : level0, level1, nLCoverScene, processMatrix, domainMeta, &
            LC_year_start, LC_year_end
    use mo_common_constants, only : landCoverPeriodsVarName
    use mo_mrm_constants, only : nRoutingStates
    use mo_mpr_global_variables, only : L0_slope
    use mo_mrm_global_variables, only : L0_fdir, L0_fAcc, L0_streamnet, &
                                        L1_L11_Id, &
                                        L11_C1, L11_C2, L11_K, L11_L1_Id, L11_Qmod, &
                                        L11_TSrout, L11_aFloodPlain, L11_colOut, L11_colOut, L11_fCol, L11_fDir, &
                                        L11_fAcc, L11_fRow, L11_fromN, L11_label, L11_length, L11_nLinkFracFPimp, &
                                        L11_netPerm, L11_qOUT, L11_qTIN, L11_qTR, L11_rOrder, L11_rowOut, L11_rowOut, &
                                        L11_sink, L11_slope, L11_tCol, L11_tRow, L11_toN, L11_xi, L11_celerity, &
                                        level11, domain_mrm
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of domain
    integer(i4), intent(in) :: iDomain

    ! domain
    integer(i4), intent(in) :: domainID

    ! list of Output paths per Domain
    character(256), dimension(:), intent(in) :: OutFile

    character(256) :: Fname

    integer(i4) :: ii

    ! number of outlets at Level 0
    integer(i4) :: Noutlet

    ! start index at level 0
    integer(i4) :: s0

    ! end index at level 0
    integer(i4) :: e0

    ! mask at level 0
    logical, dimension(:, :), allocatable :: mask0

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

    ! dummy variable for writing L11_sink
    integer(i4), allocatable :: dummy(:)

    ! mask at level 11
    logical, dimension(:, :), allocatable :: mask11

    ! dummy variable
    real(dp), dimension(:, :, :), allocatable :: dummy_d3
    real(dp), dimension(:), allocatable :: dummy_d1

    type(NcDataset) :: nc

    type(NcDimension) :: rows0, cols0, rows1, cols1, rows11, cols11, it11, lcscenes, nout

    type(NcDimension) :: links, nts, nproc

    type(NcVariable) :: var


    ! get Level1 and Level11 information about the Domain
    noutlet = domain_mrm(domainMeta%L0DataFrom(iDomain))%L0_noutlet
    s0 = level0(domainMeta%L0DataFrom(iDomain))%iStart
    e0 = level0(domainMeta%L0DataFrom(iDomain))%iEnd
    mask0 = level0(domainMeta%L0DataFrom(iDomain))%mask
    s1 = level1(iDomain)%iStart
    e1 = level1(iDomain)%iEnd
    mask1 = level1(iDomain)%mask
    s11 = level11(iDomain)%iStart
    e11 = level11(iDomain)%iEnd
    mask11 = level11(iDomain)%mask
    ncols11 = level11(iDomain)%ncols
    nrows11 = level11(iDomain)%nrows
    allocate(dummy_d3(nrows11, ncols11, nRoutingStates))

    ! set restart file name
    Fname = trim(OutFile(iDomain))

    call message('    Writing mRM restart file to ' // trim(Fname) // ' ...')

    nc = NcDataset(fname, "w")

    call write_grid_info(level0(domainMeta%L0DataFrom(iDomain)), "0", nc)
    call write_grid_info(level1(iDomain), "1", nc)
    call write_grid_info(level11(iDomain), "11", nc)

    nout = nc%setDimension("Noutlet", Noutlet)
    rows0 = nc%getDimension("nrows0")
    cols0 = nc%getDimension("ncols0")
    rows1 = nc%getDimension("nrows1")
    cols1 = nc%getDimension("ncols1")
    rows11 = nc%getDimension("nrows11")
    cols11 = nc%getDimension("ncols11")

    it11 = nc%setDimension("nIT", nRoutingStates)
    links = nc%setDimension("nLinks", size(L11_length(s11 : e11)))
    nts = nc%setDimension("TS", 1)
    nproc = nc%setDimension("Nprocesses", size(processMatrix, dim = 1))
    allocate(dummy_d1(nLCoverScene+1))
    dummy_d1(1:nLCoverScene) = LC_year_start(:)
    ! this is done because bounds are always stored as real so e.g.
    ! 1981-1990,1991-2000 is thus saved as 1981.0-1991.0,1991.0-2001.0
    ! it is translated back into ints correctly during reading
    dummy_d1(nLCoverScene+1) = LC_year_end(nLCoverScene) + 1
    lcscenes = nc%setCoordinate(trim(landCoverPeriodsVarName), nLCoverScene, dummy_d1, 0_i4)
    deallocate(dummy_d1)


    ! add processMatrix
    var = nc%setVariable("ProcessMatrix", "i32", (/nproc/))
    call var%setFillValue(nodata_i4)
    call var%setData(processMatrix(:, 1))
    call var%setAttribute("long_name", "Process Matrix")

    ! add L0 variables if processmatrix is equal to 3
    if (processMatrix(8, 1) .eq. 3_i4) then
      ! add L0_fdir, L0_fAcc, L0_slope, L0_streamnet
      var = nc%setVariable("L0_fDir", "i32", (/rows0, cols0/))
      call var%setFillValue(nodata_i4)
      call var%setData(unpack(L0_fdir(s0:e0), mask0, nodata_i4))
      call var%setAttribute("long_name", "flow direction at level 0")

      var = nc%setVariable("L0_fAcc", "i32", (/rows0, cols0/))
      call var%setFillValue(nodata_i4)
      call var%setData(unpack(L0_fAcc(s0:e0), mask0, nodata_i4))
      call var%setAttribute("long_name", "flow accumulation at level 0")

      var = nc%setVariable("L0_slope", "f64", (/rows0, cols0/))
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L0_slope(s0:e0), mask0, nodata_dp))
      call var%setAttribute("long_name", "slope at level 0")

      var = nc%setVariable("L0_streamnet", "i32", (/rows0, cols0/))
      call var%setFillValue(nodata_i4)
      call var%setData(unpack(L0_streamnet(s0:e0), mask0, nodata_i4))
      call var%setAttribute("long_name", "streamnet at level 0")
    end if

    var = nc%setVariable("L1_Id", "i32", (/rows1, cols1/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(level1(iDomain)%Id(1:e1-s1+1), mask1, nodata_i4))
    call var%setAttribute("long_name", "cell IDs at level 1")

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
    var = nc%setVariable("L11_domain_Mask", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(merge(1_i4, 0_i4,  mask11))
    call var%setAttribute("long_name", "Mask at Level 11")

    var = nc%setVariable("L11_TSrout", "i32", (/nts/))
    call var%setFillValue(nodata_i4)
    call var%setData(L11_TSrout(iDomain))
    call var%setAttribute("long_name", "routing resolution at Level 11")
    call var%setAttribute("units", "s")

    var = nc%setVariable("L11_Id", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(level11(iDomain)%Id(1:e11-s11+1), mask11, nodata_i4))
    call var%setAttribute("long_name", "cell Ids at Level 11")

    var = nc%setVariable("L11_fDir", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_fDir(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "flow Direction at Level 11")

    var = nc%setVariable("L11_fAcc", "f64", (/rows11, cols11/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(L11_fAcc(s11:e11), mask11, nodata_dp))
    call var%setAttribute("long_name", "flow accumulation at Level 11")

    var = nc%setVariable("L11_rowOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_rowOut(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Grid vertical location of the Outlet at Level 11")

    var = nc%setVariable("L11_colOut", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_colOut(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Grid horizontal location of the Outlet at Level 11")

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
    allocate(dummy(e11 - s11 + 1))
    dummy = 0_i4
    where(L11_sink(s11 : e11))
      dummy = 1_i4
    end where
    ! call var%setData(merge(1_i4, 0_i4, L11_sink(s11 : e11)))
    call var%setData(dummy)
    deallocate(dummy)
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

    var = nc%setVariable("L11_L1_Id", "i32", (/rows11, cols11/))
    call var%setFillValue(nodata_i4)
    call var%setData(unpack(L11_L1_Id(s11 : e11), mask11, nodata_i4))
    call var%setAttribute("long_name", "Mapping of L1 Id on L11")

    var = nc%setVariable("gaugeNodeList", "i32", &
            (/nc%setDimension("Ngauges", size(domain_mrm(iDomain)%gaugeNodeList(:)))/) &
            )
    call var%setFillValue(nodata_i4)
    call var%setData(domain_mrm(iDomain)%gaugeNodeList(:))
    call var%setAttribute("long_name", "cell ID of gauges")

    var = nc%setVariable("InflowGaugeNodeList", "i32", &
            (/nc%setDimension("nInflowGauges", size(domain_mrm(iDomain)%InflowGaugeNodeList(:)))/) &
            )
    call var%setFillValue(nodata_i4)
    call var%setData(domain_mrm(iDomain)%InflowGaugeNodeList(:))
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
  !>       mRM_states_<domain_id>.nc that has to be in the given
  !>       path directory. This subroutine has to be called directly
  !>       each time the mHM_eval or mRM_eval is called such that the
  !>       the states are always the same at the first simulation time
  !>       step, crucial for optimization.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"    number of domains
  !>       \param[in] "character(256) :: InFile" Input Path including trailing slash

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Sep 2015

  ! Modifications:
  ! David Schaefer Mar 2016 - mo_netcdf
  ! Stephan Thober May 2016 - split L0_OutletCoord into L0_rowOutlet & L0_colOutlet because multiple outlets could exist
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_restart_states(iDomain, domainID, InFile)

    use mo_common_variables, only : nLCoverScene
    use mo_mrm_constants, only : nRoutingStates
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, L11_Qmod, L11_Qout, L11_nLinkFracFPimp, L11_qTIN, L11_qTR, &
                                        L11_xi, level11
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of Domain
    integer(i4), intent(in) :: iDomain

    integer(i4), intent(in) :: domainID

    ! Input Path including trailing slash
    character(256), intent(in) :: InFile

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

    !TODO-RIV-TEMP: read/write restart for riv-temp process

    ! set file name
    fname = trim(InFile)

    ! get Domain info at L11 including ncells, start, end, and mask
    s11 = level11(iDomain)%iStart
    e11 = level11(iDomain)%iEnd
    mask11 = level11(iDomain)%mask

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
  !>       \param[in] "integer(i4) :: iDomain"    number of Domain
  !>       \param[in] "character(256) :: InFile" Input Path including trailing slash

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
  ! Stephan Thober May 2019 - added L0 info required for Process 3

  subroutine mrm_read_restart_config(iDomain, domainID, InFile)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level0, level1, domainMeta, processMatrix, domainMeta
    use mo_kind, only : dp, i4
    use mo_mpr_global_variables, only : L0_slope
    use mo_mrm_global_variables, only : L0_fdir, L0_fAcc, L0_streamnet, &
                                        L11_L1_Id, L11_TSrout, L11_aFloodPlain, L11_colOut, L11_fCol, &
                                        L11_fDir, L11_fAcc, L11_fRow, L11_fromN, L11_label, L11_length, L11_nOutlets, L11_netPerm, &
                                        L11_rOrder, L11_rowOut, L11_sink, L11_slope, L11_tCol, L11_tRow, L11_toN, &
                                        L1_L11_Id, domain_mrm, level11
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! number of Domain
    integer(i4), intent(in) :: iDomain

    ! domain
    integer(i4), intent(in) :: domainID

    ! Input Path including trailing slash
    character(256), intent(in) :: InFile

    character(256) :: fname

    ! Mask at Level 0
    logical, allocatable, dimension(:, :) :: mask0

    ! Mask at Level 1
    logical, allocatable, dimension(:, :) :: mask1

    ! Mask at Level 11
    logical, allocatable, dimension(:, :) :: mask11

    ! dummy, 1 dimension I4
    integer(i4), allocatable, dimension(:) :: dummyI1

    ! dummy, 2 dimension I4
    integer(i4), allocatable, dimension(:, :) :: dummyI2

    ! dummy, DP
    real(dp), allocatable, dimension(:) :: dummyD1
    real(dp), allocatable, dimension(:, :) :: dummyD2
    real(dp), allocatable :: dummyD0

    type(NcDataset) :: nc

    type(NcVariable) :: var


    ! set file name
    fname = trim(InFile)
    call message('        Reading mRM restart file:  ', trim(adjustl(Fname)), ' ...')

    ! get Domain info at L11 mask
    mask0 = level0(domainMeta%L0DataFrom(iDomain))%mask
    mask1 = level1(iDomain)%mask
    mask11 = level11(iDomain)%mask

    nc = NcDataset(fname, "r")

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read Process Matrix for check <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    var = nc%getVariable("ProcessMatrix")
    allocate(dummyI1(size(processMatrix, dim = 1)))
    call var%getData(dummyI1)
    if (dummyI1(8) .ne. processMatrix(8, 1)) then
      call error_message('***ERROR: process description for routing', raise=.false.)
      call error_message('***ERROR: given in restart file does not match', raise=.false.)
      call error_message('***ERROR: that in namelist', raise=.false.)
      call error_message('***ERROR: restart file value:. ' // num2str(dummyI1(8), '(i2)'), raise=.false.)
      call error_message('***ERROR: namelist value:..... ' // num2str(processMatrix(8, 1), '(i2)'), raise=.false.)
      call error_message('ERROR: mrm_read_restart_config')
    end if
    deallocate(dummyI1)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L0 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (processMatrix(8, 1) .eq. 3_i4) then
      ! add L0_fdir, L0_fAcc, L0_slope, L0_streamnet
      var = nc%getVariable("L0_fDir")
      call var%getData(dummyI2)
      call append(L0_fdir, pack(dummyI2, mask0))

      var = nc%getVariable("L0_fAcc")
      call var%getData(dummyI2)
      call append(L0_fAcc, pack(dummyI2, mask0))

      var = nc%getVariable("L0_slope")
      call var%getData(dummyD2)
      call append(L0_Slope, pack(dummyD2, mask0))

      var = nc%getVariable("L0_streamnet")
      call var%getData(dummyI2)
      call append(L0_streamnet, pack(dummyI2, mask0))
    end if

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
    if (iDomain .eq. 1) then
      allocate(L11_TSrout(domainMeta%nDomains))
      L11_tsRout = nodata_dp
    end if
    var = nc%getVariable("L11_TSrout")
    call var%getData(dummyD0)
    L11_TSrout(iDomain) = dummyD0

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

    ! Flow accumulation
    var = nc%getVariable("L11_fAcc")
    call var%getData(dummyD2)
    call append(L11_fAcc, pack(dummyD2, mask11))

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
    domain_mrm(iDomain)%gaugeNodeList(:) = dummyI1

    ! read InflowGaugeNodelist
    if (domain_mrm(iDomain)%nInflowGauges > 0) then
      var = nc%getVariable("InflowGaugeNodeList")
      call var%getData(dummyI1)
      domain_mrm(iDomain)%InflowgaugeNodeList(:) = dummyI1
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
