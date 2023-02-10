!> \file mo_mrm_routing.f90
!> \brief   \copybrief mo_mrm_routing
!> \details \copydetails mo_mrm_routing

!> \brief Performs runoff routing for mHM at level L11.
!> \details This module performs flood routing at a given time step through the stream network at level L11 to the sink cell.
!! The Muskingum flood routing algorithm is used.
!> \changelog
!! - Stephan Thober Aug 2015
!!   - adapted to mRM
!! - Sebastian Mueller Jun 2020
!!   - outsourcing helper functions
!> \authors Luis Samaniego
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
MODULE mo_mrm_routing

  ! This module performs runoff flood routing for mHM.

  ! Written Luis Samaniego, Dec 2012

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mRM_routing

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mRM_routing

  !    PURPOSE
  !>       \brief route water given runoff

  !>       \details This routine first performs mpr for the routing variables
  !>       if required, then accumulates the runoff to the routing resolution
  !>       and eventually routes the water in a third step. The last step is
  !>       repeated multiple times if the routing timestep is smaller than
  !>       the timestep of the hydrological timestep

  !    INTENT(IN)
  !>       \param[in] "logical :: read_states"                            whether states are derived from restart file
  !>       \param[in] "integer(i4) :: processCase"                        Process switch for routing
  !>       \param[in] "real(dp), dimension(:) :: global_routing_param"    routing parameters
  !>       \param[in] "real(dp), dimension(:) :: L1_total_runoff"         total runoff from L1 grid cells
  !>       \param[in] "real(dp), dimension(:) :: L1_areaCell"             L1 cell area
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id"            L1 cell ids on L11
  !>       \param[in] "real(dp), dimension(:) :: L11_areaCell"            L11 cell area
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id"            L11 cell ids on L1
  !>       \param[in] "integer(i4), dimension(:) :: L11_netPerm"          L11 routing order
  !>       \param[in] "integer(i4), dimension(:) :: L11_fromN"            L11 source grid cell order
  !>       \param[in] "integer(i4), dimension(:) :: L11_toN"              L11 target grid cell order
  !>       \param[in] "integer(i4) :: L11_nOutlets"                       L11 number of outlets/sinks
  !>       \param[in] "integer(i4) :: timestep"                           simulation timestep in [h]
  !>       \param[in] "real(dp) :: tsRoutFactor"                          factor between routing timestep and
  !>       hydrological timestep
  !>       \param[in] "integer(i4) :: nNodes"                             number of nodes
  !>       \param[in] "integer(i4) :: nInflowGauges"                      number of inflow gauges
  !>       \param[in] "integer(i4), dimension(:) :: InflowGaugeIndexList" index list of inflow gauges
  !>       \param[in] "logical, dimension(:) :: InflowGaugeHeadwater"     flag for headwater cell of inflow gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowGaugeNodeList"  gauge node list at L11
  !>       \param[in] "real(dp), dimension(:) :: InflowDischarge"         inflowing discharge at discharge gauge at
  !>       current day
  !>       \param[in] "integer(i4) :: nGauges"                            number of recording gauges
  !>       \param[in] "integer(i4), dimension(:) :: gaugeIndexList"       index list for outflow gauges
  !>       \param[in] "integer(i4), dimension(:) :: gaugeNodeList"        gauge node list at L11
  !>       \param[in] "logical :: map_flag"                               flag indicating whether routing resolution
  !>       iscoarser than hydrologic resolution
  !>       \param[in] "real(dp), dimension(:) :: L11_length"              L11 link length
  !>       \param[in] "real(dp), dimension(:) :: L11_slope"               L11 slope
  !>       \param[in] "real(dp), dimension(:) :: L11_FracFPimp"           L11 fraction of flood plain with impervios
  !>       cover

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:) :: L11_C1"         L11 muskingum parameter 1
  !>       \param[inout] "real(dp), dimension(:) :: L11_C2"         L11 muskingum parameter 2
  !>       \param[inout] "real(dp), dimension(:) :: L11_qOut"       total runoff from L11 grid cells
  !>       \param[inout] "real(dp), dimension(:, :) :: L11_qTIN"    L11 inflow to the reach
  !>       \param[inout] "real(dp), dimension(:, :) :: L11_qTR"     L11 routed outflow
  !>       \param[inout] "real(dp), dimension(:) :: L11_qMod"       modelled discharge at each grid cell
  !>       \param[inout] "real(dp), dimension(:) :: GaugeDischarge" modelled discharge at each gauge

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober Sep 2015 - using arguments instead of global variables
  ! Stephan Thober Sep 2015 - added variables for routing resolution higher than hydrologic resolution
  ! Stephan Thober May 2016 - added check whether gauge is actually inside modelling domain before copying simulated runoff
  ! Stephan Thober Nov 2016 - implemented second routing process i.e. adaptive timestep
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mRM_routing( &
    read_states, processCase, global_routing_param, L1_total_runoff, L1_areaCell, L1_L11_Id, &
    L11_areaCell, L11_L1_Id, L11_netPerm, L11_fromN, L11_toN, L11_nOutlets, timestep, tsRoutFactor, &
    nNodes, nInflowGauges, InflowGaugeIndexList, InflowGaugeHeadwater, InflowGaugeNodeList, &
    InflowDischarge, nGauges, gaugeIndexList, gaugeNodeList, map_flag, L11_length, L11_slope, &
    L11_FracFPimp, L11_C1, L11_C2, L11_qOut, L11_qTIN, L11_qTR, L11_qMod, GaugeDischarge &
  )

    use mo_constants, only : T0_dp
    use mo_mrm_global_variables, only : riv_temp_pcs, is_start
    use mo_mrm_mpr, only : reg_rout
    use mo_mrm_pre_routing, only : L11_runoff_acc, add_inflow
    ! use mo_mrm_riv_temp_class, only : riv_temp_type

    implicit none

    ! whether states are derived from restart file
    logical, intent(in) :: read_states
    ! Process switch for routing
    integer(i4), intent(in) :: processCase
    ! routing parameters
    real(dp), dimension(:), intent(in) :: global_routing_param
    ! total runoff from L1 grid cells
    real(dp), dimension(:), intent(in) :: L1_total_runoff
    ! L1 cell area
    real(dp), dimension(:), intent(in) :: L1_areaCell
    ! L1 cell ids on L11
    integer(i4), dimension(:), intent(in) :: L1_L11_Id
    ! L11 cell area
    real(dp), dimension(:), intent(in) :: L11_areaCell
    ! L11 cell ids on L1
    integer(i4), dimension(:), intent(in) :: L11_L1_Id
    ! L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_netPerm
    ! L11 source grid cell order
    integer(i4), dimension(:), intent(in) :: L11_fromN
    ! L11 target grid cell order
    integer(i4), dimension(:), intent(in) :: L11_toN
    ! L11 number of outlets/sinks
    integer(i4), intent(in) :: L11_nOutlets
    ! simulation timestep in [h]
    integer(i4), intent(in) :: timestep
    ! factor between routing timestep and hydrological timestep
    real(dp), intent(in) :: tsRoutFactor
    ! number of nodes
    integer(i4), intent(in) :: nNodes
    ! number of inflow gauges
    integer(i4), intent(in) :: nInflowGauges
    ! index list of inflow gauges
    integer(i4), dimension(:), intent(in) :: InflowGaugeIndexList
    ! flag for headwater cell of inflow gauge
    logical, dimension(:), intent(in) :: InflowGaugeHeadwater
    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: InflowGaugeNodeList
    ! inflowing discharge at discharge gauge at current day
    real(dp), dimension(:), intent(in) :: InflowDischarge
    ! number of recording gauges
    integer(i4), intent(in) :: nGauges
    ! index list for outflow gauges
    integer(i4), dimension(:), intent(in) :: gaugeIndexList
    ! gauge node list at L11
    integer(i4), dimension(:), intent(in) :: gaugeNodeList
    ! flag indicating whether routing resolution iscoarser than hydrologic resolution
    logical, intent(in) :: map_flag
    ! L11 link length
    real(dp), dimension(:), intent(in) :: L11_length
    ! L11 slope
    real(dp), dimension(:), intent(in) :: L11_slope
    ! L11 fraction of flood plain with impervios cover
    real(dp), dimension(:), intent(in) :: L11_FracFPimp
    ! L11 muskingum parameter 1
    real(dp), dimension(:), intent(inout) :: L11_C1
    ! L11 muskingum parameter 2
    real(dp), dimension(:), intent(inout) :: L11_C2
    ! total runoff from L11 grid cells
    real(dp), dimension(:), intent(inout) :: L11_qOut
    ! L11 inflow to the reach
    real(dp), dimension(:, :), intent(inout) :: L11_qTIN
    ! L11 routed outflow
    real(dp), dimension(:, :), intent(inout) :: L11_qTR
    ! modelled discharge at each grid cell
    real(dp), dimension(:), intent(inout) :: L11_qMod
    ! modelled discharge at each gauge
    real(dp), dimension(:), intent(inout) :: GaugeDischarge

    integer(i4) :: s11, e11 ! only for riv temp routing
    integer(i4) :: gg
    integer(i4) :: tt
    ! number of routing loops
    integer(i4) :: rout_loop
    ! variable for accumulation
    real(dp), dimension(size(L11_qMod, dim = 1)) :: L11_qAcc
    real(dp), dimension(:), allocatable :: L11_E_Acc

    if ( riv_temp_pcs%active ) then
      ! allocate accumulated temperature energy
      allocate(L11_E_Acc(size(L11_qMod, dim = 1)))
      L11_E_Acc = 0._dp ! init to zero
      ! get shortcuts for start end ending of current L11-domain
      s11 = riv_temp_pcs%s11
      e11 = riv_temp_pcs%e11
    end if

    if (is_start) is_start = .false.

    ! this is using the sealed fraction for determining the routing parameters
    ! MPR has already been done
    if (processCase .eq. 1_i4 .AND. (.not. read_states)) then
      ! for a single node model run
      if (nNodes .GT. 1) then
        call reg_rout(global_routing_param, &
                L11_length, L11_slope, L11_FracFPimp(: nNodes - L11_nOutlets), &
                real(timeStep, dp), L11_C1(: nNodes - L11_nOutlets), L11_C2(: nNodes - L11_nOutlets))
      end if
    end if

    ! =====================================================================
    ! NOW, EXECUTE ROUTING
    ! ====================================================================
    ! calculate number of routing loops
    rout_loop = max(1_i4, nint(1._dp / tsRoutFactor))

    ! runoff accumulation from L1 to L11 level
    call L11_runoff_acc( &
      L1_total_runoff, L1_areaCell, L1_L11_Id, &
      L11_areaCell, L11_L1_Id, timeStep, & ! Intent IN
      map_flag, & ! Intent IN
      L11_qOut & ! Intent OUT
    )
    ! add inflow
    call add_inflow( &
      nInflowGauges, &
      InflowGaugeIndexList, &
      InflowGaugeHeadwater, &
      InflowGaugeNodeList, &
      InflowDischarge, & ! Intent IN
      L11_qOUT & ! Intent INOUT
    )
    ! for a single node model run
    if(nNodes .GT. 1) then
      L11_qAcc = 0._dp
      ! routing multiple times if timestep is smaller than 1
      do tt = 1, rout_loop
        ! routing of water within river reaches
        call L11_routing( &
          nNodes, &
          nNodes - L11_nOutlets, &
          L11_netPerm, &
          L11_fromN, & ! Intent IN
          L11_toN, & ! Intent IN
          L11_C1, & ! Intent IN
          L11_C2, & ! Intent IN
          L11_qOut, & ! Intent IN
          nInflowGauges, & ! Intent IN
          InflowGaugeHeadwater, & ! Intent IN
          InflowGaugeNodeList, & ! Intent IN
          L11_qTIN, & ! Intent INOUT
          L11_qTR, & ! Intent INOUT
          L11_Qmod & ! Intent OUT
        )
        ! accumulate values of individual subtimesteps
        L11_qAcc = L11_qAcc + L11_qMod
        ! do the temperature routing
        if ( riv_temp_pcs%active ) then
          call riv_temp_pcs%L11_routing_E( &
            nNodes - L11_nOutlets, &
            L11_netPerm, &
            L11_fromN, & ! Intent IN
            L11_toN, & ! Intent IN
            L11_C1, & ! Intent IN
            L11_C2, & ! Intent IN
            nInflowGauges, & ! Intent IN
            InflowGaugeHeadwater, & ! Intent IN
            InflowGaugeNodeList, & ! Intent IN
            L11_qTR(:, 1), & ! Intent IN
            L11_Qmod & ! Intent IN
          )
        end if
      end do
      ! calculate mean over routing period (timestep)
      L11_qMod = L11_qAcc / real(rout_loop, dp)
    else
      L11_Qmod = L11_qOUT
      if ( riv_temp_pcs%active ) riv_temp_pcs%river_temp(s11 : e11) = &
        max(riv_temp_pcs%delta_T, riv_temp_pcs%netNode_E_out(s11 : e11) / L11_qOUT - T0_dp)
    end if

    !----------------------------------------------------------------------
    ! FOR STORING the optional arguments
    !
    ! FOR RUNOFF
    ! NOTE:: Node ID for a given gauging station is stored at gaugeindex's
    !        index in runoff. In consequence the gauges in runoff are
    !        ordered corresponing to gauge%Q(:,:)
    !----------------------------------------------------------------------
    do gg = 1, nGauges
      GaugeDischarge(gaugeIndexList(gg)) = L11_Qmod(gaugeNodeList(gg))
    end do

  end subroutine mRM_routing

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_routing

  !    PURPOSE
  !>       \brief Performs runoff routing for mHM at L11 upscaled network
  !>       (\ref fig_routing "Routing Network").
  !>       \details
  !>       Hydrograph routing is carried out with the Muskingum algorithm
  !>       \cite CMM1988.  This simplification of the St. Venant
  !>       equations is justified in mHM because the potential areas of
  !>       application of this model would hardly exhibit abruptly
  !>       changing hydrographs with supercritical flows.  The discharge
  !>       leaving the river reach located on cell \f$ i \f$ \f$
  !>       Q_{i}^{1}(t) \f$ at time step \f$ t \f$ can be determined by
  !>       \f[ Q_{i}^{1}(t) =  Q_{i}^{1}(t-1)
  !>       + c_{1} \left( Q_{i}^{0}(t-1) - Q_{i}^{1}(t-1) \right)
  !>       + c_{2} \left( Q_{i}^{0}(t)   - Q_{i}^{0}(t-1) \right) \f]
  !>       with
  !>       \f[  Q_{i}^{0}(t) = Q_{i'}(t) + Q_{i'}^{1}(t) \f]
  !>       \f[ c_{1}= \frac{\Delta t} { \kappa (1- \xi ) + \frac{\Delta t}{2} } \f]
  !>       \f[ c_{2}= \frac{ \frac{\Delta t}{2} - \kappa \xi} { \kappa (1- \xi)
  !>       + \frac{\Delta t}{2} } \f]
  !>       where
  !>       \f$ Q_{i}^{0} \f$ and \f$ Q_{i}^{1} \f$ denote the discharge
  !>       entering and leaving the river reach located on cell \f$ i \f$
  !>       respectively.
  !>       \f$ Q_{i'} \f$ is the contribution from the upstream cell \f$
  !>       i'\f$.
  !>       \f$ \kappa \f$ Muskingum travel time parameter.
  !>       \f$ \xi \f$ Muskingum attenuation parameter.
  !>       \f$ \Delta t \f$ time interval in hours.
  !>       \f$ t \f$ Time index for each \f$ \Delta t \f$ interval.
  !>       To improve performance, a routing sequence "netPerm" is
  !>       required. This permutation is determined in the mo_init_mrm
  !>       routine.

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: nNodes"                       number of network nodes = nCells1
  !>       \param[in] "integer(i4) :: nLinks"                       number of stream segment (reaches)
  !>       \param[in] "integer(i4), dimension(:) :: netPerm"        routing order of a given domain (permutation)
  !>       \param[in] "integer(i4), dimension(:) :: netLink_fromN"  from node
  !>       \param[in] "integer(i4), dimension(:) :: netLink_toN"    to node
  !>       \param[in] "real(dp), dimension(:) :: netLink_C1"        routing parameter  C1 (\cite CMM1988 p. 25-41)
  !>       \param[in] "real(dp), dimension(:) :: netLink_C2"        routing parameters C2 (id)
  !>       \param[in] "real(dp), dimension(:) :: netNode_qOUT"      Total outflow from cells (given domain) L11 at time
  !>       tt in [m3 s-1]
  !>       \param[in] "integer(i4) :: nInflowGauges"                [-]      number of inflow points
  !>       \param[in] "logical, dimension(:) :: InflowHeadwater"    [-]      if to consider headwater cells of inflow
  !>       gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList" [-]      L11 ID of inflow points

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :) :: netNode_qTIN" [m3 s-1] Total inputs at t-1 and t
  !>       \param[inout] "real(dp), dimension(:, :) :: netNode_qTR"  [m3 s-1] Transformed outflow leaving
  !>       node I (Muskingum)

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(nNodes) :: netNode_Qmod" [m3 s-1] Simulated routed discharge

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Feb 2008 - routing module (cells)
  ! Rohini Kumar   Aug 2011 - vector version of mHM-UFZ
  !                Nov 2011 - parallel version
  ! Luis Samaniego Jan 2013 - modularization, documentation
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_routing(nNodes, nLinks, netPerm, netLink_fromN, netLink_toN, netLink_C1, netLink_C2, netNode_qOUT, &
                        nInflowGauges, InflowHeadwater, InflowNodeList, netNode_qTIN, netNode_qTR, netNode_Qmod)
    implicit none

    ! number of network nodes = nCells1
    integer(i4), intent(in) :: nNodes
    ! number of stream segment (reaches)
    integer(i4), intent(in) :: nLinks
    ! routing order of a given domain (permutation)
    integer(i4), dimension(:), intent(in) :: netPerm
    ! from node
    integer(i4), dimension(:), intent(in) :: netLink_fromN
    ! to node
    integer(i4), dimension(:), intent(in) :: netLink_toN
    ! routing parameter  C1 (\cite CMM1988 p. 25-41)
    real(dp), dimension(:), intent(in) :: netLink_C1
    ! routing parameters C2 (id)
    real(dp), dimension(:), intent(in) :: netLink_C2
    ! Total outflow from cells (given domain) L11 at time tt in [m3 s-1]
    real(dp), dimension(:), intent(in) :: netNode_qOUT
    ! [-]      number of inflow points
    integer(i4), intent(in) :: nInflowGauges
    ! [-]      if to consider headwater cells of inflow gauge
    logical, dimension(:), intent(in) :: InflowHeadwater
    ! [-]      L11 ID of inflow points
    integer(i4), dimension(:), intent(in) :: InflowNodeList
    ! [m3 s-1] Total inputs at t-1 and t
    real(dp), dimension(:, :), intent(inout) :: netNode_qTIN
    ! [m3 s-1] Transformed outflow leaving node I (Muskingum)
    real(dp), dimension(:, :), intent(inout) :: netNode_qTR
    ! [m3 s-1] Simulated routed discharge
    real(dp), dimension(nNodes), intent(out) :: netNode_Qmod

    integer(i4) :: g, i, k, iNode, tNode
    ! current routing state (2)
    integer(i4), parameter :: IT = 2
    ! past routing state (1)
    integer(i4), parameter :: IT1 = 1

    ! Entry value for the auxiliary vectors
    !   netNode_qTIN(iNode,:)
    !   netNode_qTR(iNode,:)
    ! which store current and past states of
    ! incoming and outgoing of discharge at iNode
    !--------------------------------------------------------------------------
    !                             Muskingum Flood Routing
    !--------------------------------------------------------------------------
    ! initialize total input at point time IT in all nodes
    netNode_qTIN(:, IT) = 0.0_dp
    !--------------------------------------------------------------------------
    ! Links in sequential mode .... with single node
    !--------------------------------------------------------------------------
    ! ST - decent parallelization has to be done!!!
    !!$OMP parallel
    !!$OMP do private(g, i, inode, tnode)
    do k = 1, nLinks
      ! get LINK routing order -> i
      i = netPerm(k)
      iNode = netLink_fromN(i)
      tNode = netLink_toN(i)

      ! accumulate all inputs in iNode
      netNode_qTIN(iNode, IT) = netNode_qTIN(iNode, IT) + netNode_qOUT(iNode)

      ! routing iNode
      netNode_qTR(iNode, IT) = netNode_qTR(iNode, IT1)                               &
              + netLink_C1(i) * (netNode_qTIN(iNode, IT1) - netNode_qTR (iNode, IT1)) &
              + netLink_C2(i) * (netNode_qTIN(iNode, IT) - netNode_qTIN(iNode, IT1))

      ! check if the inflow from upstream cells should be deactivated
      if (nInflowGauges .GT. 0) then
        do g = 1, nInflowGauges
          ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
          if ((tNode == InflowNodeList(g)) .AND. (.NOT. InflowHeadwater(g))) netNode_qTR(iNode, IT) = 0.0_dp
        end do
      end if

      ! add routed water to downstream node
      netNode_qTIN(tNode, IT) = netNode_qTIN(tNode, IT) + netNode_qTR(iNode, IT)
    end do
    !!$OMP end do
    !!$OMP end parallel

    !--------------------------------------------------------------------------
    ! Accumulate all inputs in tNode (netNode_qOUT) ONLY for last link
    !--------------------------------------------------------------------------
    tNode = netLink_toN(netPerm(nLinks))
    netNode_qTIN(tNode, IT) = netNode_qTIN(tNode, IT) + netNode_qOUT(tNode)

    !--------------------------------------------------------------------------
    ! save modeled discharge at time step tt then shift flow storages
    ! (NOTE aggregation to daily values to be done outside)
    !--------------------------------------------------------------------------
    ! !!$OMP parallel
    ! store generated discharge
    netNode_Qmod(1 : nNodes) = netNode_qTIN(1 : nNodes, IT)
    ! backflow t-> t-1
    netNode_qTR(1 : nNodes, IT1) = netNode_qTR(1 : nNodes, IT)
    netNode_qTIN(1 : nNodes, IT1) = netNode_qTIN(1 : nNodes, IT)
    ! !!$OMP end parallel

  end subroutine L11_routing

END MODULE mo_mrm_routing
