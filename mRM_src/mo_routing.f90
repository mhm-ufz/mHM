!> \file mo_routing.f90

!> \brief Performs runoff routing for mHM at level L11.

!> \details This module performs flood routing at a given time step
!>          through the stream network at level L11 to the sink cell. 
!>          The Muskingum flood routing algorithm is used.

!> \author Luis Samaniego
!> \date Dec 2012

MODULE mo_routing

  ! This module performs runoff flood routing for mHM.

  ! Written Luis Samaniego, Dec 2012

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mRM_routing
  PUBLIC :: L11_routing      ! route water
  PUBLIC :: L11_runoff_acc   ! upscale runoff and add inflow gauge

  ! ------------------------------------------------------------------

CONTAINS

  subroutine mRM_routing(global_routing_param, iBasin, runoff, iTS, LCyearID, do_mpr_routing, discharge)
    use mo_global_variables_routing, only: &
         timeStep, &
         basin_mrm, &
         is_start, &
         L0_s, L0_e, & ! Level0 help variables
         L1_s, L1_e, & ! Level0 help variables
         L11_nNodes, L11_s, L11_e, & ! Level11 help variables
         L110_s, L110_e, & ! Level110 help variables
         L0_areaCell, &
         L0_floodPlain, & ! flood plains at L0 level
         L1_areaCell, &
         L11_aFloodPlain, & ! flood plains at L11 level
         L11_FracFPimp, & ! fraction of impervious layer at L11 scale
         L11_length, & ! link length
         L11_slope, &
         L11_C1, & ! first muskingum parameter
         L11_C2, & ! second muskigum parameter
         L1_L11_Id, &
         InflowGauge, &
         L11_qOUT, & ! routed runoff flowing out of L11 cell
         L11_netPerm, & ! routing order at L11
         L11_fromN, & ! link source at L11
         L11_toN, & ! link target at L11
         L11_qTIN, & ! inflow water into the reach at L11
         L11_qTR, & !
         L11_qMod ! final variable containing routed water
    use mo_init_mrm, only: L11_fraction_sealed_floodplain
    use mo_mpr_routing, only: reg_rout
    !ST The following dependency has to be removed
    use mo_global_variables, only: L0_LCover
         
    !
    implicit none
    ! input variables
    real(dp), dimension(5), intent(in) :: global_routing_param
    integer(i4), intent(in) :: iBasin
    real(dp), dimension(:), intent(in) :: runoff ! generated runoff for this timestep
    integer(i4), intent(in) :: LCyearID ! current land cover
    integer(i4), intent(in) :: iTS ! current timestep
    logical, optional, intent(in) :: do_mpr_routing
    ! output variables
    real(dp), optional, intent(out) :: discharge
    ! local variables
    logical :: do_mpr
    integer(i4) :: nNodes
    integer(i4) :: s11
    integer(i4) :: e11
    integer(i4) :: s110
    integer(i4) :: e110
    integer(i4) :: s1
    integer(i4) :: e1
    integer(i4) :: s0
    integer(i4) :: e0

    nNodes = L11_nNodes(iBasin)
    s11 = L11_s(iBasin)
    e11 = L11_e(iBasin)
    s110 = L110_s(iBasin)
    e110 = L110_e(iBasin)
    s1 = L1_s(iBasin)
    e1 = L1_e(iBasin)
    s0 = L0_s(iBasin)
    e0 = L0_e(iBasin)
    
    ! evaluate whether mpr should be executed
    do_mpr = .false.
    if (is_start) then
       do_mpr = .true.
       is_start = .false.
    end if
    if (present(do_mpr_routing)) do_mpr = do_mpr_routing

    ! execute mpr for routing
    if (do_mpr) then
       !-------------------------------------------------------------------
       ! estimate fraction of impervious cover in flood plains
       ! --> time independent variable: to be initalized every time
       !     with landcover change 
       ! --> Note: L11_fraction_sealed_floodplain routine is called
       !           only in case when routing process is ON
       !-------------------------------------------------------------------
       CALL L11_fraction_sealed_floodplain( nNodes-1, &
            L0_LCover(s0:e0, LCyearID), &
            ! L0_floodPlain(s110:e110),       &
            L0_floodPlain(s0:e0),       &
            L0_areaCell(s0:e0), &
            L11_aFloodPlain(s11:e11), &
            2, &
            L11_FracFPimp(s11:e11))
       ! for a single node model run
       if( nNodes .GT. 1) then
          call reg_rout( global_routing_param, &
               L11_length(s11:e11 - 1), L11_slope(s11:e11 - 1), L11_FracFPimp(s11:e11 - 1), &
               real(timeStep,dp), L11_C1(s11:e11 - 1), L11_C2(s11:e11 -1 ))
       end if 
    end if

    ! execute routing
    !-------------------------------------------------------------------
    ! routing at L11 level
    !-------------------------------------------------------------------
    ! runoff accumulation at L11 from L1 level
    call L11_runoff_acc(runoff, L1_areaCell(s1:e1), L1_L11_Id(s1:e1), timeStep, & ! Intent IN
         basin_mrm%nInflowGauges(iBasin), &
         basin_mrm%InflowGaugeIndexList(iBasin,:), &
         basin_mrm%InflowGaugeHeadwater(iBasin,:), &
         basin_mrm%InflowGaugeNodeList(iBasin,:), &
         InflowGauge%Q(iTS,:), & ! Intent IN
         L11_qOUT(s11:e11) )                                                                         ! Intent OUT
    ! for a single node model run
    if( nNodes .GT. 1) then
       ! routing of water within river reaches
       call L11_routing( nNodes, nNodes-1, &
            L11_netPerm(s11:e11), &
            L11_fromN(s11:e11), & ! Intent IN
            L11_toN(s11:e11), & ! Intent IN
            L11_C1(s11:e11), & ! Intent IN
            L11_C2(s11:e11), & ! Intent IN
            L11_qOUT(s11:e11), & ! Intent IN
            basin_mrm%nInflowGauges(iBasin), & ! Intent IN
            basin_mrm%InflowGaugeHeadwater(iBasin,:), & ! Intent IN
            basin_mrm%InflowGaugeNodeList(iBasin,:), & ! Intent IN
            L11_qTIN(s11:e11,:), & ! Intent INOUT
            L11_qTR(s11:e11,:), & ! Intent INOUT
            L11_Qmod(s11:e11) ) ! Intent OUT
    else
       L11_Qmod(s11:e11) = L11_qOUT(s11:e11) 
    end if

    ! set output variable discharge if given
    if (present(discharge)) discharge = L11_Qmod(s11)
    
  end subroutine mRM_routing

  ! ------------------------------------------------------------------

  !     NAME
  !         L11_routing

  !     PURPOSE 

  !>    \brief Performs runoff routing for mHM at L11 upscaled network
  !>     (\ref fig_routing "Routing Network").

  !>    \details

  !>    Hydrograph routing is carried out with the Muskingum algorithm
  !>    \cite CMM1988.  This simplification of the St. Venant
  !>    equations is justified in mHM because the potential areas of
  !>    application of this model would hardly exhibit abruptly
  !>    changing hydrographs with supercritical flows.  The discharge
  !>    leaving the river reach located on cell \f$ i \f$ \f$
  !>    Q_{i}^{1}(t) \f$ at time step \f$ t \f$ can be determined by

  !>     \f[ Q_{i}^{1}(t) =  Q_{i}^{1}(t-1) 
  !>             + c_{1} \left( Q_{i}^{0}(t-1) - Q_{i}^{1}(t-1) \right) 
  !>             + c_{2} \left( Q_{i}^{0}(t)   - Q_{i}^{0}(t-1) \right) \f]
  !>    with
  !>     \f[  Q_{i}^{0}(t) = Q_{i'}(t) + Q_{i'}^{1}(t) \f]
  !>     \f[ c_{1}= \frac{\Delta t} { \kappa (1- \xi ) + \frac{\Delta t}{2} } \f]
  !>     \f[ c_{2}= \frac{ \frac{\Delta t}{2} - \kappa \xi} { \kappa (1- \xi)
  !>        + \frac{\Delta t}{2} } \f]

  !>    where \n
  !>    \f$ Q_{i}^{0} \f$ and \f$ Q_{i}^{1} \f$ denote the discharge
  !>    entering and leaving the river reach located on cell \f$ i \f$
  !>    respectively. \n 
  !>    \f$ Q_{i'} \f$ is the contribution from the upstream cell \f$
  !>    i'\f$. \n 
  !>    \f$ \kappa \f$ Muskingum travel time parameter. \n
  !>    \f$ \xi \f$ Muskingum attenuation parameter. \n 
  !>    \f$ \Delta t \f$ time interval in hours. \n 
  !>    \f$ t \f$ Time index for each \f$ \Delta t \f$ interval. \n
  !>    To improve performance, a routing sequence "netPerm" is
  !>    required. This permutation is determined in the mo_init_mrm
  !>    routine.

  !     INTENT(IN)
  !>      \param[in] "integer(i4)  ::  nNodes"  number of network nodes = nCells1
  !>      \param[in] "integer(i4)  ::  nLinks"  number of stream segment (reaches)
  !>      \param[in] "integer(i4)  ::  netPerm" routing order of a given basin (permutation)
  !>      \param[in] "integer(i4)  ::  netLink_fromN" from node 
  !>      \param[in] "integer(i4)  ::  netLink_toN"   to node
  !>      \param[in] "real(dp)     ::  netLink_C1"    routing parameter  C1 (\cite CMM1988 p. 25-41)
  !>      \param[in] "real(dp)     ::  netLink_C2"    routing parameters C2 (id)
  !>      \param[in] "real(dp)     ::  netNode_qOUT"  Total outflow from cells (given basin) L11 at time tt in [m3 s-1]

  !     INTENT(INOUT)
  !>      \param[in,out] "real(dp) ::  netNode_qTIN(nNodes,2)"  Total discharge inputs at t-1 and t
  !>      \param[in,out] "real(dp) ::  netNode_qTR(nNodes,2)"   Routed outflow leaving a node (\cite TS2006 )
 
  !     INTENT(OUT)
  !>     \param[out]  "real(dp)    :: netNode_Qmod(nNodes)"   Simulated discharge [m3 s-1] 

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN


  !     RESTRICTIONS
  !>       \note A basin outlet should be identified.

  !     EXAMPLE
  !        none

  !     LITERATURE

  !     HISTORY
  !>        \author  Luis Samaniego
  !>        \date    Dec 2005

  !         Modified Luis Samaniego   Feb 2008 - routing module (cells)
  !                  Rohini Kumar     Aug 2011 - vector version of mHM-UFZ
  !                                   Nov 2011 - parallel version
  !                  Luis Samaniego   Jan 2013 - modularization, documentation
  ! ------------------------------------------------------------------
  subroutine L11_routing( &
       nNodes, nLinks, netPerm, netLink_fromN, netLink_toN,                                                   & 
       netLink_C1, netLink_C2, netNode_qOUT, nInflowGauges, InflowHeadwater, InflowNodeList, &
       netNode_qTIN, netNode_qTR, netNode_Qmod )

    implicit none

    ! Input
    integer(i4),                    intent(in)    :: nNodes        ! number of network nodes = nCells1
    integer(i4),                    intent(in)    :: nLinks        ! number of stream segment (reaches)
    ! Stream link description   network topology  ==>    netLink
    integer(i4), dimension(:),      intent(in)    :: netPerm       ! basin routing order (permutation)
    integer(i4), dimension(:),      intent(in)    :: netLink_fromN ! from node 
    integer(i4), dimension(:),      intent(in)    :: netLink_toN   ! to node
    real(dp),    dimension(:),      intent(in)    :: netLink_C1    ! [1]   routing parameter  C1 (Chow, 25-41)
    real(dp),    dimension(:),      intent(in)    :: netLink_C2    ! [1]   routing parameters C2 (")
    ! State variables
    real(dp),    dimension(:),      intent(in)    :: netNode_qOUT  ! [m3 s-1] Total outflow, all cells, basin, 
    !                                                                !          level L11 at time tt
    integer(i4),                    intent(in)    :: nInflowGauges   ! [-]      number of inflow points
    logical,     dimension(:),      intent(in)    :: InflowHeadwater ! [-]      if to consider headwater cells of inflow gauge
    integer(i4), dimension(:),      intent(in)    :: InflowNodeList  ! [-]      L11 ID of inflow points

    ! Input - Output
    real(dp),    dimension(:,:),    intent(inout) :: netNode_qTIN  ! [m3 s-1] Total inputs at t-1 and t
    real(dp),    dimension(:,:),    intent(inout) :: netNode_qTR   ! [m3 s-1] Transformed outflow leaving 
    !                                                                !          node I (Muskingum)
    ! Output
    real(dp),    dimension(nNodes), intent(out)   :: netNode_Qmod  ! [m3 s-1] Simulated routed discharge  


    ! local
    integer(i4)                                   :: i, k, iNode, tNode
    integer(i4), parameter                        :: IT  = 2              ! current routing state (2)
    integer(i4), parameter                        :: IT1 = 1              ! past routing state (1)
                                                                          ! Entry value for the auxiliary vectors
                                                                          !   netNode_qTIN(iNode,:)
                                                                          !   netNode_qTR(iNode,:)
                                                                          ! which store current and past states of
                                                                          ! incoming and outgoing of discharge at iNode 



    !--------------------------------------------------------------------------
    !                             Muskingum Flood Routing
    !--------------------------------------------------------------------------

    ! initialize total input at point time IT in all nodes
    netNode_qTIN(:,IT) = 0.0_dp
    !--------------------------------------------------------------------------
    ! Links in sequential mode .... with single node
    !--------------------------------------------------------------------------
    ! ST - decent parallelization has to be done!!!
    !!$OMP parallel
    !!$OMP do private( i, inode, tnode)
    do k = 1 , nLinks
       ! get LINK routing order -> i
       i = netPerm(k)
       iNode = netLink_fromN(i)
       tNode = netLink_toN(i)

       ! accumulate all inputs in iNode         
       netNode_qTIN(iNode,IT) = netNode_qTIN(iNode,IT) + netNode_qOUT(iNode)

       ! routing iNode
       netNode_qTR(iNode,IT) = netNode_qTR(iNode,IT1)                               &
            + netLink_C1(i) * ( netNode_qTIN(iNode,IT1) - netNode_qTR (iNode,IT1) ) &
            + netLink_C2(i) * ( netNode_qTIN(iNode,IT)  - netNode_qTIN(iNode,IT1) )

       ! check if the inflow from upstream cells should be deactivated
       if (nInflowGauges .GT. 0) then
          do i = 1, nInflowGauges
             ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
             if ( (tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) netNode_qTR(iNode,IT) = 0.0_dp
          end do
       end if
       
       ! add routed water to downstream node
       netNode_qTIN(tNode,IT) = netNode_qTIN(tNode,IT) + netNode_qTR(iNode,IT)
    end do
    !!$OMP end do
    !!$OMP end parallel

    !--------------------------------------------------------------------------
    ! Accumulate all inputs in tNode (netNode_qOUT) ONLY for last link
    !--------------------------------------------------------------------------
    tNode = netLink_toN( netPerm(nLinks) )
    netNode_qTIN(tNode,IT) = netNode_qTIN(tNode,IT) + netNode_qOUT(tNode)

    !--------------------------------------------------------------------------
    ! save modeled discharge at time step tt then shift flow storages
    ! (NOTE aggregation to daily values to be done outside)
    !--------------------------------------------------------------------------
    ! ST - decent parallelization has to be done!!!
    !!$OMP parallel
    !!$OMP do schedule( static )
    do i = 1, nNodes
       ! store generated discharge
       netNode_Qmod(i) = netNode_qTIN(i,IT)
       ! backflow t-> t-1
       netNode_qTR(i,IT1) = netNode_qTR(i,IT)
       netNode_qTIN(i,IT1)= netNode_qTIN(i,IT)
    end do
    !!$OMP end do
    !!$OMP end parallel

  end subroutine L11_routing
  ! ------------------------------------------------------------------

  !     NAME
  !         L11_runoff_acc

  !     PURPOSE
  !>        \brief total runoff accumulation at L11.

  !>        \details Accumulates runoff in space from L1 to L11.
  !>  

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "real(dp)    ::  qAll"              total runoff L1 [mm tst-1]
  !>        \param[in] "real(dp)    ::  efecArea"          effective area in [km2] 
  !>        \param[in] "integer(i4) ::  L11id"             L11 mapped on L1   
  !>        \param[in] "integer(i4) ::  TS"                time step in [s]
  !>        \param[in] "integer(i4) ::  nInflowGauges"     number of inflow points
  !>        \param[in] "integer(i4  ::  InflowIndexList"   index of inflow points
  !>        \param[in] "integer(i4) ::  InflowNodeList"    L11 ID of inflow points
  !>        \param[in] "real(dp)    ::  QInflow"           inflowing water 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: qOUT"                 aggregated runoff at L11 [m3 s-1]

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
  !         None

  !     HISTORY
  !>        \author Luis Samaniego
  !>        \date Jan 2013
  !         Modified  Matthias Zink , Mar 2014 - added inflow from upstream areas
  !                   Matthias Zink,  Dec 2014 - adopted inflow gauges to ignore headwater cells

  ! ------------------------------------------------------------------
  SUBROUTINE L11_runoff_acc(qAll, efecArea, L11id, TS, nInflowGauges, InflowIndexList, &
                            InflowHeadwater, InflowNodeList, QInflow, qOUT)

    use mo_mrm_constants, only:   HourSecs

    IMPLICIT NONE

    real(dp),    dimension(:), intent(in)  :: qall            ! [mm tst-1] total runoff l1 
    real(dp),    dimension(:), intent(in)  :: efecarea        ! [km2]      efective area at l1 
    integer(i4), dimension(:), intent(in)  :: l11id           ! l11        mapped on l1   
    integer(i4),               intent(in)  :: ts              ! [h]        time step 
    integer(i4),               intent(in)  :: nInflowGauges   ! [-]        number of inflow points
    integer(i4), dimension(:), intent(in)  :: InflowIndexList ! [-]        index of inflow points
    logical,     dimension(:), intent(in)  :: InflowHeadwater ! [-]        if to consider headwater cells of inflow gauge
    integer(i4), dimension(:), intent(in)  :: InflowNodeList  ! [-]        L11 ID of inflow points
    real(dp),    dimension(:), intent(in)  :: QInflow         ! [m3 s-1]   inflowing water 
    real(dp),    dimension(:), intent(out) :: qout            ! [m3 s-1]   aggregated runoff at l11 

                                                              ! local variables
    integer(i4)                            :: k
    REAL(dp)                               :: TST             ! [s]        time step  

    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ****  NO TUH :: OFF AT A MOMENT BUT EFFECTS ARE
    !                 INTEGRATED IN RIVER ROUTING PROCESS
    !          TUH    Triangular Unit Hydrograph
    ! ------------------------------------------------------------------

    TST = HourSecs * real( TS, dp)   ! in [s]

    !$OMP PARALLEL
    !$OMP DO SCHEDULE( STATIC )
    do k = 1, size(qOUT)
      !  Estimate specific runoff at  L11
      !  NOTE:
      !  1) Total discharge depth aggregated at L11 level [mm/TST]
      !  2) Transform  depth [mm/TST] to discharge [m3/s]
      !  Total runoff should be divided by total_area to get 
      !  specific discharge at L11. Then, to transform specific
      !  discharge from [mm/TST] to [m3/s], it should be multiplied by
      !  total_area [km2]*10^3 and divided by TST.
      !  Therefore, in this operation total_area cancels out. 
      !  The simplified equation is then:         
      qOUT(k) = sum( qAll(:) * efecArea(:),  L11id(:) .eq. k ) * 1000.0_dp / TST
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
    if (nInflowGauges .gt. 0) then
       do k = 1, nInflowGauges
          if (InflowHeadwater(k)) then 
             ! add inflowing water to water produced by upstream/headwater cells
             qOUT(InflowNodeList(k)) = qOUT(InflowNodeList(k)) + QInflow(InflowIndexList(k))
          else
             ! put only timeseries and cut upstream/headwater cells produced water for routing
             qOUT(InflowNodeList(k)) = QInflow(InflowIndexList(k))
          end if
       end do
    end if

  END SUBROUTINE L11_runoff_acc

END MODULE mo_routing
