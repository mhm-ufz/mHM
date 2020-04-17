module mo_mrm_river_temp
   use mo_common_variables, only: level0, domainMeta
   use mo_mrm_global_variables, only: L0_L11_remap, L11_bankfull_runoff_in, &
                                      L0_slope, L0_channel_elevation
   use mo_kind, only: i4, dp
   use mo_common_constants, only: nodata_dp
   use mo_append, only: append
   use mo_netcdf, only: NcVariable

   implicit none

   private

   ! calculates the river temperature
   ! public :: calc_river_temp

contains

   !    NAME
   !        L11_temperature_acc

   !    PURPOSE
   !>       \brief total runoff accumulation at L11.

   !>       \details Upscales runoff in space from L1 to L11 if routing resolution
   !>       is higher than hydrology resolution (map_flag equals .true.) or
   !>       downscales runoff from L1 to L11 if routing resolution is lower
   !>       than hydrology resolution.

   !    INTENT(IN)
   !>       \param[in] "real(dp), dimension(:) :: qall"         total runoff L1 [mm tst-1]
   !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
   !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
   !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
   !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
   !>       \param[in] "integer(i4) :: TS"                      time step in [s]
   !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
   !>       hydrologic one

   !    INTENT(OUT)
   !>       \param[out] "real(dp), dimension(:) :: qAcc" aggregated runoff at L11 [m3 s-1]

   !    HISTORY
   !>       \authors Sebastian Mueller

   !>       \date Apr 2020

   ! Modifications:
   ! Sebastian Mueller  Apr 2020 - Creation

   SUBROUTINE L11_temperature_acc(TqAll, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, TqAcc)

      use mo_common_constants, only: HourSecs, nodata_dp

      implicit none

      ! total runoff L1 [mm tst-1]
      real(dp), intent(in), dimension(:) :: Tqall
      ! effective area in [km2] at Level 1
      real(dp), intent(in), dimension(:) :: efecarea
      ! L11 Ids mapped on L1
      integer(i4), intent(in), dimension(:) :: L1_L11_Id
      ! effective area in [km2] at Level 11
      real(dp), intent(in), dimension(:) :: L11_areacell
      ! L1 Ids mapped on L11
      integer(i4), intent(in), dimension(:) :: L11_L1_Id
      ! time step in [h]
      integer(i4), intent(in) :: TS
      ! Flag indicating whether routing resolution is higher than hydrologic one
      logical, intent(in) :: map_flag
      ! aggregated runoff at L11 [m3 s-1]
      real(dp), intent(out), dimension(:) :: TqAcc

      integer(i4) :: k
      ! [s] time step
      real(dp) :: TST

      ! ------------------------------------------------------------------
      ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
      ! ------------------------------------------------------------------
      ! Hydrologic timestep in seconds
      ! TST = HourSecs*TS

      ! if (map_flag) then
      !    TqAcc = 0._dp
      !    ! loop over high-resolution cells (L1) and add discharge to
      !    ! corresponding low-resolution cells (L11)
      !    do k = 1, size(qAll, 1)
      !       TqAcc(L1_L11_Id(k)) = TqAcc(L1_L11_Id(k)) + TqAll(k)*efecArea(k)
      !    end do
      !    TqAcc = TqAcc*1000.0_dp/TST
      !    !
      ! else
      !    ! initialize qout
      !    TqAcc = nodata_dp
      !    do k = 1, size(qAcc, 1)
      !       ! map temp-energy flux from coarse L1 resolution to fine L11 resolution
      !       TqAcc(k) = TqAll(L11_L1_Id(k))
      !    end do
      !    ! adjust temp-energy flux by area cell
      !    TqAcc(:) = TqAcc(:)*L11_areaCell(:)*1000.0_dp/TST
      ! end if

   END SUBROUTINE L11_temperature_acc

   ! ------------------------------------------------------------------

   !    NAME
   !        add_inflow_E

   !    PURPOSE
   !>       \brief
   !>          Adds inflow energy flux (T*q) to the energy-flux produced at the
   !>          cell where the inflow is occurring.
   !>       \details
   !>          If a inflow gauge is given, then this routine is adding the
   !>          values to the energy-flux produced at the grid cell where the
   !>          inflow is happening. The values are not directly added to the
   !>          river network. If this cell is not a headwater then the streamflow
   !>          produced upstream will be neglected.

   !    INTENT(IN)
   !>       \param[in] "integer(i4) :: nInflowGauges"                 [-]      number of inflow points
   !>       \param[in] "integer(i4), dimension(:) :: InflowIndexList" [-]      index of inflow points
   !>       \param[in] "logical, dimension(:) :: InflowHeadwater"     [-]      if to consider headwater cells of inflow gauge
   !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList"  [-]      L11 ID of inflow points
   !>       \param[in] "real(dp), dimension(:) :: QInflow"            [m3 s-1] inflowing water
   !>       \param[in] "real(dp), dimension(:) :: TInflow"            [K]      inflowing water temperature

   !    INTENT(INOUT)
   !>       \param[inout] "real(dp), dimension(:) :: EOut" [K m3 s-1] Series of attenuated energy-flux

   !    HISTORY
   !>       \author Sebastian Mueller

   !>       \date Apr 2020

   subroutine add_inflow_E(nInflowGauges, InflowIndexList, InflowHeadwater, InflowNodeList, QInflow, TInflow, EOut)

      use mo_kind, only: dp, i4

      implicit none

      ! [-] number of inflow points
      integer(i4), intent(in) :: nInflowGauges
      ! [-] index of inflow points
      integer(i4), intent(in), dimension(:) :: InflowIndexList
      ! [-] if to consider headwater cells of inflow gauge
      logical, intent(in), dimension(:) :: InflowHeadwater
      ! [-]        L11 ID of inflow points
      integer(i4), intent(in), dimension(:) :: InflowNodeList
      ! [m3 s-1]   inflowing water
      real(dp), intent(in), dimension(:) :: QInflow
      ! [K]   inflowing water temperature
      real(dp), intent(in), dimension(:) :: TInflow
      ! [m3 s-1] Series of attenuated runoff
      real(dp), intent(inout), dimension(:) :: EOut
      ! looping var
      integer(i4) :: ii

      ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
      ! should be put after UH attenuation because it is measured runoff at this cell
      if (nInflowGauges .gt. 0) then
         do ii = 1, nInflowGauges
            if (InflowHeadwater(ii)) then
               ! add inflowing energy flux to energy flux produced by upstream/headwater cells
               EOut(InflowNodeList(ii)) = EOut(InflowNodeList(ii)) + QInflow(InflowIndexList(ii))*TInflow(InflowIndexList(ii))
            else
               ! put only timeseries and cut upstream/headwater cells produced energy flux for routing
               EOut(InflowNodeList(ii)) = QInflow(InflowIndexList(ii))*TInflow(InflowIndexList(ii))
            end if
         end do
      end if
   end subroutine add_inflow_E

   ! ------------------------------------------------------------------

   !    NAME
   !        L11_routing_E

   !    PURPOSE
   !>       \brief Performs temperature routing for mHM at L11 upscaled network
   !>       (\ref fig_routing "Routing Network").
   !>       \details
   !>       Hydrograph routing is carried out with the Muskingum algorithm
   !>       \cite CMM1988.  This simplification of the St. Venant
   !>       equations is justified in mHM because the potential areas of
   !>       application of this model would hardly exhibit abruptly
   !>       changing hydrographs with supercritical flows.

   !    INTENT(IN)
   !>       \param[in] "integer(i4) :: nNodes"                       number of network nodes = nCells1
   !>       \param[in] "integer(i4) :: nLinks"                       number of stream segment (reaches)
   !>       \param[in] "integer(i4), dimension(:) :: netPerm"        routing order of a given domain (permutation)
   !>       \param[in] "integer(i4), dimension(:) :: netLink_fromN"  from node
   !>       \param[in] "integer(i4), dimension(:) :: netLink_toN"    to node
   !>       \param[in] "real(dp), dimension(:) :: netLink_C1"        routing parameter  C1 (\cite CMM1988 p. 25-41)
   !>       \param[in] "real(dp), dimension(:) :: netLink_C2"        routing parameters C2 (id)
   !>       \param[in] "real(dp), dimension(:) :: netNode_TqOUT"      Total outflow from cells (given domain) L11 at time
   !>       tt in [m3 s-1]
   !>       \param[in] "integer(i4) :: nInflowGauges"                [-]      number of inflow points
   !>       \param[in] "logical, dimension(:) :: InflowHeadwater"    [-]      if to consider headwater cells of inflow
   !>       gauge
   !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList" [-]      L11 ID of inflow points

   !    INTENT(INOUT)
   !>       \param[inout] "real(dp), dimension(:, :) :: netNode_TqTIN" [m3 s-1] Total inputs at t-1 and t
   !>       \param[inout] "real(dp), dimension(:, :) :: netNode_TqTR"  [m3 s-1] Transformed outflow leaving node I (Muskingum)

   !    INTENT(OUT)
   !>       \param[out] "real(dp), dimension(nNodes) :: netNode_temp" [m3 s-1] Simulated routed discharge

   !    HISTORY
   !>       \authors Sebastian Mueller

   !>       \date Apr 2020

   ! Modifications:
   ! Sebastian Mueller Apr 2020 - creation

   subroutine L11_routing_E(nNodes, nLinks, netPerm, netLink_fromN, netLink_toN, netLink_C1, netLink_C2, netNode_TqOUT, &
                            nInflowGauges, InflowHeadwater, InflowNodeList, netNode_TqTIN, netNode_TqTR, netNode_temp)
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
      real(dp), dimension(:), intent(in) :: netNode_TqOUT
      ! [-]      number of inflow points
      integer(i4), intent(in) :: nInflowGauges
      ! [-]      if to consider headwater cells of inflow gauge
      logical, dimension(:), intent(in) :: InflowHeadwater
      ! [-]      L11 ID of inflow points
      integer(i4), dimension(:), intent(in) :: InflowNodeList
      ! [m3 s-1] Total inputs at t-1 and t
      real(dp), dimension(:, :), intent(inout) :: netNode_TqTIN
      ! [m3 s-1] Transformed outflow leaving
      ! node I (Muskingum)
      real(dp), dimension(:, :), intent(inout) :: netNode_TqTR
      ! [m3 s-1] Simulated routed discharge
      real(dp), dimension(nNodes), intent(out) :: netNode_temp

      integer(i4) :: i, k, iNode, tNode
      ! current routing state (2)
      integer(i4), parameter :: IT = 2
      ! past routing state (1)
      integer(i4), parameter :: IT1 = 1

      ! Entry value for the auxiliary vectors
      !   netNode_TqTIN(iNode,:)
      !   netNode_TqTR(iNode,:)
      ! which store current and past states of
      ! incoming and outgoing of discharge at iNode
      !--------------------------------------------------------------------------
      !                             Muskingum Flood Routing
      !--------------------------------------------------------------------------
      ! initialize total input at point time IT in all nodes
      ! netNode_TqTIN(:, IT) = 0.0_dp
      ! !--------------------------------------------------------------------------
      ! ! Links in sequential mode .... with single node
      ! !--------------------------------------------------------------------------
      ! do k = 1, nLinks
      !    ! get LINK routing order -> i
      !    i = netPerm(k)
      !    iNode = netLink_fromN(i)
      !    tNode = netLink_toN(i)

      !    ! accumulate all inputs in iNode
      !    netNode_TqTIN(iNode, IT) = netNode_TqTIN(iNode, IT) + netNode_TqOUT(iNode)

      !    ! routing iNode
      !    netNode_TqTR(iNode, IT) = netNode_TqTR(iNode, IT1) &
      !                              + netLink_C1(i)*(netNode_TqTIN(iNode, IT1) - netNode_TqTR(iNode, IT1)) &
      !                              + netLink_C2(i)*(netNode_TqTIN(iNode, IT) - netNode_TqTIN(iNode, IT1))

      !    ! check if the inflow from upstream cells should be deactivated
      !    if (nInflowGauges .GT. 0) then
      !       do i = 1, nInflowGauges
      !          ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
      !          if ((tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) netNode_TqTR(iNode, IT) = 0.0_dp
      !       end do
      !    end if

      !    ! add routed water to downstream node
      !    netNode_TqTIN(tNode, IT) = netNode_TqTIN(tNode, IT) + netNode_TqTR(iNode, IT)
      ! end do

      ! ! --------------------------------------------------------------------------
      ! !  Accumulate all inputs in tNode (netNode_TqOUT) ONLY for last link
      ! ! --------------------------------------------------------------------------
      ! tNode = netLink_toN(netPerm(nLinks))
      ! netNode_TqTIN(tNode, IT) = netNode_TqTIN(tNode, IT) + netNode_TqOUT(tNode)

      ! ! --------------------------------------------------------------------------
      ! !  save modeled discharge at time step tt then shift flow storages
      ! !  (NOTE aggregation to daily values to be done outside)
      ! ! --------------------------------------------------------------------------
      ! !  store generated discharge
      ! netNode_temp(1:nNodes) = netNode_TqTIN(1:nNodes, IT)
      ! !  backflow t-> t-1
      ! netNode_TqTR(1:nNodes, IT1) = netNode_TqTR(1:nNodes, IT)
      ! netNode_TqTIN(1:nNodes, IT1) = netNode_TqTIN(1:nNodes, IT)

   end subroutine L11_routing_E

end module mo_mrm_river_temp
