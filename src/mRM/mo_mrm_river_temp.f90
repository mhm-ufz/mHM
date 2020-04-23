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

   !    intent(IN)
   !>       \param[in] "real(dp), dimension(:) :: qall"         total runoff L1 [mm tst-1]
   !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
   !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
   !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
   !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
   !>       \param[in] "integer(i4) :: TS"                      time step in [s]
   !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
   !>       hydrologic one

   !    intent(OUT)
   !>       \param[out] "real(dp), dimension(:) :: qAcc" aggregated runoff at L11 [m3 s-1]

   !    HISTORY
   !>       \authors Sebastian Mueller

   !>       \date Apr 2020

   ! Modifications:
   ! Sebastian Mueller  Apr 2020 - Creation

   subroutine L11_acc_E(E_all, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, TqAcc)

      use mo_common_constants, only: HourSecs, nodata_dp

      implicit none

      ! total energy flux L1 [mm tst-1]
      real(dp), intent(in), dimension(:) :: E_all
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
      !       TqAcc(L1_L11_Id(k)) = TqAcc(L1_L11_Id(k)) + E_all(k)*efecArea(k)
      !    end do
      !    TqAcc = TqAcc*1000.0_dp/TST
      !    !
      ! else
      !    ! initialize qout
      !    TqAcc = nodata_dp
      !    do k = 1, size(qAcc, 1)
      !       ! map temp-energy flux from coarse L1 resolution to fine L11 resolution
      !       TqAcc(k) = E_all(L11_L1_Id(k))
      !    end do
      !    ! adjust temp-energy flux by area cell
      !    TqAcc(:) = TqAcc(:)*L11_areaCell(:)*1000.0_dp/TST
      ! end if

   end subroutine L11_acc_E

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

   !    intent(IN)
   !>       \param[in] "integer(i4) :: nInflowGauges"                 [-]      number of inflow points
   !>       \param[in] "integer(i4), dimension(:) :: InflowIndexList" [-]      index of inflow points
   !>       \param[in] "logical, dimension(:) :: InflowHeadwater"     [-]      if to consider headwater cells of inflow gauge
   !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList"  [-]      L11 ID of inflow points
   !>       \param[in] "real(dp), dimension(:) :: QInflow"            [m3 s-1] inflowing water
   !>       \param[in] "real(dp), dimension(:) :: TInflow"            [K]      inflowing water temperature

   !    intent(INOUT)
   !>       \param[inout] "real(dp), dimension(:) :: E_out" [K m3 s-1] Series of attenuated energy-flux

   !    HISTORY
   !>       \author Sebastian Mueller

   !>       \date Apr 2020

   subroutine add_inflow_E(nInflowGauges, InflowIndexList, InflowHeadwater, InflowNodeList, QInflow, TInflow, E_out)

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
      real(dp), intent(inout), dimension(:) :: E_out
      ! looping var
      integer(i4) :: ii

      ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
      ! should be put after UH attenuation because it is measured runoff at this cell
      if (nInflowGauges .gt. 0) then
         do ii = 1, nInflowGauges
            if (InflowHeadwater(ii)) then
               ! add inflowing energy flux to energy flux produced by upstream/headwater cells
               E_out(InflowNodeList(ii)) = E_out(InflowNodeList(ii)) + QInflow(InflowIndexList(ii))*TInflow(InflowIndexList(ii))
            else
               ! put only timeseries and cut upstream/headwater cells produced energy flux for routing
               E_out(InflowNodeList(ii)) = QInflow(InflowIndexList(ii))*TInflow(InflowIndexList(ii))
            end if
         end do
      end if
   end subroutine add_inflow_E

   ! ------------------------------------------------------------------

   !    NAME
   !        L11_routing_E

   !    PURPOSE
   !>       \brief Performs temperature energy routing for mHM at L11 upscaled network
   !>       (\ref fig_routing "Routing Network").
   !>       \details
   !>       Hydrograph routing is carried out with the Muskingum algorithm
   !>       \cite CMM1988.  This simplification of the St. Venant
   !>       equations is justified in mHM because the potential areas of
   !>       application of this model would hardly exhibit abruptly
   !>       changing hydrographs with supercritical flows.

   !    intent(IN)
   !>       \param[in] "integer(i4) :: nNodes"                       number of network nodes = nCells1
   !>       \param[in] "integer(i4) :: nLinks"                       number of stream segment (reaches)
   !>       \param[in] "integer(i4), dimension(:) :: netPerm"        routing order of a given domain (permutation)
   !>       \param[in] "integer(i4), dimension(:) :: netLink_fromN"  from node
   !>       \param[in] "integer(i4), dimension(:) :: netLink_toN"    to node
   !>       \param[in] "real(dp), dimension(:) :: netLink_C1"        routing parameter  C1 (\cite CMM1988 p. 25-41)
   !>       \param[in] "real(dp), dimension(:) :: netLink_C2"        routing parameters C2 (id)
   !>       \param[in] "real(dp), dimension(:) :: netNode_E_lat"      Total outflow from cells (given domain) L11 at time
   !>       tt in [m3 s-1]
   !>       \param[in] "integer(i4) :: nInflowGauges"                [-]      number of inflow points
   !>       \param[in] "logical, dimension(:) :: InflowHeadwater"    [-]      if to consider headwater cells of inflow
   !>       gauge
   !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList" [-]      L11 ID of inflow points

   !    intent(INOUT)
   !>       \param[inout] "real(dp), dimension(:, :) :: netNode_E_in" [m3 s-1] Total inputs at t-1 and t
   !>       \param[inout] "real(dp), dimension(:, :) :: netNode_E_rout"  [m3 s-1] Transformed energy flux leaving node I (Muskingum)

   !    intent(OUT)
   !>       \param[out] "real(dp), dimension(nNodes) :: netNode_E_gen" [m3 s-1] Generated routed energy flux

   !    HISTORY
   !>       \authors Sebastian Mueller

   !>       \date Apr 2020

   ! Modifications:
   ! Sebastian Mueller Apr 2020 - creation

   subroutine L11_routing_E(nNodes, nLinks, netPerm, netLink_fromN, netLink_toN, netLink_C1, netLink_C2, netNode_E_lat, &
                            nInflowGauges, InflowHeadwater, InflowNodeList, netNode_E_in, netNode_E_rout, netNode_E_gen)
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
      real(dp), dimension(:), intent(in) :: netNode_E_lat
      ! [-]      number of inflow points
      integer(i4), intent(in) :: nInflowGauges
      ! [-]      if to consider headwater cells of inflow gauge
      logical, dimension(:), intent(in) :: InflowHeadwater
      ! [-]      L11 ID of inflow points
      integer(i4), dimension(:), intent(in) :: InflowNodeList
      ! [m3 s-1] Total inputs at t-1 and t
      real(dp), dimension(:, :), intent(inout) :: netNode_E_in
      ! [m3 s-1] Transformed outflow leaving
      ! node I (Muskingum)
      real(dp), dimension(:, :), intent(inout) :: netNode_E_rout
      ! [m3 s-1] Simulated routed discharge
      real(dp), dimension(nNodes), intent(out) :: netNode_E_gen

      integer(i4) :: i, k, iNode, tNode
      ! current routing state (2)
      integer(i4), parameter :: IT = 2
      ! past routing state (1)
      integer(i4), parameter :: IT1 = 1

      ! Entry value for the auxiliary vectors
      !   netNode_E_in(iNode,:)
      !   netNode_E_rout(iNode,:)
      ! which store current and past states of
      ! incoming and outgoing of discharge at iNode
      !--------------------------------------------------------------------------
      !                             Muskingum Flood Routing
      !--------------------------------------------------------------------------
      ! initialize total input at point time IT in all nodes
      netNode_E_in(:, IT) = 0.0_dp
      !--------------------------------------------------------------------------
      ! Links in sequential mode .... with single node
      !--------------------------------------------------------------------------
      do k = 1, nLinks
         ! get LINK routing order -> i
         i = netPerm(k)
         iNode = netLink_fromN(i)
         tNode = netLink_toN(i)

         ! accumulate all inputs in iNode
         netNode_E_in(iNode, IT) = netNode_E_in(iNode, IT) + netNode_E_lat(iNode)

         ! routing iNode
         netNode_E_rout(iNode, IT) = netNode_E_rout(iNode, IT1) &
                                   + netLink_C1(i)*(netNode_E_in(iNode, IT1) - netNode_E_rout(iNode, IT1)) &
                                   + netLink_C2(i)*(netNode_E_in(iNode, IT) - netNode_E_in(iNode, IT1))

         ! check if the inflow from upstream cells should be deactivated
         if (nInflowGauges .GT. 0) then
            do i = 1, nInflowGauges
               ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
               if ((tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) netNode_E_rout(iNode, IT) = 0.0_dp
            end do
         end if

         ! add routed water to downstream node
         netNode_E_in(tNode, IT) = netNode_E_in(tNode, IT) + netNode_E_rout(iNode, IT)
      end do

      ! --------------------------------------------------------------------------
      !  Accumulate all inputs in tNode (netNode_E_lat) ONLY for last link
      ! --------------------------------------------------------------------------
      tNode = netLink_toN(netPerm(nLinks))
      netNode_E_in(tNode, IT) = netNode_E_in(tNode, IT) + netNode_E_lat(tNode)

      ! --------------------------------------------------------------------------
      !  save modeled discharge at time step tt then shift flow storages
      !  (NOTE aggregation to daily values to be done outside)
      ! --------------------------------------------------------------------------
      !  store generated energy fluxes
      netNode_E_gen(1:nNodes) = netNode_E_in(1:nNodes, IT)
      !  backflow t-> t-1
      netNode_E_rout(1:nNodes, IT1) = netNode_E_rout(1:nNodes, IT)
      netNode_E_in(1:nNodes, IT1) = netNode_E_in(1:nNodes, IT)

   end subroutine L11_routing_E

   ! ------------------------------------------------------------------

   !    NAME
   !        L1_total_runoff

   !    PURPOSE
   !>       \brief total runoff accumulation at level 1

   !>       \details Accumulates runoff.
   !>       \f[ q_{T} = ( q_0 + q_1 + q_2 ) * (1-fSealed) + q_{D} * fSealed \f],
   !>       where fSealed is the fraction of sealed area.

   !    intent(IN)
   !>       \param[in] "real(dp) :: fSealed_area_fraction" sealed area fraction [1]
   !>       \param[in] "real(dp) :: fast_interflow"        \f$ q_0 \f$ Fast runoff component [mm tst-1]
   !>       \param[in] "real(dp) :: slow_interflow"        \f$ q_1 \f$ Slow runoff component [mm tst-1]
   !>       \param[in] "real(dp) :: baseflow"              \f$ q_2 \f$ Baseflow [mm tsts-1]
   !>       \param[in] "real(dp) :: direct_runoff"         \f$ q_D \f$ Direct runoff from impervious areas  [mm tst-1]

   !    intent(OUT)
   !>       \param[out] "real(dp) :: total_runoff" \f$ q_T \f$ Generated runoff [mm tst-1]

   !    HISTORY
   !>       \authors Sebastian Mueller

   !>       \date Apr 2020
   subroutine L11_lateral_E( &
      fSealed_area_fraction, &
      fast_interflow, &
      slow_interflow, &
      baseflow, &
      direct_runoff, &
      temp_interflow, &
      temp_baseflow, &
      temp_direct_runoff, &
      total_runoff)

      implicit none

      ! sealed area fraction [1]
      real(dp), intent(IN) :: fSealed_area_fraction
      ! \f$ q_0 \f$ Fast runoff component [mm tst-1]
      real(dp), intent(IN) :: fast_interflow
      ! \f$ q_1 \f$ Slow runoff component [mm tst-1]
      real(dp), intent(IN) :: slow_interflow
      ! \f$ q_2 \f$ Baseflow [mm tsts-1]
      real(dp), intent(IN) :: baseflow
      ! \f$ q_D \f$ Direct runoff from impervious areas  [mm tst-1]
      real(dp), intent(IN) :: direct_runoff
      ! \f$ q_D \f$ temperature
      real(dp), intent(IN) :: temp_interflow
      ! \f$ q_D \f$ temperature
      real(dp), intent(IN) :: temp_baseflow
      ! \f$ q_D \f$ temperature
      real(dp), intent(IN) :: temp_direct_runoff
      ! \f$ q_T \f$ Generated runoff [mm tst-1]
      real(dp), intent(OUT) :: total_runoff

      total_runoff = ((baseflow + slow_interflow + fast_interflow)*(1.0_dp - fSealed_area_fraction)) + &
                     (direct_runoff*fSealed_area_fraction)

   end subroutine L11_lateral_E

end module mo_mrm_river_temp
