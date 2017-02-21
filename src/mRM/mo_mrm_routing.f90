!> \file mo_mrm_routing.f90

!> \brief Performs runoff routing for mHM at level L11.

!> \details This module performs flood routing at a given time step
!>          through the stream network at level L11 to the sink cell. 
!>          The Muskingum flood routing algorithm is used.

!> \author Luis Samaniego
!> \date Dec 2012
!  Modified
!       Stephan Thober, Aug 2015 - adapted to mRM

MODULE mo_mrm_routing

  ! This module performs runoff flood routing for mHM.

  ! Written Luis Samaniego, Dec 2012

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mRM_routing

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         mrm_routing

  !     PURPOSE
  !>        \brief route water given runoff
  !
  !>        \details This routine first performs mpr for the routing variables
  !>        if required, then accumulates the runoff to the routing resolution
  !>        and eventually routes the water in a third step. The last step is
  !>        repeated multiple times if the routing timestep is smaller than
  !>        the timestep of the hydrological timestep

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: processCase" Process switch for routing
  !>        \param[in] "real(dp), dimension(5)    :: global_routing_param" routing parameters
  !>        \param[in] "real(dp), dimension(:)    :: L1_total_runoff" total runoff from L1 grid cells
  !>        \param[in] "real(dp), dimension(:)    :: L1_areaCell" L1 cell area
  !>        \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L1 cell ids on L11
  !>        \param[in] "real(dp), dimension(:)    :: L11_areaCell" L11 cell area
  !>        \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L11 cell ids on L1
  !>        \param[in] "integer(i4), dimension(:) :: L11_netPerm" L11 routing order
  !>        \param[in] "integer(i4), dimension(:) :: L11_fromN" L11 source grid cell order
  !>        \param[in] "integer(i4), dimension(:) :: L11_toN" L11 target grid cell order
  !>        \param[in] "integer(i4)               :: timestep" simulation timestep in [h]
  !>        \param[in] "real(dp)                  :: tsRoutFactor" factor between routing timestep and hydrological timestep
  !>        \param[in] "integer(i4)               :: nNodes" number of nodes
  !>        \param[in] "integer(i4)               :: nInflowGauges" number of inflow gauges
  !>        \param[in] "integer(i4), dimension(:) :: InflowGaugeIndexList" index list of inflow gauges
  !>        \param[in] "logical, dimension(:)     :: InflowGaugeHeadwater" flag for headwater cell of inflow gauge
  !>        \param[in] "integer(i4), dimension(:) :: InflowGaugeNodeList" gauge node list at L11
  !>        \param[in] "real(dp), dimension(:)    :: InflowDischarge" inflowing discharge at discharge gauge at current day
  !>        \param[in] "integer(i4)               :: nGauges" number of recording gauges
  !>        \param[in] "integer(i4), dimension(:) :: gaugeIndexList" index list for outflow gauges
  !>        \param[in] "integer(i4), dimension(:) :: gaugeNodeList" gauge node list at L11
  !>        \param[in] "logical                   :: map_flag" flag indicating whether routing resolution is
  !>                                                 coarser than hydrologic resolution
  !>        \param[in] "integer(i4), dimension(:) :: L0_LCover" L0 land cover
  !>        \param[in] "integer(i4), dimension(:) :: L0_floodPlain" L0 fraction of flood plains
  !>        \param[in] "real(dp), dimension(:)    :: L0_areaCell" L0 cell area
  !>        \param[in] "real(dp), dimension(:)    :: L11_aFloodPlain" L11 area of flood plain
  !>        \param[in] "real(dp), dimension(:)    :: L11_length" L11 link length
  !>        \param[in] "real(dp), dimension(:)    :: L11_slope" L11 slope
  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:) :: L11_C1" L11 muskingum parameter 1
  !>        \param[inout] "real(dp), dimension(:) :: L11_C2" L11 muskingum parameter 2
  !>        \param[inout] "real(dp), dimension(:) :: L11_qOut" total runoff from L11 grid cells
  !>        \param[inout] "real(dp), dimension(:,:) :: L11_qTIN" L11 inflow to the reach
  !>        \param[inout] "real(dp), dimension(:,:) :: L11_qTR" L11 routed outflow
  !>        \param[inout] "real(dp), dimension(:) :: L11_qMod" modelled discharge at each grid cell
  !>        \param[inout] "real(dp), dimension(:) :: GaugeDischarge" modelled discharge at each gauge
  !>        \param[inout] "real(dp), dimension(:) :: L11_FracFPimp" L11 fraction of flood plain with impervios cover

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, optional :: do_mpr_routing" indicate whether routing is to be performed"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         Modified, Sep 2015, Stephan Thober - using arguments instead of global variables
  !                   Sep 2015, Stephan Thober - added variables for routing resolution higher than hydrologic resolution
  !                   May 2016, Stephan Thober - added check whether gauge is actually inside modelling domain
  !                                              before copying simulated runoff
  !                   Nov 2016, Stephan Thober - implemented second routing process i.e. adaptive timestep

  subroutine mRM_routing( &
       !
       ! general input variables =================================================
       processCase, & ! processCase for Routing
       global_routing_param, & ! routing parameters
       L1_total_runoff, & ! total runoff from L1 grid cells
       L1_areaCell, & ! L1 cell area
       L1_L11_Id, & ! L11 cell ids on Level 1
       L11_areaCell, & ! area cell [km2] at Level11
       L11_L1_Id, & ! L1 cell ids on Level 11 
       L11_netPerm, & ! L11 routing order
       L11_fromN, & ! L11 source grid cell order
       L11_toN, & ! L11 target grid cell order
       L11_nOutlets, & ! number of outlets
       timestep, & ! simulation timestep in [h]
       tsRoutFactor, & ! factor from routing temporal resolution to hydrology temporal resolution
       nNodes, & ! number of nodes
       nInflowGauges, & ! number of inflow gauges
       InflowGaugeIndexList, & ! index list of inflow gauges
       InflowGaugeHeadwater, & ! flag for headwater cell of inflow gauge
       InflowGaugeNodeList, & ! gauge node list at L11
       InflowDischarge, & ! inflowing discharge at discharge gauge at current day
       nGauges, & ! number of recording gauges
       gaugeIndexList, & ! index list for outflow gauges
       gaugeNodeList, & ! gauge node list at L11s
       map_flag, & ! flag indicating whether routing resolution is larger than hydrologic resolution
       !
       ! original routing specific input variables ===============================
       L0_LCover, & ! L0 land cover
       L0_floodPlain, & ! L0 fraction of flood plains
       L0_areaCell, & ! L0 cell area
       L11_aFloodPlain, & ! L11 area of flood plain
       L11_length, & ! L11 link length
       L11_slope, & ! L11 slope
       !
       ! general input/output variables ==========================================
       L11_C1, & ! L11 muskingum parameter 1
       L11_C2, & ! L11 muskingum parameter 2
       L11_qOut, & ! total runoff from L11 grid cells
       L11_qTIN, & ! L11 inflow to the reach
       L11_qTR, & ! L11 routed outflow
       L11_qMod, &
       GaugeDischarge, &
       !
       ! original routing specific input/output variables ========================
       L11_FracFPimp, & ! L11 fraction of flood plain with impervios cover
       !
       ! optional input variables ================================================
       do_mpr_routing &
       )
    use mo_mrm_global_variables, only: is_start
    use mo_mrm_mpr, only: reg_rout
    use mo_mrm_net_startup, only: L11_fraction_sealed_floodplain

    implicit none
    ! general input variables ====================================================
    integer(i4),               intent(in) :: processCase
    real(dp), dimension(:),    intent(in) :: global_routing_param ! routing parameters
    real(dp), dimension(:),    intent(in) :: L1_total_runoff ! total runoff from L1 grid cells
    real(dp), dimension(:),    intent(in) :: L1_areaCell ! L1 cell area
    integer(i4), dimension(:), intent(in) :: L1_L11_Id ! L11 cell ids at L1
    real(dp), dimension(:),    intent(in) :: L11_areaCell ! L11 area cell [km2]
    integer(i4), dimension(:), intent(in) :: L11_L1_Id ! L1 cell ids at L11 
    integer(i4), dimension(:), intent(in) :: L11_netPerm ! L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_fromN ! L11 source grid cell order
    integer(i4), dimension(:), intent(in) :: L11_toN ! L11 target grid cell order
    integer(i4),               intent(in) :: L11_nOutlets ! L11 number of outlets/sinks
    integer(i4),               intent(in) :: timestep ! simulation timestep in [h]
    real(dp),                  intent(in) :: tsRoutFactor ! factor from routing temporal resolution to hydrology temporal resolution
    integer(i4),               intent(in) :: nNodes ! number of nodes
    integer(i4),               intent(in) :: nInflowGauges ! number of inflow gauges
    integer(i4), dimension(:), intent(in) :: InflowGaugeIndexList ! index list of inflow gauges
    logical, dimension(:),     intent(in) :: InflowGaugeHeadwater ! flag for headwater cell of inflow gauge
    integer(i4), dimension(:), intent(in) :: InflowGaugeNodeList ! gauge node list at L11
    real(dp), dimension(:),    intent(in) :: InflowDischarge ! inflowing discharge at discharge gauge at current day
    integer(i4),               intent(in) :: nGauges ! number of recording gauges
    integer(i4), dimension(:), intent(in) :: gaugeIndexList ! index list for outflow gauges
    integer(i4), dimension(:), intent(in) :: gaugeNodeList ! gauge node list at L11
    logical,                   intent(in) :: map_flag ! flag indicating whether routing resolution > hydrology resolution
    !
    ! original routing specific input variables ===============================
    integer(i4), dimension(:), intent(in) :: L0_LCover ! L0 land cover
    integer(i4), dimension(:), intent(in) :: L0_floodPlain ! L0 fraction of flood plains
    real(dp),    dimension(:), intent(in) :: L0_areaCell ! L0 cell area
    real(dp),    dimension(:), intent(in) :: L11_aFloodPlain ! L11 area of flood plain
    real(dp),    dimension(:), intent(in) :: L11_length ! L11 link length
    real(dp),    dimension(:), intent(in) :: L11_slope ! L11 slope
    !
    ! input/output variables ==================================================
    real(dp), dimension(:), intent(inout) :: L11_C1 ! L11 muskingum parameter 1
    real(dp), dimension(:), intent(inout) :: L11_C2 ! L11 muskingum parameter 2
    real(dp), dimension(:), intent(inout) :: L11_qOut ! total runoff from L11 grid cells
    real(dp), dimension(:,:), intent(inout) :: L11_qTIN ! L11 inflow to the reach
    real(dp), dimension(:,:), intent(inout) :: L11_qTR ! L11 routed outflow
    real(dp), dimension(:), intent(inout) :: L11_qMod ! modelled discharge at each grid cell
    real(dp), dimension(:), intent(inout) :: GaugeDischarge ! modelled discharge at each gauge
    !
    ! original routing specific input/output variables ========================
    real(dp), dimension(:), intent(inout) :: L11_FracFPimp ! L11 fraction of flood plain with impervios cover
    !
    ! optional input variables ================================================
    logical, optional, intent(in) :: do_mpr_routing ! flag for performing mpr
    ! local variables =========================================================
    logical     :: do_mpr
    integer(i4) :: gg
    integer(i4) :: tt
    integer(i4) :: rout_loop ! number of routing loops
    real(dp)    :: L11_qAcc(size(L11_qMod, dim=1))     ! variable for accumulation

    ! ====================================================================
    ! FIRST, EXECUTE MPR
    ! ====================================================================
    ! evaluate whether mpr should be executed
    do_mpr = .false.
    if (is_start) then
       do_mpr = .true.
       is_start = .false.
    end if
    if (present(do_mpr_routing)) do_mpr = do_mpr_routing
       
    if (processCase .eq. 1_i4) then
       ! execute mpr for routing
       if (do_mpr) then
          !-------------------------------------------------------------------
          ! estimate fraction of impervious cover in flood plains
          ! --> time independent variable: to be initalized every time
          !     with landcover change 
          ! --> Note: L11_fraction_sealed_floodplain routine is called
          !           only in case when routing process is ON
          !-------------------------------------------------------------------
          CALL L11_fraction_sealed_floodplain( nNodes - L11_nOutlets, &
               L0_LCover, &
               L0_floodPlain,       &
               L0_areaCell, &
               L11_aFloodPlain, &
               2, &
               L11_FracFPimp) ! intent out
          ! for a single node model run
          if (nNodes .GT. 1) then
             call reg_rout( global_routing_param, &
                  L11_length, L11_slope, L11_FracFPimp(:nNodes - L11_nOutlets), &
                  real(timeStep,dp), L11_C1(:nNodes - L11_nOutlets), L11_C2(:nNodes - L11_nOutlets))
          end if
       end if
    end if

    ! =====================================================================
    ! SECOND, EXECUTE ROUTING
    ! ====================================================================
    ! calculate number of routing loops
    rout_loop = max(1_i4, nint(1._dp / tsRoutFactor))
    

    ! runoff accumulation from L1 to L11 level
    call L11_runoff_acc(L1_total_runoff, L1_areaCell, L1_L11_Id, &
         L11_areaCell, L11_L1_Id, timeStep, & ! Intent IN
         map_flag, & ! Intent IN
         L11_qOut) ! Intent OUT

    ! add inflow
    call add_inflow( nInflowGauges, &
         InflowGaugeIndexList, &
         InflowGaugeHeadwater, &
         InflowGaugeNodeList, &
         InflowDischarge, & ! Intent IN
         L11_qOUT) ! Intent INOUT

    ! for a single node model run
    if( nNodes .GT. 1) then
       ! routing multiple times if timestep is smaller than 1
       !
       L11_qAcc = 0._dp
       do tt = 1, rout_loop
          ! routing of water within river reaches
          call L11_routing( nNodes, nNodes - L11_nOutlets, &
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
               L11_Qmod) ! Intent OUT
          ! accumulate values of individual subtimesteps
          L11_qAcc = L11_qAcc + L11_qMod
       end do
       ! calculate mean over routing period (timestep)
       L11_qMod = L11_qAcc / real(rout_loop, dp)
    else
       L11_Qmod = L11_qOUT 
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

  !     NAME
  !         L11_runoff_acc

  !     PURPOSE
  !>        \brief total runoff accumulation at L11.

  !>        \details Upscales runoff in space from L1 to L11 if routing resolution
  !>                 is higher than hydrology resolution (map_flag equals .true.) or
  !>                 downscales runoff from L1 to L11 if routing resolution is lower
  !>                 than hydrology resolution.

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "real(dp)    ::  qAll"              total runoff L1 [mm tst-1]
  !>        \param[in] "real(dp)    ::  efecArea"          effective area in [km2] at Level 1
  !>        \param[in] "integer(i4) ::  L1_L11_Id"         L11 Ids mapped on L1
  !>        \param[in] "real(dp)    ::  L11_areaCell"      effective area in [km2] at Level 11
  !>        \param[in] "integer(i4) ::  L11_L1_Id"         L1 Ids mapped on L11
  !>        \param[in] "integer(i4) ::  TS"                time step in [s]
  !>        \param[in] "logical     ::  map_flag"          Flag indicating whether routing resolution is higher than hydrologic one

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: qAcc"                 aggregated runoff at L11 [m3 s-1]

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
  !                   Stephan Thober, Sep 2015 - included downscaling of runoff
  !                   Stephan Thober, Feb 2016 - refactored upscaling of discharge from L1 to L11
  !                   Stephan Thober, Feb 2016 - refactored downscaling of discharge from L1 to L11

  ! ------------------------------------------------------------------
  SUBROUTINE L11_runoff_acc(qAll, efecArea, L1_L11_Id, &
       L11_areaCell, L11_L1_Id, &
       TS, map_flag, qAcc)

    use mo_mrm_constants, only: HourSecs, nodata_dp

    IMPLICIT NONE

    real(dp),    intent(in) :: qall(:) ! [mm tst-1] total runoff l1 
    real(dp),    intent(in) :: efecarea(:) ! [km2]      efective area at l1 
    integer(i4), intent(in) :: L1_L11_Id(:) ! l1 ids     mapped on l11
    real(dp),    intent(in) :: L11_areacell(:) ! [km2]      efective area at l11
    integer(i4), intent(in) :: L11_L1_Id(:) ! l11 ids    mapped on l1
    integer(i4), intent(in) :: TS ! [h] time step 
    logical,     intent(in) :: map_flag ! true when routing resolution is larger than hydrologic resolution
    real(dp),    intent(out) :: qAcc(:) ! [m3 s-1]   aggregated runoff at l11 

                                                              ! local variables
    integer(i4)             :: k
    real(dp)                :: TST             ! [s] time step

    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ------------------------------------------------------------------

    ! Hydrologic timestep in seconds
    TST = HourSecs * TS 

    if (map_flag) then
       ! Estimate specific runoff at  L11
       ! NOTE:
       ! 1) Total discharge depth aggregated at L11 level [mm/TST]
       ! 2) Transform  depth [mm/TST] to discharge [m3/s]
       ! Total runoff should be divided by total_area to get 
       ! specific discharge at L11. Then, to transform specific
       ! discharge from [mm/TST] to [m3/s], it should be multiplied by
       ! total_area [km2]*10^3 and divided by TST.
       ! Therefore, in this operation total_area cancels out. 
       qAcc = 0._dp
       ! loop over high-resolution cells (L1) and add discharge to
       ! corresponding low-resolution cells (L11)
       do k = 1, size(qAll, 1)
          qAcc(L1_L11_Id(k)) = qAcc(L1_L11_Id(k)) + qAll(k) * efecArea(k) 
       end do
       qAcc = qAcc * 1000.0_dp / TST
       !
    else
       ! initialize qout
       qAcc = nodata_dp
       do k = 1, size(qAcc, 1)
          ! map flux from coarse L1 resolution to fine L11 resolution
          qAcc(k) = qAll(L11_L1_Id(k))
       end do
       ! adjust flux by area cell
       qAcc(:) = qAcc(:) * L11_areaCell(:) * 1000.0_dp / TST
    end if

  END SUBROUTINE L11_runoff_acc

  ! ------------------------------------------------------------------

  !     NAME
  !         add_inflow

  !     PURPOSE 

  !>    \brief
  !>        Adds inflow discharge to the runoff produced at the
  !>        cell where the inflow is occurring.

  !>    \details

  !>        If a inflow gauge is given, then this routine is adding the
  !>        values to the runoff produced at the grid cell where the
  !>        inflow is happening. The values are not directly added to the
  !>        river network. If this cell is not a headwater then the streamflow
  !>        produced upstream will be neglected.
  
  !     INTENT(IN)
  !>      \param[in] "integer(i4)  :: nInflowGauges"   number of inflow gauges
  !>      \param[in] "integer(i4)  :: InflowIndexList" index of inflow points
  !>      \param[in] "logical      :: InflowHeadwater" flag to consider headwater cells of inflow gauge
  !>      \param[in] "integer(i4)  :: InflowNodeList"  L11 ID of inflow points
  !>      \param[in] "real(dp)     :: QInflow"         [m3 s-1] inflowing water

  !     INTENT(INOUT)
  !>      \param[in,out] "real(dp) :: qOut"            [m3 s-1] Series of attenuated runoff 
 
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
  !        none

  !     EXAMPLE
  !        none

  !     LITERATURE

  !     HISTORY
  !>        \author  Stephan Thober & Matthias Zink
  !>        \date    Jul 2016

  !         Modified 
  ! ------------------------------------------------------------------
  subroutine add_inflow(nInflowGauges, InflowIndexList, InflowHeadwater, &
       InflowNodeList, QInflow, & ! Intent IN
       qOut) ! Intent INOUT

    use mo_kind, only: i4, dp
    
    implicit none

    ! input variables
    integer(i4), intent(in)    :: nInflowGauges      ! [-] number of inflow points
    integer(i4), intent(in)    :: InflowIndexList(:) ! [-] index of inflow points
    logical,     intent(in)    :: InflowHeadwater(:) ! [-] if to consider headwater cells of inflow gauge
    integer(i4), intent(in)    :: InflowNodeList(:)  ! [-]        L11 ID of inflow points
    real(dp),    intent(in)    :: QInflow(:)         ! [m3 s-1]   inflowing water
    ! output variables
    real(dp),    intent(inout) :: qOut(:)            ! [m3 s-1] Series of attenuated runoff 

    ! local variables
    integer(i4) :: ii
    
    ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
    ! should be put after UH attenuation because it is measured runoff at this cell 
    if (nInflowGauges .gt. 0) then
       do ii = 1, nInflowGauges
          if (InflowHeadwater(ii)) then 
             ! add inflowing water to water produced by upstream/headwater cells
             qOut(InflowNodeList(ii)) = qOut(InflowNodeList(ii)) + QInflow(InflowIndexList(ii))
          else
             ! put only timeseries and cut upstream/headwater cells produced water for routing
             qOut(InflowNodeList(ii)) = QInflow(InflowIndexList(ii))
          end if
       end do
    end if
  end subroutine add_inflow
    
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
       nNodes, nLinks, netPerm, netLink_fromN, netLink_toN, &
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
    !                                                              !          node I (Muskingum)
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
    ! !!$OMP parallel
    ! store generated discharge
    netNode_Qmod(1:nNodes) = netNode_qTIN(1:nNodes,IT)
    ! backflow t-> t-1
    netNode_qTR(1:nNodes,IT1) = netNode_qTR(1:nNodes,IT)
    netNode_qTIN(1:nNodes,IT1)= netNode_qTIN(1:nNodes,IT)
    ! !!$OMP end parallel

  end subroutine L11_routing

END MODULE mo_mrm_routing
