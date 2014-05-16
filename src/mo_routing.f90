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

  PUBLIC :: L11_routing      ! single routine

  ! ------------------------------------------------------------------

CONTAINS

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
  !>    required. This permutation is determined in the mo_net_startup
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
       netLink_C1, netLink_C2, netNode_qOUT, netNode_qTIN, netNode_qTR, netNode_Qmod )

    implicit none

    ! Input
    integer(i4),                      intent(in)    :: nNodes        ! number of network nodes = nCells1
    integer(i4),                      intent(in)    :: nLinks        ! number of stream segment (reaches)
    ! Stream link description   network topology  ==>    netLink
    integer(i4), dimension(:),        intent(in)    :: netPerm       ! basin routing order (permutation)
    integer(i4), dimension(:),        intent(in)    :: netLink_fromN ! from node 
    integer(i4), dimension(:),        intent(in)    :: netLink_toN   ! to node
    real(dp),    dimension(:),        intent(in)    :: netLink_C1    ! [1]   routing parameter  C1 (Chow, 25-41)
    real(dp),    dimension(:),        intent(in)    :: netLink_C2    ! [1]   routing parameters C2 (")
    ! State variables
    real(dp),    dimension(:),        intent(in)    :: netNode_qOUT  ! [m3 s-1] Total outflow, all cells, basin, 
    !                                                                !          level L11 at time tt 
    ! Input - Output
    real(dp),    dimension(:,:), intent(inout)      :: netNode_qTIN  ! [m3 s-1] Total inputs at t-1 and t
    real(dp),    dimension(:,:), intent(inout)      :: netNode_qTR   ! [m3 s-1] Transformed outflow leaving 
    !                                                                !          node I (Muskingum)
    ! Output
    real(dp),    dimension(nNodes),   intent(out)   :: netNode_Qmod  ! [m3 s-1] Simulated routed discharge  


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

       ! add routing to tNode
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

END MODULE mo_routing
