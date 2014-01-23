!> \file mo_runoff.f90

!> \brief Runoff generation for the  unsaturated zone, saturated zone (or groundwater zone),
!>        and runoff accumulation.

!> \details This module generates the runoff for the unsaturated and saturated zones and provides
!>          runoff accumulation.

!> \authors Vladyslav Prykhodko
!> \date Dec 2012

MODULE mo_runoff
  
  USE mo_kind, ONLY: i4, dp
  USE mo_constants, ONLY: eps_dp

  IMPLICIT NONE

  PUBLIC ::   runoff_unsat_zone
  PUBLIC ::   runoff_sat_zone
  PUBLIC ::   L1_total_runoff
  PUBLIC ::   L11_runoff_acc

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         runoff_unsat_zone

  !  PURPOSE
  !>        \brief Runoff generation for the saturated zone.

  !>        \details Calculates the runoff generation for the unsaturated zone.\n
  !>                 Calculates percolation, interflow and baseflow.\n
  !>                 Updates upper soil and groundwater storages.\n

  !     CALLING SEQUENCE
  !         runoff_unsat_zone(k1, kp, coeff_up, alpha, karst_loss, pefec_soil, &
  !                           unsat_thresh, slow_interflow, perc, fast_interflow, &
  !                           sat_storage, unsat_storage )

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: pefec_soil"         Input to the soil layer [mm]
  !>        \param[in] "real(dp) :: k0"                 Recession coefficient of the upper
  !>                                                    reservoir, upper outlet [d-1]
  !>        \param[in] "real(dp) :: k1"                 Recession coefficient of the upper reservoir,
  !>                                                    lower outlet [d-1]
  !>        \param[in] "real(dp) :: kp"                 Percolation coefficient [d-1]
  !>        \param[in] "real(dp) :: alpha"              Exponent for the upper reservoir [-]
  !>        \param[in] "real(dp) :: unsat_thresh"       Threshold water depth in upper reservoir
  !>                                                    (for Runoff contribution) [mm]
  !>        \param[in] "real(dp) :: karst_loss"         Karstic percolation loss [-]

  !     INTENT(INOUT)
  !>        \param[in,out] "real(dp) :: unsat_storage   Upper soil storage [mm]
  !>        \param[in,out] "real(dp) :: sat_storage     Groundwater storage [mm]

  !     INTENT(OUT)
  !>        \param[out] "real(dp) ::  perc              Percolation [mm d-1]
  !>        \param[out] "real(dp) ::  fast_interflow    Fast runoff component [mm d-1]
  !>        \param[out] "real(dp) ::  slow_interflow    Slow runoff component [mm d-1]

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
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Vladyslav Prykhodko
  !>        \date    Dec 2012

  !         Created  LS, Dec 2005
  !         Modified LS, Feb 2006 - fast response
  !                  LS, Feb 2007 - MaxInter
  !                  RK, May 2007 - fracArea, errors in Qmod
  !                  LS, Dec 2012 - variable names and process sat. zone
  !                  LS, Jan 2013 - total runoff accumulation L11
  !                  JM, Aug 2013  - ordering of arguments changed
  ! ------------------------------------------------------------------

  SUBROUTINE runoff_unsat_zone(k1, kp, k0, alpha,            &   ! Intent IN
       karst_loss, pefec_soil, unsat_thresh,                 &   ! Intent IN
       sat_storage, unsat_storage,                           &   ! Intent INOUT
       slow_interflow, fast_interflow, perc)                     ! Intent OUT

    IMPLICIT NONE

    REAL(dp), INTENT(IN)    ::  k1
    REAL(dp), INTENT(IN)    ::  kp
    REAL(dp), INTENT(IN)    ::  k0
    REAL(dp), INTENT(IN)    ::  alpha
    REAL(dp), INTENT(IN)    ::  karst_loss
    REAL(dp), INTENT(IN)    ::  pefec_soil
    REAL(dp), INTENT(IN)    ::  unsat_thresh
    REAL(dp), INTENT(INOUT) ::  sat_storage
    REAL(dp), INTENT(INOUT) ::  unsat_storage
    REAL(dp), INTENT(OUT)   ::  slow_interflow
    REAL(dp), INTENT(OUT)   ::  fast_interflow
    REAL(dp), INTENT(OUT)   ::  perc

    !---------------------------------------------------------------
    ! SOIL LAYER BETWEEN UNSATURATED AND SATURATED ZONE
    !---------------------------------------------------------------
    ! HERE input is from last soil Horizon...
    !pefec_soil = Cell1_soil(k, nHorizons_mHM)%infil

    unsat_storage = unsat_storage + pefec_soil
   
    ! FAST INTERFLOW WITH THRESHOLD BEHAVIOUR
    fast_interflow   = 0.0_dp
    if(unsat_storage > unsat_thresh) then
       fast_interflow   = MIN( (k0*(unsat_storage - unsat_thresh)), &
            (unsat_storage - eps_dp))
    end if
    unsat_storage = unsat_storage -  fast_interflow
   
    ! SLOW PERMANENT INTERFLOW
    slow_interflow= 0.0_dp
   
    if( unsat_storage > eps_dp ) then
       slow_interflow= min( (k1*(unsat_storage**(1.0_dp + alpha) )), &
            (unsat_storage - eps_dp))
    end if
    unsat_storage = unsat_storage - slow_interflow
   
    !--------------------------------------------------------
    ! PERCOLATION FROM SOIL LAYER TO THE SATURATED ZONE
    !--------------------------------------------------------
    perc  = kp * unsat_storage
   
    ! Taking into account for the KARSTIC aquifers
    !*** karstic loss gain or loss if Karstic aquifer is present...
    if(unsat_storage > perc ) then
       unsat_storage = unsat_storage - perc
       sat_storage = sat_storage + perc * karst_loss
    else
       sat_storage = sat_storage + unsat_storage * karst_loss
       unsat_storage = 0.0_dp
    end if 
   
  END SUBROUTINE runoff_unsat_zone

  ! ------------------------------------------------------------------

  !     NAME
  !         runoff_sat_zone

  !     PURPOSE
  !>        \brief Runoff generation for the saturated zone.

  !>        \details Calculates the runoff generation for the saturated zone.
  !>                 If the level of the ground water reservoir is zero, then
  !>                 the baseflow is also zero.
  !>                 If the level of the ground water reservoir is greater than zero, then
  !>                 the baseflow is equal to baseflow recession coefficient times the level
  !>                 of the ground water reservoir, which
  !>                 will be then reduced by the value of baseflow.

  !     CALLING SEQUENCE
  !         runoff_sat_zone(k2, baseflow, sat_storage)

  !     INTENT(IN)
  !>        \param[in] "real(dp) ::  k2"                   Baseflow recession coefficient [d-1]

  !     INTENT(INOUT)
  !>        \param[in,out] "real(dp) ::  sat_storage"      Groundwater storage [mm]

  !     INTENT(OUT)
  !>        \param[out] "real(dp) ::  baseflow"            Baseflow [mm d-1]

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
  !         K=0.1
  !         B=5.0
  !         S=10.0
  !         runoff_sat_zone(K,B,S)

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Vladyslav Prykhodko
  !>        \date Dec 2012
  !         Modified JM, Aug 2013  - ordering of arguments changed

  SUBROUTINE runoff_sat_zone(k2, sat_storage, baseflow)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)    :: k2
    REAL(dp), INTENT(INOUT) :: sat_storage
    REAL(dp), INTENT(OUT)   :: baseflow

    if (sat_storage > 0.0_dp) then
       baseflow   = k2 * sat_storage
       sat_storage = sat_storage - baseflow
    else
       baseflow   = 0.0_dp
       sat_storage = 0.0_dp
    end if

  END SUBROUTINE runoff_sat_zone


  ! ------------------------------------------------------------------

  !     NAME
  !        L1_total_runoff

  !     PURPOSE
  !>        \brief total runoff accumulation at level 1

  !>        \details Accumulates runoff.
  !>        \f[ q_{T} = q_0 + q_1 + q_2 + q_{D} \f]

  !     CALLING SEQUENCE
  !         runoff_accum(fSealed_area_fraction, fast_interflow, slow_interflow, baseflow,  direct_runoff, total_runoff)

  !     INTENT(IN)
  !>        \param[in] "real(dp)  ::  fSealed_area_fraction      sealed area fraction [1]
  !>        \param[in] "real(dp)  ::  fast_interflow            \f$ q_0 \f$ Fast runoff component [mm tst-1]
  !>        \param[in] "real(dp)  ::  slow_interflow            \f$ q_1 \f$ Slow runoff component [mm tst-1] 
  !>        \param[in] "real(dp)  ::  baseflow"                 \f$ q_2 \f$ Baseflow [mm tsts-1]   
  !>        \param[in] "real(dp)  ::  direct_runoff"            \f$ q_D \f$ Direct runoff from impervious areas  [mm tst-1]

  !     INTENT(INOUT)
  !        None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: total_runoff"              \f$ q_T \f$ Generated runoff [mm tst-1]

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
  !>        \author Vladyslav Prykhodko
  !>        \date Dec 2012
  !         Modified  RK, Jul 2013  - A Mosiac approach is implemented for processes accounted
  !                                   within the permeamble & impervious area. 

  SUBROUTINE L1_total_runoff(fSealed_area_fraction, fast_interflow, slow_interflow, &
                             baseflow, direct_runoff, total_runoff)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)    :: fSealed_area_fraction
    REAL(dp), INTENT(IN)    :: fast_interflow
    REAL(dp), INTENT(IN)    :: slow_interflow
    REAL(dp), INTENT(IN)    :: baseflow
    REAL(dp), INTENT(IN)    :: direct_runoff
    REAL(dp), INTENT(OUT)   :: total_runoff

    total_runoff = ( (baseflow + slow_interflow + fast_interflow)*(1.0_dp-fSealed_area_fraction) ) + &
                   ( direct_runoff*fSealed_area_fraction )

  END SUBROUTINE L1_total_runoff
  
 
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


  ! ------------------------------------------------------------------
  SUBROUTINE L11_runoff_acc(qAll,efecArea, L11id, TS, qOUT)
    use mo_mhm_constants, only:   HourSecs

    IMPLICIT NONE

    real(dp),    dimension(:), intent(in)  :: qall       ! [mm tst-1] total runoff l1 
    real(dp),    dimension(:), intent(in)  :: efecarea   ! [km2]      efective area at l1 
    integer(i4), dimension(:), intent(in)  :: l11id      ! l11        mapped on l1   
    integer(i4),               intent(in)  :: ts         ! [h]        time step 
    real(dp),    dimension(:), intent(out) :: qout       ! [m3 s-1]   aggregated runoff at l11 

    ! local variables
    integer(i4)                            :: k
    REAL(dp)                               :: TST        ! [s]        time step  

    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ****  NO TUH :: OFF AT A MOMENT BUT EFFECTS ARE
    !                 INTEGRATED IN RIVER ROUTING PROCESS
    !          TUH    Triangular Unit Hydrograph
    ! ------------------------------------------------------------------

    TST = HourSecs * real( TS, dp)   ! in [s]

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
      qOUT(k) = sum( qAll(:) * efecArea(:),  L11id(:) == k ) * 1000.0_dp / TST
    end do

  END SUBROUTINE L11_runoff_acc
  
END MODULE mo_runoff
