!> \file mo_runoff.f90
!> \brief \copybrief mo_runoff
!> \details \copydetails mo_runoff

!> \brief Runoff generation.
!> \details This module generates the runoff for the unsaturated and saturated zones and provides runoff accumulation.
!> \authors Vladyslav Prykhodko
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_runoff

  USE mo_kind, ONLY : dp
  USE mo_common_constants, ONLY : eps_dp

  IMPLICIT NONE

  PUBLIC :: runoff_unsat_zone
  PUBLIC :: runoff_sat_zone
  PUBLIC :: L1_total_runoff

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        runoff_unsat_zone

  !    PURPOSE
  !>       \brief Runoff generation for the saturated zone.

  !>       \details Calculates the runoff generation for the unsaturated zone.
  !>       Calculates percolation, interflow and baseflow.
  !>       Updates upper soil and groundwater storages.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: k1"           Recession coefficient of the upper reservoir,lower outlet [TS-1]
  !>       \param[in] "REAL(dp) :: kp"           Percolation coefficient [TS-1]
  !>       \param[in] "REAL(dp) :: k0"           Recession coefficient of the upperreservoir, upper outlet [TS-1]
  !>       \param[in] "REAL(dp) :: alpha"        Exponent for the upper reservoir [-]
  !>       \param[in] "REAL(dp) :: karst_loss"   Karstic percolation loss [-]
  !>       \param[in] "REAL(dp) :: pefec_soil"   Input to the soil layer [mm]
  !>       \param[in] "REAL(dp) :: unsat_thresh" Threshold water depth in upper reservoir(for Runoff contribution) [mm]

  !    INTENT(INOUT)
  !>       \param[inout] "REAL(dp) :: sat_storage"   Groundwater storage [mm]
  !>       \param[inout] "REAL(dp) :: unsat_storage" Upper soil storage [mm]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: slow_interflow" Slow runoff component [mm TS-1]
  !>       \param[out] "REAL(dp) :: fast_interflow" Fast runoff component [mm TS-1]
  !>       \param[out] "REAL(dp) :: perc"           Percolation [mm TS-1]

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! LS Feb 2006 - fast response
  ! LS Feb 2007 - MaxInter
  ! RK May 2007 - fracArea, errors in Qmod
  ! LS Dec 2012 - variable names and process sat. zone
  ! LS Jan 2013 - total runoff accumulation L11
  ! JM Aug 2013 - ordering of arguments changed
  ! Robert Schweppe Jun 2018 - refactoring and reformatting


  SUBROUTINE runoff_unsat_zone(k1, kp, k0, alpha, karst_loss, pefec_soil, unsat_thresh, sat_storage, unsat_storage, &
                              slow_interflow, fast_interflow, perc)
    implicit none

    ! Recession coefficient of the upper reservoir,lower outlet [TS-1]
    REAL(dp), INTENT(IN) :: k1

    ! Percolation coefficient [TS-1]
    REAL(dp), INTENT(IN) :: kp

    ! Recession coefficient of the upperreservoir, upper outlet [TS-1]
    REAL(dp), INTENT(IN) :: k0

    ! Exponent for the upper reservoir [-]
    REAL(dp), INTENT(IN) :: alpha

    ! Karstic percolation loss [-]
    REAL(dp), INTENT(IN) :: karst_loss

    ! Input to the soil layer [mm]
    REAL(dp), INTENT(IN) :: pefec_soil

    ! Threshold water depth in upper reservoir(for Runoff contribution) [mm]
    REAL(dp), INTENT(IN) :: unsat_thresh

    ! Groundwater storage [mm]
    REAL(dp), INTENT(INOUT) :: sat_storage

    ! Upper soil storage [mm]
    REAL(dp), INTENT(INOUT) :: unsat_storage

    ! Slow runoff component [mm TS-1]
    REAL(dp), INTENT(OUT) :: slow_interflow

    ! Fast runoff component [mm TS-1]
    REAL(dp), INTENT(OUT) :: fast_interflow

    ! Percolation [mm TS-1]
    REAL(dp), INTENT(OUT) :: perc


    !---------------------------------------------------------------
    ! SOIL LAYER BETWEEN UNSATURATED AND SATURATED ZONE
    !---------------------------------------------------------------
    ! HERE input is from last soil Horizon...
    !pefec_soil = Cell1_soil(k, nHorizons_mHM)%infil
    unsat_storage = unsat_storage + pefec_soil

    ! FAST INTERFLOW WITH THRESHOLD BEHAVIOUR
    fast_interflow = 0.0_dp
    if(unsat_storage > unsat_thresh) then
      fast_interflow = MIN((k0 * (unsat_storage - unsat_thresh)), &
              (unsat_storage - eps_dp))
    end if
    unsat_storage = unsat_storage - fast_interflow

    ! SLOW PERMANENT INTERFLOW
    slow_interflow = 0.0_dp

    if(unsat_storage > eps_dp) then
      slow_interflow = min((k1 * (unsat_storage**(1.0_dp + alpha))), &
              (unsat_storage - eps_dp))
    end if
    unsat_storage = unsat_storage - slow_interflow

    !--------------------------------------------------------
    ! PERCOLATION FROM SOIL LAYER TO THE SATURATED ZONE
    !--------------------------------------------------------
    perc = kp * unsat_storage

    ! Taking into account for the KARSTIC aquifers
    !*** karstic loss gain or loss if Karstic aquifer is present...
    if(unsat_storage > perc) then
      unsat_storage = unsat_storage - perc
      sat_storage = sat_storage + perc * karst_loss
    else
      sat_storage = sat_storage + unsat_storage * karst_loss
      unsat_storage = 0.0_dp
    end if

  END SUBROUTINE runoff_unsat_zone

  ! ------------------------------------------------------------------

  !    NAME
  !        runoff_sat_zone

  !    PURPOSE
  !>       \brief Runoff generation for the saturated zone.

  !>       \details Calculates the runoff generation for the saturated zone.
  !>       If the level of the ground water reservoir is zero, then
  !>       the baseflow is also zero.
  !>       If the level of the ground water reservoir is greater than zero, then
  !>       the baseflow is equal to baseflow recession coefficient times the level
  !>       of the ground water reservoir, which
  !>       will be then reduced by the value of baseflow.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: k2" Baseflow recession coefficient [TS-1]

  !    INTENT(INOUT)
  !>       \param[inout] "REAL(dp) :: sat_storage" Groundwater storage [mm]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: baseflow" Baseflow [mm TS-1]

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! JM Aug 2013 - ordering of arguments changed
  ! Robert Schweppe Jun 2018 - refactoring and reformatting


  SUBROUTINE runoff_sat_zone(k2, sat_storage, baseflow)
    implicit none

    ! Baseflow recession coefficient [TS-1]
    REAL(dp), INTENT(IN) :: k2

    ! Groundwater storage [mm]
    REAL(dp), INTENT(INOUT) :: sat_storage

    ! Baseflow [mm TS-1]
    REAL(dp), INTENT(OUT) :: baseflow


    if (sat_storage > 0.0_dp) then
      baseflow = k2 * sat_storage
      sat_storage = sat_storage - baseflow
    else
      baseflow = 0.0_dp
      sat_storage = 0.0_dp
    end if

  END SUBROUTINE runoff_sat_zone


  ! ------------------------------------------------------------------

  !    NAME
  !        L1_total_runoff

  !    PURPOSE
  !>       \brief total runoff accumulation at level 1

  !>       \details Accumulates runoff.
  !>       \f[ q_{T} = ( q_0 + q_1 + q_2 ) * (1-fSealed) + q_{D} * fSealed \f],
  !>       where fSealed is the fraction of sealed area.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: fSealed_area_fraction" sealed area fraction [1]
  !>       \param[in] "REAL(dp) :: fast_interflow"        \f$ q_0 \f$ Fast runoff component [mm TS-1]
  !>       \param[in] "REAL(dp) :: slow_interflow"        \f$ q_1 \f$ Slow runoff component [mm TS-1]
  !>       \param[in] "REAL(dp) :: baseflow"              \f$ q_2 \f$ Baseflow [mm TS-1]
  !>       \param[in] "REAL(dp) :: direct_runoff"         \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: total_runoff" \f$ q_T \f$ Generated runoff [mm TS-1]

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! RK Jul 2013 - A Mosiac approach is implemented for processes accounted within the permeamble & impervious area.
  ! ST May 2015 - updated equation in the documentation
  ! Robert Schweppe Jun 2018 - refactoring and reformatting


  SUBROUTINE L1_total_runoff(fSealed_area_fraction, fast_interflow, slow_interflow, baseflow, direct_runoff, &
                            total_runoff)
    implicit none

    ! sealed area fraction [1]
    REAL(dp), INTENT(IN) :: fSealed_area_fraction

    ! \f$ q_0 \f$ Fast runoff component [mm TS-1]
    REAL(dp), INTENT(IN) :: fast_interflow

    ! \f$ q_1 \f$ Slow runoff component [mm TS-1]
    REAL(dp), INTENT(IN) :: slow_interflow

    ! \f$ q_2 \f$ Baseflow [mm TS-1]
    REAL(dp), INTENT(IN) :: baseflow

    ! \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]
    REAL(dp), INTENT(IN) :: direct_runoff

    ! \f$ q_T \f$ Generated runoff [mm TS-1]
    REAL(dp), INTENT(OUT) :: total_runoff


    total_runoff = ((baseflow + slow_interflow + fast_interflow) * (1.0_dp - fSealed_area_fraction)) + &
            (direct_runoff * fSealed_area_fraction)

  END SUBROUTINE L1_total_runoff

END MODULE mo_runoff
