!> \file mo_runoff.f90
!> \copydoc mo_runoff

!> \brief Runoff generation.
!> \details This module generates the runoff for the unsaturated and saturated zones and provides runoff accumulation.
!> \authors Vladyslav Prykhodko
!> \date Dec 2012
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_runoff

  USE mo_kind, ONLY : dp
  USE mo_constants, ONLY : eps_dp

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
    REAL(dp), INTENT(IN) :: k1 !< Recession coefficient of the upper reservoir,lower outlet [TS-1]
    REAL(dp), INTENT(IN) :: kp !< Percolation coefficient [TS-1]
    REAL(dp), INTENT(IN) :: k0 !< Recession coefficient of the upperreservoir, upper outlet [TS-1]
    REAL(dp), INTENT(IN) :: alpha !< Exponent for the upper reservoir [-]
    REAL(dp), INTENT(IN) :: karst_loss !< Karstic percolation loss [-]
    REAL(dp), INTENT(IN) :: pefec_soil !< Input to the soil layer [mm]
    REAL(dp), INTENT(IN) :: unsat_thresh !< Threshold water depth in upper reservoir(for Runoff contribution) [mm]
    REAL(dp), INTENT(INOUT) :: sat_storage !< Groundwater storage [mm]
    REAL(dp), INTENT(INOUT) :: unsat_storage !< Upper soil storage [mm]
    REAL(dp), INTENT(OUT) :: slow_interflow !< Slow runoff component [mm TS-1]
    REAL(dp), INTENT(OUT) :: fast_interflow !< Fast runoff component [mm TS-1]
    REAL(dp), INTENT(OUT) :: perc !< Percolation [mm TS-1]

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

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! JM Aug 2013 - ordering of arguments changed
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE runoff_sat_zone(k2, sat_storage, baseflow)
    implicit none
    REAL(dp), INTENT(IN) :: k2 !< Baseflow recession coefficient [TS-1]
    REAL(dp), INTENT(INOUT) :: sat_storage !< Groundwater storage [mm]
    REAL(dp), INTENT(OUT) :: baseflow !< Baseflow [mm TS-1]

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
    REAL(dp), INTENT(IN) :: fSealed_area_fraction !< sealed area fraction [1]
    REAL(dp), INTENT(IN) :: fast_interflow !< \f$ q_0 \f$ Fast runoff component [mm TS-1]
    REAL(dp), INTENT(IN) :: slow_interflow !< \f$ q_1 \f$ Slow runoff component [mm TS-1]
    REAL(dp), INTENT(IN) :: baseflow !< \f$ q_2 \f$ Baseflow [mm TS-1]
    REAL(dp), INTENT(IN) :: direct_runoff !< \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]
    REAL(dp), INTENT(OUT) :: total_runoff !< \f$ q_T \f$ Generated runoff [mm TS-1]

    total_runoff = ((baseflow + slow_interflow + fast_interflow) * (1.0_dp - fSealed_area_fraction)) + &
            (direct_runoff * fSealed_area_fraction)

  END SUBROUTINE L1_total_runoff

END MODULE mo_runoff
