!> \file mo_snow_accum_melt.f90
!> \copydoc mo_snow_accum_melt

!> \brief Snow melting and accumulation.
!> \details This module calculates snow melting and accumulation.
!> \authors Vladyslav Prykhodko
!> \date Dec 2012
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_snow_accum_melt

  USE mo_kind, ONLY : dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: snow_accum_melt

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        snow_accum_melt

  !    PURPOSE
  !>       \brief Snow melting and accumulation.

  !>       \details Separates throughfall into rain and snow by comparing the temperature with the treshhold.
  !>       by comparing the temperature with the treshhold.
  !>       Calculates degree daily factor.
  !>       Calculates snow melting rates.
  !>       Calculates snow, rain and effective precipitation depth
  !>       and snow pack.

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! JM Aug 2013 - ordering of arguments changed
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE snow_accum_melt(deg_day_incr, deg_day_max, deg_day_noprec, prec, temperature, temperature_thresh, thrfall, &
                            snow_pack, deg_day, melt, prec_effect, rain, snow)
    implicit none
    REAL(dp), INTENT(IN) :: deg_day_incr !< Increase of the Degree-day factor per mm of increase in precipitation [s-1 degreeC-1]
    REAL(dp), INTENT(IN) :: deg_day_max !< Maximum Degree-day factor [m-1 degreeC-1]
    REAL(dp), INTENT(IN) :: deg_day_noprec !< Degree-day factor with no precipitation [m-1 degreeC-1]
    REAL(dp), INTENT(IN) :: prec !< Daily mean precipitation [m]
    REAL(dp), INTENT(IN) :: temperature !< Daily mean temperature [degreeC]
    REAL(dp), INTENT(IN) :: temperature_thresh !< Threshold temperature for snow/rain [degreeC]
    REAL(dp), INTENT(IN) :: thrfall !< Throughfall [m TS-1]
    REAL(dp), INTENT(INOUT) :: snow_pack !< Snow pack [m]
    REAL(dp), INTENT(OUT) :: deg_day !< Degree-day factor  [m s-1 degreeC-1]
    REAL(dp), INTENT(OUT) :: melt !< Melting snow depth [m TS-1]
    REAL(dp), INTENT(OUT) :: prec_effect !< Effective precipitation depth (snow melt + rain) [m]
    REAL(dp), INTENT(OUT) :: rain !< Rain precipitation depth [m]
    REAL(dp), INTENT(OUT) :: snow !< Snow precipitation depth [m]

    ! Auxiliary helping variable [-]
    REAL(dp) :: aux_help

    !separate throughfall into rain and snow
    if(temperature >  temperature_thresh) then
      snow = 0.0_dp
      rain = thrfall
    else
      snow = thrfall
      rain = 0.0_dp
    end if

    ! calculate degree daily factor
    if (prec <= (deg_day_max - deg_day_noprec) / deg_day_incr) then
      deg_day = deg_day_noprec + deg_day_incr * prec
    else
      deg_day = deg_day_max
    end if

    ! melting/snow accumulation
    if (temperature > temperature_thresh) then
      ! melting
      if (snow_pack > 0.0_dp) then
        aux_help = deg_day * (temperature - temperature_thresh)
        if (aux_help > snow_pack) then
          melt = snow_pack
          snow_pack = 0.0_dp
        else
          melt = aux_help
          snow_pack = snow_pack - aux_help
        end if
      else
        melt = 0.0_dp
        snow_pack = 0.0_dp
      end if
    else
      ! snow accumulation
      melt = 0.0_dp
      snow_pack = snow_pack + snow
    end if

    ! effective precipitation
    prec_effect = melt + rain

  END SUBROUTINE snow_accum_melt

END MODULE mo_snow_accum_melt
