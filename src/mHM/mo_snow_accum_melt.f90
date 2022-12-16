!> \file mo_snow_accum_melt.f90
!> \brief \copybrief mo_snow_accum_melt
!> \details \copydetails mo_snow_accum_melt

!> \brief Snow melting and accumulation.
!> \details This module calculates snow melting and accumulation.
!> \authors Vladyslav Prykhodko
!> \date Dec 2012
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

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: deg_day_incr"       Increase of the Degree-day factor per mm of increasein
  !>       precipitation [s-1 degreeC-1]
  !>       \param[in] "REAL(dp) :: deg_day_max"        Maximum Degree-day factor [m-1 degreeC-1]
  !>       \param[in] "REAL(dp) :: deg_day_noprec"     Degree-day factor with no precipitation [m-1 degreeC-1]
  !>       \param[in] "REAL(dp) :: prec"               Daily mean precipitation [m]
  !>       \param[in] "REAL(dp) :: temperature"        Daily mean temperature [degreeC]
  !>       \param[in] "REAL(dp) :: temperature_thresh" Threshold temperature for snow/rain [degreeC]
  !>       \param[in] "REAL(dp) :: thrfall"            Throughfall [m TS-1]

  !    INTENT(INOUT)
  !>       \param[inout] "REAL(dp) :: snow_pack" Snow pack [m]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: deg_day"     Degree-day factor  [m s-1 degreeC-1]
  !>       \param[out] "REAL(dp) :: melt"        Melting snow depth [m TS-1]
  !>       \param[out] "REAL(dp) :: prec_effect" Effective precipitation depth (snow melt + rain) [m]
  !>       \param[out] "REAL(dp) :: rain"        Rain precipitation depth [m]
  !>       \param[out] "REAL(dp) :: snow"        Snow precipitation depth [m]

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! JM Aug 2013 - ordering of arguments changed
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE snow_accum_melt(deg_day_incr, deg_day_max, deg_day_noprec, prec, temperature, temperature_thresh, thrfall, &
                            snow_pack, deg_day, melt, prec_effect, rain, snow)
    implicit none

    ! Increase of the Degree-day factor per mm of increasein precipitation [s-1 degreeC-1]
    REAL(dp), INTENT(IN) :: deg_day_incr

    ! Maximum Degree-day factor [m-1 degreeC-1]
    REAL(dp), INTENT(IN) :: deg_day_max

    ! Degree-day factor with no precipitation [m-1 degreeC-1]
    REAL(dp), INTENT(IN) :: deg_day_noprec

    ! Daily mean precipitation [m]
    REAL(dp), INTENT(IN) :: prec

    ! Daily mean temperature [degreeC]
    REAL(dp), INTENT(IN) :: temperature

    ! Threshold temperature for snow/rain [degreeC]
    REAL(dp), INTENT(IN) :: temperature_thresh

    ! Throughfall [m TS-1]
    REAL(dp), INTENT(IN) :: thrfall

    ! Snow pack [m]
    REAL(dp), INTENT(INOUT) :: snow_pack

    ! Degree-day factor  [m s-1 degreeC-1]
    REAL(dp), INTENT(OUT) :: deg_day

    ! Melting snow depth [m TS-1]
    REAL(dp), INTENT(OUT) :: melt

    ! Effective precipitation depth (snow melt + rain) [m]
    REAL(dp), INTENT(OUT) :: prec_effect

    ! Rain precipitation depth [m]
    REAL(dp), INTENT(OUT) :: rain

    ! Snow precipitation depth [m]
    REAL(dp), INTENT(OUT) :: snow

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
