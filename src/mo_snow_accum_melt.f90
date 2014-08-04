!> \file mo_snow_accum_melt.f90

!> \brief Snow melting and accumulation.

!> \details This module calculates snow melting and accumulation.

!> \authors Vladyslav Prykhodko
!> \date Dec 2012

MODULE mo_snow_accum_melt

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::   snow_accum_melt

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         snow_accum_melt

  !     PURPOSE
  !>        \brief Snow melting and accumulation.

  !>        \details Separates throughfall into rain and snow by comparing the temperature with the treshhold.
  !>                  by comparing the temperature with the treshhold.
  !>                  Calculates degree daily factor.
  !>                  Calculates snow melting rates.
  !>                  Calculates snow, rain and effective precipitation depth
  !>                  and snow pack.

  !     CALLING SEQUENCE
  !         snow_accum_melt(deg_day_incr, deg_day_max, deg_day_noprec, prec, temperature, temperature_tresh, thrfall, &
  !                         deg_day, snow_pack, melt,  prec_effect, rain, snow)

  !     INTENT(IN)

  !>        \param[in] "real(dp) ::  deg_day_incr"         Increase of the Degree-day factor per mm of increase
  !>                                                       in precipitation [s-1 degreeC-1]
  !>        \param[in] "real(dp) ::  deg_day_max"          Maximum Degree-day factor [m-1 degreeC-1]
  !>        \param[in] "real(dp) ::  deg_day_noprec"       Degree-day factor with no precipitation [m-1 degreeC-1]
  !>        \param[in] "real(dp) ::  prec"                 Daily mean precipitation [m]
  !>        \param[in] "real(dp) ::  temperature"               Daily mean temperature [degreeC]
  !>        \param[in] "real(dp) ::  temperature_thresh"        Threshold temperature for snow/rain [degreeC]
  !>        \param[in] "real(dp) ::  thrfall"              Throughfall [m s-1]

  !     INTENT(INOUT)
  !>        \param[in,out] "real(dp) ::  snow_pack"        Snow pack [m]

  !     INTENT(OUT)
  !>        \param[out] "real(dp) ::  deg_day"             Degree-day factor  [m s-1 degreeC-1]
  !>        \param[out] "real(dp) ::  melt"                Melting snow depth [m s-1]
  !>        \param[out] "real(dp) ::  prec_effect"         Effective precipitation depth (snow melt + rain) [m]
  !>        \param[out] "real(dp) ::  rain"                Rain precipitation depth [m]
  !>        \param[out] "real(dp) ::  snow"                Snow precipitation depth [m]

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
  !>        \author Vladyslav Prykhodko
  !>        \date Dec 2012
  !         Modified JM, Aug 2013 - ordering of arguments changed

  SUBROUTINE snow_accum_melt(deg_day_incr, deg_day_max, deg_day_noprec, prec, &
       temperature, temperature_thresh, thrfall, &
       snow_pack, deg_day, melt,  prec_effect, rain, snow)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)     :: deg_day_incr
    REAL(dp), INTENT(IN)     :: deg_day_max
    REAL(dp), INTENT(IN)     :: deg_day_noprec
    REAL(dp), INTENT(IN)     :: prec
    REAL(dp), INTENT(IN)     :: temperature
    REAL(dp), INTENT(IN)     :: temperature_thresh
    REAL(dp), INTENT(IN)     :: thrfall
    REAL(dp), INTENT(INOUT)  :: snow_pack
    REAL(dp), INTENT(OUT)    :: deg_day
    REAL(dp), INTENT(OUT)    :: melt
    REAL(dp), INTENT(OUT)    :: prec_effect
    REAL(dp), INTENT(OUT)    :: rain
    REAL(dp), INTENT(OUT)    :: snow

    ! local variables
    REAL(dp)                 :: aux_help          !Auxiliary helping variable [-]

    !separate throughfall into rain and snow
    if(temperature >  temperature_thresh) then
       snow = 0.0_dp
       rain = thrfall
    else
       snow = thrfall
       rain = 0.0_dp
    end if

    ! calculate degree daily factor
    if ( prec <= (deg_day_max - deg_day_noprec) / deg_day_incr ) then
       deg_day = deg_day_noprec + deg_day_incr * prec
    else
       deg_day = deg_day_max
    end if

    ! melting/snow accumulation
    if (temperature > temperature_thresh) then
       ! melting
       if ( snow_pack > 0.0_dp ) then
          aux_help = deg_day * ( temperature - temperature_thresh )
          if ( aux_help > snow_pack ) then
             melt      = snow_pack
             snow_pack = 0.0_dp
          else
             melt      = aux_help
             snow_pack = snow_pack - aux_help
          end if
       else
          melt      = 0.0_dp
          snow_pack = 0.0_dp
       end if
    else
       ! snow accumulation
       melt      = 0.0_dp
       snow_pack = snow_pack + snow
    end if

    ! effective precipitation
    prec_effect = melt + rain

  END SUBROUTINE snow_accum_melt

END MODULE mo_snow_accum_melt
