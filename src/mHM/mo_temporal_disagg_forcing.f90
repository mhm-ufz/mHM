!> \file mo_temporal_disagg_forcing.f90

!> \brief Temporal disaggregation of daily input values

!> \details Calculate actual values for precipitation, PET and temperature from daily mean inputs

!>  \note There is not PET correction for aspect in this routine. Use pet * fasp before or after the routine.

!> \authors Matthias Cuntz
!> \date Dec 2012

MODULE mo_temporal_disagg_forcing

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: temporal_disagg_forcing ! Temporally distribute forcings onto time step

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         temporal_disagg_forcing

  !     PURPOSE
  !>        \brief Temporally distribute daily mean forcings onto time step

  !>        \details Calculates actual precipitation, PET and temperature from daily mean inputs.\n
  !>        Precipitation and PET are distributed with predefined factors onto the day.\n
  !>        Temperature gets a predefined amplitude added on day and substracted at night.

  !     CALLING SEQUENCE
  !         elemental pure subroutine temporal_disagg_forcing(isday, ntimesteps_day, prec_day, pet_day, temp_day, &
  !             fday_prec, fday_pet, fday_temp, fnight_prec, fnight_pet, fnight_temp, prec, &
  !             pet, temp)

  !     INTENT(IN)
  !>        \param[in] "logical  :: isday"            is day or night
  !>        \param[in] "real(dp) :: ntimesteps_day"   # of time steps per day
  !>        \param[in] "real(dp) :: prec_day"         Daily mean precipitation [mm/s]
  !>        \param[in] "real(dp) :: pet_day"          Daily mean ET [mm/s]
  !>        \param[in] "real(dp) :: temp_day"         Daily mean air temperature [K]
  !>        \param[in] "real(dp) :: fday_prec"        Daytime fraction of precipitation
  !>        \param[in] "real(dp) :: fday_pet"         Daytime fraction of PET
  !>        \param[in] "real(dp) :: fday_temp"        Daytime air temparture increase
  !>        \param[in] "real(dp) :: fnight_prec"      Daytime fraction of precipitation
  !>        \param[in] "real(dp) :: fnight_pet"       Daytime fraction of PET
  !>        \param[in] "real(dp) :: fnight_temp"      Daytime air temparture increase

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: prec"           Actual precipitation [mm/s]
  !>        \param[out] "real(dp) :: pet"            Reference ET [mm/s]
  !>        \param[out] "real(dp) :: temp"           Air temperature [K]

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !>        \note There is no PET correction for aspect in this routine. Use pet * fasp before or after the routine.

  !     EXAMPLE
  !>        \n Example:
  !>        \code
  !>        call temporal_disagg_forcing((tt <= 6) .or. ((tt>=19) .and. (tt<=24)), &
  !>                 real(NTSTEPDAY, dp), prec_day(k), pet_day(k), temp_day(k), &
  !>                 fpre(month,2), fpet(month,2), ftem(month,2), &
  !>                 fpre(month,1), fpet(month,1), ftem(month,1), &
  !>                 prec(k), pet(k), temp(k))
  !>        pet(k) = pet(k) * fasp(k)
  !>        \endcode

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Dec 2012
  !
  !         Modifications
  !         R. Kumar, July 2013 --> added default initalization of meteorological values
  !                                 also added the if loop around the temporal dissagg. step
  elemental pure subroutine temporal_disagg_forcing(isday, ntimesteps_day, prec_day, pet_day, temp_day, &
       fday_prec, fday_pet, fday_temp, fnight_prec, fnight_pet, fnight_temp, prec, &
       pet, temp)

    implicit none

    ! Intent variables
    logical,  intent(in)    :: isday          ! is day or night
    real(dp), intent(in)    :: ntimesteps_day ! # of time steps per day
    real(dp), intent(in)    :: prec_day       ! Daily mean precipitation [mm/s]
    real(dp), intent(in)    :: pet_day        ! Daily mean ET [mm/s]
    real(dp), intent(in)    :: temp_day       ! Daily mean air temperature [K]
    real(dp), intent(in)    :: fday_prec      ! Daytime fraction of precipitation
    real(dp), intent(in)    :: fday_pet       ! Daytime fraction of PET
    real(dp), intent(in)    :: fday_temp      ! Daytime air temparture increase
    real(dp), intent(in)    :: fnight_prec    ! Daytime fraction of precipitation
    real(dp), intent(in)    :: fnight_pet     ! Daytime fraction of PET
    real(dp), intent(in)    :: fnight_temp    ! Daytime air temparture increase
    real(dp), intent(out)   :: prec           ! actual precipitation [mm/s]
    real(dp), intent(out)   :: pet            ! Reference ET [mm/s]
    real(dp), intent(out)   :: temp           ! Air temperature [K]

    ! default vaule used if ntimesteps_day = 1 (i.e., e.g. daily values)
    prec =  prec_day
    pet  =  pet_day 
    temp =  temp_day

    ! Distribute Prec, PET and Temp into time steps night/day
    if(ntimesteps_day .gt. 1.0_dp) then 
      if (isday) then ! DAY-TIME
         prec = 2.0_dp * prec_day * fday_prec / ntimesteps_day
         pet  = 2.0_dp * pet_day  * fday_pet  / ntimesteps_day
         temp =          temp_day + fday_temp
      else            ! NIGHT-TIME
         prec = 2.0_dp * prec_day * fnight_prec / ntimesteps_day
         pet  = 2.0_dp * pet_day  * fnight_pet  / ntimesteps_day
         temp =          temp_day + fnight_temp
      end if
    end if

  end subroutine temporal_disagg_forcing

END MODULE mo_temporal_disagg_forcing
