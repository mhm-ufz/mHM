!>       \file mo_temporal_disagg_forcing.f90

!>       \brief Temporal disaggregation of daily input values

!>       \details Calculate actual values for precipitation, PET and temperature from daily mean inputs
!>       ote There is not PET correction for aspect in this routine. Use pet * fasp before or after the routine.

!>       \authors Matthias Cuntz

!>       \date Dec 2012

! Modifications:

MODULE mo_temporal_disagg_forcing

  USE mo_kind, ONLY : dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: temporal_disagg_forcing ! Temporally distribute forcings onto time step
  PUBLIC :: temporal_disagg_meteo_weights ! Temporally distribute meteo value by weights correction
  PUBLIC :: temporal_disagg_flux_daynight ! Temporally distribute meteo flux by day/night correction
  PUBLIC :: temporal_disagg_state_daynight ! Temporally distribute meteo state variable by day/night correction

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        temporal_disagg_forcing

  !    PURPOSE
  !>       \brief Temporally distribute daily mean forcings onto time step

  !>       \details Calculates actual precipitation, PET and temperature from daily mean inputs.
  !>       Precipitation and PET are distributed with predefined factors onto the day.
  !>       Temperature gets a predefined amplitude added on day and substracted at night.
  !>       Alternatively, weights for each hour and month can be given and disaggregation is
  !>       using these as factors for PET and temperature. Precipitation is distributed uniformly.

  !    INTENT(IN)
  !>       \param[in] "logical :: isday"              is day or night
  !>       \param[in] "real(dp) :: ntimesteps_day"    # of time steps per day
  !>       \param[in] "real(dp) :: prec_day"          Daily mean precipitation [mm/d]
  !>       \param[in] "real(dp) :: pet_day"           Daily mean ET [mm/d]
  !>       \param[in] "real(dp) :: temp_day"          Daily mean air temperature [K]
  !>       \param[in] "real(dp) :: fday_prec"         Daytime fraction of precipitation
  !>       \param[in] "real(dp) :: fday_pet"          Daytime fraction of PET
  !>       \param[in] "real(dp) :: fday_temp"         Daytime air temparture increase
  !>       \param[in] "real(dp) :: fnight_prec"       Daytime fraction of precipitation
  !>       \param[in] "real(dp) :: fnight_pet"        Daytime fraction of PET
  !>       \param[in] "real(dp) :: fnight_temp"       Daytime air temparture increase
  !>       \param[in] "real(dp) :: temp_weights"      weights for average temperature
  !>       \param[in] "real(dp) :: pet_weights"       weights for PET
  !>       \param[in] "real(dp) :: pre_weights"       weights for precipitation
  !>       \param[in] "logical :: read_meteo_weights" flag indicating that weights should be used

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: prec" Actual precipitation [mm/d]
  !>       \param[out] "real(dp) :: pet"  Reference ET [mm/d]
  !>       \param[out] "real(dp) :: temp" Air temperature [K]

  !    HISTORY
  !>       \authors Matthias Cuntz

  !>       \date Dec 2012

  ! Modifications:
  ! S. Thober Jan 2017 - > added disaggregation based on weights given in nc file
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Sebastian Mueller Jul 2020 - separate routines for weights and daynight disaggregation

  elemental pure subroutine temporal_disagg_forcing( &
      isday, ntimesteps_day, &
      prec_day, pet_day, temp_day, &
      fday_prec, fday_pet, fday_temp, &
      fnight_prec, fnight_pet, fnight_temp, &
      temp_weights, pet_weights, pre_weights, read_meteo_weights, &
      prec, pet, temp &
  )
    use mo_constants, only : T0_dp  ! 273.15 - Celcius <-> Kelvin [K]

    implicit none

    ! is day or night
    logical, intent(in) :: isday
    ! # of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    ! Daily mean precipitation [mm/d]
    real(dp), intent(in) :: prec_day
    ! Daily mean ET [mm/d]
    real(dp), intent(in) :: pet_day
    ! Daily mean air temperature [K]
    real(dp), intent(in) :: temp_day
    ! Daytime fraction of precipitation
    real(dp), intent(in) :: fday_prec
    ! Daytime fraction of PET
    real(dp), intent(in) :: fday_pet
    ! Daytime air temparture increase
    real(dp), intent(in) :: fday_temp
    ! Daytime fraction of precipitation
    real(dp), intent(in) :: fnight_prec
    ! Daytime fraction of PET
    real(dp), intent(in) :: fnight_pet
    ! Daytime air temparture increase
    real(dp), intent(in) :: fnight_temp
    ! weights for average temperature
    real(dp), intent(in) :: temp_weights
    ! weights for PET
    real(dp), intent(in) :: pet_weights
    ! weights for precipitation
    real(dp), intent(in) :: pre_weights
    ! flag indicating that weights should be used
    logical, intent(in) :: read_meteo_weights
    ! Actual precipitation [mm/d]
    real(dp), intent(out) :: prec
    ! Reference ET [mm/d]
    real(dp), intent(out) :: pet
    ! Air temperature [K]
    real(dp), intent(out) :: temp

    if (read_meteo_weights) then
      ! all meteo forcings are disaggregated with given weights
      call temporal_disagg_meteo_weights(pet_day, pet_weights, pet)
      call temporal_disagg_meteo_weights(prec_day, pre_weights, prec)
      call temporal_disagg_meteo_weights(temp_day, temp_weights, temp, weights_correction=T0_dp)
    else
      ! all meteo forcings are disaggregated with day-night correction values
      call temporal_disagg_flux_daynight(isday, ntimesteps_day, pet_day, fday_pet, fnight_pet, pet)
      call temporal_disagg_flux_daynight(isday, ntimesteps_day, prec_day, fday_prec, fnight_prec, prec)
      call temporal_disagg_state_daynight(isday, ntimesteps_day, temp_day, fday_temp, fnight_temp, temp, add_correction=.true.)
    end if

  end subroutine temporal_disagg_forcing

  ! ------------------------------------------------------------------

  !    NAME
  !        temporal_disagg_meteo_weights

  !    PURPOSE
  !>       \brief Temporally distribute daily mean forcings onto time step

  !>       \details Calculates actual meteo forcings from daily mean inputs.
  !>       They are distributed with predefined weights onto the day.

  !    INTENT(IN)
  !>       \param[in] "logical :: isday"              is day or night
  !>       \param[in] "real(dp) :: ntimesteps_day"    # of time steps per day
  !>       \param[in] "real(dp) :: meteo_val_day"     Daily mean meteo value

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: meteo_val"        Actual meteo value

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: weights_correction" Additive correction value before weights are applied

  !    HISTORY
  !>       \authors Sebastian Mueller

  !>       \date Jun 2020

  elemental subroutine temporal_disagg_meteo_weights( &
    meteo_val_day, &
    meteo_val_weights, &
    meteo_val, &
    weights_correction &
  )
    implicit none

    ! Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    ! weights for meteo_val
    real(dp), intent(in) :: meteo_val_weights
    ! Actual meteo_val
    real(dp), intent(out) :: meteo_val
    ! Additive correction before applying weights (e.g. for temperature conversion)
    real(dp), intent(in), optional :: weights_correction

    real(dp) :: weights_correction_

    ! default values (can't be initialized directly in a pure function)
    weights_correction_ = 0.0_dp
    ! set potential given optional values
    if ( present(weights_correction) ) weights_correction_ = weights_correction

    ! use weights to distribute meteo values
    meteo_val = (meteo_val_day + weights_correction_) * meteo_val_weights - weights_correction_

  end subroutine temporal_disagg_meteo_weights

  ! ------------------------------------------------------------------

  !    NAME
  !        temporal_disagg_flux_daynight

  !    PURPOSE
  !>       \brief Temporally distribute daily mean forcings onto time step

  !>       \details Calculates actual meteo forcings from daily mean inputs.
  !>       They are distributed with predefined factors/summands for day and night.

  !    INTENT(IN)
  !>       \param[in] "logical :: isday"              is day or night
  !>       \param[in] "real(dp) :: ntimesteps_day"    # of time steps per day
  !>       \param[in] "real(dp) :: meteo_val_day"     Daily mean meteo value
  !>       \param[in] "real(dp) :: fday_meteo_val"    Daytime fraction of meteo value
  !>       \param[in] "real(dp) :: fnight_meteo_val"  Nighttime fraction of meteo value

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: meteo_val"        Actual meteo value

  !    HISTORY
  !>       \authors Sebastian Mueller

  !>       \date Jun 2020

  elemental subroutine temporal_disagg_flux_daynight( &
    isday, &
    ntimesteps_day, &
    meteo_val_day, &
    fday_meteo_val, &
    fnight_meteo_val, &
    meteo_val &
  )
    implicit none

    ! is day (False for night)
    logical, intent(in) :: isday
    ! # of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    ! Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    ! Daytime fraction of meteo_val
    real(dp), intent(in) :: fday_meteo_val
    ! Nighttime fraction of meteo_val
    real(dp), intent(in) :: fnight_meteo_val
    ! Actual meteo_val
    real(dp), intent(out) :: meteo_val

    ! Distribute into time steps night/day
    if(ntimesteps_day .gt. 1.0_dp) then
      if ( isday ) then ! DAY-TIME
        meteo_val = 2.0_dp * meteo_val_day * fday_meteo_val / ntimesteps_day
      else            ! NIGHT-TIME
        meteo_val = 2.0_dp * meteo_val_day * fnight_meteo_val / ntimesteps_day
      end if
    else
      ! default vaule used if ntimesteps_day = 1 (i.e., e.g. daily values)
      meteo_val = meteo_val_day
    end if

  end subroutine temporal_disagg_flux_daynight

  ! ------------------------------------------------------------------

  !    NAME
  !        temporal_disagg_state_daynight

  !    PURPOSE
  !>       \brief Temporally distribute daily mean state forcings onto time step

  !>       \details Calculates meteo forcings from daily mean inputs of a state variable.
  !>       They are distributed with predefined factors/summands for day and night.

  !    INTENT(IN)
  !>       \param[in] "logical :: isday"              is day or night
  !>       \param[in] "real(dp) :: ntimesteps_day"    # of time steps per day
  !>       \param[in] "real(dp) :: meteo_val_day"     Daily mean meteo value
  !>       \param[in] "real(dp) :: fday_meteo_val"    Daytime fraction/addition of meteo value
  !>       \param[in] "real(dp) :: fnight_meteo_val"  Nighttime fraction/addition of meteo value

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: meteo_val"        Actual meteo value

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: add_correction"   State if correction should be added instead of multiplied

  !    HISTORY
  !>       \authors Sebastian Mueller

  !>       \date Jun 2020

  elemental subroutine temporal_disagg_state_daynight( &
    isday, &
    ntimesteps_day, &
    meteo_val_day, &
    fday_meteo_val, &
    fnight_meteo_val, &
    meteo_val, &
    add_correction &
  )
    implicit none

    ! is day (False for night)
    logical, intent(in) :: isday
    ! # of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    ! Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    ! Daytime fraction of meteo_val
    real(dp), intent(in) :: fday_meteo_val
    ! Nighttime fraction of meteo_val
    real(dp), intent(in) :: fnight_meteo_val
    ! Actual meteo_val
    real(dp), intent(out) :: meteo_val
    ! if True, correcting values will be added (e.g. for temperature), otherwise used as percentage
    logical, intent(in), optional :: add_correction

    logical :: add_correction_

    ! default values (can't be initialized directly in a pure function)
    add_correction_ = .False.
    ! set potential given optional values
    if ( present(add_correction) ) add_correction_ = add_correction

    ! Distribute into time steps night/day
    if(ntimesteps_day .gt. 1.0_dp) then
      if ( add_correction_ ) then  ! for e.g. temperature
        if ( isday ) then ! DAY-TIME
          meteo_val = meteo_val_day + fday_meteo_val
        else            ! NIGHT-TIME
          meteo_val = meteo_val_day + fnight_meteo_val
        end if
      else
        if ( isday ) then ! DAY-TIME
          meteo_val = 2.0_dp * meteo_val_day * fday_meteo_val
        else            ! NIGHT-TIME
          meteo_val = 2.0_dp * meteo_val_day * fnight_meteo_val
        end if
      end if
    else
      ! default vaule used if ntimesteps_day = 1 (i.e., e.g. daily values)
      meteo_val = meteo_val_day
    end if

  end subroutine temporal_disagg_state_daynight

END MODULE mo_temporal_disagg_forcing
