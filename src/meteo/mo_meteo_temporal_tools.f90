!> \file mo_meteo_temporal_tools.f90
!> \brief \copybrief mo_meteo_temporal_tools
!> \details \copydetails mo_meteo_temporal_tools

!> \brief Temporal disaggregation of daily input values
!> \details Calculate actual values for precipitation, PET and temperature from daily mean inputs
!> \note There is not PET correction for aspect in this routine. Use pet * fasp before or after the routine.
!> \authors Matthias Cuntz
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_meteo
MODULE mo_meteo_temporal_tools

  USE mo_kind, ONLY : dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: temporal_disagg_meteo_weights ! Temporally distribute meteo value by weights correction
  PUBLIC :: temporal_disagg_flux_daynight ! Temporally distribute meteo flux by day/night correction
  PUBLIC :: temporal_disagg_state_daynight ! Temporally distribute meteo state variable by day/night correction

CONTAINS

  !> \brief Temporally distribute daily mean forcings onto time step
  !> \details Calculates actual meteo forcings from daily mean inputs.
  !! They are distributed with predefined weights onto the day.
  !> \authors Sebastian Mueller
  !> \date Jun 2020
  elemental subroutine temporal_disagg_meteo_weights( &
    meteo_val_day, &
    meteo_val_weights, &
    meteo_val, &
    weights_correction &
  )
    implicit none

    !> Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    !> weights for meteo_val
    real(dp), intent(in) :: meteo_val_weights
    !> Actual meteo_val
    real(dp), intent(out) :: meteo_val
    !> Additive correction before applying weights (e.g. for temperature conversion)
    real(dp), intent(in), optional :: weights_correction

    real(dp) :: weights_correction_

    ! default values (can't be initialized directly in a pure function)
    weights_correction_ = 0.0_dp
    ! set potential given optional values
    if ( present(weights_correction) ) weights_correction_ = weights_correction

    ! use weights to distribute meteo values
    meteo_val = (meteo_val_day + weights_correction_) * meteo_val_weights - weights_correction_

  end subroutine temporal_disagg_meteo_weights

  !> \brief Temporally distribute daily mean forcings onto time step
  !> \details Calculates actual meteo forcings from daily mean inputs.
  !! They are distributed with predefined factors/summands for day and night.
  !> \authors Sebastian Mueller
  !> \date Jun 2020
  elemental subroutine temporal_disagg_flux_daynight( &
    isday, &
    ntimesteps_day, &
    meteo_val_day, &
    fday_meteo_val, &
    fnight_meteo_val, &
    meteo_val &
  )
    implicit none

    !> is day (False for night)
    logical, intent(in) :: isday
    !> number of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    !> Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    !> Daytime fraction of meteo_val
    real(dp), intent(in) :: fday_meteo_val
    !> Nighttime fraction of meteo_val
    real(dp), intent(in) :: fnight_meteo_val
    !> Actual meteo_val
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

  !> \brief Temporally distribute daily mean state forcings onto time step
  !> \details Calculates meteo forcings from daily mean inputs of a state variable.
  !> They are distributed with predefined factors/summands for day and night.
  !> \authors Sebastian Mueller
  !> \date Jun 2020
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

    !> is day (False for night)
    logical, intent(in) :: isday
    !> number of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    !> Daily meteo_val
    real(dp), intent(in) :: meteo_val_day
    !> Daytime fraction of meteo_val
    real(dp), intent(in) :: fday_meteo_val
    !> Nighttime fraction of meteo_val
    real(dp), intent(in) :: fnight_meteo_val
    !> Actual meteo_val
    real(dp), intent(out) :: meteo_val
    !> if True, correcting values will be added (e.g. for temperature), otherwise used as percentage
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

END MODULE mo_meteo_temporal_tools
