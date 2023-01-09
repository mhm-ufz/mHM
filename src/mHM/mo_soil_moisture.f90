!> \file mo_soil_moisture.f90
!> \brief \copybrief mo_soil_moisture
!> \details \copydetails mo_soil_moisture

!> \brief Soil moisture of the different layers
!> \details Soil moisture in the different layers is calculated with
!! infiltration as \f$ (\theta / \theta_{sat})^\beta \f$
!! Then evapotranspiration is calculated from PET with a soil water stress  factor \f$ f_{SM} \f$
!! either using  the Feddes equation - processCase(3) = 1:
!! \f[ f_{SM} = \frac{\theta - \theta_\mathit{pwp}}{\theta_\mathit{fc} - \theta_\mathit{pwp}} \f]
!! or using the Jarvis equation - processCase(3) = 2:
!! \f[ f_{SM} = \frac{1}{\theta_\mathit{stress-index-C1}}
!! \frac{\theta - \theta_\mathit{pwp}}{\theta_\mathit{sat} - \theta_\mathit{pwp}} \f]
!> \authors Matthias Cuntz, Luis Samaniego
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_soil_moisture

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PUBLIC :: feddes_et_reduction
  PUBLIC :: jarvis_et_reduction

  PUBLIC :: soil_moisture  ! Soil moisture in different soil horizons

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        soil_moisture

  !    PURPOSE
  !>       \brief Soil moisture in different soil horizons

  !>       \details Infiltration \f$I\f$ from one layer \f$k-1\f$ to the next \f$k\f$ on
  !>       pervious areas is calculated as (omit \f$t\f$)
  !>       \f[ I[k] = I[k-1] (\theta[k] / \theta_{sat}[k])^{\beta[k]} \f]
  !>       Then soil moisture can be calculated as (omit \f$k\f$)
  !>       \f[ \theta[t] = \theta[t-1] + I[t] - \mathit{ET}[t] \f]
  !>       with \f$ \mathit{ET} \f$ (omit \f$[k,t]\f$) being
  !>       \f[ \mathit{ET} = f_\mathrm{roots} \cdot f_{SM} \cdot \mathit{PET} \f].

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: processCase"                   1 - Feddes equation for PET reduction2 - Jarvis
  !>       equation for PET reduction3 - Jarvis equation for PET reduction and FC dependency on root fraction coefficient
  !>       \param[in] "real(dp) :: frac_sealed"                      Fraction of sealed area
  !>       \param[in] "real(dp) :: water_thresh_sealed"              Threshhold water depth in impervious areas [mm TS-1]
  !>       \param[in] "real(dp) :: pet"                              Reference evapotranspiration [mm TS-1]
  !>       \param[in] "real(dp) :: evap_coeff"                       Evaporation coefficent for free-water surface of
  !>       that current month
  !>       \param[in] "real(dp), dimension(:) :: soil_moist_sat"     Saturation soil moisture for each horizon [mm]
  !>       \param[in] "real(dp), dimension(:) :: frac_roots"         Fraction of Roots in soil horizon
  !>       \param[in] "real(dp), dimension(:) :: soil_moist_FC"      Soil moisture below which actual ET is reduced [mm]
  !>       \param[in] "real(dp), dimension(:) :: wilting_point"      Permanent wilting point for each horizon [mm]
  !>       \param[in] "real(dp), dimension(:) :: soil_moist_exponen" Exponential parameter to how non-linear is the soil
  !>       water retention
  !>       \param[in] "real(dp) :: jarvis_thresh_c1"                 Jarvis critical value for normalized soil water
  !>       content
  !>       \param[in] "real(dp) :: aet_canopy"                       Actual ET from canopy [mm TS-1]

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp) :: prec_effec"                                       Effective precipitation (rain +
  !>       snow melt) [mm]
  !>       \param[inout] "real(dp) :: runoff_sealed"                                    Direct runoff from impervious
  !>       areas
  !>       \param[inout] "real(dp) :: storage_sealed"                                   Retention storage of impervious
  !>       areas
  !>       \param[inout] "real(dp), dimension(size(soil_moist_sat, 1)) :: infiltration" Recharge, infiltration intensity
  !>       oreffective precipitation of each horizon [mm TS-1]
  !>       \param[inout] "real(dp), dimension(size(soil_moist_sat, 1)) :: soil_moist"   Soil moisture of each horizon
  !>       [mm]

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(size(soil_moist_sat, 1)) :: aet" actual ET [mm TS-1]
  !>       \param[out] "real(dp) :: aet_sealed"                              actual ET from free-water surfaces,i.e
  !>       impervious cover [mm TS-1]

  !    HISTORY
  !>       \authors Matthias Cuntz

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! M. Cuneyd Demirel, Simon Stisen Jun 2020 - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4

  subroutine soil_moisture(processCase, frac_sealed, water_thresh_sealed, pet, evap_coeff, soil_moist_sat, frac_roots, &
                          soil_moist_FC, wilting_point, soil_moist_exponen, jarvis_thresh_c1, aet_canopy, prec_effec, &
                          runoff_sealed, storage_sealed, infiltration, soil_moist, aet, aet_sealed)

    use mo_common_constants, only : eps_dp

    implicit none

    ! 1 - Feddes equation for PET reduction2 - Jarvis equation for PET reduction3 - Jarvis equation for PET reduction
    ! and FC dependency on root fraction coefficient
    integer(i4), intent(in) :: processCase

    ! Fraction of sealed area
    real(dp), intent(in) :: frac_sealed

    ! Threshhold water depth in impervious areas [mm TS-1]
    real(dp), intent(in) :: water_thresh_sealed

    ! Reference evapotranspiration [mm TS-1]
    real(dp), intent(in) :: pet

    ! Evaporation coefficent for free-water surface of that current month
    real(dp), intent(in) :: evap_coeff

    ! Saturation soil moisture for each horizon [mm]
    real(dp), dimension(:), intent(in) :: soil_moist_sat

    ! Fraction of Roots in soil horizon
    real(dp), dimension(:), intent(in) :: frac_roots

    ! Soil moisture below which actual ET is reduced [mm]
    real(dp), dimension(:), intent(in) :: soil_moist_FC

    ! Permanent wilting point for each horizon [mm]
    real(dp), dimension(:), intent(in) :: wilting_point

    ! Exponential parameter to how non-linear is the soil water retention
    real(dp), dimension(:), intent(in) :: soil_moist_exponen

    ! Jarvis critical value for normalized soil water content
    real(dp), intent(in) :: jarvis_thresh_c1

    ! Actual ET from canopy [mm TS-1]
    real(dp), intent(in) :: aet_canopy

    ! Effective precipitation (rain + snow melt) [mm]
    real(dp), intent(in) :: prec_effec

    ! Direct runoff from impervious areas
    real(dp), intent(inout) :: runoff_sealed

    ! Retention storage of impervious areas
    real(dp), intent(inout) :: storage_sealed

    ! Recharge, infiltration intensity oreffective precipitation of each horizon [mm TS-1]
    real(dp), dimension(size(soil_moist_sat, 1)), intent(inout) :: infiltration

    ! Soil moisture of each horizon [mm]
    real(dp), dimension(size(soil_moist_sat, 1)), intent(inout) :: soil_moist

    ! actual ET [mm TS-1]
    real(dp), dimension(size(soil_moist_sat, 1)), intent(out) :: aet

    ! actual ET from free-water surfaces,i.e impervious cover [mm TS-1]
    real(dp), intent(out) :: aet_sealed

    ! counter
    integer(i4) :: hh

    ! Effective Prec or infiltration from above
    real(dp) :: prec_effec_soil

    ! Runoof fraction
    real(dp) :: frac_runoff

    ! PET reduction factor according to actual soil moisture
    real(dp) :: soil_stress_factor

    ! temporary variable for misc use
    real(dp) :: tmp


    ! ----------------------------------------------------------------
    ! IMPERVIOUS COVER PROCESS
    ! ----------------------------------------------------------------
    runoff_sealed = 0.0_dp
    aet_sealed = 0.0_dp

    if (frac_sealed > 0.0_dp) then
      tmp = storage_sealed + prec_effec

      if (tmp > water_thresh_sealed) then
        runoff_sealed = tmp - water_thresh_sealed
        storage_sealed = water_thresh_sealed
      else
        runoff_sealed = 0.0_dp
        storage_sealed = tmp
      end if

      ! aET from sealed area is propotional to the available water content
      if(water_thresh_sealed .gt. eps_dp) then
        aet_sealed = (pet / evap_coeff - aet_canopy) * (storage_sealed / water_thresh_sealed)
        ! numerical problem
        if (aet_sealed .lt. 0.0_dp) aet_sealed = 0.0_dp
      else
        aet_sealed = huge(1.0_dp)
      end if

      ! sealed storage updata
      if (storage_sealed .gt. aet_sealed) then
        storage_sealed = storage_sealed - aet_sealed
      else
        aet_sealed = storage_sealed
        storage_sealed = 0.0_dp
      end if

    end if
    ! ----------------------------------------------------------------
    ! N-LAYER SOIL MODULE
    ! ----------------------------------------------------------------
    aet(:) = 0.0_dp
    infiltration(:) = 0.0_dp

    ! for 1st layer input is prec_effec
    prec_effec_soil = prec_effec

    do hh = 1, size(soil_moist_sat, 1) ! nHorizons
      ! input for other layers is the infiltration from its immediate upper layer will be input
      if (hh .NE. 1) prec_effec_soil = infiltration(hh - 1)

      !  start processing for soil moisture process
      !  BASED ON SMs as its upper LIMIT

      if (soil_moist(hh) > soil_moist_sat(hh)) then
        infiltration(hh) = prec_effec_soil
      else
        ! to avoid underflow -- or numerical errors
        if(soil_moist(hh) > eps_dp) then
          !frac_runoff = (soil_moist(hh) / soil_moist_sat(hh))**soil_moist_exponen(hh)
          frac_runoff = exp(soil_moist_exponen(hh) * log(soil_moist(hh) / soil_moist_sat(hh)))
        else
          frac_runoff = 0.0_dp
        end if
        tmp = prec_effec_soil * (1.0_dp - frac_runoff)

        if ((soil_moist(hh) + tmp) > soil_moist_sat(hh)) then
          infiltration(hh) = prec_effec_soil + (soil_moist(hh) - soil_moist_sat(hh))
          soil_moist(hh) = soil_moist_sat(hh)
        else
          infiltration(hh) = prec_effec_soil - tmp
          soil_moist(hh) = soil_moist(hh) + tmp
        end if
      end if

      !             aET calculations

      !  Satisfying ET demand sequentially from top to the bottom layer
      !  Note that the potential ET for the first soil layer is reduced after
      !  satisfying ET demands of the canopy surface

      aet(hh) = pet - aet_canopy                                                     ! First layer
      if (hh /= 1) aet(hh) = aet(hh) - sum(aet(1 : hh - 1), mask = (aet(1 : hh - 1) > 0.0_dp)) ! remaining layers

      ! estimate fraction of ET demand based on root fraction and SM status
      select case(processCase)
        ! FEDDES EQUATION: https://doi.org/10.1016/0022-1694(76)90017-2
      case(1 , 4)
        soil_stress_factor = feddes_et_reduction(soil_moist(hh), soil_moist_FC(hh), wilting_point(hh), &
                             frac_roots(hh))
        ! JARVIS EQUATION: https://doi.org/10.1016/0022-1694(89)90050-4
      case(2 , 3)
        !!!!!!!!! INTRODUCING STRESS FACTOR FOR SOIL MOISTURE ET REDUCTION !!!!!!!!!!!!!!!!!
        soil_stress_factor = jarvis_et_reduction(soil_moist(hh), soil_moist_sat(hh), wilting_point(hh), &
                             frac_roots(hh), jarvis_thresh_c1)
      end select

      aet(hh) = aet(hh) * soil_stress_factor

      ! avoid numerical error
      if(aet(hh) < 0.0_dp) aet(hh) = 0.0_dp

      ! reduce SM state
      if(soil_moist(hh) > aet(hh)) then
        soil_moist(hh) = soil_moist(hh) - aet(hh)
      else
        aet(hh) = soil_moist(hh) - eps_dp
        soil_moist(hh) = eps_dp
      end if

      ! avoid numerical error of underflow
      if(soil_moist(hh) < eps_dp) soil_moist(hh) = eps_dp

    end do ! hh



  end subroutine soil_moisture


  ! ------------------------------------------------------------------

  !    NAME
  !        feddes_et_reduction

  !    PURPOSE
  !>       \brief stress factor for reducing evapotranspiration based on actual soil moisture

  !>       \details Potential evapotranspiration is reduced to 0 if SM is lower PWP. PET is equal
  !>       fraction of roots if soil moisture is exceeding field capacity. If soil moisture is
  !>       in between PWP and FC PET is reduced by fraction of roots times a stress factor.

  !>       The ET reduction factor \f$ f \f$ is estimated as
  !>       \f[ f = \left\{
  !>       \begin{array}{lr}
  !>       f_{roots}  & if \theta \ge \theta_{fc}\\
  !>       f_{roots} \cdot \frac{\theta - \theta_\mathit{pwp}}{\theta_\mathit{fc} - \theta_\mathit{pwp}} &
  !>       if \theta < \theta_{fc} \\
  !>       0 & if \theta < \theta_{pwp}
  !>       \end{array}
  !>       \right. \f]

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: soil_moist"    Soil moisture of each horizon [mm]
  !>       \param[in] "real(dp) :: soil_moist_FC" Soil moisture below which actual ET is reduced [mm]
  !>       \param[in] "real(dp) :: wilting_point" Permanent wilting point
  !>       \param[in] "real(dp) :: frac_roots"    Fraction of Roots in soil horizon is reduced [mm]

  !    RETURN
  !>       \return real(dp) :: feddes_et_reduction; et reduction factor

  !    HISTORY
  !>       \authors Matthias Cuntz, Cueneyd Demirel, Matthias Zink

  !>       \date March 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! M. Cuneyd Demirel, Simon Stisen Jun 2020 - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4

  elemental pure FUNCTION feddes_et_reduction(soil_moist, soil_moist_FC, wilting_point, frac_roots)
    implicit none

    ! Soil moisture of each horizon [mm]
    real(dp), intent(in) :: soil_moist

    ! Soil moisture below which actual ET is reduced [mm]
    real(dp), intent(in) :: soil_moist_FC

    ! Permanent wilting point
    real(dp), intent(in) :: wilting_point

    ! Fraction of Roots in soil horizon is reduced [mm]
    real(dp), intent(in) :: frac_roots

    ! reference evapotranspiration in [mm s-1]
    real(dp) :: feddes_et_reduction


    !    SM >= FC
    if (soil_moist >= soil_moist_FC) then
      feddes_et_reduction = frac_roots
      ! PW < SM < FC
    else if (soil_moist > wilting_point) then
      feddes_et_reduction = frac_roots * (soil_moist - wilting_point) / (soil_moist_FC - wilting_point)
      ! SM <= PW
    else
      feddes_et_reduction = 0.0_dp
    end if

  END FUNCTION feddes_et_reduction

  ! ------------------------------------------------------------------

  !    NAME
  !        jarvis_et_reduction

  !    PURPOSE
  !>       \brief stress factor for reducing evapotranspiration based on actual soil moisture

  !>       \details The soil moisture stress factor is estimated based on the normalized soil water
  !>       content. The normalized soil water content \f$ \theta_{norm} \f$ is estimated as:
  !>       \f[ \theta_{norm} =  \frac{\theta - \theta_\mathit{pwp}}
  !>       {\theta_{sat} - \theta_{pwp}}  \f]
  !>       The ET reduction factor \f$ f \f$ is estimated as
  !>       \f[ f = \left\{
  !>       \begin{array}{lr}
  !>       f_{roots}  & if \theta_{norm} \ge jarvis\_sm\_threshold\_c1 \\
  !>       f_{roots}\frac{\theta_{norm}}{jarvis\_sm\_threshold\_c1}  &
  !>       if  \theta_{norm} < jarvis\_sm\_threshold\_c1 \\
  !>       \end{array}
  !>       \right. \f]

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: soil_moist"       Soil moisture of each horizon [mm]
  !>       \param[in] "real(dp) :: soil_moist_sat"   saturated Soil moisture content [mm]
  !>       \param[in] "real(dp) :: wilting_point"    Permanent wilting point
  !>       \param[in] "real(dp) :: frac_roots"       Fraction of Roots in soil horizon is reduced [mm]
  !>       \param[in] "real(dp) :: jarvis_thresh_c1" parameter C1 from Jarvis formulation

  !    RETURN
  !>       \return real(dp) :: jarvis_et_reduction; et reduction factor

  !    HISTORY
  !>       \authors Cueneyd Demirel, Matthias Zink

  !>       \date March 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! M. Cuneyd Demirel, Simon Stisen Jun 2020 - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4

  elemental pure FUNCTION jarvis_et_reduction(soil_moist, soil_moist_sat, wilting_point, frac_roots, jarvis_thresh_c1)
    implicit none

    ! Soil moisture of each horizon [mm]
    real(dp), intent(in) :: soil_moist

    ! saturated Soil moisture content [mm]
    real(dp), intent(in) :: soil_moist_sat

    ! Permanent wilting point
    real(dp), intent(in) :: wilting_point

    ! Fraction of Roots in soil horizon is reduced [mm]
    real(dp), intent(in) :: frac_roots

    ! parameter C1 from Jarvis formulation
    real(dp), intent(in) :: jarvis_thresh_c1

    ! reference evapotranspiration in [mm]
    real(dp) :: jarvis_et_reduction

    ! normalized soil water content
    real(dp) :: theta_inorm


    ! Calculating normalized Soil Water Content
    theta_inorm = (soil_moist - wilting_point) / (soil_moist_sat - wilting_point)

    ! correct for numerical unaccuracies
    if (theta_inorm .LT. 0.0_dp)    theta_inorm = 0.0_dp
    if (theta_inorm .GT. 1.0_dp)    theta_inorm = 1.0_dp

    ! estimate fraction of ET demand based on root fraction and SM status using theta_inorm
    ! theta_inorm >= jarvis_thresh_c1
    if (theta_inorm .GE. jarvis_thresh_c1) then
      jarvis_et_reduction = frac_roots
      ! 0 < theta_inorm < jarvis_thresh_c1
    else if (theta_inorm .LT. jarvis_thresh_c1) then
      jarvis_et_reduction = frac_roots * (theta_inorm / jarvis_thresh_c1)
    end if

  END FUNCTION jarvis_et_reduction

END MODULE mo_soil_moisture
