!> \file mo_soil_moisture.f90

!> \brief Soil moisture of the different layers

!> \details Soil moisture in the different layers is calculated with
!> infiltration as \f$ (\theta / \theta_{sat})^\beta \f$ \n
!> Then evapotranspiration is calculated from PET with a soil water reduction factor
!> \f$ \frac{\theta - \theta_\mathit{pwp}}{\theta_\mathit{fc} - \theta_\mathit{pwp}} \f$.

!> \authors Matthias Cuntz, Luis Samaniego
!> \date Dec 2012

MODULE mo_soil_moisture

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: soil_moisture  ! Soil moisture in different soil horizons

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         soil_moisture

  !     PURPOSE
  !>        \brief Soil moisture in different soil horizons

  !>        \details Infiltration \f$I\f$ from one layer \f$k-1\f$ to the next \f$k\f$ on
  !>        pervious areas is calculated as (omit \f$t\f$)
  !>        \f[ I[k] = I[k-1] (\theta[k] / \theta_{sat}[k])^{\beta[k]} \f]
  !>        Then soil moisture can be calculated as (omit \f$k\f$)
  !>        \f[ \theta[t] = \theta[t-1] + I[t] - \mathit{ET}[t] \f]
  !>        with \f$ \mathit{ET} \f$ (omit \f$[k,t]\f$) being \f$ \mathit{PET} \f$ if
  !>        \f$ \theta > \theta_\mathit{thresh} \f$ else
  !>        \f[ \mathit{ET} = f_\mathrm{roots} \mathit{PET} \f]
  !>        and
  !>        \f[ f_\mathrm{roots} = \frac{\theta - \theta_\mathit{pwp}}{\theta_\mathit{fc} - \theta_\mathit{pwp}} \f]

  !     CALLING SEQUENCE
  !         subroutine soil_moisture(frac_sealed, water_thresh_sealed, pet, &
  !             evap_coeff, soil_moist_sat, frac_roots, soil_moist_FC, wilting_point, &
  !             soil_moist_exponen, aet_canopy, prec_effec, runoff_sealed, storage_sealed, &
  !             infiltration, soil_moist, aet, aet_sealed)

  !     INTENT(IN)

  !>        \param[in] "real(dp)               :: frac_sealed"
  !>                                              Fraction of sealed area
  !>        \param[in] "real(dp)               :: water_thresh_sealed"
  !>                                              Threshhold water depth in impervious areas [mm/s]
  !>        \param[in] "real(dp)               :: pet"
  !>                                              Reference evapotranspiration [mm/s]
  !>        \param[in] "real(dp)               :: evap_coeff"
  !>                                              Evaporation coefficent for free-water surface of that current month
  !>        \param[in] "real(dp), dimension(:) :: soil_moist_sat"
  !>                                              Saturation soil moisture for each horizon [mm]
  !>        \param[in] "real(dp), dimension(:) :: frac_roots"
  !>                                              Fraction of Roots in soil horizon
  !>        \param[in] "real(dp), dimension(:) :: soil_moist_FC"
  !>                                              Soil moisture below which actual ET is reduced [mm]
  !>        \param[in] "real(dp), dimension(:) :: wilting_point"
  !>                                              Permanent wilting point for each horizon [mm]
  !>        \param[in] "real(dp), dimension(:) :: soil_moist_exponen"
  !>                                              Exponential parameter to how non-linear is the soil water retention
  !>        \param[in] "real(dp)               :: aet_canopy"
  !>                                              Actual ET from canopy [mm/s]

  !     INTENT(INOUT)
  !>        \param[in,out] "real(dp)               :: prec_effec"       Effective precipitation (rain + snow melt) [mm]
  !>        \param[in,out] "real(dp)               :: runoff_sealed"    Direct runoff from impervious areas
  !>        \param[in,out] "real(dp)               :: storage_sealed"   Retention storage of impervious areas
  !>        \param[in,out] "real(dp), dimension(:) :: infiltration"     Recharge, infiltration intensity or
  !>                                                                    effective precipitation of each horizon [mm/s]
  !>        \param[in,out] "real(dp), dimension(:) :: soil_moist"       Soil moisture of each horizon [mm]

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:) :: aet"                 actual ET [mm/s]
  !>        \param[out] "real(dp)               :: aet_sealed"          actual ET from free-water surfaces,
  !>                                                                    i.e impervious cover [mm/s]

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Dec 2012
  !         Modified RK, July 2013 - A Mosiac apporach is implemented for processes accounted
  !                                  within the permeamble & impervious area. Precipitation and 
  !                                  effective PET intensity are same for both areas.
  !                                - changes made for variables "water_thresh_sealed" when it becomes
  !                                  zero

  subroutine soil_moisture(frac_sealed, water_thresh_sealed, pet, &
       evap_coeff, soil_moist_sat, frac_roots, soil_moist_FC, wilting_point, &
       soil_moist_exponen, aet_canopy, prec_effec, runoff_sealed, storage_sealed, &
       infiltration, soil_moist, aet, aet_sealed)

    use mo_constants, only: eps_dp

    implicit none

    ! Intent variables
    real(dp),                                    intent(in)    :: frac_sealed         ! fraction of sealed area
    real(dp),                                    intent(in)    :: water_thresh_sealed ! Threshhold water depth in impervious 
    !                                                                                 ! areas [mm/s]
    real(dp),                                    intent(in)    :: pet                 ! Reference evapotranspiration [mm/s]
    real(dp),                                    intent(in)    :: evap_coeff          ! Evaporation coefficent for free-water 
    !                                                                                 ! surface of that current month
    real(dp), dimension(:),                      intent(in)    :: soil_moist_sat      ! Saturation soil moisture 
    !                                                                                 ! for each horizon [mm]
    real(dp), dimension(:),                      intent(in)    :: frac_roots          ! Fraction of Roots in soil horizon
    real(dp), dimension(:),                      intent(in)    :: soil_moist_FC       ! Soil moisture below which actual ET 
    !                                                                                 ! is reduced [mm]
    real(dp), dimension(:),                      intent(in)    :: wilting_point       ! Permanent wilting point 
    !                                                                                 ! for each horizon [mm]
    real(dp), dimension(:),                      intent(in)    :: soil_moist_exponen  ! Exponential parameter to how non-linear 
    !                                                                                 ! is the soil water retention
    real(dp),                                    intent(in)    :: aet_canopy          ! actual ET from canopy [mm/s]
    real(dp),                                    intent(inout) :: prec_effec          ! Effective precipitation 
    !                                                                                 ! (rain + snow melt) [mm]
    real(dp),                                    intent(inout) :: runoff_sealed       ! Direct runoff from impervious areas
    real(dp),                                    intent(inout) :: storage_sealed      ! Retention storage of impervious areas
    real(dp), dimension(size(soil_moist_sat,1)), intent(inout) :: infiltration        ! Recharge, infiltration intensity or
    !                                                                                 ! effective precipitation 
    !                                                                                 ! for each horizon [mm/s]
    real(dp), dimension(size(soil_moist_sat,1)), intent(inout) :: soil_moist          ! Soil moisture of each horizon [mm]
    real(dp), dimension(size(soil_moist_sat,1)), intent(out)   :: aet                 ! actual ET [mm/s]
    real(dp),                                    intent(out)   :: aet_sealed          ! actual ET from free-water surfaces,
    !                                                                                 ! i.e impervious cover [mm/s]

    ! Local variables
    integer(i4) :: hh              ! counter
    real(dp)    :: prec_effec_soil ! Effective Prec or infiltration from above
    real(dp)    :: frac_runoff     ! Runoof fraction
    real(dp)    :: tmp             ! temporary variable for misc use

    ! ----------------------------------------------------------------
    ! IMPERVIOUS COVER PROCESS
    ! ----------------------------------------------------------------
    runoff_sealed = 0.0_dp
    aet_sealed    = 0.0_dp

    if (frac_sealed > 0.0_dp) then
       tmp = storage_sealed + prec_effec

       if (tmp > water_thresh_sealed) then
          runoff_sealed  = tmp - water_thresh_sealed
          storage_sealed = water_thresh_sealed
       else
          runoff_sealed  = 0.0_dp
          storage_sealed = tmp
       end if

       ! aET from sealed area is propotional to the available water content
       if(water_thresh_sealed .gt. eps_dp ) then 
          aet_sealed = ( pet / evap_coeff - aet_canopy) * (storage_sealed / water_thresh_sealed)
          ! numerical problem
          if (aet_sealed .lt. 0.0_dp) aet_sealed = 0.0_dp
       else
          aet_sealed = huge(1.0_dp)
       end if

       ! sealed storage updata
       if (storage_sealed .gt. aet_sealed) then
          storage_sealed = storage_sealed - aet_sealed
       else
          aet_sealed     = storage_sealed
          storage_sealed = 0.0_dp
       end if

    end if

    ! ----------------------------------------------------------------
    ! N-LAYER SOIL MODULE
    ! ----------------------------------------------------------------
    aet(:)          = 0.0_dp
    infiltration(:) = 0.0_dp

    ! for 1st layer input is prec_effec
    prec_effec_soil = prec_effec
    
    do hh = 1, size(soil_moist_sat,1) ! nHorizons
       ! input for other layers is the infiltration from its immediate upper layer will be input
       if (hh .NE. 1) prec_effec_soil = infiltration(hh-1)

       !  start processing for soil moisture process
       !  BASED ON SMs as its upper LIMIT

       if (soil_moist(hh) > soil_moist_sat(hh)) then
          infiltration(hh) = prec_effec_soil
       else
          ! to avoid underflow -- or numerical errors
          if(soil_moist(hh) > eps_dp) then
             !frac_runoff = (soil_moist(hh) / soil_moist_sat(hh))**soil_moist_exponen(hh)
             frac_runoff = exp(soil_moist_exponen(hh)*log(soil_moist(hh)/soil_moist_sat(hh)))
          else
             frac_runoff = 0.0_dp
          end if
          tmp = prec_effec_soil * (1.0_dp - frac_runoff)

          if ( (soil_moist(hh) + tmp) > soil_moist_sat(hh) ) then
             infiltration(hh) = prec_effec_soil + ( soil_moist(hh) - soil_moist_sat(hh) )
             soil_moist(hh)   = soil_moist_sat(hh)
          else
             infiltration(hh) = prec_effec_soil - tmp
             soil_moist(hh)   = soil_moist(hh)  + tmp
          end if
       end if

       !             aET calculations

       !  Satisfying ET demand sequentially from top to the bottom layer
       !  Note that the potential ET for the first soil layer is reduced after
       !  satisfying ET demands of the canopy surface

       aet(hh) = pet - aet_canopy                                                     ! First layer
       if (hh /= 1) aet(hh) = aet(hh) - sum(aet(1:hh-1), mask=(aet(1:hh-1) > 0.0_dp)) ! remaining layers

       ! estimate fraction of ET demand based on root fraction and SM status

       !    SM >= FC
       if ( soil_moist(hh) >= soil_moist_FC(hh) ) then
          tmp = frac_roots(hh)
          ! PW < SM < FC
       else if ( (soil_moist(hh) < soil_moist_FC(hh)) .AND.  &
            (soil_moist(hh) > wilting_point(hh))        ) then
          tmp = frac_roots(hh) * (soil_moist(hh) - wilting_point(hh)) / (soil_moist_FC(hh) - wilting_point(hh))
          ! SM <= PW
       else if ( soil_moist(hh) <= wilting_point(hh) ) then
          tmp = 0.0_dp
       else
          stop 'Error soil_moisture: tmp used uninitialised.'
       end if
       aet(hh) = aet(hh) * tmp
       ! avoid numerical error
       if(aet(hh) < 0.0_dp) aet(hh) = 0.0_dp

       ! reduce SM state
       if(soil_moist(hh) > aet(hh)) then
          soil_moist(hh) = soil_moist(hh) - aet(hh)
       else
          aet(hh)        = soil_moist(hh) - eps_dp
          soil_moist(hh) = eps_dp
       end if

       ! avoid numerical error of underflow
       if(soil_moist(hh) < eps_dp) soil_moist(hh) = eps_dp

    end do ! hh

  
  end subroutine soil_moisture

END MODULE mo_soil_moisture
