!> \file mo_lake_evaporation.f90

!> \brief lake evaporation module

!> \details  lake evaporation module for calculating the evaporative loss
!>            from the surface of the lake. All new methods of estimating
!>            lake evaporation are to be placed here.

!> \authors Pallav K Shrestha, Stephan Thober, Luis Samaniego

!> \date March 2019

! Modifications:
! Pallav Shrestha, Oct 2026 - v6 migration (minimal natural lake, Penman + Priestley-Taylor):
!     * ported from v5.13 mo_mlm_evaporation - lake_evap_penman1948 and
!       lake_evap_priestleytaylor1972 plus their shared helpers; the method is selected by
!       the evaporation_method argument (1 = Penman, 2 = Priestley-Taylor). The Harbeck,
!       Linacre and de Bruin methods and the evaporation-predictor update are out of minimal
!       scope and not ported. (Priestley-Taylor needs only temperature and downward radiation,
!       so it runs on the test-domain meteo without windspeed/dew-point fixtures.)
!     * SI units kept throughout (the Penman kernel was already pure SI in v5).
!     * daily rate, disaggregated by the caller: the Penman combination equation
!       (Shuttleworth 1993, Handbook of Hydrology Eq. 4.2.30) is a DAILY formulation -
!       its empirical wind function 6.43(1 + 0.536 u) is calibrated for daily open-water
!       evaporation. So the routine keeps the v5 daily integration (DaySecs) and returns
!       mm d-1; the daily -> model-step temporal disaggregation is done by the caller with
!       the v6 meteo tools, matching how mo_meteo_container::meteo_update_pet treats the
!       Penman-Monteith and the other PET methods (all computed daily, then disaggregated).
!     * no global-variable connections: the run-configuration switches (varying_albedo,
!       albedo_method, varying_lambda, varying_psychro) and the input windspeed height
!       are threaded as arguments instead of read from mo_mlm_global_variables. The
!       empirical Penman/albedo/windspeed constants continue to come from mo_mlm_constants
!       (recreated in the v6 tree), as in v5; generic physical constants come from
!       mo_constants (forces).
! Pallav Shrestha, Oct 2026 - added a calibration-free Priestley-Taylor method
!     (evaporation_method = 3) following Su and Singh (2023, WRR), adapted to the open-water
!     Penman formulation. The constant Priestley-Taylor coefficient (method 2) is replaced by an
!     analytical coefficient alpha = 1 + (gamma/Delta) 6.43(1 + 0.536 u2)(es - e)/(Rn + Ah),
!     evaluated with a CLIMATOLOGICAL windspeed u2 (not the instantaneous wind) and applied to the
!     live equilibrium term. Decoupling alpha from the per-step wind keeps evaporation robust to
!     noisy instantaneous wind, while alpha still varies with radiation, temperature and humidity
!     (so it spans a range about ~1.26 instead of being fixed). Like Penman (method 1) it needs
!     temperature, downward short- and long-wave radiation, windspeed and dew-point temperature.

MODULE mo_lake_evaporation

  USE mo_kind,   only : i4, i8, dp
  IMPLICIT NONE
  private


  PUBLIC  :: lake_evaporation         ! wrapper for lake evaporation process (Penman or Priestley-Taylor)
  PRIVATE :: lake_evap_penman1948     ! Penman, using temperatures, downward radiations and windspeed
  PRIVATE :: lake_evap_priestleytaylor1972  ! Priestley-Taylor, using temperatures and downward radiations
  PRIVATE :: lake_evap_cfpriestleytaylor2023 ! calibration-free Priestley-Taylor, Su and Singh (2023) open-water adaptation
  PRIVATE :: windspeed_correction     ! corrects windspeed for height of recording
  PUBLIC  :: albedo_water             ! estimate of albedo for water for given latitude and doy, at noon
  PRIVATE :: hourangle                ! estimate of solar hour angle in degrees
  PUBLIC  :: spc_heat_vap_water       ! estimate of specific heat of vaporization of water for given surface water temperature
  PUBLIC  :: psychro_constant         ! estimate of psychrometric constant for given elevation and surface water temperature
  PUBLIC  :: atmospheric_pressure     ! estimate of atmospheric pressure for given elevation

  ! ------------------------------------------------------------------

CONTAINS



  ! ------------------------------------------------------------------

  !     NAME
  !         lake_evaporation

  !     PURPOSE
  !>        \brief wrapper for lake evaporation process

  !>        \details wrapper for lake evaporation process. Minimal natural-lake scope
  !>                 evaluates the Penman (1948) method only; the water albedo is either
  !>                 held constant or estimated for the current hour and day of year.

  !    CALLING SEQUENCE

  !     INTENT(IN)

  !     INTENT(IN), OPTIONAL

  !     RETURN

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE


  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Jan 2022

  ! Modifications:
  ! Pallav Shrestha, Oct 2026 - v6 migration: reduced to the Penman (1948) path; the
  !     run-configuration switches and the input windspeed height are arguments. The
  !     daily evaporation is returned (the caller disaggregates to the model step, as
  !     mo_meteo_container does for the land PET methods).



  recursive SUBROUTINE lake_evaporation( evaporation_method, hour, doy, temp, tdew, windspeed, ssrd, strd, &
                                tin, tout, tpre, lake_level, qin, qout, pre, latitude, longitude, &
                                input_windspeed_height, climatological_windspeed, varying_albedo, albedo_method, &
                                varying_lambda, varying_psychro, evap, evap_components)

    use mo_meteo_constants, only : albedo_w_const

    IMPLICIT NONE

    ! lake evaporation method (1 - Penman 1948, 2 - Priestley-Taylor 1972)
    integer(i4), intent(in) :: evaporation_method

    ! hour of the day
    integer(i4), intent(in) :: hour

    ! day of year
    integer(i4), intent(in) :: doy

    ! [degC] average temperature
    real(dp), intent(in) :: temp

    ! [degC] dew point temperature (used by Penman only)
    real(dp), intent(in) :: tdew

    ! [m s-1] average windspeed (used by Penman only)
    real(dp), intent(in) :: windspeed

    ! [W m-2] downward short-wave radiation
    real(dp), intent(in) :: ssrd

    ! [W m-2] downward long-wave radiation
    real(dp), intent(in) :: strd

    ! [degC] temperature of lake inflow
    real(dp), intent(in) :: tin

    ! [degC] temperature of lake outflow
    real(dp), intent(in) :: tout

    ! [degC] temperature of lake precipitation
    real(dp), intent(in) :: tpre

    ! [m asl] lake surface elevation
    real(dp), intent(in) :: lake_level

    ! [mm d-1] rate of lake inflow
    real(dp), intent(in) :: qin

    ! [mm d-1] rate of lake outflow
    real(dp), intent(in) :: qout

    ! [mm d-1] rate of lake precipitation
    real(dp), intent(in) :: pre

    ! [degree] average latitude of the lake
    real(dp), intent(in) :: latitude

    ! [degree] average longitude of the lake
    real(dp), intent(in) :: longitude

    ! [m] height for which the input windspeed data is provided
    real(dp), intent(in) :: input_windspeed_height

    ! [m s-1] climatological (long-term mean) 2 m windspeed, used by the calibration-free Priestley-Taylor
    ! method (3) in its coefficient instead of the instantaneous wind
    real(dp), intent(in) :: climatological_windspeed

    ! use a varying (vs constant) water albedo
    logical, intent(in) :: varying_albedo

    ! method for calculating water albedo (1 - Fresnel, 2 - Grishchenko)
    integer(i4), intent(in) :: albedo_method

    ! use a varying (vs constant) latent heat of vaporisation of water
    logical, intent(in) :: varying_lambda

    ! use a varying (vs constant) psychrometric constant
    logical, intent(in) :: varying_psychro

    ! [mm d-1] lake evaporation (daily; the caller disaggregates to the model step)
    real(dp), intent(out) :: evap

    ! [mm d-1] lake evaporation components from Penman (daily)
    real(dp), dimension(:), intent(out) :: evap_components

    !! local vars

    real(dp) :: alb_w ! albedo of water
    real(dp) :: sha   ! solar hour angle


    ! Water albedo (shared by both methods): constant, or estimated for the current hour/doy.
    if (varying_albedo) then
      sha   = hourangle(longitude, doy, hour)
      alb_w = albedo_water(latitude, doy, sha, albedo_method)
    else
      alb_w = albedo_w_const
    end if

    select case (evaporation_method)

    case (1) ! Penman (1948) - temperatures, downward radiations and windspeed

      evap_components(1:3) = lake_evap_penman1948(lake_level, alb_w, ssrd, strd, temp, tdew, windspeed, &
                                                  tin, tout, tpre, qin, qout, pre, &
                                                  input_windspeed_height, varying_lambda, varying_psychro)
      evap = max( sum(evap_components(1:3)), 0._dp ) ! no condensation: Penman daily evaporation floored at zero

    case (2) ! Priestley-Taylor (1972) - temperatures and downward radiations (no windspeed/VPD)

      evap = lake_evap_priestleytaylor1972(lake_level, alb_w, ssrd, strd, temp, &
                                           tin, tout, tpre, qin, qout, pre, &
                                           varying_lambda, varying_psychro)
      evap_components(1:3) = 0.0_dp ! the 3-term breakdown is a Penman concept; unused for Priestley-Taylor

    case (3) ! calibration-free Priestley-Taylor - analytical coefficient, climatological wind

      evap = lake_evap_cfpriestleytaylor2023(lake_level, alb_w, ssrd, strd, temp, tdew, &
                                              climatological_windspeed, tin, tout, tpre, qin, qout, pre, &
                                              varying_lambda, varying_psychro)
      evap_components(1:3) = 0.0_dp ! the 3-term breakdown is a Penman concept; unused for Priestley-Taylor

    case default

      evap = 0.0_dp
      evap_components(1:3) = 0.0_dp

    end select


  END SUBROUTINE lake_evaporation


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_evap_penman1948

  !     PURPOSE
  !>        \brief calculates the evaporative loss from the surface of the lake
  !>                 using Penman (1948)

  !>        \details calculates the evaporative loss from the surface of the lake
  !>                 using Penman (1948), given temperatures, downward radiations
  !>                 and windspeed

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: lake_level"            [m asl]   lake surface elevation
  !>       \param[in] "real(dp) :: alb"                   [-]       water albedo
  !>       \param[in] "real(dp) :: ssrd"                  [W m-2]   downward short-wave radiation
  !>       \param[in] "real(dp) :: strd"                  [W m-2]   downward long-wave radiation
  !>       \param[in] "real(dp) :: Temp"                  [degC]    average temperature
  !>       \param[in] "real(dp) :: Tdew"                  [degC]    dew point temperature
  !>       \param[in] "real(dp) :: windspeed"            [m s-1]   windspeed
  !>       \param[in] "real(dp) :: Tin"                   [degC]    temperature of lake inflow
  !>       \param[in] "real(dp) :: Tout"                  [degC]    temperature of lake outflow
  !>       \param[in] "real(dp) :: Tpre"                  [degC]    temperature of lake precipitation
  !>       \param[in] "real(dp) :: qin"                   [mm d-1]  rate of lake inflow
  !>       \param[in] "real(dp) :: qout"                  [mm d-1]  rate of lake outflow
  !>       \param[in] "real(dp) :: pre"                   [mm d-1]  rate of lake precipitation
  !>       \param[in] "real(dp) :: input_windspeed_height" [m]     height of the input windspeed
  !>       \param[in] "logical  :: varying_lambda"        constant / varying latent heat of vaporisation
  !>       \param[in] "logical  :: varying_psychro"       constant / varying psychrometric constant

  !     INTENT(IN), OPTIONAL

  !     RETURN
  !         \return real(dp), dimension(3) :: lake_evap_penman1948 &mdash; Penman (1948) lake
  !                 evaporation components [mm d-1] (1 - radiation, 2 - advection, 3 - VPD)

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !
  !           Penman, H. L. (1948). Natural evaporation from open water, bare soil
  !           and grass. Proceedings of the Royal Society of London. Series A,
  !           Mathematical and Physical Sciences, 193(1032), 120-145.
  !
  !           Shuttleworth, W. J. (1993) Evaporation. In: Maidment, D.R., (ed.)
  !           Handbook of Hydrology, New York: McGraw-Hill, 4.1-4.53.

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date February 2022
  !
  ! Modifications:
  ! Pallav Shrestha, Oct 2026 - v6 migration: kept the daily DaySecs integration (the
  !     equation and its 6.43/0.536 wind function are daily, Shuttleworth 1993 Eq. 4.2.30);
  !     varying_lambda, varying_psychro and the input windspeed height passed as arguments
  !     (no globals). The daily result is disaggregated to the model step by the caller.


  recursive FUNCTION lake_evap_penman1948(lake_level, alb, ssrd, strd, Temp, Tdew, windspeed, Tin, Tout, Tpre, &
                                qin, qout, pre, input_windspeed_height, varying_lambda, varying_psychro)

    use mo_constants,     only : Psychro_dp, SpecHeatET_dp, cp_w_dp, DaySecs, sigma_dp, T0_dp
    use mo_meteo_constants, only : penman1948_1, penman1948_2, boltzmann_e, Emissivity_w
    use mo_pet,           only : sat_vap_pressure, slope_satpressure


    IMPLICIT NONE

    real(dp), intent(IN) ::   lake_level, &  ! [m asl]   lake surface elevation
                              alb, &         ! [-]       water albedo
                              ssrd, &        ! [W m-2]   downward short-wave radiation
                              strd, &        ! [W m-2]   downward long-wave radiation
                              Temp, &        ! [degC]    average temperature
                              Tdew, &        ! [degC]    dew point temperature
                              windspeed, &   ! [m s-1]   windspeed
                              Tin, &         ! [degC]    temperature of lake inflow
                              Tout, &        ! [degC]    temperature of lake outflow
                              Tpre, &        ! [degC]    temperature of lake precipitation
                              qin, &         ! [mm d-1]  rate of lake inflow
                              qout, &        ! [mm d-1]  rate of lake outflow
                              pre            ! [mm d-1]  rate of lake precipitation

    ! [m] height for which the input windspeed data is provided
    real(dp), intent(IN) ::   input_windspeed_height

    ! constant / varying latent heat of vaporisation and psychrometric constant
    logical, intent(IN) ::    varying_lambda, varying_psychro


    !! local vars
    real(dp) :: windspeed_2m, satVapPress, ambVapPress, Ah, ssr, str, stro, netrad, SpecHeat_w, Psychro

    ! lake evaporation in [mm per day] and the components (1 - radiation, 2 - advection, 3 - VPD)
    real(dp), dimension(3) :: lake_evap_penman1948



    !----------------------------------------------------------------------------
    ! Lake evaporation by Penman 1948
    !
    !        d                    r      6.43 ( 1 + 0.536 u )
    ! E = -------- (Rn + Ah) + ------- . -------------------- . (es - e)
    !      d + r                d + r              ^
    !
    ! where, E is lake evaporation in mm d-1
    !        d is slope of saturation vapor pressure vs temperature curve, in kPa degC-1
    !        r is Psychrometric constant, in kPa degC-1
    !        ^ is latent heat of vaporization of water, in MJ kg-1 or MJ m-2 mm-1
    !        Rn is the net radiation, in mm d-1
    !        Ah is the advection rate, in mm d-1
    !        u is the wind speed, in m s-1
    !        es is saturation vapor pressure at temperature T, kPa
    !        e is ambient vapor pressure found using dew point temperature Tdew, in kPa
    !
    !        Rn = Sn + Ln
    !
    !        where, Sn (or "ssr") is net shortwave radiation
    !               Ln (or "str") is net longwave radiation
    !
    !               Sn = St (1 - a)
    !
    !               where, St (or "ssrd") is downward shortwave radiation
    !                      a is albedo (for open water)
    !
    !               Ln = Li - Lo
    !
    !               where, Li (or "strd") is downward longwave radiation
    !                      Lo is outgoing longwave radiation
    !
    !                      Boltzmann equation gives
    !
    !                      Lo = e * sigma * T^4
    !
    !                      where, e is emissivity of water, 0.96
    !                             sigma is Stefan-Boltzmann constant, in W m-2 K-4
    !                             T is water temperature, in Kelvin
    !
    !
    !        Ah = rho_w * cp_w * ( qin * Tin - qout * Tout + pre * Tpre )
    !
    !        where, Ah is advection rate, in J m-2 d-1
    !               rho_w is density of water, in kg m-3
    !               cp_w is specific heat capacity of water, in J kg-1 C-1
    !               qin is lake inflow rate, mm d-1
    !               qout is lake outflow rate, mm d-1
    !               pre is lake precipitation rate, mm d-1
    !               Tin is temperature of lake inflow, C
    !               Tout is temperature of lake outflow, C
    !               Tpre is temperature of lake precipitation, C
    !
    ! Note (v6 migration): the equation and its empirical wind function 6.43(1 + 0.536 u)
    ! are daily (Shuttleworth 1993, Eq. 4.2.30), so the daily DaySecs integration is kept
    ! and the routine returns mm d-1. The caller disaggregates this daily value to the
    ! model step, as mo_meteo_container::meteo_update_pet does for the land PET methods.
    !----------------------------------------------------------------------------

    ! Calculate net shortwave radiation (W m-2)
    ssr = ssrd * (1 - alb)

    ! Calculate outgoing longwave radiation from Boltzmann equation (W m-2)
    stro = Emissivity_w * sigma_dp * (Temp + T0_dp) ** boltzmann_e
      ! PS - using air temperature as water surface temperature until
      !      lake temperature procedure is added and connected with river
      !      temperature routing

    ! Calculate net longwave radiation (W m-2)
    str = strd - stro

    ! Calculate net radiation (W m-2)
    netrad = ssr + str

    ! Calculate advection rate (J m-2 d-1)
    Ah = cp_w_dp * 1000._dp * ( qin * Tin - qout * Tout + pre * Tpre )

    ! Converting input windspeed to standard 2m windspeed
    windspeed_2m = windspeed_correction (windspeed, input_windspeed_height)


    ! Specific heat of vaporisation of water
    if (varying_lambda) then
      SpecHeat_w = spc_heat_vap_water(Temp) ! considering water surface temperature to be equal to air temperature
    else
      SpecHeat_w = SpecHeatET_dp
    end if

    ! Estimate Psychrometric constant
    if (varying_psychro) then
      Psychro = psychro_constant(lake_level, Temp) ! considering water surface temperature to be equal to air temperature
    else
      Psychro = Psychro_dp
    end if

    ! Calculate es and e (kPa)
    satVapPress = sat_vap_pressure(Temp) ! svp at water surface temperature (assumed to be tavg)
    ambVapPress = sat_vap_pressure(Tdew) ! avp obtained using svp formula with tdew


    ! "Evaporation but no condensation"
    ! We only apply Penman for non-negative values of netrad
    if ( netrad >= 0._dp ) then

      if ( netrad * DaySecs + Ah  >= 0._dp ) then ! In case Ah is lower than netrad, evap should be cutoff at 0.
                                                  ! Theoretically, could happen for large inflow with cold temperature, but does this happen?

        ! Calculate E components (Penman, 1948) in mm d-1
        lake_evap_penman1948(1) = slope_satpressure(Temp) / ( slope_satpressure(Temp) + Psychro ) &
                                * ( netrad * DaySecs ) / SpecHeat_w
        lake_evap_penman1948(2) = slope_satpressure(Temp) / ( slope_satpressure(Temp) + Psychro ) &
                                * ( Ah ) / SpecHeat_w
        lake_evap_penman1948(3) = ( Psychro / ( slope_satpressure(Temp) + Psychro) ) &
                                * penman1948_1 * ( 1._dp + penman1948_2 * windspeed_2m ) / (SpecHeat_w / 1000000._dp) &
                                * ( satVapPress - ambVapPress )
        ! No condensation: when the air is wetter than the water surface (es - e < 0) the
        ! vapour-transfer term turns negative; clamp it to zero so the aerodynamic path never
        ! adds water to the lake (consistent with the evaporation-only design of the leaves).
        lake_evap_penman1948(3) = max( lake_evap_penman1948(3), 0._dp )
      else
        lake_evap_penman1948(1:3) = 0._dp
      end if

    else
      lake_evap_penman1948(1:3) = 0._dp
    end if


  END FUNCTION lake_evap_penman1948


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_evap_priestleytaylor1972

  !     PURPOSE
  !>        \brief calculates the evaporative loss from the surface of the lake
  !>                 using Priestley and Taylor (1972)

  !>        \details calculates the evaporative loss from the surface of the lake
  !>                 using Priestley and Taylor (1972), given temperatures and
  !>                 downward radiations (no windspeed / vapour-pressure-deficit term)

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: lake_level"            [m asl]   lake surface elevation
  !>       \param[in] "real(dp) :: alb"                   [-]       water albedo
  !>       \param[in] "real(dp) :: ssrd"                  [W m-2]   downward short-wave radiation
  !>       \param[in] "real(dp) :: strd"                  [W m-2]   downward long-wave radiation
  !>       \param[in] "real(dp) :: Temp"                  [degC]    average temperature
  !>       \param[in] "real(dp) :: Tin"                   [degC]    temperature of lake inflow
  !>       \param[in] "real(dp) :: Tout"                  [degC]    temperature of lake outflow
  !>       \param[in] "real(dp) :: Tpre"                  [degC]    temperature of lake precipitation
  !>       \param[in] "real(dp) :: qin"                   [mm d-1]  rate of lake inflow
  !>       \param[in] "real(dp) :: qout"                  [mm d-1]  rate of lake outflow
  !>       \param[in] "real(dp) :: pre"                   [mm d-1]  rate of lake precipitation
  !>       \param[in] "logical  :: varying_lambda"        constant / varying latent heat of vaporisation
  !>       \param[in] "logical  :: varying_psychro"       constant / varying psychrometric constant

  !     INTENT(IN), OPTIONAL

  !     RETURN
  !         \return real(dp) :: lake_evap_priestleytaylor1972 &mdash; lake evaporation from
  !                 Priestley and Taylor (1972) [mm d-1]

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !
  !           Priestley, C. H. B. and Taylor, R. J. (1972) On the assessment of surface heat flux
  !           and evaporation using large-scale parameters. Monthly Weather Review, 100, 81-92.

  !     HISTORY
  !>        \author Pallav Shrestha, Rohini Kumar
  !>        \date March 2022
  !
  ! Modifications:
  ! Pallav Shrestha, Oct 2026 - v6 migration: kept the daily DaySecs integration (returns mm d-1,
  !     the caller disaggregates to the model step); varying_lambda and varying_psychro passed as
  !     arguments (no globals).


  recursive FUNCTION lake_evap_priestleytaylor1972(lake_level, alb, ssrd, strd, Temp, Tin, Tout, Tpre, &
                                         qin, qout, pre, varying_lambda, varying_psychro)

    use mo_constants,     only : Psychro_dp, SpecHeatET_dp, cp_w_dp, DaySecs, sigma_dp, T0_dp
    use mo_meteo_constants, only : ptcoeff, boltzmann_e, Emissivity_w
    use mo_pet,           only : slope_satpressure


    IMPLICIT NONE

    real(dp), intent(IN) ::   lake_level, &  ! [m asl]   lake surface elevation
                              alb, &         ! [-]       water albedo
                              ssrd, &        ! [W m-2]   downward short-wave radiation
                              strd, &        ! [W m-2]   downward long-wave radiation
                              Temp, &        ! [degC]    average temperature
                              Tin, &         ! [degC]    temperature of lake inflow
                              Tout, &        ! [degC]    temperature of lake outflow
                              Tpre, &        ! [degC]    temperature of lake precipitation
                              qin, &         ! [mm d-1]  rate of lake inflow
                              qout, &        ! [mm d-1]  rate of lake outflow
                              pre            ! [mm d-1]  rate of lake precipitation

    ! constant / varying latent heat of vaporisation and psychrometric constant
    logical, intent(IN) ::    varying_lambda, varying_psychro


    !! local vars
    real(dp) :: Ah, ssr, str, stro, netrad, SpecHeat_w, Psychro

    ! lake evaporation in [mm per day]
    real(dp) :: lake_evap_priestleytaylor1972



    !----------------------------------------------------------------------------
    ! Lake evaporation by Priestley and Taylor 1972
    !
    !          d
    ! E = a -------- (Rn + Ah)
    !        d + r
    !
    ! where, E is lake evaporation in mm d-1
    !        a is the Priestley-Taylor coefficient
    !        d is slope of saturation vapor pressure vs temperature curve, in kPa degC-1
    !        r is Psychrometric constant, in kPa degC-1
    !        Rn is the net radiation, in mm d-1
    !        Ah is the advection rate, in mm d-1
    !
    !        Rn = Sn + Ln  (net short-wave + net long-wave radiation), with the outgoing
    !        long-wave Lo = e * sigma * T^4 (Boltzmann), as in the Penman routine.
    !        Ah = rho_w * cp_w * ( qin * Tin - qout * Tout + pre * Tpre ).
    !
    ! Note (v6 migration): the equation is daily (returns mm d-1); the caller disaggregates it
    ! to the model step, matching how mo_meteo_container treats the land PET methods.
    !----------------------------------------------------------------------------

    ! Calculate net shortwave radiation (W m-2)
    ssr = ssrd * (1 - alb)

    ! Calculate outgoing longwave radiation from Boltzmann equation (W m-2)
    stro = Emissivity_w * sigma_dp * (Temp + T0_dp) ** boltzmann_e
      ! PS - using air temperature as water surface temperature until
      !      lake temperature procedure is added and connected with river
      !      temperature routing

    ! Calculate net longwave radiation (W m-2)
    str = strd - stro

    ! Calculate net radiation (W m-2)
    netrad = ssr + str

    ! Calculate advection rate (J m-2 d-1)
    Ah = cp_w_dp * 1000._dp * ( qin * Tin - qout * Tout + pre * Tpre )


    ! Specific heat of vaporisation of water
    if (varying_lambda) then
      SpecHeat_w = spc_heat_vap_water(Temp) ! considering water surface temperature to be equal to air temperature
    else
      SpecHeat_w = SpecHeatET_dp
    end if

    ! Estimate Psychrometric constant
    if (varying_psychro) then
      Psychro = psychro_constant(lake_level, Temp) ! considering water surface temperature to be equal to air temperature
    else
      Psychro = Psychro_dp
    end if


    ! "Evaporation but no condensation"
    ! We only apply Priestley-Taylor for non-negative values of netrad
    if ( netrad >= 0._dp ) then

      if ( netrad * DaySecs + Ah  >= 0._dp ) then ! In case Ah is lower than netrad, evap should be cutoff at 0.

        ! Calculate E by Priestley-Taylor (1972) in mm d-1
        lake_evap_priestleytaylor1972 = ptcoeff * slope_satpressure(Temp) / ( slope_satpressure(Temp) + Psychro ) &
                                * ( netrad * DaySecs + Ah ) / SpecHeat_w
      else
        lake_evap_priestleytaylor1972 = 0._dp
      end if

    else
      lake_evap_priestleytaylor1972 = 0._dp
    end if


  END FUNCTION lake_evap_priestleytaylor1972


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_evap_cfpriestleytaylor2023

  !     PURPOSE
  !>        \brief calibration-free Priestley-Taylor lake evaporation

  !>        \details Estimates open-water (lake) evaporation with the Priestley-Taylor functional
  !>                 form but an analytically derived, calibration-free coefficient, following
  !>                 Su and Singh (2023) adapted to the open-water Penman formulation. The constant
  !>                 Priestley-Taylor coefficient is replaced by
  !>
  !>                     alpha = 1 + (gamma / Delta) * f(u2bar) * (es - e) / (Rn + Ah),
  !>
  !>                 where f(u2bar) = 6.43 (1 + 0.536 u2bar) is the daily open-water wind function
  !>                 evaluated with a CLIMATOLOGICAL 2 m windspeed u2bar rather than the instantaneous
  !>                 wind. Decoupling the coefficient from the per-step wind keeps the evaporation
  !>                 estimate robust to noisy instantaneous wind, while alpha still varies with the
  !>                 available energy (Rn + Ah), the vapour-pressure deficit (es - e), the curve
  !>                 slope Delta and the psychrometric constant gamma. The coefficient is then applied
  !>                 to the equilibrium (radiation + advection) term:
  !>
  !>                     E = alpha * Delta/(Delta + gamma) * (Rn + Ah) / lambda.
  !>
  !>                 alpha is bounded to [1, cfpt_coeff_max]: the lower bound of 1 matches the
  !>                 no-condensation design (the es - e < 0 condensation regime is zeroed as in the
  !>                 Penman and Priestley-Taylor routines), and the upper cap guards the limit
  !>                 (Rn + Ah) -> 0, where the coefficient would otherwise diverge; the available
  !>                 energy used in the alpha denominator is floored at cfpt_avail_floor.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: lake_level"      [m asl]   lake surface elevation
  !>       \param[in] "real(dp) :: alb"             [-]       water albedo
  !>       \param[in] "real(dp) :: ssrd"            [W m-2]   downward short-wave radiation
  !>       \param[in] "real(dp) :: strd"            [W m-2]   downward long-wave radiation
  !>       \param[in] "real(dp) :: Temp"            [degC]    average temperature
  !>       \param[in] "real(dp) :: Tdew"            [degC]    dew point temperature
  !>       \param[in] "real(dp) :: windspeed_clim"  [m s-1]   climatological 2 m windspeed
  !>       \param[in] "real(dp) :: Tin"             [degC]    temperature of lake inflow
  !>       \param[in] "real(dp) :: Tout"            [degC]    temperature of lake outflow
  !>       \param[in] "real(dp) :: Tpre"            [degC]    temperature of lake precipitation
  !>       \param[in] "real(dp) :: qin"             [mm d-1]  rate of lake inflow
  !>       \param[in] "real(dp) :: qout"            [mm d-1]  rate of lake outflow
  !>       \param[in] "real(dp) :: pre"             [mm d-1]  rate of lake precipitation
  !>       \param[in] "logical  :: varying_lambda"            constant / varying latent heat of vaporisation
  !>       \param[in] "logical  :: varying_psychro"           constant / varying psychrometric constant

  !     RETURN
  !>        \return real(dp) :: lake_evap_cfpriestleytaylor2023 &mdash; lake evaporation [mm d-1]

  !     LITERATURE
  !>        Su, Q. and Singh, V. P. (2023) Calibration-free Priestley-Taylor method for reference
  !>        evapotranspiration estimation. Water Resources Research, 59, e2022WR033198.
  !>        Priestley, C. H. B. and Taylor, R. J. (1972) On the assessment of surface heat flux
  !>        and evaporation using large-scale parameters. Mon. Weather Rev., 100, 81-92.

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Oct 2026

  recursive FUNCTION lake_evap_cfpriestleytaylor2023(lake_level, alb, ssrd, strd, Temp, Tdew, windspeed_clim, &
                                            Tin, Tout, Tpre, qin, qout, pre, varying_lambda, varying_psychro)

    use mo_constants,     only : Psychro_dp, SpecHeatET_dp, cp_w_dp, DaySecs, sigma_dp, T0_dp
    use mo_meteo_constants, only : penman1948_1, penman1948_2, boltzmann_e, Emissivity_w, &
                                 cfpt_avail_floor, cfpt_coeff_max
    use mo_pet,           only : sat_vap_pressure, slope_satpressure


    IMPLICIT NONE

    real(dp), intent(IN) ::   lake_level, &     ! [m asl]   lake surface elevation
                              alb, &            ! [-]       water albedo
                              ssrd, &           ! [W m-2]   downward short-wave radiation
                              strd, &           ! [W m-2]   downward long-wave radiation
                              Temp, &           ! [degC]    average temperature
                              Tdew, &           ! [degC]    dew point temperature
                              windspeed_clim, & ! [m s-1]   climatological 2 m windspeed
                              Tin, &            ! [degC]    temperature of lake inflow
                              Tout, &           ! [degC]    temperature of lake outflow
                              Tpre, &           ! [degC]    temperature of lake precipitation
                              qin, &            ! [mm d-1]  rate of lake inflow
                              qout, &           ! [mm d-1]  rate of lake outflow
                              pre               ! [mm d-1]  rate of lake precipitation

    ! constant / varying latent heat of vaporisation and psychrometric constant
    logical, intent(IN) ::    varying_lambda, varying_psychro


    !! local vars
    real(dp) :: satVapPress, ambVapPress, Ah, ssr, str, stro, netrad, SpecHeat_w, Psychro
    real(dp) :: avail_energy, aero_group, alpha, slope

    ! lake evaporation in [mm per day]
    real(dp) :: lake_evap_cfpriestleytaylor2023



    !----------------------------------------------------------------------------
    ! Calibration-free Priestley-Taylor, Su and Singh (2023), open-water adaptation
    !
    !              d
    ! E = alpha -------- (Rn + Ah)
    !            d + r
    !
    !                  r   6.43 ( 1 + 0.536 u2bar ) ( es - e )
    ! alpha = 1  +  --- . --------------------------------------  , clamped to [1, cfpt_coeff_max]
    !                d                  Rn + Ah
    !
    ! where, E is lake evaporation in mm d-1
    !        alpha is the calibration-free Priestley-Taylor coefficient
    !        d is slope of saturation vapour pressure vs temperature curve, in kPa degC-1
    !        r is psychrometric constant, in kPa degC-1
    !        u2bar is the CLIMATOLOGICAL 2 m windspeed (not the instantaneous wind)
    !        es - e is the vapour-pressure deficit, in kPa
    !        Rn is the net radiation and Ah the advection, as in the Penman and Priestley-Taylor
    !        routines (Rn = Sn + Ln, outgoing Lo = e sigma T^4; Ah = rho_w cp_w (qin Tin - qout Tout + pre Tpre)).
    !
    ! The coefficient is formed from the SAME quantities as the Penman routine, so that it reproduces
    ! the Penman combination equation exactly if evaluated with the instantaneous wind; using the
    ! climatological wind in alpha is the deliberate decoupling that makes the method distinct from,
    ! and more robust than, Penman.
    !
    ! Note (v6 migration): the equation is daily (returns mm d-1); the caller disaggregates it to the
    ! model step, matching how mo_meteo_container treats the land PET methods.
    !----------------------------------------------------------------------------

    ! Calculate net shortwave radiation (W m-2)
    ssr = ssrd * (1 - alb)

    ! Calculate outgoing longwave radiation from Boltzmann equation (W m-2)
    stro = Emissivity_w * sigma_dp * (Temp + T0_dp) ** boltzmann_e
      ! PS - using air temperature as water surface temperature until
      !      lake temperature procedure is added and connected with river
      !      temperature routing

    ! Calculate net longwave radiation (W m-2)
    str = strd - stro

    ! Calculate net radiation (W m-2)
    netrad = ssr + str

    ! Calculate advection rate (J m-2 d-1)
    Ah = cp_w_dp * 1000._dp * ( qin * Tin - qout * Tout + pre * Tpre )


    ! Specific heat of vaporisation of water
    if (varying_lambda) then
      SpecHeat_w = spc_heat_vap_water(Temp) ! considering water surface temperature to be equal to air temperature
    else
      SpecHeat_w = SpecHeatET_dp
    end if

    ! Estimate psychrometric constant
    if (varying_psychro) then
      Psychro = psychro_constant(lake_level, Temp) ! considering water surface temperature to be equal to air temperature
    else
      Psychro = Psychro_dp
    end if

    ! Calculate es and e (kPa)
    satVapPress = sat_vap_pressure(Temp) ! svp at water surface temperature (assumed to be tavg)
    ambVapPress = sat_vap_pressure(Tdew) ! avp obtained using svp formula with tdew

    slope = slope_satpressure(Temp)


    ! "Evaporation but no condensation"
    ! We only apply the method for non-negative net radiation and non-negative available energy
    if ( netrad >= 0._dp ) then

      if ( netrad * DaySecs + Ah  >= 0._dp ) then ! In case Ah is lower than netrad, evap should be cutoff at 0.

        ! Available energy A = Rn + Ah (J m-2 d-1), floored to keep the coefficient finite as A -> 0
        avail_energy = netrad * DaySecs + Ah
        avail_energy = max( avail_energy, cfpt_avail_floor )

        ! Aerodynamic group of the coefficient, using the CLIMATOLOGICAL wind. The 6.43(1+0.536 u)
        ! wind function carries the kPa -> MJ m-2 mm-1 scaling (1e6) that the Penman VPD term uses,
        ! so that, with the instantaneous wind, alpha would reproduce the full Penman evaporation.
        aero_group = ( Psychro / slope ) &
                   * penman1948_1 * ( 1._dp + penman1948_2 * windspeed_clim ) * 1000000._dp &
                   * ( satVapPress - ambVapPress ) / avail_energy

        ! Calibration-free coefficient, clamped to [1, cfpt_coeff_max]
        alpha = 1._dp + aero_group
        alpha = min( max( alpha, 1._dp ), cfpt_coeff_max )

        ! Apply alpha to the equilibrium (radiation + advection) term, in mm d-1
        lake_evap_cfpriestleytaylor2023 = alpha * slope / ( slope + Psychro ) &
                                * ( netrad * DaySecs + Ah ) / SpecHeat_w
      else
        lake_evap_cfpriestleytaylor2023 = 0._dp
      end if

    else
      lake_evap_cfpriestleytaylor2023 = 0._dp
    end if


  END FUNCTION lake_evap_cfpriestleytaylor2023


  ! ------------------------------------------------------------------

  !     NAME
  !         windspeed_correction

  !     PURPOSE
  !>        \brief corrects windspeed for the height of recording

  !>        \details corrects windspeed to the standard 2 m height using
  !>                 Hellmann's exponential law

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: windspeed_in"            [m s-1]  input windspeed
  !>       \param[in] "real(dp) :: input_windspeed_height"  [m]      height of the input windspeed

  !     RETURN
  !         \return real(dp) :: windspeed_correction &mdash; windspeed corrected to 2 m [m s-1]

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date February 2022
  !
  ! Modifications:
  ! Pallav Shrestha, Oct 2026 - v6 migration: the input windspeed height is an argument
  !     (was mo_mlm_global_variables::in_windspeed_ht).

  elemental pure FUNCTION windspeed_correction(windspeed_in, input_windspeed_height)

    USE mo_meteo_constants,             ONLY: std_windspeed_ht, Hellmann_exp

    IMPLICIT NONE

    real(dp), intent(IN) ::   windspeed_in            ! [m s-1]  input windspeed
    real(dp), intent(IN) ::   input_windspeed_height  ! [m]      height of the input windspeed

    ! corrected windspeed
    real(dp) :: windspeed_correction


    ! Converting input windspeed to standard 2m windspeed [Hellmann exponential law]

    windspeed_correction = windspeed_in * ( std_windspeed_ht / input_windspeed_height) ** Hellmann_exp


  END FUNCTION windspeed_correction


  ! ------------------------------------------------------------------

  !     NAME
  !         albedo_water

  !     PURPOSE
  !>        \brief estimates albedo of still open water surface

  !>        \details estimates albedo of still open water surface in fraction of 1

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp)    ::  lat"                       [degree]  average lake latitude
  !>       \param[in] "integer(i4) ::  doy"                       [-]       day of year
  !>       \param[in] "real(dp)    ::  sha"                       [degree]  solar hour angle
  !>       \param[in] "integer(i4) ::  method"                    [-]       method for albedo estimation

  !     RETURN
  !         \return real(dp) :: albedo_water &mdash; albedo of still open water

  !     LITERATURE
  !
  !            Cogley, J. G. (1979). The albedo of water as a function of latitude.
  !            Monthly Weather Review, 107(6), 775-781.
  !            https://doi.org/10.1175/1520-0493(1979)107<0775:TAOWAA>2.0.CO;2
  !
  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date February 2022
  !
  ! Modifications:
  !

  elemental pure FUNCTION albedo_water(lat, doy, sha, method)

    USE mo_constants,       ONLY: YearDays, PIO2_dp, TWOPI_D, deg2rad_dp, rad2deg_dp
    USE mo_meteo_constants, ONLY: DuffieDelta1, DuffieDelta2
    USE mo_meteo_constants,   ONLY: RefracIndex_w, fresnel_c, albedo_w_low, albedo_w_high, &
                                  grishchenko_k, grishchenko_l, grishchenko_s, grishchenko_c


    IMPLICIT NONE

    real(dp),     intent(IN) ::   lat            ! [degree]  average lake latitude
    integer(i4),  intent(IN) ::   doy            ! [-]       day of year
    real(dp),     intent(IN) ::   sha            ! [degree]  solar hour angle
    integer(i4),  intent(IN) ::   method         ! [-]       method for albedo estimation

    ! estimated water albedo
    real(dp) :: albedo_water

    ! local vars
    real(dp) :: delta, Z, r, E



    ! declination of the sun above the celestial equator (radians)
    delta = DuffieDelta1 * sin(TWOPI_D * doy / YearDays - DuffieDelta2)

    ! Solar elevation angle (radians) (complementary angle to solar zenith angle)
    E = asin( sin(lat*deg2rad_dp) * sin(delta) + cos(lat*deg2rad_dp) * cos(delta) * cos(sha*deg2rad_dp) )

    ! Adjust E if sun has already set
    if (E < 0._dp) E = 0._dp

    ! Solar Zenith angle (radians)
    Z = PIO2_dp - E


    ! Select method
    select case (method)
    case (1) ! Fresnel

      ! Get angle of refraction of water from Snell's Law (radians)
      r = asin(sin(Z) / RefracIndex_w)

      ! Estimate water albedo from Fresnel's formula
      albedo_water = fresnel_c * ( ( sin(Z - r)**2 / sin(Z + r)**2 ) + &
                                   ( tan(Z - r)**2 / tan(Z + r)**2 ) )

      ! Water albedo in fraction
      albedo_water = albedo_water / 100._dp

      ! Constrain the estimated albedo
      if (albedo_water .lt. albedo_w_low)   albedo_water = albedo_w_low
      if (albedo_water .gt. albedo_w_high)  albedo_water = albedo_w_high


    case (2) ! Grishchenko

      ! Solar elevation angle (degrees)
      E = E * rad2deg_dp

      ! Estimate water albedo from Weibul Fit on Grishchenko (1959)
      albedo_water = grishchenko_s * grishchenko_k / grishchenko_l * (E/10._dp/grishchenko_l) ** (grishchenko_k - 1._dp)
      albedo_water = albedo_water * exp(-(E/10._dp/grishchenko_l) ** (grishchenko_k))
      albedo_water = albedo_water + grishchenko_c
      albedo_water = albedo_water / 100._dp

    end select


  END FUNCTION albedo_water


  ! ------------------------------------------------------------------

  !     NAME
  !         hourangle

  !     PURPOSE
  !>        \brief estimates solar hour angle in degrees

  !>        \details estimates solar hour angle in degrees

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp)    ::  lon"      [degree]  average lake longitude
  !>       \param[in] "integer(i4) ::  doy"      [-]       day of year
  !>       \param[in] "integer(i4) ::  hour"     [-]       hour of day

  !     RETURN
  !         \return real(dp) :: hourangle &mdash; solar hour angle

  !     LITERATURE
  !
  !                Meeus, J. (1991). Astronomical algorithms. Richmond.
  !
  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date December 2022
  !
  ! Modifications:
  !

  elemental pure FUNCTION hourangle(lon, doy, hour)

    USE mo_constants,                 ONLY: YearDays, DayHours, TWOPI_D
    USE mo_meteo_constants,             ONLY: meeus_1, meeus_2, meeus_3, meeus_4, meeus_5, meeus_6

    IMPLICIT NONE

    real(dp),     intent(IN) ::   lon            ! [degree]  average lake longitude
    integer(i4),  intent(IN) ::   doy            ! [-]       day of year
    integer(i4),  intent(IN) ::   hour           ! [-]       hour of the day

    ! estimated solar hour angle
    real(dp) :: hourangle

    ! local vars
    real(dp) :: fy, eot, toff, tz, tst


    ! Get fractional year (radians)
    fy    = TWOPI_D / YearDays * (doy - 1 + (hour - 0.5_dp * DayHours)/DayHours )
    ! PS - might need check for leap years and use 366 instead of YearDays.

    ! Calculate equation of time (minutes)
    eot   = meeus_1 * (meeus_2 + meeus_3 * cos(fy) + meeus_4 * sin(fy) + meeus_5 * cos(2*fy) + meeus_6 * sin(2*fy))

    ! Estimate time zone
    tz    = nint(lon / 15)

    ! Calculate time offset (minutes)
    toff  = eot + 4 * lon - 60 * tz

    ! True solar time (minutes)
    tst   = hour * 60 + toff

    ! Solar hour angle (degrees)
    hourangle = tst / 4 - 180


  END FUNCTION hourangle


  ! ------------------------------------------------------------------

  !     NAME
  !         spc_heat_vap_water

  !     PURPOSE
  !>        \brief estimates specific heat of vaporization of water

  !>        \details estimates specific heat of vaporization of water

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp) ::  T"                          [degC]    water surface temperature

  !     RETURN
  !         \return real(dp) :: spc_heat_vap_water &mdash; specific heat of vaporization of water

  !     LITERATURE
  !
  !           Shuttleworth, W. J. (1993) Evaporation. In: Maidment, D.R., (ed.)
  !           Handbook of Hydrology, New York: McGraw-Hill, 4.1-4.53.
  !
  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date March 2022
  !
  ! Modifications:
  !

  elemental pure FUNCTION spc_heat_vap_water(T)

    USE mo_meteo_constants,   ONLY: lat_heat_vap_c1, lat_heat_vap_c2

    IMPLICIT NONE

    real(dp),     intent(IN) ::   T        ! [degC]    water surface temperature

    ! estimated specific heat of vaporisation of water
    real(dp) :: spc_heat_vap_water


    !----------------------------------------------------------------------------
    !
    ! Latent heat of vaporization of water has virtually linear relationship
    ! with Surface temperature for the range of [0, 100] deg C
    !
    ! ^ = 2.501 - 0.002361 T
    !
    ! where, ^ is Latent heat of vaporization of water, in MJ kg-1
    !
    ! s = ^ * rho_w
    !
    ! where, s is the specific heat of vaporization of water, in MJ m-2 mm-1
    !        rho_w is the density of water, in kg m-3
    !
    !        rho_w = 1000 kg m-3 = 1 kg m-2 mm-1
    !
    !----------------------------------------------------------------------------


    ! Estimate Latent heat of vaporization or specific heat of vaporization of water
    spc_heat_vap_water = (lat_heat_vap_c1 - lat_heat_vap_c2 * T) * 1000000.0_dp


  END FUNCTION spc_heat_vap_water



  ! ------------------------------------------------------------------

  !     NAME
  !         psychro_constant

  !     PURPOSE
  !>        \brief estimates psychrometric constant

  !>        \details estimates psychrometric constant

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp)    ::  elev"         [m asl]    elevation
  !>       \param[in] "real(dp)    ::  T"            [degC]     surface water temperature

  !     RETURN
  !         \return real(dp) :: psychro_constant &mdash; psychrometric constant

  !     LITERATURE
  !
  !           Shuttleworth, W. J. (1993) Evaporation. In: Maidment, D.R., (ed.)
  !           Handbook of Hydrology, New York: McGraw-Hill, 4.1-4.53.
  !
  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date March 2022
  !
  ! Modifications:
  !

  elemental pure FUNCTION psychro_constant(elev, T)

    USE mo_meteo_constants,  ONLY: psychro_constant_c1

    IMPLICIT NONE

    real(dp),     intent(IN) ::   elev         ! [m asl]    elevation
    real(dp),     intent(IN) ::   T            ! [degC]     surface water temperature

    ! estimated psychrometric constant
    real(dp) :: psychro_constant

    !----------------------------------------------------------------------------
    !
    ! Psychrometric constant is not constant and varies strongly with the elevation
    ! of the location through atmospheric pressure and water surface temperature
    ! through latent heat of vaporisation of water as follows -
    !
    !                P
    ! r = 0.0016286 ---
    !                ^
    !
    ! where, r is the psychrometric constant, kPa degC-1
    !        P is the atmospheric pressure, in kPa
    !        ^ is Latent heat of vaporization of water, in MJ kg-1
    !
    !----------------------------------------------------------------------------


    ! Estimate psychrometric constant
    psychro_constant = psychro_constant_c1 * atmospheric_pressure(elev) / (spc_heat_vap_water(T) / 1000000.0_dp)


  END FUNCTION psychro_constant



  ! ------------------------------------------------------------------

  !     NAME
  !         atmospheric_pressure

  !     PURPOSE
  !>        \brief estimates atmospheric pressure from elevation

  !>        \details estimates atmospheric pressure from elevation

  !    CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp)    ::  elev"         [m asl]    elevation

  !     RETURN
  !         \return real(dp) :: atmospheric_pressure &mdash; atmospheric pressure

  !     LITERATURE
  !
  !           Shuttleworth, W. J. (1993) Evaporation. In: Maidment, D.R., (ed.)
  !           Handbook of Hydrology, New York: McGraw-Hill, 4.1-4.53.
  !
  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date March 2022
  !
  ! Modifications:
  !

  elemental pure FUNCTION atmospheric_pressure(elev)

    USE mo_meteo_constants,  ONLY: atm_pressure_c1, atm_pressure_c2, &
                                 atm_pressure_c3, atm_pressure_c4

    IMPLICIT NONE

    real(dp),     intent(IN) ::   elev         ! [m asl]    elevation

    ! estimated atmospheric pressure
    real(dp) :: atmospheric_pressure

    !----------------------------------------------------------------------------
    !
    ! Atmospheric pressure varies with elevation, which can be estimated as -
    !
    !                   ( 293 - 0.0065 z ) ^ 5.256
    !        P = 101.3 (------------------)
    !                   (      293       )
    !
    !        where, z is the elevation in m a.s.l.
    !
    !----------------------------------------------------------------------------


    ! Estimate atmospheric pressure
    atmospheric_pressure = atm_pressure_c1 * ( (atm_pressure_c2 - atm_pressure_c3 * elev) / atm_pressure_c2) ** atm_pressure_c4


  END FUNCTION atmospheric_pressure


END MODULE mo_lake_evaporation
