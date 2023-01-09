!> \file mo_pet.f90
!> \brief   \copybrief mo_pet
!> \details \copydetails mo_pet

!> \brief Module for calculating reference/potential evapotranspiration  [mm d-1]
!> \details This module calculates PET [mm/d] based on one of the methods
!!       - Hargreaves-Samani (1982)
!!       - Priestly-Taylor (1972)
!!       - Penman-Monteith FAO (1998)
!> \authors Matthias Zink, Christoph Schneider, Matthias Cuntz
!> \date Apr 2014
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_pet

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PUBLIC :: pet_hargreaves ! Hargreaves-Samani
  PUBLIC :: pet_priestly   ! Priestley-Taylor
  PUBLIC :: pet_penman     ! Penman-Monteith
  PUBLIC :: slope_satpressure
  PUBLIC :: extraterr_rad_approx
  PUBLIC :: sat_vap_pressure


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        pet_hargreaves

  !    PURPOSE
  !>       \brief Reference Evapotranspiration after Hargreaves

  !>       \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the Hargreaves-Samani
  !>       (1982)
  !>       model for a given cell by applying the equation
  !>       \f[ PET = HarSamCoeff * R_a * (T_{avg} +  HarSamConst) * \sqrt{ T_{max} - T_{min}} \f]
  !>       where \f$ R_a\;[W\;m^{-2}]\f$ is the incoming solar radiation and
  !>       \f$ T_{avg}, T_{max} \f$ and \f$ T_{min}\f$  \f$ [ ^{\circ}C]\f$ are the mean, maximum,
  !>       and minimum daily temperatures at the given day, respectively.
  !>       \note Hargreaves, G.H., and Samani, Z.A. (1982). "Estimating potential evapotranspiration."
  !>       Tech. Note, J. Irrig. and drain. Engrg., ASCE, 108(3):225-230.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: HarSamCoeff" coefficient of Hargreaves-Samani equation [-]
  !>       \param[in] "real(dp) :: HarSamConst" constant    of Hargreaves-Samani equation [-]
  !>       \param[in] "real(dp) :: tavg"        daily men temperature \f$[^{\circ}C]\f$
  !>       \param[in] "real(dp) :: tmax"        daily maximum of temp \f$[^{\circ}C]\f$
  !>       \param[in] "real(dp) :: tmin"        daily minimum of temp \f$[^{\circ}C]\f$
  !>       \param[in] "real(dp) :: latitude"    latitude of the cell for Ra estimation \f$[radians]\f$
  !>       \param[in] "integer(i4) :: doy"      day of year for Ra estimation

  !    RETURN
  !>       \return real(dp) :: pet_hargreaves &mdash; Hargreaves-Samani pot. evapotranspiration [mm d-1]

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION pet_hargreaves(HarSamCoeff, HarSamConst, tavg, tmax, tmin, latitude, doy)

    use mo_constants, only : deg2rad_dp
    use mo_utils, only : LE

    implicit none

    ! coefficient of Hargreaves-Samani equation [-]
    real(dp), intent(in) :: HarSamCoeff

    ! constant    of Hargreaves-Samani equation [-]
    real(dp), intent(in) :: HarSamConst

    ! daily men temperature \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tavg

    ! daily maximum of temp \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tmax

    ! daily minimum of temp \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tmin

    ! latitude of the cell for Ra estimation \f$[radians]\f$
    real(dp), intent(in) :: latitude

    ! day of year for Ra estimation
    integer(i4), intent(in) :: doy

    ! reference evapotranspiration in [mm d-1]
    real(dp) :: pet_hargreaves

    ! tmax-Tmin
    ! correction for shity input data (tmax<tmin) and to avoid numerical errors ! MZMZMZMZ
    real(dp) :: delta_temp


    delta_temp = tmax - tmin
    if(LE(delta_temp, 0.0_dp) .or. LE(tavg, -HarSamConst)) then
      pet_hargreaves = 0.0_dp
    else
      pet_hargreaves = HarSamCoeff * extraterr_rad_approx(doy, deg2rad_dp * latitude) * (tavg + HarSamConst) * sqrt(delta_temp)
    end if

  END FUNCTION pet_hargreaves


  ! ------------------------------------------------------------------

  !    NAME
  !        pet_priestly

  !    PURPOSE
  !>       \brief Reference Evapotranspiration after Priestly-Taylor

  !>       \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the
  !>       Priestly-Taylor (1972) model for every given cell by applying the equation
  !>       \f[ PET = \alpha * \frac{\Delta}{(\gamma + \Delta)} * R_n \f]
  !>       where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation \f$\Delta =  f(T_{avg})\f$ is the slope
  !>       of the saturation-vapour pressure curve and \f$\alpha\f$ is a emperical coefficient.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: PrieTayParam" Priestley-Taylor coefficient \f$ \alpha [-] \f$
  !>       \param[in] "real(dp) :: Rn"           net solar radiation \f$ [W\;m^{-2}] \f$
  !>       \param[in] "real(dp) :: Tavg"

  !    RETURN
  !>       \return real(dp) :: pet_priestly &mdash; Priestley-Taylor pot. evapotranspiration [mm d-1]

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION pet_priestly(PrieTayParam, Rn, tavg)

    use mo_constants, only : DaySecs, Psychro_dp, SpecHeatET_dp

    implicit none

    ! Priestley-Taylor coefficient \f$ \alpha [-] \f$
    real(dp), intent(in) :: PrieTayParam

    ! net solar radiation \f$ [W\;m^{-2}] \f$
    real(dp), intent(in) :: Rn

    real(dp), intent(in) :: Tavg

    ! reference evapotranspiration in [mm d-1]
    real(dp) :: pet_priestly

    ! save slope of saturation vapor pressure curve
    real(dp) :: delta


    delta = slope_satpressure(Tavg) ! slope of saturation vapor pressure curve
    ! in [mm d-1]
    pet_priestly = PrieTayParam * delta / (Psychro_dp + delta) * (Rn * DaySecs / SpecHeatET_dp)

  END FUNCTION pet_priestly


  ! ------------------------------------------------------------------

  !    NAME
  !        pet_penman

  !    PURPOSE
  !>       \brief Reference Evapotranspiration after Penman-Monteith

  !>       \details Calculates the reference evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the
  !>       Penman-Monteith model for every given cell by applying the equation
  !>       \f[ PET = \frac{1}{\lambda}  \cdot
  !>       \frac{\Delta \cdot R_n + \rho \cdot c_p \cdot (e_s-e) \cdot \frac{a_sh}{r_a}}
  !>       {\Delta + \gamma \cdot \frac{a_sh}{a_s} \cdot \left( 1 + \frac{r_s}{r_a} \right) }         \f]
  !>       where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation,
  !>       \f$\Delta\;[kPa\;K^{-1}]\f$ is the slope of the saturation-vapour pressure curve,
  !>       \f$ \lambda\;[MJ\;kg^{-1}] \f$ is the latent heat of vaporization,
  !>       \f$ (e_s-e)\;[kPa] \f$ is the vapour pressure deficit of the air,
  !>       \f$ \rho\;[kg\;m^{-3}] \f$ is the mean atmospheric density,
  !>       \f$ c_p=1005.0\;J\;kg^{-1}\;K^{-1} \f$ is the specific heat of the air,
  !>       \f$ \gamma [kPa\;K^{-1}] \f$ is the psychrometric constant,
  !>       \f$ r_s [s m^{-1}] \f$ is the bulk canopy resistance,
  !>       \f$ r_a [s m^{-1}] \f$ is the aerodynamic resistance,
  !>       \f$ a_s [1] \f$ is the fraction of one-sided leaf area covered by stomata
  !>       (1 if stomata are on one side only, 2 if they are on both sides) and
  !>       \f$ a_{sh} [-] \f$ is the fraction of projected area exchanging sensible heat with the air (2)
  !>       Implementation refers to the so-called Penman-Montheith equation for transpiration.
  !>       Adjusting the arguments \f$ a_{sh} \f$ and \f$ a_s \f$ we obtain the corrected MU equation (for details
  !>       see Schymanski and Or, 2017). If \f$ a_{sh} = 1 = a_s \f$ Penman-Montheith equation for transpiration
  !>       is preserved. For reproducing characteristics of symmetrical amphistomatous leaves use
  !>       \f$ a_{sh} = 2 = a_s \f$, in which case the classic PM equation is only missing a factor
  !>       of 2 in the nominator, as pointed out by Jarvis and McNaughton (1986, Eq. A9).
  !>       These analytical solutions eliminated the non-linearity problem of the saturation vapour pressure curve,
  !>       but they do not consider the dependency of the long-wave component of the soil surface or leaf energy balance
  !>       (\f$ R_l \f$) on soil or leaf temperature (\f$ T_l \f$). We assume that net radiation
  !>       equals the absorbed short-wave radiation, i.e. \f$ R_N = R_s \f$ (p.79 in Monteith and Unsworth, 2013).

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: net_rad"                net radiation \f$[W m^{-2}]\f$
  !>       \param[in] "real(dp) :: tavg"                   average daily temperature \f$[^{\circ}C]\f$
  !>       \param[in] "real(dp) :: act_vap_pressure"       actual vapur pressure \f$[kPa]\f$
  !>       \param[in] "real(dp) :: aerodyn_resistance"     aerodynmaical resistance \f$s\;m^{-1}\f$
  !>       \param[in] "real(dp) :: bulksurface_resistance" bulk surface resistance  \f$s\;m^{-1}\f$
  !>       \param[in] "real(dp) :: a_s"                    fraction of one-sided leaf area covered by stomata \f$1\f$
  !>       \param[in] "real(dp) :: a_sh"                   fraction of projected area exchanging sensible heat with the
  !>       air \f$1\f$

  !    RETURN
  !>       \return real(dp) :: pet_penman &mdash; Reference Evapotranspiration [mm d-1]

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2014

  ! Modifications:
  ! Johannes Brenner Nov 2017 - include arguments a_s and a_sh to enable corrected MU approach
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION pet_penman(net_rad, tavg, act_vap_pressure, aerodyn_resistance, bulksurface_resistance, a_s, &
                                    a_sh)

    use mo_constants, only : DaySecs, Psychro_dp, SpecHeatET_dp, cp0_dp, rho0_dp

    implicit none

    ! net radiation \f$[W m^{-2}]\f$
    real(dp), intent(in) :: net_rad

    ! average daily temperature \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tavg

    ! actual vapur pressure \f$[kPa]\f$
    real(dp), intent(in) :: act_vap_pressure

    ! aerodynmaical resistance \f$s\;m^{-1}\f$
    real(dp), intent(in) :: aerodyn_resistance

    ! bulk surface resistance  \f$s\;m^{-1}\f$
    real(dp), intent(in) :: bulksurface_resistance

    ! fraction of one-sided leaf area covered by stomata \f$1\f$
    real(dp), intent(in) :: a_s

    ! fraction of projected area exchanging sensible heat with the air \f$1\f$
    real(dp), intent(in) :: a_sh

    ! reference evapotranspiration in [mm d-1]
    real(dp) :: pet_penman


    pet_penman = DaySecs / SpecHeatET_dp * & ! conversion factor [W m-2] to [mm d-1]
            (slope_satpressure(tavg) * net_rad + &
                    rho0_dp * cp0_dp * (sat_vap_pressure(tavg) - act_vap_pressure) * a_sh / aerodyn_resistance) / &
            (slope_satpressure(tavg) + Psychro_dp * a_sh / a_s * (1.0_dp + bulksurface_resistance / aerodyn_resistance))

  END FUNCTION pet_penman

  ! ------------------------------------------------------------------

  !    NAME
  !        extraterr_rad_approx

  !    PURPOSE
  !>       \brief Approximation of extraterrestrial radiation

  !>       \details Approximation of extraterrestrial radiation at the top of the atmosphere \f$ R_a \f$
  !>       after Duffie and Beckman (1980).
  !>       \f$ R_a \f$ is converted from \f$ [J\;m^{-2}\;d^{-1}] \f$ in \f$ [mm\;d^{-1}]\f$ .
  !>       \f[ R_a   = \frac{86400}{ \pi \cdot \lambda} \cdot E_0 \cdot
  !>       d_r \cdot (\omega \cdot \sin(latitude) \cdot \sin(\delta) + \cos(latitude) \cdot \cos(\delta) \cdot
  !>       \sin(\omega) \f]
  !>       where \f$ E_0=1367\;J\;m^{-2}\;s^{-1} \f$ is the solar constant and
  !>       It is dependent on the following sub equations:
  !>       The relative distance Earth-Sun:
  !>       \f[ d_r =  1 + 0.033 \cdot \cos \left( \frac{2 \cdot \pi \cdot doy}{365} \right) \f]
  !>       in which doy is the day of the year.
  !>       The solar declination [radians] defined by
  !>       \f[ \delta  = 0.4093 \cdot \sin\left( \frac{2 \cdot \pi \cdot doy}{365} - 1.405 \right) \f]
  !>       The sunset hour angle [radians]:
  !>       \f[ \omega  = \arccos( - \tan(latitude) * \tan(\delta) )  \f]
  !>       <                 \f$ \lambda = 2.45 \cdot 10^6\;J\;m^{-2}\;mm^{-1} \f$ is the latent heat of vaporization.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: doy"   day of year [-]
  !>       \param[in] "real(dp) :: latitude" latitude [rad]

  !    RETURN
  !>       \return real(dp) :: extraterr_rad_approx &mdash; extraterrestrial radiation approximation \f$[W\;m^{-2}]\f$

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION extraterr_rad_approx(doy, latitude)

    use mo_constants, only : DaySecs, YearDays, PI_D, SolarConst_dp, SpecHeatET_dp, TWOPI_D
    use mo_mhm_constants, only : DuffieDelta1, DuffieDelta2, DuffieDr

    implicit none

    ! day of year [-]
    integer(i4), intent(in) :: doy

    ! latitude [rad]
    real(dp), intent(in) :: latitude

    ! extraterrestrial radiation
    real(dp) :: extraterr_rad_approx

    real(dp) :: dr, delta

    real(dp) :: omega

    real(dp) :: arg


    ! inverse relative distance Earth-Sun - correction for eccentricity of Earths orbit around the sun
    dr = 1.0_dp + DuffieDr * cos(TWOPI_D * doy / YearDays)
    ! declination of the sun above the celestial equator in radians
    delta = DuffieDelta1 * sin(TWOPI_D * doy / YearDays - DuffieDelta2)

    ! arccos(x) is only defined between PI and 0 (for x between -1 and 1)
    ! check limits
    arg = - tan(latitude) * tan(delta)
    if(arg .lt. -1.0_dp) arg = -1.0_dp
    if(arg .gt.  1.0_dp) arg = 1.0_dp

    ! sunrise hour angle in radians
    omega = acos(arg)

    ! Ra - converted from [J m-2 d-1] in [mm d-1]
    extraterr_rad_approx = DaySecs / PI_D / SpecHeatET_dp * SolarConst_dp * &
            dr * (omega * sin(latitude) * sin(delta) + cos(latitude) * cos(delta) * sin(omega))

  end FUNCTION extraterr_rad_approx


  ! ------------------------------------------------------------------

  !    NAME
  !        slope_satpressure

  !    PURPOSE
  !>       \brief slope of saturation vapour pressure curve

  !>       \details slope of saturation vapour pressure curve after Tetens
  !>       \f[ \Delta = \frac{0.6108 * e_s(T_a)}{e^(2 \cdot \log(T_a + 237.3))} \f]

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: tavg" average daily temperature \f$[^{\circ}C]\f$

  !    RETURN
  !>       \return real(dp) :: slope_satpressure &mdash;  slope of saturation vapour pressure curve
  !>       \f$[kPa\;K{-1}]\f$

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION slope_satpressure(tavg)

    use mo_mhm_constants, only : satpressureslope1, tetens_c3

    implicit none

    ! average daily temperature \f$[^{\circ}C]\f$
    real(dp), intent(in) :: tavg

    ! slope of saturation vapour pressure curve
    real(dp) :: slope_satpressure


    slope_satpressure = satpressureslope1 * sat_vap_pressure(tavg) / exp(2.0_dp * log(Tavg + tetens_c3))

  END FUNCTION slope_satpressure

  ! ------------------------------------------------------------------

  !    NAME
  !        sat_vap_pressure

  !    PURPOSE
  !>       \brief calculation of the saturation vapour pressure

  !>       \details Calculation of the saturation vapour pressure
  !>       \f[ e_s(T_a) = 0.6108 \cdot \exp \left( \frac{17.27 \cdot T_a}{T_a + 237.3} \right)  \f]

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: tavg" temperature [degC]

  !    RETURN
  !>       \return real(dp) :: sat_vap_pressure &mdash; saturation vapour pressure [kPa]

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure FUNCTION sat_vap_pressure(tavg)

    use mo_mhm_constants, only : tetens_c1, tetens_c2, tetens_c3

    implicit none

    ! temperature [degC]
    real(dp), intent(in) :: tavg

    ! saturation vapour pressure [kPa]
    real(dp) :: sat_vap_pressure


    sat_vap_pressure = tetens_c1 * exp(tetens_c2 * tavg / (tavg + tetens_c3))

  END FUNCTION sat_vap_pressure
  !
END MODULE mo_pet
