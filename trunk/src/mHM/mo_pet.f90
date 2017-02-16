!> \file mo_pet.f90

!> \brief Module for calculating reference/potential evapotranspiration  [mm s-1]

!> \details This module calculates PET [mm/s] based on one of the methods \n
!>          - Hargreaves-Samani (1982) \n
!>          - Priestly-Taylor (1972) \n
!>          - Penman-Monteith FAO (1998) \n

!> \author Matthias Zink, Christoph Schneider, Matthias Cuntz
!> \date   Apr 2014

MODULE mo_pet

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind,      ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pet_hargreaves ! Hargreaves-Samani
  PUBLIC :: pet_priestly   ! Priestley-Taylor
  PUBLIC :: pet_penman     ! Penman-Monteith

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         pet_hargreaves

  !>        \brief Reference Evapotranspiration after Hargreaves

  !>        \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the Hargreaves-Samani (1982) 
  !>                 model for a given cell by applying the equation
  !>                 \f[ PET = HarSamCoeff * R_a * (T_{avg} +  HarSamConst) * \sqrt{ T_{max} - T_{min}} \f]
  !>                 where \f$ R_a\;[W\;m^{-2}]\f$ is the incoming solar radiation and 
  !>                 \f$ T_{avg}, T_{max} \f$ and \f$ T_{min}\f$  \f$ [ ^{\circ}C]\f$ are the mean, maximum,
  !>                 and minimum daily temperatures at the given day, respectively.

  !     INTENT(IN)
  !>        \param[in] "real(dp),    intent(in) :: HarSamCoeff" coefficient of Hargreaves-Samani equation [-]
  !>        \param[in] "real(dp),    intent(in) :: HarSamConst" constant    of Hargreaves-Samani equation [-]
  !>        \param[in] "real(dp),    intent(in) :: tavg"        daily men temperature \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp),    intent(in) :: tmax"        daily maximum of temp \f$[^{\circ}C]\f$ 
  !>        \param[in] "real(dp),    intent(in) :: tmin"        daily minimum of temp \f$[^{\circ}C]\f$
  !>        \param[in] "real(dp),    intent(in) :: latitude"   latitude of the cell for Ra estimation \f$[radians]\f$

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: pet_hargreaves &mdash; Hargreaves-Samani pot. evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         \note Hargreaves, G.H., and Samani, Z.A. (1982). "Estimating potential evapotranspiration."
  !             Tech. Note, J. Irrig. and drain. Engrg., ASCE, 108(3):225-230.

  !     HISTORY
  !>        \author   Matthias Zink
  !>        \date     Dec 2012

  elemental pure FUNCTION pet_hargreaves(HarSamCoeff, HarSamConst, tavg, tmax, tmin, latitude, doy)

    use mo_constants,     only: deg2rad_dp
    use mo_utils,         only: LE

    implicit none

    real(dp),    intent(in) :: HarSamCoeff           !  coefficient of Hargreaves-Samani equation
    real(dp),    intent(in) :: HarSamConst           !  constatnt   of Hargreaves-Samani equation
    real(dp),    intent(in) :: tavg                  !  daily men temperature
    real(dp),    intent(in) :: tmax                  !  daily maximum of temp.
    real(dp),    intent(in) :: tmin                  !  daily minimum of temp.
    real(dp),    intent(in) :: latitude              ! latitude of the cell for Ra estimation
    integer(i4), intent(in) :: doy                   ! day of year for Ra estimation

    real(dp)                :: pet_hargreaves        ! reference evapotranspiration in [mm s-1]

    ! local
    real(dp)                :: delta_temp            ! tmax-Tmin

    ! correction for shity input data (tmax<tmin) and to avoid numerical errors ! MZMZMZMZ
    delta_temp = tmax - tmin
    if( LE(delta_temp, 0.0_dp) .or. LE(tavg, -HarSamConst) ) then
       pet_hargreaves = 0.0_dp
    else
       pet_hargreaves = HarSamCoeff * extraterr_rad_approx(doy, deg2rad_dp * latitude) * (tavg + HarSamConst) * sqrt(delta_temp) 
    end if
    
  END FUNCTION pet_hargreaves


  ! ------------------------------------------------------------------

  !     NAME
  !         pet_priestly

  !>        \brief Reference Evapotranspiration after Priestly-Taylor

  !>        \details Calculates the Reference Evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the 
  !>                 Priestly-Taylor (1972) model for every given cell by applying the equation
  !>                 \f[ PET = \alpha * \frac{\Delta}{(\gamma + \Delta)} * R_n \f]
  !>                 where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation \f$\Delta =  f(T_{avg})\f$ is the slope 
  !>                 of the saturation-vapour pressure curve and \f$\alpha\f$ is a emperical coefficient.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: PrieTayParam" Priestley-Taylor coefficient \f$ \alpha [-] \f$
  !>        \param[in] "real(dp) :: Rn"           net solar radiation \f$ [W\;m^{-2}] \f$
  !>        \param[in] "real(dp) :: Tavg"         daily mean air temperature \f$ [ ^{\circ}C]\f$  

  !     INTENT(INOUT) 
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>         \return real(dp) :: pet_priestly &mdash; Priestley-Taylor pot. evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Priestley, C.H.B., and R.J. Taylor. 1972. On the assessment of surface heat flux and evaporation using
  !>                   large-scale parameters. Mon. Weather Rev., 100:81-82.
  !>         \note ASAE Standards. 1998. EP406.2: heating, cooling, and ventilating greenhouses. St. Joseph, MI, USA.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION pet_priestly(PrieTayParam, Rn, tavg)
!  FUNCTION pet_priestly(PrieTayParam, Rn, tavg)

    use mo_mhm_constants, only: DaySecs
    use mo_constants,     only: Psychro_dp, SpecHeatET_dp 

    implicit none

    real(dp), intent(in) :: PrieTayParam       ! Priestley-Taylor coefficient
    real(dp), intent(in) :: Rn
    real(dp), intent(in) :: Tavg
    real(dp)             :: pet_priestly       ! reference evapotranspiration in [mm s-1]

    real(dp)             :: delta              ! save slope of saturation vapor pressure curve

    delta        = slope_satpressure(Tavg) ! slope of saturation vapor pressure curve
    ! in [mm d-1] 
    pet_priestly = PrieTayParam * delta / (Psychro_dp + delta) * ( Rn * DaySecs / SpecHeatET_dp ) 
   
  END FUNCTION pet_priestly



  ! ------------------------------------------------------------------

  !     NAME
  !         pet_penman

  !>        \brief Reference Evapotranspiration after Penman-Monteith

  !>        \details Calculates the reference evapotranspiration \f$ [mm\;d^{-1}] \f$ based on the 
  !>                 Penman-Monteith model for every given cell by applying the equation
  !>                 \f[ PET = \frac{1}{\lambda}  \cdot
  !>                           \frac{\Delta \cdot R_n + \frac{\rho \cdot c_p \cdot (e_s-e)}{r_a}}
  !>                           {\Delta + \gamma \left( 1 + \frac{r_s}{r_a} \right) }         \f]
  !>                 where \f$R_n\;[W\;m^{-2}]\f$ is the net solar radiation,
  !>                 \f$\Delta\;[kPa\;K^{-1}]\f$ is the slope of the saturation-vapour pressure curve, 
  !>                 \f$ \lambda\;[MJ\;kg^{-1}] \f$ is the latent heat of vaporization, 
  !>                 \f$ (e_s-e)\;[kPa] \f$ is the vapour pressure deficit of the air, 
  !>                 \f$ \rho\;[kg\;m^{-3}] \f$ is the mean atmospheric density, 
  !>                 \f$ c_p=1005.0\;J\;kg^{-1}\;K^{-1} \f$ is the specific heat of the air, 
  !>                 \f$ \gamma [kPa\;K^{-1}] \f$ is the psychrometric constant,
  !>                 \f$ r_s [s m^{-1}] \f$ is the bulk canopy resistance and
  !>                 \f$ r_a [s m^{-1}] \f$ is the aerodynamic resistance.

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: net_rad"                net radiation \f$[W m^{-2}]\f$
  !>        \param[in] "real(dp), intent(in) :: tavg"                   average daily temperature \f$[^{\circ}C]\f$ 
  !>        \param[in] "real(dp), intent(in) :: act_vap_pressure"       actual vapur pressure \f$[kPa]\f$ 
  !>        \param[in] "real(dp), intent(in) :: aerodyn_resistance"     aerodynmaical resistance \f$s\;m^{-1}\f$
  !>        \param[in] "real(dp), intent(in) :: bulksurface_resistance" bulk surface resistance  \f$s\;m^{-1}\f$
  !>        \param[in] "real(dp)             :: pet_penman"             reference evapotranspiration \f$[mm\;s-1]\f$

  !     INTENT(INOUT) 
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>         \return real(dp) :: pet_penman &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>        \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.  

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION pet_penman(net_rad, tavg, act_vap_pressure, aerodyn_resistance, bulksurface_resistance)

    use mo_mhm_constants, only: DaySecs
    use mo_constants,     only: Psychro_dp, SpecHeatET_dp, rho0_dp, cp0_dp 

    implicit none

    real(dp), intent(in) :: net_rad                ! net radiation
    real(dp), intent(in) :: tavg                   ! average daily temperature
    real(dp), intent(in) :: act_vap_pressure       ! actual vapur pressure
    real(dp), intent(in) :: aerodyn_resistance     ! aerodynmaical resistance
    real(dp), intent(in) :: bulksurface_resistance ! bulk surface resistance
    real(dp)             :: pet_penman             ! reference evapotranspiration in [mm s-1]

    pet_penman =  DaySecs / SpecHeatET_dp  *           & ! conversion factor [W m-2] to [mm d-1]
                  (slope_satpressure(tavg) * net_rad + &
                  rho0_dp * cp0_dp * (sat_vap_pressure(tavg) - act_vap_pressure ) / aerodyn_resistance) / &
                  (slope_satpressure(tavg) + Psychro_dp * (1.0_dp + bulksurface_resistance/aerodyn_resistance))
    
  END FUNCTION pet_penman

  ! ------------------------------------------------------------------

  !     NAME
  !         extraterr_rad_approx

  !>        \brief Approximation of extraterrestrial radiation 

  !>        \details Approximation of extraterrestrial radiation at the top of the atmosphere \f$ R_a \f$
  !>                 after Duffie and Beckman (1980). 
  !>                 \f$ R_a \f$ is converted from \f$ [J\;m^{-2}\;d^{-1}] \f$ in \f$ [mm\;d^{-1}]\f$ . 
  !>                 \f[ R_a   = \frac{86400}{ \pi \cdot \lambda} \cdot E_0 \cdot  
  !>                  d_r \cdot (\omega \cdot \sin(latitude) \cdot \sin(\delta) + \cos(latitude) \cdot \cos(\delta) \cdot
  !>                 \sin(\omega) \f] 
  !>                 where \f$ E_0=1367\;J\;m^{-2}\;s^{-1} \f$ is the solar constant and 
  !<                 \f$ \lambda = 2.45 \cdot 10^6\;J\;m^{-2}\;mm^{-1} \f$ is the latent heat of vaporization. \n
  !>                 It is dependent on the following sub equations:\n 
  !>                 The relative distance Earth-Sun:
  !>                 \f[ d_r =  1 + 0.033 \cdot \cos \left( \frac{2 \cdot \pi \cdot doy}{365} \right) \f]            
  !>                 in which doy is the day of the year.\n 
  !>                 The solar declination [radians] defined by
  !>                 \f[ \delta  = 0.4093 \cdot \sin\left( \frac{2 \cdot \pi \cdot doy}{365} - 1.405 \right) \f]
  !>                 The sunset hour angle [radians]:
  !>                 \f[ \omega  = \arccos( - \tan(latitude) * \tan(\delta) )  \f]                 

  !     INTENT(IN)
  !>        \param[in] "integer(i4), intent(in) :: doy" day of year [-]
  !>        \param[in] "real(dp),    intent(in) :: latitude" latitude [rad]

  !     INTENT(INOUT) 
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>         \return real(dp) :: extraterr_rad_approx &mdash; extraterrestrial radiation approximation \f$[W\;m^{-2}]\f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>        \note Duffie, J.A. and W.A. Beckman. 1980. Solar engineering of thermal processes.
  !>            John Wiley and Sons, New York. pp. 1-109.

  !     HISTORY
  !>        \author   Matthias Zink
  !>        \date     Apr 2014
  !         Modified  R. Kumar and M. Zink,    June 2016 - correction to include NaN in acos(arg) 

  elemental pure FUNCTION extraterr_rad_approx(doy, latitude) 

    use mo_constants,     only: SolarConst_dp, SpecHeatET_dp, PI_D, TWOPI_D
    use mo_mhm_constants, only: DuffieDr, DuffieDelta1, DuffieDelta2, YearDays, DaySecs

    implicit none    

    integer(i4), intent(in)             :: doy
    real(dp),    intent(in)             :: latitude             ! latitude [rad]
    real(dp)                            :: extraterr_rad_approx ! extraterrestrial radiation

    ! local
    real(dp)                            :: dr, delta
    real(dp)                            :: omega
    real(dp)                            :: arg
    
    ! inverse relative distance Earth-Sun - correction for eccentricity of Earths orbit around the sun
    dr     =  1.0_dp + DuffieDr * cos( TWOPI_D * doy / YearDays )            
    ! declination of the sun above the celestial equator in radians
    delta  =       DuffieDelta1 * sin( TWOPI_D * doy / YearDays - DuffieDelta2 ) 

    ! arccos(x) is only defined between PI and 0 (for x between -1 and 1)
    ! check limits
    arg = - tan(latitude) * tan(delta)
    if( arg .lt. -1.0_dp ) arg = -1.0_dp
    if( arg .gt.  1.0_dp ) arg =  1.0_dp

    ! sunrise hour angle in radians
    omega  = acos( arg )
       
    ! Ra - converted from [J m-2 d-1] in [mm d-1]
    extraterr_rad_approx   = DaySecs / PI_D / SpecHeatET_dp * SolarConst_dp *  &
         dr * (omega * sin(latitude) * sin(delta) + cos(latitude) * cos(delta) * sin(omega))

  end FUNCTION extraterr_rad_approx


  ! ------------------------------------------------------------------

  !     NAME
  !         slope_satpressure

  !>        \brief slope of saturation vapour pressure curve

  !>        \details slope of saturation vapour pressure curve after Tetens
  !>                 \f[ \Delta = \frac{0.6108 * e_s(T_a)}{e^(2 \cdot \log(T_a + 237.3))} \f]

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" average daily temperature \f$[^{\circ}C]\f$ 

  !     INTENT(INOUT) 
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>         \return real(dp) :: slope_satpressure &mdash;  slope of saturation vapour pressure curve 
  !>                             \f$[kPa\;K{-1}]\f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Tetens, O., 1930. Ueber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         \note Murray, F.W. 1967. On the computation of saturation vapor pressure. J. Appl. Meteor. 6: 203-204.
  !>         \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.  

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  ! 
  elemental pure FUNCTION slope_satpressure(tavg)

    use mo_mhm_constants, only: satpressureslope1, tetens_c3

    implicit none

    real(dp), intent(in) :: tavg                    ! average daily temperature
    real(dp)             :: slope_satpressure       ! slope of saturation vapour pressure curve


    slope_satpressure = satpressureslope1 * sat_vap_pressure(tavg) / exp(2.0_dp*log(Tavg + tetens_c3))
    
  END FUNCTION slope_satpressure

  ! ------------------------------------------------------------------

  !     NAME
  !         sat_vap_pressure

  !>        \brief calculation of the saturation vapour pressure

  !>        \details Calculation of the saturation vapour pressure
  !>                 \f[ e_s(T_a) = 0.6108 \cdot \exp \left( \frac{17.27 \cdot T_a}{T_a + 237.3} \right)  \f]

  !     INTENT(IN)
  !>        \param[in] "real(dp), intent(in) :: tavg" temperature [degC]

  !     INTENT(INOUT) 
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: sat_vap_pressure &mdash; saturation vapour pressure [kPa]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         \note Tetens, O., 1930. Ueber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         \note Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.  

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION sat_vap_pressure(tavg)

    use mo_mhm_constants, only: tetens_c1, tetens_c2, tetens_c3 

    implicit none

    real(dp), intent(in) :: tavg                      ! temperature [degC]
    real(dp)             :: sat_vap_pressure          ! saturation vapour pressure [kPa]

    sat_vap_pressure = tetens_c1 * exp(tetens_c2 * tavg / (tavg + tetens_c3))

  END FUNCTION sat_vap_pressure
  !
END MODULE mo_pet
