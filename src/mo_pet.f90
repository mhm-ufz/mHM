!> \file mo_refer_et.f90

!> \brief Module for calculating Reference ET [mm/s]

!> \details This module calculates ReferET [mm/s] based on one of the methods
!>  - Hargreaves-Samani (1982)
!>  - Priestly-Taylor   (1972)

!> \author Matthias Zink, Christoph Schneider, Matthias Cuntz
!> \date   Apr 2014

MODULE mo_pet

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind,      ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pet_hargreaves ! Hargreaves-Samani
  PUBLIC :: pet_priestly   ! Priestley-Taylor
  PUBLIC :: pet_penman     ! Penman Monteith

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         pet_hargreaves

  !>        \brief Reference Evapotranspiration after Hargreaves

  !>        \details Calculates the Reference Evapotranspiration [mm/d] based on the Hargreaves-Samani (1982) model for a given
  !>         cell by applying the equation
  !>         \f[ RET = \beta * R_a * (T_{avg} +  17.8) * \sqrt{ T_{max} - T_{min}} \f]
  !>         with \f$R_a\,[W m^{-2}]\f$ the extraterrestrial radiation (prescribed by a given yearly cycle resolved
  !>         at daily basis as given in \ref TopOfAtmosphereRadiation and
  !>         input data \f$T_{avg}, T_{max} \f$ and \f$ T_{min}\f$ the mean, maximum,
  !>         and minimum daily temperatures \f$ [ ^0C]\f$ at a given day and \f$\beta\f$ a conversion factor including site specifics.

  !     INTENT(IN)
  !>        XXX\param[in] "real(dp) :: Ra"          top-of-the-atmosphere-radiation \f$ [W m^{-2}]\f$
  !>        XXX\param[in] "real(dp) :: Tavg"        daily mean air temperature \f$ [ ^0C]\f$
  !>        XXX\param[in] "real(dp) :: Tmax"        maximum daily temperature \f$ [ ^0C]\f$
  !>        XXX\param[in] "real(dp) :: Tmin"        minimum daily temperature \f$ [ ^0C]\f$
  !>        XXXX\param[in] "real(dp) :: gamma"       parameter

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
  !>         \return real(dp) :: ReferET &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Hargreaves, G.H., and Samani, Z.A. (1982). "Estimating potential evapotranspiration."
  !             Tech. Note, J. Irrig. and drain. Engrg., ASCE, 108(3):225-230.

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Dec 2012

  elemental pure FUNCTION pet_hargreaves(HarSamCoeff, HarSamConst, tavg, tmax, tmin, latitude, doy)

    use mo_constants,     only: deg2rad_dp

    implicit none

    real(dp),    intent(in) :: HarSamCoeff           !  coefficient of Hargreaves-Samani equation
    real(dp),    intent(in) :: HarSamConst           !  constatnt   of Hargreaves-Samani equation
    real(dp),    intent(in) :: tavg                  !  daily men temperature
    real(dp),    intent(in) :: tmax                  !  daily maximum of temp.
    real(dp),    intent(in) :: tmin                  !  daily minimum of temp.
    real(dp),    intent(in) :: latitude              ! latitude of the cell for Ra estimation
    integer(i4), intent(in) :: doy                   ! day of year for Ra estimation

    real(dp)                :: pet_hargreaves        ! reference evapotranspiration in [mm s-1]

    real(dp)                :: delta_temp        ! tmax-Tmin

    ! correction for shity input data ! MZMZMZMZ
    if (tmax - tmin .GT. 0.0_dp) then 
      delta_temp     = tmax - tmin ! MZMZMZMZ
    else
      delta_temp     =  1.0e-05_dp  ! MZMZMZMZ
    end if
    ! in [mm s-1]   
  ! to avoid numerical errors
  if (tavg .lt. -17.8_dp) then
     pet_hargreaves = 1.0E-5_dp
  else
    pet_hargreaves = HarSamCoeff * extraterr_rad_approx(doy, deg2rad_dp * latitude) * (tavg + HarSamConst) * sqrt(delta_temp) 
  end if


  END FUNCTION pet_hargreaves


  ! ------------------------------------------------------------------

  !     NAME
  !         pet_priestly

  !>        \brief Reference Evapotranspiration after Priestly-Taylor

  !>        \details Calculates the Reference Evapotranspiration [mm/d] based on the Priestly-Taylor (1972) model
  !>        for every given cell by applying the equation
  !>        \f[ RET = \alpha * \frac{\Delta}{(\gamma + \Delta)} * R_n \f]
  !>        where \f$R_n\,[W m^{-2}]\f$ is the net solar radiation \f$\Delta =  f(T_{avg})\f$ is the slope 
  !>        of the saturation-vapour pressure curve, which is a function of the mean daily temperature 
  !>        \f$ [ ^0C]\f$ and \f$\alpha\f$ is a land cover dependent coefficient.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: Rn"           net solar radiation \f$ [W m^{-2}]\f$
  !>        \param[in] "real(dp) :: Tavg"         daily mean air temperature \f$ [ ^0C]\f$  
  !>        \param[in] "real(dp) :: PrieTaycoeff" Priestley-Taylor coefficient

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
  !>         \return real(dp) :: ReferET &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         Priestley, C.H.B., and R.J. Taylor. 1972. On the assessment of surface heat flux and evaporation using
  !>                   large-scale parameters. Mon. Weather Rev., 100:81-82.
  !>         ASAE Standards. 1998. EP406.2: heating, cooling, and ventilating greenhouses. St. Joseph, MI, USA.

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

!  elemental pure FUNCTION pet_priestly(PrieTayParam, Rn, Tavg)
  FUNCTION pet_priestly(PrieTayParam, Rn, Tavg)

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

  !>        \brief XXX

  !>        \details XXX

  !     INTENT(IN)
  !>        XXX\param[in] "real(dp) :: Rn"           net solar radiation \f$ [W m^{-2}]\f$
  !>        XXX\param[in] "real(dp) :: Tavg"         daily mean air temperature \f$ [ ^0C]\f$  
  !>        XXX\param[in] "real(dp) :: PrieTaycoeff" Priestley-Taylor coefficient

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
  !>         \return real(dp) :: ReferET &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
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

    pet_penman =  DaySecs / SpecHeatET_dp  *           &                          ! conversion factor [W m-2] to [mm d-1]
                  (slope_satpressure(tavg) * net_rad + &
                  rho0_dp * cp0_dp * (sat_vap_pressure(tavg) - act_vap_pressure ) / aerodyn_resistance) / &
                  (slope_satpressure(tavg) + Psychro_dp * (1 + bulksurface_resistance/aerodyn_resistance))
    
  END FUNCTION pet_penman





  ! Ra as calculated by Duffie and Beckman (1980)
  elemental pure FUNCTION extraterr_rad_approx(doy, latitude) 


  !     LITERATURE
  !>        Duffie, J.A. and W.A. Beckman. 1980. Solar engineering of thermal processes.
  !>            John Wiley and Sons, New York. pp. 1-109.

    use mo_constants,     only: SolarConst_dp, SpecHeatET_dp, PI_D, TWOPI_D
    use mo_mhm_constants, only: DuffieDr, DuffieDelta1, DuffieDelta2, YearDays, DaySecs

    implicit none    

    integer(i4), intent(in)             :: doy
    real(dp),    intent(in)             :: latitude             ! latitude [rad]
    real(dp)                            :: extraterr_rad_approx ! extraterrestrial radiation
    real(dp)                            :: dr, delta
    real(dp)                            :: omega
    
    ! inverse relative distance Earth-Sun - correction for eccentricity of Earths orbit around the sun
    dr     =  1.0_dp + DuffieDr * cos( TWOPI_D * doy / YearDays )            
    ! declination of the sun above the celestial equator in radians
    delta  =       DuffieDelta1 * sin( TWOPI_D * doy / YearDays - DuffieDelta2 ) 
    ! sunrise hour angle in radians
    omega  = acos( - tan(latitude) * tan(delta) )                  
    
    ! Ra - converted from [J m-2 d-1] in [mm d-1]
    extraterr_rad_approx   = DaySecs / PI_D / SpecHeatET_dp * SolarConst_dp *  &
         dr * (omega * sin(latitude) * sin(delta) + cos(latitude) * cos(delta) * sin(omega))

  end FUNCTION extraterr_rad_approx


  ! ------------------------------------------------------------------

  !     NAME
  !         slope_satpressure

  !>        \brief slope of saturation vapour pressure curve

  !>        \details XXX

  !     INTENT(IN)
  !>        XXX
  !>        XXX
  !>        XXX

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
  !>         XXX

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         Tetens, O., 1930. Uber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         Murray, F.W. 1967. On the computation of saturation vapor pressure. J. Appl. Meteor. 6: 203-204.
  !>         Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.  

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  ! 
  elemental pure FUNCTION slope_satpressure(Tavg)

    use mo_mhm_constants, only: satpressureslope1, tetens_c3

    implicit none

    real(dp), intent(in) :: tavg
    real(dp)             :: slope_satpressure       ! slope of saturation vapour pressure curve


    slope_satpressure = satpressureslope1 * sat_vap_pressure(tavg) / exp(2*log(Tavg + tetens_c3))
    
  END FUNCTION slope_satpressure

  ! ------------------------------------------------------------------

  !     NAME
  !         sat_vap_pressure

  !>        \brief XXX

  !>        \details XXX

  !     INTENT(IN)
  !>        XXX\param[in] "real(dp) :: Rn"           net solar radiation \f$ [W m^{-2}]\f$
  !>        XXX\param[in] "real(dp) :: Tavg"         daily mean air temperature \f$ [ ^0C]\f$  
  !>        XXX\param[in] "real(dp) :: PrieTaycoeff" Priestley-Taylor coefficient

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
  !>         XXX\return real(dp) :: ReferET &mdash; Reference Evapotranspiration [mm s-1]

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>         Tetens, O., 1930. Uber einige meteorologische Begriffe. z. Geophys. 6:297-309.
  !>         Allen, R. G. R., Pereira, L., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for 
  !>             computing crop water requirements - FAO Irrigation and drainage paper 56. Rome.  

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    Apr 2014

  elemental pure FUNCTION sat_vap_pressure(temp)

    use mo_mhm_constants, only:tetens_c1, tetens_c2, tetens_c3 

    implicit none

    real(dp), intent(in) :: temp                      ! temperature [degC]
    real(dp)             :: sat_vap_pressure          ! saturation vapour pressure [kPa]

    sat_vap_pressure = tetens_c1 * exp(tetens_c2 * temp / (temp + tetens_c3))

  END FUNCTION sat_vap_pressure
  !
END MODULE mo_pet
