!> \file mo_refer_et.f90

!> \brief Module for calculating Reference ET [mm/s]

!> \details This module calculates ReferET [mm/s] based on one of the methods
!>  - Hargreaves-Samani (1982)
!>  - Priestly-Taylor   (1972)

!> \authors Christoph Schneider, Matthias Cuntz
!> \date Dec 2012

MODULE mo_pet

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind,      ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pet_hargreaves ! Hargreaves
  PUBLIC :: pet_priestly   ! Priestley-Taylor

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
  !>        \param[in] "real(dp) :: Ra"          top-of-the-atmosphere-radiation \f$ [W m^{-2}]\f$
  !>        \param[in] "real(dp) :: Tavg"        daily mean air temperature \f$ [ ^0C]\f$
  !>        \param[in] "real(dp) :: Tmax"        maximum daily temperature \f$ [ ^0C]\f$
  !>        \param[in] "real(dp) :: Tmin"        minimum daily temperature \f$ [ ^0C]\f$
  !>        \param[in] "real(dp) :: gamma"       parameter

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
  !>        \author Christoph Schneider
  !>        \date Dec 2012

  elemental pure FUNCTION pet_hargreaves(coefficient, constant, tavg, tmax, tmin, latitude, doy)

    use mo_constants,     only: deg2rad_dp

    implicit none

    real(dp),    intent(in) :: coefficient           !  coefficient of Hargreaves-Samani equation
    real(dp),    intent(in) :: constant              !  constatnt   of Hargreaves-Samani equation
    real(dp),    intent(in) :: tavg                  !  daily men temperature
    real(dp),    intent(in) :: tmax                  !  daily maximum of temp.
    real(dp),    intent(in) :: tmin                  !  daily minimum of temp.
    real(dp),    intent(in) :: latitude              ! latitude of the cell for Ra estimation
    integer(i4), intent(in) :: doy                   ! day of year for Ra estimation

    real(dp)                :: pet_hargreaves        ! reference evapotranspiration in [mm s-1]

    ! local
    real(dp)                :: Ra                    ! extraterrestrial radiation in [mm s-1]


    Ra             = extraterr_rad_approx(doy, deg2rad_dp * latitude)
    pet_hargreaves = coefficient * Ra * (Tavg + constant) * sqrt(Tmax - Tmin) ! in [mm s-1]

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
  !>        \param[in] "real(dp) :: Rn"          net solar radiation \f$ [W m^{-2}]\f$
  !>        \param[in] "real(dp) :: Tavg"        daily mean air temperature \f$ [ ^0C]\f$  
  !>        \param[in] "real(dp) :: gamma"       parameter

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
  !         None

  !     HISTORY
  !>        \authors Christoph Schneider
  !>        \date Dec 2012, Sep 2013
  !>    Worked on the documentation
  elemental pure FUNCTION pet_priestly(Rn, Tavg, gamma)

    use mo_mhm_constants, only: DeltaPriestly1, DeltaPriestly2
    use mo_constants,     only: Psychro_dp, T0_dp, SpecHeatET_dp !, PI_D 

    implicit none

    real(dp), intent(in) :: Rn
    real(dp), intent(in) :: Tavg
    real(dp), intent(in) :: gamma              ! parameter
    real(dp)             :: pet_priestly       ! reference evapotranspiration in [mm s-1]

    real(dp)             :: delta              ! slope of saturation-to-vapor-pressure

    delta = DeltaPriestly1 * exp(DeltaPriestly2 * (Tavg - T0_dp)) ! slope of saturation-to-vapor pressure curve

    pet_priestly = gamma * delta / (Psychro_dp + delta) * (1.0_dp/SpecHeatET_dp * Rn) ! in [mm s-1]

  END FUNCTION pet_priestly

  ! Ra as calculated by Duffie and Beckman (1980)
  elemental pure FUNCTION extraterr_rad_approx(doy, latitude) 

    use mo_constants,     only: SolarConst_dp, SpecHeatET_dp, PI_D, TWOPI_D
    use mo_mhm_constants, only: DaySecs

    implicit none    

    integer(i4), intent(in)             :: doy
    real(dp),    intent(in)             :: latitude             ! latitude [rad]
    real(dp)                            :: extraterr_rad_approx ! extraterrestrial radiation
    real(dp)                            :: dr, delta
    real(dp)                            :: omega
    
    ! correction for eccentricity of Earths orbit around the sun
    dr     =  1.0_dp + 0.0330_dp * cos( TWOPI_D * doy / 365.0_dp )            
    ! declination of the sun above the celestial equator in radians
    delta  =           0.4093_dp * sin( TWOPI_D * doy / 365.0_dp - 1.39_dp ) 
    ! sunrise hour angle in radians
    omega  = acos( - tan(latitude) * tan(delta) )                  
    
    ! Ra - converted from [J m-2 d-1] in [mm d-1]
    extraterr_rad_approx   = DaySecs / PI_D / SpecHeatET_dp * SolarConst_dp *  &
         dr * (omega * sin(latitude) * sin(delta) + cos(latitude) * cos(delta) * sin(omega))
  end FUNCTION extraterr_rad_approx
  !
END MODULE mo_pet
