!> \file mo_radiation.f90

!> \brief Extraterrestrial Radiation

!> \details This module provides the calculation of Extraterrestrial radiation at a given day of year.\n
!> It also prvodes a routine to calculate net radiation from the individual components.
!> But this will be removed in a later version.

!> \authors Christoph Schneider, Matthias Cuntz
!> \date Dec 2012

MODULE mo_radiation

  USE mo_kind,      ONLY: i4, dp
  
  PRIVATE

  PUBLIC :: SolarNetRadiation
  PUBLIC :: ExtraterrestrialRadiation

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         ExtraterrestrialRadiation

  !>        \brief Top-of-the-atmosphere radiation, i.e. extraterrestrial radiation

  !>        \details Calculates the Extraterrestrial radiation [W m-2] based on Duffie & Beckman (1980)
  !>        \f[ R_{0} = c_{solar} / \pi * d_r * \left( \omega \sin(\delta) \sin(\phi) + 
  !>                  \cos(\delta) \cos(\phi) \sin(\omega) \right) \f]
  !>        with \f$ c_{solar} \f$ the solar constant (1367 W \f$m^{-2}\f$), \f$\phi\f$ the latitude, 
  !>        \f$\delta\f$ the declination of the sun above the celestial equator,
  !>        \f[ \delta = 0.409 \sin(2\pi doy/365 - 1.39) \f]
  !>        \f$doy\f$ the day of the year, \f$\omega\f$ the angle at sunrise
  !>        \f[ \omega = \arccos(-\tan(\phi)\tan(\delta)) \f]
  !>        and \f$d_r\f$ a correction for the eccentricity of Earth's orbit
  !>        \f[ d_r = 1 + 0.033 \cos(2\pi doy/365) \f]

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: md"         day-of-the-year
  !>        \param[in] "real(dp)    :: phiLat"     latitude in radians
 

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
  !>        \return real(dp) :: ExtraterrestrialRadiation &mdash; Extraterrestrial radiation \f$[W m^{-2}]\f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Christoph Schneider
  !>        \date Dec 2012

  ! Extraterrestrial-radiation calculated by Duffie and Beckman (1980)
  elemental pure FUNCTION ExtraterrestrialRadiation(md, phiLat)
    USE mo_constants,     ONLY: SolarConst_dp, pi_dp
    USE mo_mhm_constants, ONLY: DuffieDr, DuffieDelta1, DuffieDelta2, YearDays
    
    IMPLICIT NONE

    integer(i4), intent(in) :: md     ! day of year
    real(dp),    intent(in) :: phiLat ! latitude in radians
    real(dp)                :: ExtraterrestrialRadiation

    ! Local variables
    real(dp)                :: dr, delta
    real(dp)                :: omega

    ! correction for eccentricity of Earth''s orbit around sun
    dr    = 1.0_dp + DuffieDr * cos ( 2.0_dp * pi_dp * md / YearDays )  
    ! declination of the sun above celestial equator in radians             
    delta = DuffieDelta1 * sin( 2.0_dp * pi_dp * md / YearDays - DuffieDelta2 ) 
    ! sunrise hour angle in radians    
    omega  = acos( - tan(phiLat) * tan(delta) )                
    ! in [W m-2]                     
    ExtraterrestrialRadiation = SolarConst_dp / pi_dp * dr * ( omega * sin(delta) * sin(phiLat) + &
         cos(phiLat) * cos(delta) * sin(omega) )                       

  END FUNCTION ExtraterrestrialRadiation

  ! ------------------------------------------------------------------

  !     NAME
  !         SolarNetRadiation(Albedo, Emmissivity, LW_in, SW_in, Tavg)

  !>        \brief Net radiation of the sun

  !>        \details Calculates the SolarNetRadiation [W m-2] based on long and shortwave components.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: Albedo"       Albedo of the land surface [-]
  !>        \param[in] "real(dp) :: Emmissivity"  Emmissivity of the land surface [-]
  !>        \param[in] "real(dp) :: LW_in"        Ingoing long wave radiation \f$[W m^{-2}]\f$
  !>        \param[in] "real(dp) :: LW_in"        Ingoing short wave radiation \f$[W m^{-2}]\f$  
  !>        \param[in] "real(dp) :: T_avg"        daily mean air temperature [\f$^0C]\f$

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
  !>         \return real(dp) :: SolarNetRadiation &mdash; Solar net radiation \f$[W m^{-2}]\f$

  !     RESTRICTIONS
  !>       \note Outgoing short wave radiation SW_out is given by Albedo*SW_in, outgoing long wave radiation LW_out
  !>             by StBoltzmann * Emmissivity * Tavg^4.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Christoph Schneider
  !>        \date Dec 2012
  
  elemental pure FUNCTION SolarNetRadiation(Albedo, Emmissivity, LW_in, SW_in, Tavg)
    USE mo_mhm_constants,     ONLY: StBoltzmann
    IMPLICIT NONE

    real(dp), intent(in) :: Albedo, Emmissivity, SW_in, LW_in, Tavg
    real(dp)             :: SolarNetRadiation

    !          Rn     =     SW_in       - SW_out  +  LW_in -  LW_out
    SolarNetRadiation = SW_in * (1.0_dp - Albedo) +  LW_in -  StBoltzmann * Emmissivity * Tavg**4 ! in [W m-2]

  END FUNCTION SolarNetRadiation

END MODULE mo_radiation
