!> \file mo_refer_et.f90

!> \brief Module for calculating Reference ET [mm/s]

!> \details This module calculates ReferET [mm/s] based on one of the methods
!>  - Hargreaves-Samani (1982)
!>  - Priestly-Taylor (1972)

!> \authors Christoph Schneider, Matthias Cuntz
!> \date Dec 2012

MODULE mo_refer_et

  ! This module is for the UFZ CHS mesoscale hydrologic model mHM.

  USE mo_kind,      ONLY: dp
  USE mo_constants, ONLY: Psychro_dp, T0_dp, SpecHeatET_dp 
  USE mo_mhm_constants, ONLY: HargreavesConst, DeltaPriestly1, DeltaPriestly2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: refer_et_hargreaves ! Hargreaves
  PUBLIC :: refer_et_priestly   ! Priestley-Taylor

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         refer_et_hargreaves

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

  elemental pure FUNCTION refer_et_hargreaves(Ra, Tavg, Tmax, Tmin, gamma)

    real(dp), intent(in) :: Ra
    real(dp), intent(in) :: Tavg
    real(dp), intent(in) :: Tmax
    real(dp), intent(in) :: Tmin
    real(dp), intent(in) :: gamma                 ! parameter
    real(dp)             :: refer_et_hargreaves   ! reference evapotranspiration in [mm s-1]

    refer_et_hargreaves = gamma * (1.0_dp/SpecHeatET_dp * Ra) * (Tavg + HargreavesConst) * sqrt(Tmax - Tmin) ! in [mm s-1]

  END FUNCTION refer_et_hargreaves


  ! ------------------------------------------------------------------

  !     NAME
  !         refer_et_priestly

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
  elemental pure FUNCTION refer_et_priestly(Rn, Tavg, gamma)

    real(dp), intent(in) :: Rn
    real(dp), intent(in) :: Tavg
    real(dp), intent(in) :: gamma              ! parameter
    real(dp)             :: refer_et_priestly  ! reference evapotranspiration in [mm s-1]

    real(dp)             :: delta              ! slope of saturation-to-vapor-pressure

    delta = DeltaPriestly1 * exp(DeltaPriestly2 * (Tavg - T0_dp)) ! slope of saturation-to-vapor pressure curve

    refer_et_priestly = gamma * delta / (Psychro_dp + delta) * (1.0_dp/SpecHeatET_dp * Rn) ! in [mm s-1]

  END FUNCTION refer_et_priestly

END MODULE mo_refer_et
