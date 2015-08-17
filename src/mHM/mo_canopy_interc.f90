!> \file mo_canopy_interc.f90

!> \brief Canopy interception.

!> \details This module deals with processes related to canopy interception, evaporation and throughfall.

!> \authors Vladyslav Prykhodko
!> \date Dec 2012
!   Modified RK, Sep 2013 - Documentation updated (formula and a short description added) 

MODULE mo_canopy_interc

  USE mo_kind,      ONLY: dp
  USE mo_constants, ONLY: twothird_dp, eps_dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: canopy_interc ! Canopy interception

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         canopy_interc

  !     PURPOSE
  !>        \brief Canopy interception.

  !>        \details  Calculates throughfall.
  !>                  Updates interception and evaporation intensity from canopy.\n
  
  !>         Throughfall (\f$F\f$) is estimated as a function of the incoming precipitation (\f$P\f$), 
  !>         the current status of the canopy water content (\f$C\f$), and the max. water 
  !          content(\f$C_{max}\f$) that can be intecepted by the vegetation.\n       
  !>         \f[ F = Max( (P + C - C_{max}), 0) \f]
  
  !>         Evaporation (\f$E\f$) from canopy is estimated as a fraction of the potential 
  !>         evapotranspiration(\f$E_{p}\f$) depending on the current status of the canopy  
  !>         water content (\f$C\f$) and the max. water content(\f$C_{max}\f$) that can be 
  !>         intecepted by the vegetation. \n
  !>         \f[ E = E_{p}(C/C_{max})^{2/3} \f]
  
  
  !     CALLING SEQUENCE
  !         canopy_interc(pet, interc_month_max, interc_max, precip, throughfall, evap_canopy, interc)

  !     INDENT(IN)
  !>        \param[in] "real(dp) ::  precip"               Daily mean precipitation [mm]
  !>        \param[in] "real(dp) ::  pet"                  Potential evapotranspiration [mm s-1]
  !>        \param[in] "real(dp) ::  interc_max"           Maximum interception [mm]

  !     INDENT(INOUT)
  !>        \param[in,out] "real(dp)  ::  interc"          Interception [mm]
  
  !     INDENT(OUT)
  !>        \param[out] "real(dp)  ::  evap_canopy"        Real evaporation intensity from canopy[mm s-1]
  !>        \param[out] "real(dp)  ::  throughfall"        Throughfall [mm s-1]


  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE


  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Vladyslav Prykhodko
  !>        \date Dec 2012
  !         Modified JM, Aug 2013 - ordering of arguments changed
  !                  RK, Sep 2013 - Documentation updated (formula and a short description added) 

  ELEMENTAL PURE SUBROUTINE canopy_interc(pet, interc_max, precip, interc, throughfall, evap_canopy)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)     :: pet
    REAL(dp), INTENT(IN)     :: interc_max
    REAL(dp), INTENT(IN)     :: precip
    REAL(dp), INTENT(INOUT)  :: interc
    REAL(dp), INTENT(OUT)    :: throughfall
    REAL(dp), INTENT(OUT)    :: evap_canopy
    
    ! local variables
    REAL(dp)                 :: aux_help ! Auxiliary helping variable [-]

    !===============================================
    ! Canopy Interception
    ! Canopy storage (actualize)
    ! 1st rains -> 2nd Interception -> 3rd ETP
    !===============================================

    aux_help = interc + precip
    if (aux_help >= interc_max) then
       throughfall = aux_help - interc_max
       interc      = interc_max
    else
       throughfall = 0.0_dp
       interc      = aux_help
    end if

    ! New module for evaporation from canopy surface
    ! [power (2/3) is based on the paper of Liang et al. 1994 & Deardorf, 1978]
    if (interc_max > eps_dp) then
       evap_canopy = pet * (interc/interc_max)**twothird_dp       
    else
       ! in case interc_max is 
       evap_canopy = 0.0_dp
    end if

    ! numerical problem
    if (evap_canopy < 0.0_dp) evap_canopy = 0.0_dp ! this should never appear

    if (interc > evap_canopy) then
       interc = interc - evap_canopy
    else
       evap_canopy = interc
       interc      = 0.0_dp
    end if

  END SUBROUTINE canopy_interc

END MODULE mo_canopy_interc
