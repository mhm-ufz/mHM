!> \file mo_canopy_interc.f90
!> \brief \copybrief mo_canopy_interc
!> \details \copydetails mo_canopy_interc

!> \brief Canopy interception.
!> \details This module deals with processes related to canopy interception, evaporation and throughfall.
!> \changelog
!! - RK Sep 2013
!!   - Documentation updated (formula and a short description added)
!> \authors Vladyslav Prykhodko
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_canopy_interc

  USE mo_kind, ONLY : dp
  USE mo_common_constants, ONLY : eps_dp
  USE mo_constants, ONLY : twothird_dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: canopy_interc ! Canopy interception

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        canopy_interc

  !    PURPOSE
  !>       \brief Canopy interception.

  !>       \details Calculates throughfall.
  !>       Updates interception and evaporation intensity from canopy.
  !>       Throughfall (\f$F\f$) is estimated as a function of the incoming precipitation (\f$P\f$),
  !>       the current status of the canopy water content (\f$C\f$), and the max. water
  !>       \f[ F = Max( (P + C - C_{max}), 0) \f]
  !>       Evaporation (\f$E\f$) from canopy is estimated as a fraction of the potential
  !>       evapotranspiration(\f$E_{p}\f$) depending on the current status of the canopy
  !>       water content (\f$C\f$) and the max. water content(\f$C_{max}\f$) that can be
  !>       intecepted by the vegetation.
  !>       \f[ E = E_{p}(C/C_{max})^{2/3} \f]
  !>       ADDITIONAL INFORMATION
  !>       content(\f$C_{max}\f$) that can be intecepted by the vegetation.
  !>       canopy_interc(pet, interc_month_max, interc_max, precip, throughfall, evap_canopy, interc)

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: pet"        Potential evapotranspiration [mm TS-1]
  !>       \param[in] "REAL(dp) :: interc_max" Maximum interception [mm]
  !>       \param[in] "REAL(dp) :: precip"     Daily mean precipitation [mm]

  !    INTENT(INOUT)
  !>       \param[inout] "REAL(dp) :: interc" Interception [mm]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: throughfall" Throughfall [mm TS-1]
  !>       \param[out] "REAL(dp) :: evap_canopy" Real evaporation intensity from canopy[mm TS-1]

  !    HISTORY
  !>       \authors Vladyslav Prykhodko

  !>       \date Dec 2012

  ! Modifications:
  ! JM Aug 2013 - ordering of arguments changed
  ! RK Sep 2013 - Documentation updated (formula and a short description added)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  ELEMENTAL PURE SUBROUTINE canopy_interc(pet, interc_max, precip, interc, throughfall, evap_canopy)
    implicit none

    ! Potential evapotranspiration [mm TS-1]
    REAL(dp), INTENT(IN) :: pet

    ! Maximum interception [mm]
    REAL(dp), INTENT(IN) :: interc_max

    ! Daily mean precipitation [mm]
    REAL(dp), INTENT(IN) :: precip

    ! Interception [mm]
    REAL(dp), INTENT(INOUT) :: interc

    ! Throughfall [mm TS-1]
    REAL(dp), INTENT(OUT) :: throughfall

    ! Real evaporation intensity from canopy[mm TS-1]
    REAL(dp), INTENT(OUT) :: evap_canopy

    ! Auxiliary helping variable [-]
    REAL(dp) :: aux_help


    !===============================================
    ! Canopy Interception
    ! Canopy storage (actualize)
    ! 1st rains -> 2nd Interception -> 3rd ETP
    !===============================================
    aux_help = interc + precip
    if (aux_help >= interc_max) then
      throughfall = aux_help - interc_max
      interc = interc_max
    else
      throughfall = 0.0_dp
      interc = aux_help
    end if

    ! New module for evaporation from canopy surface
    ! [power (2/3) is based on the paper of Liang et al. 1994 & Deardorf, 1978]
    if (interc_max > eps_dp) then
      evap_canopy = pet * (interc / interc_max)**twothird_dp
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
      interc = 0.0_dp
    end if

  END SUBROUTINE canopy_interc

END MODULE mo_canopy_interc
