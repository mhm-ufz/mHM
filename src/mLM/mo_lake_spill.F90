!> \file mo_lake_spill.F90
!> \copydoc mo_lake_spill

!> \brief lake spill module

!> \details  lake spill module calculates the water spill out of the lake over the dam
!>           crest, depending upon the spillway characteristics and the lake level.

!> \authors Pallav K Shrestha, Luis Samaniego

!> \date October 2018

!> \version 0.1
!> \changelog
!! - October 2026: v6 component ownership and restart integration; retained donor arithmetic.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
MODULE mo_lake_spill

  USE mo_kind,   only : i4, i8, dp
  IMPLICIT NONE
  private



  PUBLIC :: lake_spill        ! lake spill over the dam crest
  PUBLIC :: lake_spill_generalized_equation
                              ! generalized weir equation for the spill over the crest
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         lake_spill

  !     PURPOSE
  !>        \brief lake spill

  !>        \details calculates the spill based on the current lake level. Spill is
  !>                 produced only when the lake level rises above the dam crest level;
  !>                 otherwise it is zero.


  !     CALLING SEQUENCE

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Feb 2021

  !     MODIFICATIONS
  !        Pallav Shrestha, Oct 2026 - mLM v6 migration: ported from mo_mlm_spill into the
  !                                    v6 tree. The equation was already in SI units and is
  !                                    unchanged (level [m], length [m], spill [m3 s-1]).
  !                                    Consistent nomenclature (DCL -> lake_spill_reference_level,
  !                                    h -> lake_level, C_d -> spill_coefficient_of_discharge,
  !                                    L_s -> spillway_length). Dropped the unused
  !                                    mo_mlm_constants / mo_constants imports (no globals).

  recursive SUBROUTINE lake_spill(lake_level, lake_spill_reference_level, spill_coefficient_of_discharge, &
                        spillway_length, spill_exponent, spill)

    IMPLICIT NONE


    real(dp), intent(IN) :: lake_level !< [masl]  current lake level

    real(dp), intent(in) :: lake_spill_reference_level !< reference level [m]
    real(dp), intent(in) :: spill_coefficient_of_discharge !< effective SI discharge coefficient
    real(dp), intent(in) :: spillway_length !< crest length [m]
    real(dp), intent(in) :: spill_exponent !< head exponent

    real(dp), intent(OUT) :: spill !< [m3 s-1] spill flow rate



    !! local vars

    ! INITIALIZE
    ! -------------------
    spill           = 0.0_dp

    ! SPILL or NO SPILL?
    ! -----------------------------
    spillChk: if ( lake_level > lake_spill_reference_level ) then

      ! Lake spills if the water level is above the dam crest level
      call lake_spill_generalized_equation(lake_level, & ! Intent IN
        lake_spill_reference_level, spill_coefficient_of_discharge, spillway_length, spill_exponent, & ! Intent IN
        spill) ! Intent OUT

    end if spillChk


  END SUBROUTINE lake_spill




  ! ------------------------------------------------------------------

  !     NAME
  !         lake_spill_generalized_equation

  !     PURPOSE
  !>        \brief Calculates the lake spill / natural lake outflow over the crest

  !>        \details lake spill module calculates the water spill out of the lake over the
  !>                 dam crest depending upon the lake level and the spill parameters.

  !     CALLING SEQUENCE

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         USBR (1998) Design of Small Dams - A water resources technical publication

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Feb 2021

  !     MODIFICATIONS
  !        Pallav Shrestha, Oct 2026 - mLM v6 migration: consistent nomenclature
  !                                    (h0 -> reference_level, C_d ->
  !                                    spill_coefficient_of_discharge, L_s -> spillway_length)
  !                                    and dropped the unused mo_mlm_constants / mo_constants
  !                                    imports. Equation and SI units unchanged.


  ! (ver. generalized equation)
  recursive SUBROUTINE lake_spill_generalized_equation(lake_level, reference_level, &
                                             spill_coefficient_of_discharge, spillway_length, &
                                             spill_exponent, spill)

    USE mo_message,                       ONLY: message


    IMPLICIT NONE

    real(dp), intent(IN) :: lake_level !< [masl]  current lake level

    real(dp), intent(IN) :: reference_level !< [masl]  reference level (top of the conservation rule curve)

    real(dp), intent(IN) :: spill_coefficient_of_discharge !< [-]  spill coefficient model parameter
    real(dp), intent(IN) :: spillway_length !< [m]  effective spillway length
    real(dp), intent(IN) :: spill_exponent !< [-]  spill exponent model parameter

    real(dp), intent(OUT) :: spill !< [m3 s-1] spill flow rate


    !! local vars
    real(dp)              ::  deltaH           ! [m]      pressure head above reference_level

    !---------------------------------------------------------------------
    !
    ! **************** Flow over spillways ******************************
    !
    !
    !       For ogee shaped free/ uncontrolled spillways (USBR 1998, pp. 365):
    !
    !
    !             Q = Cd . L . H1^(3/2)
    !
    !
    !
    !       For ogee shaped gated/ controlled spillways (USBR 1998, pp. 376):
    !
    !             Q = Cd . D. L . sqrt(2g) . H2^(1/2)
    !
    !
    !
    !                             where, Cd : coefficient of discharge
    !                                     g : accelation due to gravity
    !                                     L : effective length of the spillways across
    !                                         the direction of flow
    !                                     D : gate opening
    !                                     H1: head above the spillway crest
    !                                     H2: head above center of D
    !
    !
    !       Generalized spill equation:
    !
    !             Q = alpha . deltaH ^ beta
    !
    !                        where,  deltaH : head above the reference_level
    !                                 alpha : spill coefficient parameter
    !                                  beta : spill exponent parameter
    !
    !---------------------------------------------------------------------

    spill = 0.0_dp

    ! head above the reference level
    deltaH = lake_level - reference_level

    ! check and communicate with user
    if (deltaH < 0) then
      call message('ERROR: Attempting to spill while water level is below reference level - lake_spill subroutine')
      stop 1
    end if

    ! discharge calculation
    spill = spill_coefficient_of_discharge * spillway_length * deltaH ** spill_exponent


  END SUBROUTINE lake_spill_generalized_equation




END MODULE mo_lake_spill
