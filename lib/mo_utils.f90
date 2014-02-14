!> \file mo_utils.f90

!> \brief General utilities for the CHS library

!> \details This module provides general utilities such as comparisons of two reals.

!> \authors Matthias Cuntz, Juliane Mai
!> \date Feb 2014
MODULE mo_utils

  ! Written  Matthias Cuntz, Juliane Mai, Feb 2014
  ! Modified Matthias Cuntz, Juliane Mai, Feb 2014 - equal, notequal

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Matthias Cuntz, Juliane Mai

  USE mo_kind, ONLY: sp, dp

  IMPLICIT NONE

  PUBLIC :: equal        ! a == b, a .eq. b
  PUBLIC :: notequal     ! a /= b, a .ne. b
  PUBLIC :: greaterequal ! a >= b, a .ge. b
  PUBLIC :: lesserequal  ! a <= b, a .le. b
  PUBLIC :: eq ! a == b, a .eq. b
  PUBLIC :: ne ! a /= b, a .ne. b
  PUBLIC :: ge ! a >= b, a .ge. b
  PUBLIC :: le ! a <= b, a .le. b

  ! ------------------------------------------------------------------

  !     NAME
  !         equal / notequal / greaterequal / lesserequal

  !     PURPOSE
  !         Elemental function returning .true. or .false. depending if the reals are equal or not.
  !
  !>        \brief Comparison of real values.
  !
  !>        \details Compares two reals if they are numerically equal or not, i.e.
  !>        equal: \f[ |\frac{a-b}{b}| < \epsilon \f]
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: a"        First number to compare
  !>        \param[in] "real(sp/dp) :: b"        Second number to compare
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp) :: equal &mdash; \f$ a == b \f$ logically true or false
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         vec1 = (/ 1., 2., 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 1., 3., -999., 10., 6. /)
  !         isequal = equal(vec1, vec2)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \authors Matthias Cuntz, Juliane Mai
  !>        \date Feb 2014
  !         Modified, Matthias Cuntz, Juliane Mai, Feb 2014 - sp, dp
  INTERFACE equal
     MODULE PROCEDURE equal_sp, equal_dp
  END INTERFACE equal

  INTERFACE notequal
     MODULE PROCEDURE notequal_sp, notequal_dp
  END INTERFACE notequal

  INTERFACE greaterequal
     MODULE PROCEDURE greaterequal_sp, greaterequal_dp
  END INTERFACE greaterequal

  INTERFACE lesserequal
     MODULE PROCEDURE lesserequal_sp, lesserequal_dp
  END INTERFACE lesserequal

  INTERFACE eq
     MODULE PROCEDURE equal_sp, equal_dp
  END INTERFACE eq

  INTERFACE ne
     MODULE PROCEDURE notequal_sp, notequal_dp
  END INTERFACE ne

  INTERFACE ge
     MODULE PROCEDURE greaterequal_sp, greaterequal_dp
  END INTERFACE ge

  INTERFACE le
     MODULE PROCEDURE lesserequal_sp, lesserequal_dp
  END INTERFACE le

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION equal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: equal_dp

    if ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) then
       equal_dp = .false.
    else
       equal_dp = .true.
    endif

  END FUNCTION equal_dp


  ELEMENTAL PURE FUNCTION equal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: equal_sp

    if ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) then
       equal_sp = .false.
    else
       equal_sp = .true.
    endif

  END FUNCTION equal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION greaterequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: greaterequal_dp

    greaterequal_dp = .true.
    ! 1st part is /=, 2nd part is the a<b
    if ( ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) .and. (a < b) ) greaterequal_dp = .false.

  END FUNCTION greaterequal_dp


  ELEMENTAL PURE FUNCTION greaterequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: greaterequal_sp

    greaterequal_sp = .true.
    ! 1st part is /=, 2nd part is the a<b
    if ( ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) .and. (a < b) ) greaterequal_sp = .false.

  END FUNCTION greaterequal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION lesserequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: lesserequal_dp

    lesserequal_dp = .true.
    ! 1st part is /=, 2nd part is the a>b
    if ( ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) .and. (a > b) ) lesserequal_dp = .false.

  END FUNCTION lesserequal_dp


  ELEMENTAL PURE FUNCTION lesserequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: lesserequal_sp

    lesserequal_sp = .true.
    ! 1st part is /=, 2nd part is the a>b
    if ( ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) .and. (a > b) ) lesserequal_sp = .false.

  END FUNCTION lesserequal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION notequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: notequal_dp

    if ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) then
       notequal_dp = .true.
    else
       notequal_dp = .false.
    endif

  END FUNCTION notequal_dp


  ELEMENTAL PURE FUNCTION notequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: notequal_sp

    if ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) then
       notequal_sp = .true.
    else
       notequal_sp = .false.
    endif

  END FUNCTION notequal_sp

END MODULE mo_utils
