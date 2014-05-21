!> \file mo_utils.f90

!> \brief General utilities for the CHS library

!> \details This module provides general utilities such as comparisons of two reals.

!> \authors Matthias Cuntz, Juliane Mai
!> \date Feb 2014
MODULE mo_utils

  ! Written  Matthias Cuntz, Juliane Mai, Feb 2014
  ! Modified Matthias Cuntz, Juliane Mai, Feb 2014 - equal, notequal
  ! Modified Matthias Cuntz,              May 2014 - swap

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

  USE mo_kind, ONLY: sp, dp, i4

  IMPLICIT NONE

  PUBLIC :: equal        ! a == b, a .eq. b
  PUBLIC :: notequal     ! a /= b, a .ne. b
  PUBLIC :: greaterequal ! a >= b, a .ge. b
  PUBLIC :: lesserequal  ! a <= b, a .le. b
  PUBLIC :: eq           ! a == b, a .eq. b
  PUBLIC :: ne           ! a /= b, a .ne. b
  PUBLIC :: ge           ! a >= b, a .ge. b
  PUBLIC :: le           ! a <= b, a .le. b
  PUBLIC :: locate       ! Find closest values in monotonic series
  PUBLIC :: swap         ! swaps arrays or elements of an array

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

  !     NAME
  !         swap

  !     PURPOSE
  !         Swap to values/arrays or two elements in 1D-array.
  !
  !>        \brief Swap to values or two elements in array.
  !
  !>        \details Swaps either two entities, i.e. scalars, vectors, matrices,
  !>                 or two elements in a vector.
  !>                 The call is either \n
  !>                   call swap(x,y) \n
  !>                 or \n
  !>                   call swap(vec,i,j)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: i"               Index of first element to be swapped with second [case swap(vec,i,j)]
  !>        \param[in] "integer(i4)    :: j"               Index of second element to be swapped with first [case swap(vec,i,j)]
  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(sp/dp/i4) :: x[(:,...)]"   First scalar or array to swap with second [case swap(x,y)]
  !>        \param[inout] "real(sp/dp/i4) :: y[(:[,:])]"   Second scalar or array to swap with first [case swap(x,y)]
  !>
  !>        \param[inout] "real(sp/dp/i4) :: x(:)"         Vector of which to elements are swapped [case swap(vec,i,j)]

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
  !     RESTRICTIONS
  !         No mask or undef.
  !
  !     EXAMPLE
  !         vec1 = (/ 1., 2., 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 1., 3., -999., 10., 6. /)
  !         call swap(vec1, vec2)
  !         call swap(vec1, 1, 3)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE swap
     MODULE PROCEDURE &
          swap_xy_dp, swap_xy_sp, swap_xy_i4, &
          swap_vec_dp, swap_vec_sp, swap_vec_i4
  END INTERFACE swap


  ! ------------------------------------------------------------------

  !     NAME
  !         locate

  !     PURPOSE
  !         Find closest values in monotonic series
  !
  !>         \brief Find closest values in monotonic series, returns the index.
  !
  !>         \details Given an array x(1:n), and given a value y,
  !>         returns a value j such that y is between
  !>         x(j) and x(j+1).\n
  !
  !>         x must be monotonically increasing.\n
  !>         j=0 or j=N is returned to indicate that x is out of range.
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp/sp)    :: x(:)"            Sorted array
  !>        \param[in] "real(dp/sp)    :: y[(:)]"          Value(s) of which the closest match in x(:) is wanted
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !>       \return     integer(i4) :: index[(:)] &mdash; index(es) of x so that y is between x(index) and x(index+1)
  !
  !     RESTRICTIONS
  !>       \note x must be monotonically increasing.\n
  !
  !     EXAMPLE
  !         x = (/ 1., 2., 3., -999., 5., 6. /)
  !         y = (/ 1.1, 5.6 /)
  !         ii = locate(x, y)
  !         -> ii == (/ 1, 5 /)
  !         y = 1.1
  !         ii = locate(x, y)
  !         -> ii == 1
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE locate
     MODULE PROCEDURE locate_0d_dp, locate_0d_sp, locate_1d_dp, locate_1d_sp
  END INTERFACE locate

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

  ! ------------------------------------------------------------------

  ! Given an array x(1:N), and given a value y, returns a value j such that y is between
  !  x(j) and x(j+1). x must be monotonically increasing.
  !  j=0 or j=N is returned to indicate that x is out of range.

  FUNCTION locate_0d_dp(x,y)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp),               INTENT(IN) :: y
    INTEGER(i4)                        :: locate_0d_dp

    INTEGER(i4), dimension(1) :: c

    c = minloc(abs(x-y))
    if (le(x(c(1)),y)) then
       locate_0d_dp = c(1)
    else
       locate_0d_dp = c(1)-1
    endif

  END FUNCTION locate_0d_dp

  FUNCTION locate_0d_sp(x,y)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    REAL(sp),               INTENT(IN) :: y
    INTEGER(i4)                        :: locate_0d_sp

    INTEGER(i4), dimension(1) :: c

    c = minloc(abs(x-y))
    if (le(x(c(1)),y)) then
       locate_0d_sp = c(1)
    else
       locate_0d_sp = c(1)-1
    endif

  END FUNCTION locate_0d_sp

  FUNCTION locate_1d_dp(x,y)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp), DIMENSION(:), INTENT(IN) :: y
    INTEGER(i4), DIMENSION(:), allocatable :: locate_1d_dp

    INTEGER(i4) :: ny, i
    INTEGER(i4), dimension(1) :: c


    ny = size(y)
    if (.not. allocated(locate_1d_dp)) allocate(locate_1d_dp(ny))

    do i=1, ny
       c = minloc(abs(x-y(i)))
       if (le(x(c(1)),y(i))) then
          locate_1d_dp(i) = c(1)
       else
          locate_1d_dp(i) = c(1)-1
       endif
    end do

  END FUNCTION locate_1d_dp

  FUNCTION locate_1d_sp(x,y)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    REAL(sp), DIMENSION(:), INTENT(IN) :: y
    INTEGER(i4), DIMENSION(:), allocatable :: locate_1d_sp

    INTEGER(i4) :: ny, i
    INTEGER(i4), dimension(1) :: c


    ny = size(y)
    if (.not. allocated(locate_1d_sp)) allocate(locate_1d_sp(ny))

    do i=1, ny
       c = minloc(abs(x-y(i)))
       if (le(x(c(1)),y(i))) then
          locate_1d_sp(i) = c(1)
       else
          locate_1d_sp(i) = c(1)-1
       endif
    end do

  END FUNCTION locate_1d_sp

  ! ------------------------------------------------------------------

  elemental pure subroutine swap_xy_dp(x,y)

    implicit none

    real(dp), intent(inout) :: x
    real(dp), intent(inout) :: y

    real(dp) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_dp

  elemental pure subroutine swap_xy_sp(x,y)

    implicit none

    real(sp), intent(inout) :: x
    real(sp), intent(inout) :: y

    real(sp) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_sp

  elemental pure subroutine swap_xy_i4(x,y)

    implicit none

    integer(i4), intent(inout) :: x
    integer(i4), intent(inout) :: y

    integer(i4) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_i4


  subroutine swap_vec_dp(x,i1,i2)

    implicit none

    real(dp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    real(dp) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_dp

  subroutine swap_vec_sp(x,i1,i2)

    implicit none

    real(sp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    real(sp) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_sp

  subroutine swap_vec_i4(x,i1,i2)

    implicit none

    integer(i4),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    integer(i4) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_i4

END MODULE mo_utils
