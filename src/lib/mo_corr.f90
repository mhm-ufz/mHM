MODULE mo_corr

  ! This module provides autocorrelation function calculations

  ! Rewritten December 2019, Sebastian Mueller

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
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2019 Sebastian Mueller


  USE mo_kind, ONLY : i4, sp, dp

  Implicit NONE

  PUBLIC :: autocorr     ! Autocorrelation coefficient at lag k = autocoeffk(k)/autocoeffk(0)

  ! ------------------------------------------------------------------

  !     NAME
  !         autocorr

  !     PURPOSE
  !         Element at lag k of autocorrelation function
  !             autocorr(x,k) = autocoeffk(x,k)/autocoeffk(x,0).
  !
  !         If an optinal mask is given, the calculations are only over those locations that correspond
  !         to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         ak = autocorr(x, k, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: x(:)        Time series
  !         integer(i4) :: k[(:)]      Lag for autocorrelation

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: ak[(:)]     Coefficient of autocorrelation function at lag k

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(vec).
  !                                    If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Last autocorrelation element
  !         acorr = autocorr(x,size(x)/2)
  !         -> see also example in test directory

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  !         Modified, Stephan Thober, Nov 2012 - added 1d version
  !         Modified, Sebastian Mueller, Dec 2019 - rewritten
  INTERFACE autocorr
    MODULE PROCEDURE autocorr_sp, autocorr_dp, autocorr_1d_sp, autocorr_1d_dp
  END INTERFACE autocorr

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION autocorr_dp(x, k, mask) result(acf)

    USE mo_moment, ONLY: moment

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), INTENT(IN) :: k
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp) :: acf

    INTEGER(i4) :: kk  ! absolute value of lag k
    INTEGER(i4) :: nn  ! number of true values in mask
    REAL(dp) :: n   ! real of nn
    REAL(dp) :: mean, var  ! moments
    REAL(dp), DIMENSION(size(x)) :: x_shift
    LOGICAL, DIMENSION(size(x)) :: maske

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error autocorr_dp: size(mask) /= size(x)'
      maske = mask
    end if

    ! calculate mean and population variance of x
    call moment(x, mean=mean, variance=var, mask=maske, sample=.false.)
    ! shift x to 0 mean
    x_shift = x - mean
    ! use absolute value of k
    kk = abs(k)
    nn = size(x)
    n = real(count(maske(1 : nn - kk).and.maske(1 + kk : nn)), dp)
    ! calculate the auto-correlation function
    acf = sum(x_shift(1 : nn - kk) * x_shift(1 + kk : nn), mask = (maske(1 : nn - kk) .and. maske(1 + kk : nn))) / n / var

  END FUNCTION autocorr_dp

  FUNCTION autocorr_sp(x, k, mask) result(acf)

    USE mo_moment, ONLY: moment

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), INTENT(IN) :: k
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp) :: acf

    INTEGER(i4) :: kk  ! absolute value of lag k
    INTEGER(i4) :: nn  ! number of true values in mask
    REAL(sp) :: n   ! real of nn
    REAL(sp) :: mean, var  ! moments
    REAL(sp), DIMENSION(size(x)) :: x_shift
    LOGICAL, DIMENSION(size(x)) :: maske

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error autocorr_sp: size(mask) /= size(x)'
      maske = mask
    end if

    ! calculate mean and population variance of x
    call moment(x, mean=mean, variance=var, mask=maske, sample=.false.)
    ! shift x to 0 mean
    x_shift = x - mean
    ! use absolute value of k
    kk = abs(k)
    nn = size(x)
    n = real(count(maske(1 : nn - kk).and.maske(1 + kk : nn)), sp)
    ! calculate the auto-correlation function
    acf = sum(x_shift(1 : nn - kk) * x_shift(1 + kk : nn), mask = (maske(1 : nn - kk) .and. maske(1 + kk : nn))) / n / var

  END FUNCTION autocorr_sp

  FUNCTION autocorr_1d_dp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: k
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4) :: i
    REAL(dp), DIMENSION(size(k)) :: acf

    if (present(mask)) then
      do i = 1, size(k)
        acf(i) = autocorr(x, k(i), mask)
      end do
    else
      do i = 1, size(k)
        acf(i) = autocorr(x, k(i))
      end do
    end if


  END FUNCTION autocorr_1d_dp

  FUNCTION autocorr_1d_sp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: k
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4) :: i
    REAL(sp), DIMENSION(size(k)) :: acf

    if (present(mask)) then
      do i = 1, size(k)
        acf(i) = autocorr(x, k(i), mask)
      end do
    else
      do i = 1, size(k)
        acf(i) = autocorr(x, k(i))
      end do
    end if

  END FUNCTION autocorr_1d_sp

END MODULE mo_corr
