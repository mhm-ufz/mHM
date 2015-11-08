MODULE mo_corr

  ! This module provides auto- and crosscorrelation function calculations

  ! Literature
  !   Corr, FFT
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  ! Written March 2011, Matthias Cuntz

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

  ! Copyright 2011 Matthias Cuntz


  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred.

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  USE mo_kind,      ONLY: i4, sp, dp, spc, dpc
  USE mo_constants, ONLY: TWOPI_sp, TWOPI_dp, PI_dp

  Implicit NONE

  PUBLIC :: autocoeffk   ! coeff_k so that autocorr = coeff_k/coeff_0
  PUBLIC :: autocorr     ! Autocorrelation coefficient at lag k = autocoeffk(k)/autocoeffk(0)
  PUBLIC :: corr         ! Correlation (function) with optional highpass filtering: covariance=correlation(1)/n
  PUBLIC :: crosscoeffk  ! coeff_k so that crosscorr = coeff_k/coeff_0, crosscoeffk(0) = covariance
  PUBLIC :: crosscorr    ! Crosscorrelation coefficient at lag k = crosscoeffk(k)/crosscoeffk(0)

  ! ------------------------------------------------------------------

  !     NAME
  !         autocoeffk

  !     PURPOSE
  !         Coefficient at lag k so that autocorrelation coefficient at lag k
  !         is autocoeffk(x,k)/autocoeffk(x,0).
  !
  !         If an optinal mask is given, the calculations are only over those locations that correspond
  !         to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         ak = autocoeffk(x, k, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: x(:)        Time series
  !         integer(i4) :: k[(:)]      Lag for autocorrelation

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: ak[(:)]     coefficient so that ak/autocoeffk(x,0) is the autocorrelation coefficient at lag k

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
  !         ! Variance = autocoeffk(0)
  !         var = autocoeffk(x,0)
  !         -> see also example in test directory

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  !         Modified, Stephan Thober, Nov 2012 - added 1d version
  INTERFACE autocoeffk
     MODULE PROCEDURE autocoeffk_sp, autocoeffk_dp, &
          autocoeffk_1d_dp, autocoeffk_1d_sp
  END INTERFACE autocoeffk

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

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  !         Modified, Stephan Thober, Nov 2012 - added 1d version
  INTERFACE autocorr
     MODULE PROCEDURE autocorr_sp, autocorr_dp, &
          autocorr_1d_sp, autocorr_1d_dp
  END INTERFACE autocorr

  ! ------------------------------------------------------------------

  !     NAME
  !         corr

  !     PURPOSE
  !         Computes the correlation of two real data sets data1 and data2 of length N (includ-
  !         ing any user-supplied zero padding) with Fast Fourier Transform (FFT). N must be an integer
  !         power of 2 for the FFT routine. Corr takes only the elements up to the last power of 2 and
  !         returns the number of elements used in nadjust.
  !         The answer is returned as the function corr, an array of length N (nadjust).
  !         The answer is stored in wrap-around order, i.e., correlations at increasingly negative lags
  !         are in corr(N) on down to corr(N/2+1), while correlations at increasingly positive lags are
  !         in correl(1) (zero lag) on up to correl(N/2). Sign convention of this routine: if data1 lags
  !         data2, i.e., is shifted to the right of it, then correl will show a peak at positive lags.
  !
  !         Optional high-pass filtering of the time series in Fourier space is implemented.
  !
  !         Note covariance(x,y) = corr(x,y)/n

  !     CALLING SEQUENCE
  !         cfunc = corr(data1,data2,nadjust=nadjust,nhigh=nhigh,nwin=nwin)

  !     INTENT(IN)
  !         real(sp/dp) :: data1(:)             1st time series
  !         real(sp/dp) :: data2(:)             2nd time series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: corr(size(data1))    Correlation function between data1 and data2 in wrap-around order

  !     INTENT(IN), OPTIONAL
  !         integer(i4) :: nhigh                If >0 then nhigh upper frequencies are filtered. nwin then defines
  !                                             the used window function for filtering. (default: 0)
  !         integer(i4) :: nwin                 Window function for highpass filtering (default: 1)
  !                                             0: no filtering
  !                                             1: ideal highpass, i.e. cut out the nhigh upper frequencies
  !                                             2: linear interpolation of 0 to 1 from highest to highest-nhigh
  !                                                frequency; similar Bartlett window

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         integer(i4) :: nadjust              Actual used number of elements.

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Covariance function: covariance(x,y) = corr(x,y)/n
  !         wT = corr(w, T, nadjust=nwT)
  !         wT(1:nwT) = wT(1:nwT) / real(nwT,dp)
  !         -> see also example in test directory

  !     LITERATURE
  !         WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  INTERFACE corr
     MODULE PROCEDURE corr_sp, corr_dp
  END INTERFACE corr

  ! ------------------------------------------------------------------

  !     NAME
  !         crosscoeffk

  !     PURPOSE
  !         Coefficient at lag k so that crosscorrelation coefficient at lag k
  !         is crosscoeffk(x,y,k)/crosscoeffk(x,y,0).
  !         -> crosscoeffk(x,y,0) = covariance(x,y)
  !
  !         If an optinal mask is given, the calculations are only over those locations that correspond
  !         to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         ck = crosscoeffk(x, y, k, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: x(:)        1st time series
  !         real(sp/dp) :: y(:)        2nd time series
  !         integer(i4) :: k           Lag for crosscorrelation

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: ck          coefficient so that ck/crosscoeffk(x,0) is the crosscorrelation coefficient at lag k

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
  !         ! covariance = crosscoeffk(0)
  !         cov = crosscoeffk(x,y,0)
  !         -> see also example in test directory

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  INTERFACE crosscoeffk
     MODULE PROCEDURE crosscoeffk_sp, crosscoeffk_dp
  END INTERFACE crosscoeffk

  ! ------------------------------------------------------------------

  !     NAME
  !         crosscorr

  !     PURPOSE
  !         Element at lag k of crosscorrelation function
  !             crosscorr(x,y,k) = crosscoeffk(x,y,k)/crosscoeffk(x,y,0).
  !
  !         If an optinal mask is given, the calculations are only over those locations that correspond
  !         to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         ck = crosscorr(x, y, k, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: x(:)        1st time series
  !         real(sp/dp) :: y(:)        2nd time series
  !         integer(i4) :: k           Lag for crosscorrelation

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: ck          Coefficient of crosscorrelation function at lag k

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
  !         ! Last crosscorrelation element
  !         ccorr = crosscorr(x,y,size(x)/2)
  !         -> see also example in test directory

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  INTERFACE crosscorr
     MODULE PROCEDURE crosscorr_sp, crosscorr_dp
  END INTERFACE crosscorr

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

  ! Private routines, mostly from numerical recipes
  INTERFACE arth
     MODULE PROCEDURE arth_sp, arth_dp, arth_i4
  END INTERFACE arth
  INTERFACE four1
     MODULE PROCEDURE four1_sp, four1_dp
  END INTERFACE four1
  INTERFACE fourrow
     MODULE PROCEDURE fourrow_sp, fourrow_dp
  END INTERFACE fourrow
  INTERFACE realft
     MODULE PROCEDURE realft_sp, realft_dp
  END INTERFACE realft
  INTERFACE swap
     MODULE PROCEDURE & ! swap_i4, &
          !swap_sp, swap_1d_sp, &
          !swap_dp, & !swap_1d_dp, &
          ! swap_spc, &
          swap_1d_spc, &
          ! swap_dpc, &
          swap_1d_dpc !, &
          ! masked_swap_sp, masked_swap_1d_sp, masked_swap_2d_sp, &
          ! masked_swap_dp, masked_swap_1d_dp, masked_swap_2d_dp, &
          ! masked_swap_spc, masked_swap_1d_spc, masked_swap_2d_spc, &
          ! masked_swap_dpc, masked_swap_1d_dpc, masked_swap_2d_dpc
  END INTERFACE swap
  !INTERFACE zroots_unity
  !   MODULE PROCEDURE zroots_unity_sp, zroots_unity_dp
  !END INTERFACE zroots_unity

  INTEGER(i4), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Returns an array of length n containing an arithmetic progression whose
  ! first value is first and whose increment is increment. If first and
  ! increment have rank greater than zero, returns an array of one larger rank,
  ! with the last subscript having size n and indexing the progressions.
  FUNCTION arth_sp(first,increment,n)

    IMPLICIT NONE

    REAL(sp),    INTENT(IN) :: first, increment
    INTEGER(i4), INTENT(IN) :: n
    REAL(sp), DIMENSION(n)  :: arth_sp

    INTEGER(i4) :: k, k2
    REAL(sp) :: temp

    if (n > 0) arth_sp(1)=first
    if (n <= NPAR_ARTH) then
       do k=2, n
          arth_sp(k) = arth_sp(k-1)+increment
       end do
    else
       do k=2, NPAR2_ARTH
          arth_sp(k) = arth_sp(k-1)+increment
       end do
       temp = increment*NPAR2_ARTH
       k = NPAR2_ARTH
       do
          if (k >= n) exit
          k2 = k+k
          arth_sp(k+1:min(k2,n)) = temp+arth_sp(1:min(k,n-k))
          temp = temp+temp
          k = k2
       end do
    end if

  END FUNCTION arth_sp


  FUNCTION arth_dp(first,increment,n)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN) :: first,increment
    INTEGER(i4), INTENT(IN) :: n
    REAL(dp), DIMENSION(n)  :: arth_dp

    INTEGER(i4) :: k, k2
    REAL(dp) :: temp

    if (n > 0) arth_dp(1)=first
    if (n <= NPAR_ARTH) then
       do k=2, n
          arth_dp(k) = arth_dp(k-1)+increment
       end do
    else
       do k=2, NPAR2_ARTH
          arth_dp(k) = arth_dp(k-1)+increment
       end do
       temp = increment*NPAR2_ARTH
       k = NPAR2_ARTH
       do
          if (k >= n) exit
          k2 = k+k
          arth_dp(k+1:min(k2,n)) = temp+arth_dp(1:min(k,n-k))
          temp = temp+temp
          k = k2
       end do
    end if

  END FUNCTION arth_dp


  FUNCTION arth_i4(first,increment,n)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)   :: first, increment, n
    INTEGER(i4), DIMENSION(n) :: arth_i4

    INTEGER(i4) :: k, k2, temp

    if (n > 0) arth_i4(1)=first
    if (n <= NPAR_ARTH) then
       do k=2, n
          arth_i4(k) = arth_i4(k-1)+increment
       end do
    else
       do k=2, NPAR2_ARTH
          arth_i4(k) = arth_i4(k-1)+increment
       end do
       temp = increment*NPAR2_ARTH
       k = NPAR2_ARTH
       do
          if (k >= n) exit
          k2 = k+k
          arth_i4(k+1:min(k2,n)) = temp+arth_i4(1:min(k,n-k))
          temp = temp+temp
          k = k2
       end do
    end if

  END FUNCTION arth_i4

  ! ------------------------------------------------------------------

  FUNCTION autocoeffk_dp(x, k, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: autocoeffk_dp

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocoeffk_dp: size(mask) /= size(x)'
       autocoeffk_dp = crosscoeffk(x, x, k, mask)
    else
       autocoeffk_dp = crosscoeffk(x, x, k)
    endif

  END FUNCTION autocoeffk_dp


  FUNCTION autocoeffk_sp(x, k, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: autocoeffk_sp

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocoeffk_sp: size(mask) /= size(x)'
       autocoeffk_sp = crosscoeffk(x, x, k, mask)
    else
       autocoeffk_sp = crosscoeffk(x, x, k)
    endif

  END FUNCTION autocoeffk_sp

  FUNCTION autocoeffk_1d_dp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4), DIMENSION(:),           INTENT(IN)  :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    INTEGER(i4)                                      :: i
    REAL(dp),    DIMENSION(size(k))                  :: acf

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocoeffk_1d_dp: size(mask) /= size(x)'
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i), mask)
       end do
    else
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i))
       end do
    endif

  END FUNCTION autocoeffk_1d_dp


  FUNCTION autocoeffk_1d_sp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4), DIMENSION(:),           INTENT(IN)  :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    INTEGER(i4)                                      :: i
    REAL(sp),    DIMENSION(size(k))                  :: acf

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocoeffk_1d_sp: size(mask) /= size(x)'
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i), mask)
       end do
    else
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i))
       end do
    endif

  END FUNCTION autocoeffk_1d_sp

  ! ------------------------------------------------------------------

  FUNCTION autocorr_dp(x, k, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: autocorr_dp

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocorr_1d_dp: size(mask) /= size(x)'
       autocorr_dp = crosscoeffk(x, x, k, mask) / crosscoeffk(x, x, 0, mask)
    else
       autocorr_dp = crosscoeffk(x, x, k) / crosscoeffk(x, x, 0)
    endif

  END FUNCTION autocorr_dp


  FUNCTION autocorr_sp(x, k, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: autocorr_sp

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocorr_1d_sp: size(mask) /= size(x)'
       autocorr_sp = crosscoeffk(x, x, k, mask) / crosscoeffk(x, x, 0, mask)
    else
       autocorr_sp = crosscoeffk(x, x, k) / crosscoeffk(x, x, 0)
    endif

  END FUNCTION autocorr_sp

  FUNCTION autocorr_1d_dp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(dp),   DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),DIMENSION(:),           INTENT(IN)  :: k
    LOGICAL,    DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    INTEGER(i4)                                     :: i
    REAL(dp),   DIMENSION(size(k))                  :: acf
    REAL(dp)                                        :: c0

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocorr_dp: size(mask) /= size(x)'
       c0 = crosscoeffk(x, x, 0, mask)
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i), mask) / c0
       end do
    else
       c0 = crosscoeffk(x, x, 0)
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i)) / c0
       end do
    endif

  END FUNCTION autocorr_1d_dp

  FUNCTION autocorr_1d_sp(x, k, mask) result(acf)

    IMPLICIT NONE

    REAL(sp),   DIMENSION(:),           INTENT(IN)  :: x
    INTEGER(i4),DIMENSION(:),           INTENT(IN)  :: k
    LOGICAL,    DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    INTEGER(i4)                                     :: i
    REAL(sp),   DIMENSION(size(k))                  :: acf
    REAL(sp)                                        :: c0

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error autocorr_sp: size(mask) /= size(x)'
       c0 = crosscoeffk(x, x, 0, mask)
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i), mask) / c0
       end do
    else
       c0 = crosscoeffk(x, x, 0)
       do i = 1, size(k)
          acf(i) = crosscoeffk(x, x, k(i)) / c0
       end do
    endif

  END FUNCTION autocorr_1d_sp

  ! ------------------------------------------------------------------

  FUNCTION corr_dp(data1,data2,nadjust,nhigh,nwin)

    IMPLICIT NONE

    REAL(dp),     DIMENSION(:),                     INTENT(IN)  :: data1, data2
    INTEGER(i4),                          OPTIONAL, INTENT(OUT) :: nadjust
    INTEGER(i4),                          OPTIONAL, INTENT(IN)  :: nhigh
    INTEGER(i4),                          OPTIONAL, INTENT(IN)  :: nwin
    REAL(dp),     DIMENSION(size(data1))                        :: corr_dp

    REAL(dp),     DIMENSION(:), ALLOCATABLE :: dat1, dat2, corrout
    COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: cdat1, cdat2
    REAL(dp),     DIMENSION(:), ALLOCATABLE :: win1
    !COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: cwin1
    INTEGER(i4) :: i, n, no2, iwin, ihigh
    REAL(dp)    :: no2_1, ihigh1

    n = size(data1)
    if (size(data2) /= n) stop 'Error corr_dp: size(data1) /= size(data2)'

    if (present(nwin)) then
       iwin = nwin
    else
       iwin = 1
    endif

    if (present(nhigh)) then
       ihigh = nhigh
    else
       ihigh = 0
    endif

    if (iand(n,n-1) /= 0) then
       if (present(nadjust)) then
          n = 2**floor(log(real(n,dp))/log(2.0_dp))
          nadjust = n
       else
          stop 'Error corr_dp: size(data1) must be a power of 2'
       endif
    else
       if (present(nadjust)) then
          nadjust = n
       endif
    endif

    allocate(dat1(n))
    allocate(dat2(n))
    dat1 = data1(1:n)
    dat2 = data2(1:n)

    no2   = n/2
    no2_1 = 1.0_dp / real(no2,dp)
    allocate(cdat1(no2))
    allocate(cdat2(no2))
    allocate(corrout(n))

    ! FFT
    call realft(dat1,1,cdat1)
    call realft(dat2,1,cdat2)

    ! Highpass
    if (ihigh > 0) then
       ! FxH
       allocate(win1(no2))
       !allocate(cwin1(no2))
       select case(iwin)
       case(0) ! no window
          win1(1:no2)  = 1.0_dp
       case(1) ! ideal high pass filter
          win1(1:ihigh)     = 0.0_dp
          win1(ihigh+1:no2) = 1.0_dp
       case(2) ! similar Bartlett window
          ihigh1 = 1.0_dp / real(ihigh,dp)
          forall(i=1:ihigh) win1(i) = real(i-1,dp) * ihigh1
          win1(ihigh+1:no2) = 1.0_dp
       case default
          stop 'Unimplemented window option in corr_dp'
       end select
       !cwin1 = cmplx(win1, win1, kind=dpc)
       ! low pass
       ! cdat1(1:no2) = cdat1(1:no2)*cwin1(1:no2)
       ! cdat2(1:no2) = cdat2(1:no2)*cwin1(1:no2)
       cdat1(1:no2) = cdat1(1:no2) * win1(1:no2)
       cdat2(1:no2) = cdat2(1:no2) * win1(1:no2)
    endif

    ! FxF*
    cdat1(1)  = cmplx(real(cdat1(1))*real(cdat2(1))*no2_1, &
         aimag(cdat1(1))*aimag(cdat2(1))*no2_1, kind=dpc)
    cdat1(2:) = cdat1(2:)*conjg(cdat2(2:))*no2_1

    ! IFFT
    call realft(corrout,-1,cdat1)
    corr_dp(1:n) = corrout(1:n)
    if (size(corr_dp) > n) corr_dp(n+1:) = 0.0_dp

    deallocate(dat1)
    deallocate(dat2)
    deallocate(cdat1)
    deallocate(cdat2)
    deallocate(corrout)
    if (ihigh > 0) then
       deallocate(win1)
       !deallocate(cwin1)
    endif

  END FUNCTION corr_dp


  FUNCTION corr_sp(data1,data2,nadjust,nhigh,nwin)

    IMPLICIT NONE

    REAL(sp),     DIMENSION(:),                     INTENT(IN)  :: data1, data2
    INTEGER(i4),                          OPTIONAL, INTENT(OUT) :: nadjust
    INTEGER(i4),                          OPTIONAL, INTENT(IN)  :: nhigh
    INTEGER(i4),                          OPTIONAL, INTENT(IN)  :: nwin
    REAL(sp),     DIMENSION(size(data1))                        :: corr_sp

    REAL(sp),     DIMENSION(:), ALLOCATABLE :: dat1, dat2, corrout
    COMPLEX(spc), DIMENSION(:), ALLOCATABLE :: cdat1, cdat2
    REAL(sp),     DIMENSION(:), ALLOCATABLE :: win1
    !COMPLEX(spc), DIMENSION(:), ALLOCATABLE :: cwin1
    INTEGER(i4) :: i, n, no2, iwin, ihigh
    REAL(sp)    :: no2_1, ihigh1

    n = size(data1)
    if (size(data2) /= n) stop 'Error corr_sp: size(data1) /= size(data2)'

    if (present(nwin)) then
       iwin = nwin
    else
       iwin = 1
    endif

    if (present(nhigh)) then
       ihigh = nhigh
    else
       ihigh = 0
    endif

    if (iand(n,n-1) /= 0) then
       if (present(nadjust)) then
          n = 2**floor(log(real(n,sp))/log(2.0_sp))
          nadjust = n
       else
          stop 'Error corr_sp: size(data1) must be a power of 2'
       endif
    else
       if (present(nadjust)) then
          nadjust = n
       endif
    endif

    allocate(dat1(n))
    allocate(dat2(n))
    dat1 = data1(1:n)
    dat2 = data2(1:n)

    no2   = n/2
    no2_1 = 1.0_sp / real(no2,sp)
    allocate(cdat1(no2))
    allocate(cdat2(no2))
    allocate(corrout(n))

    ! FFT
    call realft(dat1,1,cdat1)
    call realft(dat2,1,cdat2)

    ! Highpass
    if (ihigh > 0) then
       ! FxH
       allocate(win1(no2))
       !allocate(cwin1(no2))
       select case(iwin)
       case(0) ! no window
          win1(1:no2)  = 1.0_sp
       case(1) ! ideal high pass filter
          win1(1:ihigh)     = 0.0_sp
          win1(ihigh+1:no2) = 1.0_sp
       case(2) ! similar Bartlett window
          ihigh1 = 1.0_sp / real(ihigh,sp)
          forall(i=1:ihigh) win1(i) = real(i-1,sp) * ihigh1
          win1(ihigh+1:no2) = 1.0_sp
       case default
          stop 'Unimplemented window option in corr_sp'
       end select
       !cwin1 = cmplx(win1, win1, kind=spc)
       ! low pass
       ! cdat1(1:no2) = cdat1(1:no2)*cwin1(1:no2)
       ! cdat2(1:no2) = cdat2(1:no2)*cwin1(1:no2)
       cdat1(1:no2) = cdat1(1:no2) * win1(1:no2)
       cdat2(1:no2) = cdat2(1:no2) * win1(1:no2)
    endif

    ! FxF*
    cdat1(1)  = cmplx(real(cdat1(1))*real(cdat2(1))*no2_1, &
         aimag(cdat1(1))*aimag(cdat2(1))*no2_1, kind=spc)
    cdat1(2:) = cdat1(2:)*conjg(cdat2(2:))*no2_1

    ! IFFT
    call realft(corrout,-1,cdat1)
    corr_sp(1:n) = corrout(1:n)
    if (size(corr_sp) > n) corr_sp(n+1:) = 0.0_sp

    deallocate(dat1)
    deallocate(dat2)
    deallocate(cdat1)
    deallocate(cdat2)
    deallocate(corrout)
    if (ihigh > 0) then
       deallocate(win1)
       !deallocate(cwin1)
    endif

  END FUNCTION corr_sp

  ! ------------------------------------------------------------------

  FUNCTION crosscoeffk_dp(x, y, k, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x
    REAL(dp), DIMENSION(:),           INTENT(IN)  :: y
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: crosscoeffk_dp

    INTEGER(i4) :: nn  ! number of true values in mask
    INTEGER(i4) :: nnn ! number of true values in mask .and. shifted mask by lag k
    REAL(dp)    :: n   ! real of nn or nnn
    INTEGER(i4) :: kk  ! absolute value of lag k
    REAL(dp)    :: ave

    REAL(dp),  DIMENSION(size(x)) :: xdash
    REAL(dp),  DIMENSION(size(x)) :: ydash
    LOGICAL,   DIMENSION(size(x)) :: maske

    maske(:) = .true.
    if (present(mask)) then
       if (size(x) /= size(y)) stop 'Error crosscoeffk_dp: size(x) /= size(y)'
       if (size(mask) /= size(x)) stop 'Error crosscoeffk_dp: size(mask) /= size(x)'
       maske = mask
    endif

    nn    = count(maske)
    n     = real(nn,dp)
    ! crosscoeffk(x, y, k) = crosscoeffk(y, x, -k)
    if (k >= 0) then
       ave   = sum(x(:), mask=maske)/n
       xdash = x - ave
       ave   = sum(y(:), mask=maske)/n
       ydash = y - ave
    else
       ave   = sum(y(:), mask=maske)/n
       xdash = y - ave
       ave   = sum(x(:), mask=maske)/n
       ydash = x - ave
    endif
    kk = abs(k)
    nnn = size(x,1)
    n  = real(count(maske(1:nnn-kk).and.maske(1+kk:nnn)),dp)
    crosscoeffk_dp = sum(xdash(1:nnn-kk)*ydash(1+kk:nnn), mask=(maske(1:nnn-kk).and.maske(1+kk:nnn))) / n
    
  END FUNCTION crosscoeffk_dp

  FUNCTION crosscoeffk_sp(x, y, k, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x
    REAL(sp), DIMENSION(:),           INTENT(IN)  :: y
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: crosscoeffk_sp

    INTEGER(i4) :: nn  ! number of true values in mask
    INTEGER(i4) :: nnn ! number of true values in mask .and. shifted mask by lag k
    REAL(sp)    :: n   ! real of nn or nnn
    INTEGER(i4) :: kk  ! absolute value of lag k
    REAL(sp)    :: ave

    REAL(sp),  DIMENSION(size(x)) :: xdash
    REAL(sp),  DIMENSION(size(x)) :: ydash
    LOGICAL,   DIMENSION(size(x)) :: maske

    maske(:) = .true.
    if (present(mask)) then
       if (size(x) /= size(y)) stop 'Error crosscoeffk_sp: size(x) /= size(y)'
       if (size(mask) /= size(x)) stop 'Error crosscoeffk_sp: size(mask) /= size(x)'
       maske = mask
    endif

    nn    = count(maske)
    n     = real(nn,sp)
    ! crosscoeffk(x, y, k) = crosscoeffk(y, x, -k)
    if (k >= 0) then
       ave   = sum(x(:), mask=maske)/n
       xdash = x - ave
       ave   = sum(y(:), mask=maske)/n
       ydash = y - ave
    else
       ave   = sum(y(:), mask=maske)/n
       xdash = y - ave
       ave   = sum(x(:), mask=maske)/n
       ydash = x - ave
    endif
    kk = abs(k)
    nnn = size(x,1)
    n  = real(count(maske(1:nnn-kk).and.maske(1+kk:nnn)),sp)
    crosscoeffk_sp = sum(xdash(1:nnn-kk)*ydash(1+kk:nnn), mask=(maske(1:nnn-kk).and.maske(1+kk:nnn))) / n
    
  END FUNCTION crosscoeffk_sp

  ! ------------------------------------------------------------------

  FUNCTION crosscorr_dp(x, y, k, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x
    REAL(dp), DIMENSION(:),           INTENT(IN)  :: y
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: crosscorr_dp

    if (size(x) /= size(y)) stop 'Error crosscorr_dp: size(x) /= size(y)'
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error crosscorr_dp: size(mask) /= size(x)'
       crosscorr_dp = crosscoeffk(x, y, k, mask) / crosscoeffk(x, y, 0, mask)
    else
       crosscorr_dp = crosscoeffk(x, y, k) / crosscoeffk(x, y, 0)
    endif

  END FUNCTION crosscorr_dp

  FUNCTION crosscorr_sp(x, y, k, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x
    REAL(sp), DIMENSION(:),           INTENT(IN)  :: y
    INTEGER(i4),                      INTENT(IN)  :: k
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: crosscorr_sp

    if (size(x) /= size(y)) stop 'Error crosscorr_sp: size(x) /= size(y)'
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error crosscorr_sp: size(mask) /= size(x)'
       crosscorr_sp = crosscoeffk(x, y, k, mask) / crosscoeffk(x, y, 0, mask)
    else
       crosscorr_sp = crosscoeffk(x, y, k) / crosscoeffk(x, y, 0)
    endif

  END FUNCTION crosscorr_sp

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Replaces a complex array data by its discrete Fourier transform, if isign is input as 1;
  ! or replaces data by its inverse discrete Fourier transform times the size of data, if isign
  ! is input as -1. The size of data must be an integer power of 2. Parallelismis achieved
  ! by internally reshaping the input array to two dimensions. (Use this version if fourrow is
  ! faster than fourcol on your machine.)

  SUBROUTINE four1_sp(data,isign)

    IMPLICIT NONE

    COMPLEX(spc), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(i4),                INTENT(IN)    :: isign

    COMPLEX(spc), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(dpc), DIMENSION(:),   ALLOCATABLE :: w,wp
    REAL(dp),     DIMENSION(:),   ALLOCATABLE :: theta
    INTEGER(i4) :: n, m1, m2, j

    n = size(data)
    if (iand(n,n-1) /= 0) stop 'Error four1_sp: size(data1) must be a power of 2'

    m1 = 2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2 = n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat = reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta = arth(0,isign,m1)*TWOPI_dp/real(n,dp)
    wp    = cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w     = cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2, m2
       w = w*wp+w
       dat(:,j) = dat(:,j)*cmplx(w,kind=spc)
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)

  END SUBROUTINE four1_sp


  SUBROUTINE four1_dp(data,isign)

    IMPLICIT NONE

    COMPLEX(dpc), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(i4),                INTENT(IN)    :: isign

    COMPLEX(dpc), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(dpc), DIMENSION(:),   ALLOCATABLE :: w,wp
    REAL(dp),     DIMENSION(:),   ALLOCATABLE :: theta
    INTEGER(i4) :: n, m1, m2, j

    n=size(data)
    if (iand(n,n-1) /= 0) stop 'Error four1_dp: size(data1) must be a power of 2'

    m1 = 2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2 = n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat = reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta = arth(0,isign,m1)*TWOPI_dp/real(n,dp)
    wp    = cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w     = cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2, m2
       w = w*wp+w
       dat(:,j) = dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)

  END SUBROUTINE four1_dp

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Replaces each row (constant first index) of data(1:M,1:N) by its discrete Fourier trans-
  ! form (transform on second index), if isign is input as 1; or replaces each row of data
  ! by N times its inverse discrete Fourier transform, if isign is input as -1. N must be an
  ! integer power of 2. Parallelism is M-fold on the first index of data.

  SUBROUTINE fourrow_sp(data,isign)

    IMPLICIT NONE

    COMPLEX(spc), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(i4),                  INTENT(IN)    :: isign

    INTEGER(i4) :: n, i, istep, j, m, mmax, n2
    REAL(dp)    :: theta
    COMPLEX(spc), DIMENSION(size(data,1)) :: temp
    COMPLEX(dpc) :: w, wp
    COMPLEX(spc) :: ws

    n=size(data,2)
    if (iand(n,n-1) /= 0) stop 'Error fourrow_sp: size(data,2) must be a power of 2'
    n2 = n/2
    j = n2
    do i=1, n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m = n2
       do
          if (m < 2 .or. j < m) exit
          j = j-m
          m = m/2
       end do
       j = j+m
    end do
    mmax = 1
    do
       if (n <= mmax) exit
       istep = 2*mmax
       theta = PI_dp/real(isign*mmax,dp)
       wp = cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w  = cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1, mmax
          ws = cmplx(w,kind=spc)
          do i=m, n, istep
             j = i+mmax
             temp = ws*data(:,j)
             data(:,j) = data(:,i)-temp
             data(:,i) = data(:,i)+temp
          end do
          w = w*wp+w
       end do
       mmax = istep
    end do

  END SUBROUTINE fourrow_sp


  SUBROUTINE fourrow_dp(data,isign)

    IMPLICIT NONE

    COMPLEX(dpc), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(i4),                  INTENT(IN)    :: isign

    INTEGER(i4) :: n, i, istep, j, m, mmax, n2
    REAL(dp)    :: theta
    COMPLEX(dpc), DIMENSION(size(data,1)) :: temp
    COMPLEX(dpc) :: w, wp
    COMPLEX(dpc) :: ws

    n=size(data,2)
    if (iand(n,n-1) /= 0) stop 'Error fourrow_dp: size(data,2) must be a power of 2'
    n2 = n/2
    j  = n2
    do i=1, n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m = n2
       do
          if (m < 2 .or. j < m) exit
          j = j-m
          m = m/2
       end do
       j = j+m
    end do
    mmax = 1
    do
       if (n <= mmax) exit
       istep = 2*mmax
       theta = PI_dp/real(isign*mmax,dp)
       wp = cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w  = cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1, mmax
          ws = w
          do i=m, n, istep
             j = i+mmax
             temp = ws*data(:,j)
             data(:,j) = data(:,i)-temp
             data(:,i) = data(:,i)+temp
          end do
          w = w*wp+w
       end do
       mmax = istep
    end do

  END SUBROUTINE fourrow_dp

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! When isign=1, calculates the Fourier transform of a set of N real-valued data points,
  ! input in the array data. If the optional argument zdata is not present, the data are replaced
  ! by the positive frequency half of its complex Fourier transform. There al-valued first and
  ! last components of the complex transform are returned as elements data(1) and data(2),
  ! respectively. If the complex array zdata of lengthN/2 ispresent, data is unchanged and
  ! the transform is returned in zdata. N must be a power of 2. If isign=-1, this routine
  ! calculates the inverse transform of a complex data array if it is the transform of real data.
  ! (Result in this case must be multiplied by 2/N.) The data can be supplied either in data,
  ! with zdata absent, or inzdata.
  SUBROUTINE realft_dp(data,isign,zdata)

    IMPLICIT NONE

    REAL(dp),     DIMENSION(:), INTENT(INOUT)    :: data
    INTEGER(i4),                INTENT(IN)       :: isign
    COMPLEX(dpc), DIMENSION(:), OPTIONAL, TARGET :: zdata

    INTEGER(i4) :: n, nh, nq
    COMPLEX(dpc), DIMENSION(size(data)/4)   :: w
    COMPLEX(dpc), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(dpc), DIMENSION(:), POINTER :: cdata
    COMPLEX(dpc) :: z
    REAL(dp) :: c1=0.5_dp, c2

    n = size(data)
    if (iand(n,n-1) /= 0) stop 'Error realft_dp: size(data) must be a power of 2'
    nh = n/2
    nq = n/4
    if (present(zdata)) then
       if (n/2 /= size(zdata)) stop 'Error realft_dp: size(zdata) /= size(data)/2'
       cdata => zdata
       if (isign == 1) cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=dpc)
    else
       allocate(cdata(n/2))
       cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=dpc)
    end if
    if (isign == 1) then
       c2 = -0.5_dp
       call four1(cdata,+1)
    else
       c2 = 0.5_dp
    end if
    w  = zroots_unity_dp(sign(n,isign),n/4)
    w  = cmplx(-aimag(w),real(w),kind=dpc)
    h1 = c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2 = c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)       = h1+w(2:nq)*h2
    cdata(nh:nq+2:-1) = conjg(h1-w(2:nq)*h2)
    z = cdata(1)
    if (isign == 1) then
       cdata(1) = cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
    else
       cdata(1) = cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
       call four1(cdata,-1)
    end if
    if (present(zdata)) then
       if (isign /= 1) then
          data(1:n-1:2) = real(cdata)
          data(2:n:2)   = aimag(cdata)
       end if
    else
       data(1:n-1:2) = real(cdata)
       data(2:n:2) = aimag(cdata)
       deallocate(cdata)
    end if

  END SUBROUTINE realft_dp


  SUBROUTINE realft_sp(data,isign,zdata)

    IMPLICIT NONE

    REAL(sp),     DIMENSION(:), INTENT(INOUT)    :: data
    INTEGER(i4),                INTENT(IN)       :: isign
    COMPLEX(spc), DIMENSION(:), OPTIONAL, TARGET :: zdata

    INTEGER(i4) :: n, nh, nq
    COMPLEX(spc), DIMENSION(size(data)/4)   :: w
    COMPLEX(spc), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(spc), DIMENSION(:), POINTER :: cdata
    COMPLEX(spc) :: z
    REAL(sp) :: c1=0.5_sp, c2

    n = size(data)
    if (iand(n,n-1) /= 0) stop 'Error realft_sp: size(data) must be a power of 2'
    nh = n/2
    nq = n/4
    if (present(zdata)) then
       if (n/2 /= size(zdata)) stop 'Error realft_sp: size(zdata) /= size(data)/2'
       cdata => zdata
       if (isign == 1) cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    else
       allocate(cdata(n/2))
       cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    end if
    if (isign == 1) then
       c2 = -0.5_sp
       call four1(cdata,+1)
    else
       c2 = 0.5_sp
    end if
    w  = zroots_unity_sp(sign(n,isign),n/4)
    w  = cmplx(-aimag(w),real(w),kind=spc)
    h1 = c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2 = c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)       = h1+w(2:nq)*h2
    cdata(nh:nq+2:-1) = conjg(h1-w(2:nq)*h2)
    z = cdata(1)
    if (isign == 1) then
       cdata(1) = cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
    else
       cdata(1) = cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
       call four1(cdata,-1)
    end if
    if (present(zdata)) then
       if (isign /= 1) then
          data(1:n-1:2) = real(cdata)
          data(2:n:2)   = aimag(cdata)
       end if
    else
       data(1:n-1:2) = real(cdata)
       data(2:n:2) = aimag(cdata)
       deallocate(cdata)
    end if

  END SUBROUTINE realft_sp

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Swaps the corresponding elements of a and b. If mask is present, performs
  ! the swap only where mask is true. (Following code is the unmasked case.
  ! For speed at runtime, the masked case is implemented by overloading, not
  ! by testing for the optional argument.)
  ! SUBROUTINE swap_i4(a,b)
  !   INTEGER(i4), INTENT(INOUT) :: a,b
  !   INTEGER(i4) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_i4


  ! SUBROUTINE swap_sp(a,b)
  !   REAL(sp), INTENT(INOUT) :: a,b
  !   REAL(sp) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_sp


  ! SUBROUTINE masked_swap_sp(a,b,mask)
  !   REAL(sp), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   REAL(sp) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_sp


  ! SUBROUTINE masked_swap_1d_sp(a,b,mask)
  !   REAL(sp), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   REAL(sp), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_1d_sp


  ! SUBROUTINE masked_swap_2d_sp(a,b,mask)
  !   REAL(sp), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   REAL(sp), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_2d_sp


  ! SUBROUTINE swap_1d_sp(a,b)
  !   REAL(sp), DIMENSION(:), INTENT(INOUT) :: a,b
  !   REAL(sp), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_1d_sp


  ! SUBROUTINE swap_dp(a,b)
  !   REAL(dp), INTENT(INOUT) :: a,b
  !   REAL(dp) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_dp


  ! SUBROUTINE swap_1d_dp(a,b)
  !   REAL(dp), DIMENSION(:), INTENT(INOUT) :: a,b
  !   REAL(dp), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_1d_dp


  ! SUBROUTINE masked_swap_dp(a,b,mask)
  !   REAL(dp), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   REAL(dp) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_dp


  ! SUBROUTINE masked_swap_1d_dp(a,b,mask)
  !   REAL(dp), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   REAL(dp), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_1d_dp


  ! SUBROUTINE masked_swap_2d_dp(a,b,mask)
  !   REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   REAL(dp), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_2d_dp


  ! SUBROUTINE swap_spc(a,b)
  !   COMPLEX(spc), INTENT(INOUT) :: a,b
  !   COMPLEX(spc) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_spc


  ! SUBROUTINE masked_swap_spc(a,b,mask)
  !   COMPLEX(spc), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   COMPLEX(spc) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_spc


  ! SUBROUTINE masked_swap_1d_spc(a,b,mask)
  !   COMPLEX(spc), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   COMPLEX(spc), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_1d_spc


  ! SUBROUTINE masked_swap_2d_spc(a,b,mask)
  !   COMPLEX(spc), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   COMPLEX(spc), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_2d_spc


  SUBROUTINE swap_1d_spc(a,b)
    COMPLEX(spc), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(spc), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_spc


  ! SUBROUTINE swap_dpc(a,b)
  !   COMPLEX(dpc), INTENT(INOUT) :: a,b
  !   COMPLEX(dpc) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_dpc


  ! SUBROUTINE masked_swap_dpc(a,b,mask)
  !   COMPLEX(dpc), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   COMPLEX(dpc) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_dpc


  ! SUBROUTINE masked_swap_1d_dpc(a,b,mask)
  !   COMPLEX(dpc), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   COMPLEX(dpc), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_1d_dpc


  ! SUBROUTINE masked_swap_2d_dpc(a,b,mask)
  !   COMPLEX(dpc), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   COMPLEX(dpc), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_2d_dpc


  SUBROUTINE swap_1d_dpc(a,b)
    COMPLEX(dpc), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(dpc), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_dpc

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Returns a complex array containing nn consecutive powers of the nth
  ! complex root of unity.
  FUNCTION zroots_unity_dp(n,nn)

    INTEGER(i4), INTENT(IN)     :: n,nn
    COMPLEX(dpc), DIMENSION(nn) :: zroots_unity_dp
    INTEGER(i4) :: k
    REAL(dp)    :: theta

    zroots_unity_dp(1) = 1.0_dp
    theta = TWOPI_dp/n
    k = 1
    do
       if (k >= nn) exit
       zroots_unity_dp(k+1) = cmplx(cos(k*theta),sin(k*theta),kind=dpc)
       zroots_unity_dp(k+2:min(2*k,nn)) = zroots_unity_dp(k+1) * &
            zroots_unity_dp(2:min(k,nn-k))
       k = 2*k
    end do

  END FUNCTION zroots_unity_dp


  FUNCTION zroots_unity_sp(n,nn)

    INTEGER(i4), INTENT(IN)     :: n,nn
    COMPLEX(spc), DIMENSION(nn) :: zroots_unity_sp
    INTEGER(i4) :: k
    REAL(sp)    :: theta

    zroots_unity_sp(1) = 1.0_sp
    theta = TWOPI_sp/n
    k = 1
    do
       if (k >= nn) exit
       zroots_unity_sp(k+1) = cmplx(cos(k*theta),sin(k*theta),kind=spc)
       zroots_unity_sp(k+2:min(2*k,nn)) = zroots_unity_sp(k+1) * &
            zroots_unity_sp(2:min(k,nn-k))
       k = 2*k
    end do

  END FUNCTION zroots_unity_sp

  ! ------------------------------------------------------------------

END MODULE mo_corr
