MODULE mo_mad

  ! This module provides a median absolute deviation (MAD) test

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

  USE mo_kind,       ONLY: i4, sp, dp
  USE mo_percentile, ONLY: median

  Implicit NONE

  PUBLIC :: mad             ! Mean absolute deviation test

  ! ------------------------------------------------------------------

  !     NAME
  !         mad

  !     PURPOSE
  !         Mean absolute deviation test with optional z-value (default: 7) and mask,
  !         and 3 optional variants: 0. raw values (default), 1. 1st derivative, 2. 2nd derivative
  !         Returns mask with true everywhere except where <(median-MAD*z/0.6745) or >(md+MAD*z/0.6745)
  !
  !         If an optinal mask is given, the mad test is only performed on those locations that correspond
  !         to true values in the mask.
  !
  !         If tout is given mad returns the array with the enteries exceeding the treshold
  !         being set to the threshold. With this setting arrays are accepted. 
  !         tout accepts: u. upper values are cut at the threshold,
  !         l. lower values are cut at the threshold, b. upper and lower values are cut at the threshold
  !         With this setting only the variant 0 is available (no argument implemented).

  !     CALLING SEQUENCE
  !         out = mad(vec, z=z, mask=mask, deriv=deriv)
  !
  !         out = mad(arr, z=z, mask=mask, tout=tout)

  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers
  !
  !         real(sp/dp) :: arr        nD-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         logical :: out            mask with true everywhere except where input deviates more
  !                                   than z standard deviations from median
  !
  !         arr     :: out            Array with values exceeding the threshold being trimmed. 

  !     INTENT(IN), OPTIONAL
  !         real(sp/dp) :: z          Input is allowed to deviate maximum z standard deviations from the median (default: 7)
  !         integer(i4) :: deriv      0: Act on raw input (default: 0)
  !                                   1: Use first derivatives
  !                                   2: Use 2nd derivatives
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.
  !                                   nD-array if tout is used.
  !         character(1):: tout       u: Trim only values above mad
  !                                   l: Trim only values below mad
  !                                   l: Trim values below and above mad

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         1st derivative is
  !             d = datin[1:n]-datin[0:n-1]
  !         because mean of left and right would give 0 for spikes.

  !     EXAMPLE
  !         vec = (/ -0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62, &
  !                   2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30, &
  !                   4.64,5.34,5.42,8.01 /)
  !         mask(:) = true
  !         ! Sets last entry false
  !         mask = mask .and. mad(vec, z=4., deriv=0, mask=mask)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar 2011
  !         mad_val added, Matthias Kelbling, May 2018
  INTERFACE mad
     MODULE PROCEDURE mad_sp, mad_dp, mad_val_dp, mad_val_sp
  END INTERFACE mad

  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------

  !     NAME
  !         mad_val

  !     PURPOSE
  !         Mean absolute deviation test with optional z-value (default: 7) and mask,
  !         and 3 optional variants: 0. raw values (default), 1. 1st derivative, 2. 2nd derivative
  !         Returns mask with true everywhere except where <(median-MAD*z/0.6745) or >(md+MAD*z/0.6745)
  !
  !         If an optinal mask is given, the mad test is only performed on those locations that correspond
  !         to true values in the mask.

  !     CALLING SEQUENCE
  !         out = mad(vec, z=z, mask=mask, deriv=deriv, mval = mval)

  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         logical :: out            mask with true everywhere except where input deviates more
  !                                   than z standard deviations from median

  !     INTENT(IN), OPTIONAL
  !         real(sp/dp) :: z          Input is allowed to deviate maximum z standard deviations from the median (default: 7)
  !         integer(i4) :: deriv      0: Act on raw input (default: 0)
  !                                   1: Use first derivatives
  !                                   2: Use 2nd derivatives
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.
  !         real(sp/dp) :: mval       Mask-Value. A value to mask for calculation.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         1st derivative is
  !             d = datin[1:n]-datin[0:n-1]
  !         because mean of left and right would give 0 for spikes.

  !     EXAMPLE
  !         vec = (/ -0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62, &
  !                   2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30, &
  !                   4.64,5.34,5.42,8.01 /)
  !         mask(:) = true
  !         ! Sets last entry false
  !         mask = mask .and. mad(vec, z=4., deriv=0, mask=mask)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Kelbling, Dec 2017
  !                   Using the mad-code from Matthias Cuntz

  ! ------------------------------------------------------------------


  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION mad_dp(arr, z, mask, deriv)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arr
    REAL(dp),                  OPTIONAL, INTENT(IN) :: z
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: deriv
    LOGICAL,     DIMENSION(size(arr))               :: mad_dp

    LOGICAL,  DIMENSION(size(arr)) :: maske
    REAL(dp), DIMENSION(size(arr)) :: d
    LOGICAL,  DIMENSION(size(arr)) :: dmask
    INTEGER(i4) :: n, ideriv
    !INTEGER(i4) :: m
    REAL(dp)    :: iz, med, mabsdev, thresh

    n = size(arr)
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= n) stop 'Error mad_dp: size(mask) /= size(arr)'
       maske = mask
    endif
    if (present(z)) then
       iz = z
    else
       iz = 7.0_dp
    endif
    if (present(deriv)) then
       ideriv = deriv
    else
       ideriv = 0
    endif

    select case(ideriv)
    case(0) ! pure values
       !m       = count(maske)
       med     = median(arr,mask=maske)
       mabsdev = median(abs(arr-med),mask=maske)
       thresh  = mabsdev * iz/0.6745_dp
       mad_dp     = (arr .ge. (med-thresh)) .and. (arr .le. (med+thresh)) .and. maske
    case(1) ! 1st derivative
       ! d(1:n-2) = 0.5_dp* (arr(3:n) - arr(1:n-2)) ! does not work -> ask Clemens & Matthias M
       d(1:n-1)     = arr(2:n) - arr(1:n-1)
       dmask(1:n-1) = maske(2:n) .and. maske(1:n-1)
       !m            = count(dmask(1:n-1))
       med          = median(d(1:n-1),mask=dmask(1:n-1))
       mabsdev      = median(abs(d(1:n-1)-med),mask=dmask(1:n-1))
       thresh       = mabsdev * iz/0.6745_dp
       mad_dp(n)       = .true.
       mad_dp(1:n-1)   = (d(1:n-1) .ge. (med-thresh)) .and. (d(1:n-1) .le. (med+thresh)) .and. dmask(1:n-1)
    case(2) ! -2nd derivative
       d(1:n-2)     = arr(2:n-1) + arr(2:n-1) - arr(1:n-2) - arr(3:n)
       dmask(1:n-2) = maske(2:n-1) .and. maske(1:n-2) .and. maske(3:n)
       !m            = count(dmask(1:n-2))
       med          = median(d(1:n-2),mask=dmask(1:n-2))
       mabsdev      = median(abs(d(1:n-2)-med),mask=dmask(1:n-2))
       thresh       = mabsdev * iz/0.6745_dp
       mad_dp(1)       = .true.
       mad_dp(n)       = .true.
       mad_dp(2:n-1)   = (d(1:n-2) .ge. (med-thresh)) .and. (d(1:n-2) .le. (med+thresh)) .and. dmask(1:n-2)
    case default
       stop 'Unimplemented option in mad_dp'
    end select

  END FUNCTION mad_dp


  FUNCTION mad_sp(arr, z, mask, deriv)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arr
    REAL(sp),                  OPTIONAL, INTENT(IN) :: z
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: deriv
    LOGICAL,     DIMENSION(size(arr))               :: mad_sp

    LOGICAL,  DIMENSION(size(arr)) :: maske
    REAL(sp), DIMENSION(size(arr)) :: d
    LOGICAL,  DIMENSION(size(arr)) :: dmask
    INTEGER(i4) :: n, ideriv
    !INTEGER(i4) :: m
    REAL(sp)    :: iz, med, mabsdev, thresh

    n = size(arr)
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= n) stop 'Error mad_sp: size(mask) /= size(arr)'
       maske = mask
    endif
    if (present(z)) then
       iz = z
    else
       iz = 7.0_sp
    endif
    if (present(deriv)) then
       ideriv = deriv
    else
       ideriv = 0
    endif

    select case(ideriv)
    case(0) ! pure values
       !m       = count(maske)
       med     = median(arr,mask=maske)
       mabsdev = median(abs(arr-med),mask=maske)
       thresh  = mabsdev * iz/0.6745_sp
       mad_sp     = (arr .ge. (med-thresh)) .and. (arr .le. (med+thresh)) .and. maske
    case(1) ! 1st derivative
       ! d(1:n-2) = 0.5_sp* (arr(3:n) - arr(1:n-2)) ! does not work -> ask Clemens & Matthias M
       d(1:n-1)     = arr(2:n) - arr(1:n-1)
       dmask(1:n-1) = maske(2:n) .and. maske(1:n-1)
       !m            = count(dmask(1:n-1))
       med          = median(d(1:n-1),mask=dmask(1:n-1))
       mabsdev      = median(abs(d(1:n-1)-med),mask=dmask(1:n-1))
       thresh       = mabsdev * iz/0.6745_sp
       mad_sp(n)       = .true.
       mad_sp(1:n-1)   = (d(1:n-1) .ge. (med-thresh)) .and. (d(1:n-1) .le. (med+thresh)) .and. dmask(1:n-1)
    case(2) ! -2nd derivative
       d(1:n-2)     = arr(2:n-1) + arr(2:n-1) - arr(1:n-2) - arr(3:n)
       dmask(1:n-2) = maske(2:n-1) .and. maske(1:n-2) .and. maske(3:n)
       !m            = count(dmask(1:n-2))
       med          = median(d(1:n-2),mask=dmask(1:n-2))
       mabsdev      = median(abs(d(1:n-2)-med),mask=dmask(1:n-2))
       thresh       = mabsdev * iz/0.6745_sp
       mad_sp(1)       = .true.
       mad_sp(n)       = .true.
       mad_sp(2:n-1)   = (d(1:n-2) .ge. (med-thresh)) .and. (d(1:n-2) .le. (med+thresh)) .and. dmask(1:n-2)
    case default
       stop 'Unimplemented option in mad_sp'
    end select

  END FUNCTION mad_sp

  ! ------------------------------------------------------------------

  FUNCTION mad_val_dp(arr, z, mask, tout, mval)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arr
    REAL(dp),                  OPTIONAL, INTENT(IN) :: z, mval
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp),    DIMENSION(size(arr))               :: mad_val_dp
    CHARACTER(1)                                    :: tout ! type out
    ! u : cut upper; l : cut lower; b : cut upper and lower

    LOGICAL,  DIMENSION(size(arr)) :: maske
    INTEGER(i4) :: n
    REAL(dp)    :: iz, med, mabsdev, thresh

    n = size(arr)
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= n) stop 'Error mad_val_dp: size(mask) /= size(arr)'
       maske = mask
    endif
    if (present(z)) then
       iz = z
    else
       iz = 7.0_dp
    endif

    if (present(mval)) then
       where (abs(arr - mval) .lt. tiny(1._dp) ) maske = .false.
       ! reset if no values remain
       if (.not. any(maske)) then
          where ( abs(arr - mval) .lt. tiny(1._dp) ) maske = .true.
       end if
    endif

    med     = median(arr,mask=maske)
    mabsdev = median(abs(arr-med),mask=maske)
    thresh  = mabsdev * iz/0.6745_dp
    mad_val_dp = arr

    select case(tout)
    case("u")
      ! print *, "The threshold is set to", med, "+", thresh
      where ((mad_val_dp .gt. (med+thresh)) &
           .and. maske) mad_val_dp = med+thresh
    case("l")
      ! print *, "The threshold is set to", med, "-", thresh
      where ((mad_val_dp .lt. (med-thresh)) &
           .and. maske) mad_val_dp = med-thresh
    case("b")
      ! print *, "The threshold is set to", med, "+/-", thresh
      where ((mad_val_dp .gt. (med+thresh)) &
           .and. maske) mad_val_dp = med+thresh
      where ((mad_val_dp .lt. (med-thresh)) &
           .and. maske) mad_val_dp = med-thresh
    case default
       stop 'Unimplemented option in mad_val_dp'
    end select

  END FUNCTION mad_val_dp

  ! ------------------------------------------------------------------

  FUNCTION mad_val_sp(arr, z, mask, tout, mval)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arr
    REAL(sp),                  OPTIONAL, INTENT(IN) :: z, mval
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp),    DIMENSION(size(arr))               :: mad_val_sp
    CHARACTER(1)                                    :: tout ! type out
    ! u : cut upper; l : cut lower; b : cut upper and lower

    LOGICAL,  DIMENSION(size(arr)) :: maske
    INTEGER(i4) :: n
    REAL(sp)    :: iz, med, mabsdev, thresh

    n = size(arr)
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= n) stop 'Error mad_val_sp: size(mask) /= size(arr)'
       maske = mask
    endif
    if (present(z)) then
       iz = z
    else
       iz = 7.0_sp
    endif

    if (present(mval)) then
       where (abs(arr - mval) .lt. tiny(1._sp)) maske = .false.
       ! reset if no values remain
       if (.not. any(maske)) then
          where ( abs(arr - mval) .lt. tiny(1._dp) ) maske = .true.
       end if
    endif

    med     = median(arr,mask=maske)
    mabsdev = median(abs(arr-med),mask=maske)
    thresh  = mabsdev * iz/0.6745_sp
    mad_val_sp = arr
    select case(tout)
    case("u")
      print *, "The threshold is set to", med, "+", thresh
      where ((mad_val_sp .gt. (med+thresh)) &
           .and. maske) mad_val_sp = med+thresh
    case("l")
      print *, "The threshold is set to", med, "-", thresh
      where ((mad_val_sp .lt. (med-thresh)) &
           .and. maske) mad_val_sp = med-thresh
    case("b")
      print *, "The threshold is set to", med, "+/-", thresh
      where ((mad_val_sp .gt. (med+thresh)) &
           .and. maske) mad_val_sp = med+thresh
      where ((mad_val_sp .lt. (med-thresh)) &
           .and. maske) mad_val_sp = med-thresh
    case default
       stop 'Unimplemented option in mad_val_sp'
    end select

  END FUNCTION mad_val_sp

  ! ------------------------------------------------------------------


END MODULE mo_mad
