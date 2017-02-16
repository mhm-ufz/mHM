MODULE mo_percentile

  ! This module provides routines for median and percentiles.

  ! Written   Matthias Cuntz,              Mar 2011
  ! Modified, Juliane Mai,                 Jul 2012 - different interpolation schemes in percentiles
  !           Matthias Cuntz, Juliane Mai, Jul 2012 - uses previous of ksmallest to half execution time
  !           Matthias Cuntz,              May 2014 - removed numerical recipes

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Matthias Cuntz or Juliane Mai.

  ! Copyright 2011-2014 Matthias Cuntz, Juliane Mai, Stephan Thober

  USE mo_kind, ONLY: i4, sp, dp

  Implicit NONE

  PUBLIC :: median          ! Median
  PUBLIC :: n_element       ! The n-th smallest value in an array
  PUBLIC :: percentile      ! The value below which a certain percent of the input fall
  PUBLIC :: qmedian         ! Quick median calculation, rearranges input

  ! ------------------------------------------------------------------

  !     NAME
  !         median

  !     PURPOSE
  !         Returns the median of the values in an array.
  !         If size is even, then the mean of the size/2 and size/2+1 element is the median.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.

  !     CALLING SEQUENCE
  !         out = median(vec,mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out        median of values in input array

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 5.5
  !         out = median(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar 2011
  !         Modified, Matthias Cuntz, Juliane Mai, Jul 2012 - uses previous of ksmallest to half execution time
  INTERFACE median
     MODULE PROCEDURE median_sp, median_dp
  END INTERFACE median

  ! ------------------------------------------------------------------

  !     NAME
  !         n_element

  !     PURPOSE
  !         Returns the n-th smallest value in an array.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.

  !     CALLING SEQUENCE
  !         out = n_element(vec, n, mask=mask, before=before, previous=previous, after=after, next=next)

  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out        n-th smallest value in input array

  !     INTENT(IN), OPTIONAL
  !         integer(i4) :: n          index of sorted array
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         real(sp/dp) :: before      (n-1)-th smallest value in input array, e.g. for median/percentile calculations
  !         real(sp/dp) :: previous    same as before
  !         real(sp/dp) :: after       (n+1)-th smallest value in input array
  !         real(sp/dp) :: next        same as after

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 4
  !         out = n_element(vec,4)
  !         -> see also example in test directory

  !     LITERATURE
  !         Niklaus Wirth. "Algorithms and Data Structures". Prentice-Hall, Inc., 1985. ISBN 0-13-022005-1.

  !     HISTORY
  !         Written,  Matthias Cuntz, May 2014 - based on qmedian
  INTERFACE n_element
     MODULE PROCEDURE n_element_dp, n_element_sp
  END INTERFACE n_element

  ! ------------------------------------------------------------------

  !     NAME
  !         percentile

  !     PURPOSE
  !         Returns the value below which a certain percent of array values fall.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.
  !
  !         Different definitions can be applied to interpolate the stepwise CDF of the given data.
  !         (1) Inverse empirical CDF (no interpolation, default MATHEMATICA)
  !         (2) Linear interpolation (California method)
  !         (3) Element numbered closest
  !         (4) Linear interpolation (hydrologist method)
  !         (5) Mean-based estimate (Weibull method, default IMSL)
  !         (6) Mode-based estimate
  !         (7) Median-based estimate
  !         (8) normal distribution estimate
  !
  !         See: http://reference.wolfram.com/mathematica/tutorial/BasicStatistics.html

  !     CALLING SEQUENCE
  !         out = percentile(vec,k,mask=mask,mode_in=mode)

  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers
  !         real(sp/dp) :: k[(:)]     Percentage of percentile, can be 1 dimensional

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out[(size(k))]   k-th percentile of values in input array, can be
  !                                         1 dimensional corresponding to k

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.
  !         integer(i4) :: mode_in    Specifies the interpolation scheme applied.
  !                                   Default:
  !                                       Inverse empirical CDF (no interpolation, default Mathematica)
  !                                       mode_in = 1_i4

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 10.
  !         out = percentile(vec,95.)
  !         ! Returns (10.,8)
  !         out = percentile(vec,(/95.,80./))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar  2011
  !         Modified, Stephan Thober, Dec  2011 - added 1 dimensional version
  !                   Juliane Mai,    July 2012 - different interpolation schemes
  !         Modified, Matthias Cuntz, Juliane Mai, Jul 2012 - uses previous of ksmallest to half execution time
  INTERFACE percentile
     MODULE PROCEDURE percentile_0d_sp, percentile_0d_dp, percentile_1d_sp, percentile_1d_dp
  END INTERFACE percentile

  ! ------------------------------------------------------------------

  !     NAME
  !         qmedian

  !     PURPOSE
  !         Quick calculation of the median thereby rearranging the input array.

  !     CALLING SEQUENCE
  !         out = qmedian(vec)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers.
  !                                   Wil be rearranged on output.

  !     INTENT(OUT)
  !         real(sp/dp) :: out        median of values in input array

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 5.5
  !         out = qmedian(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         Niklaus Wirth. "Algorithms and Data Structures". Prentice-Hall, Inc., 1985. ISBN 0-13-022005-1.

  !     HISTORY
  !         Written, Filip Hroch as part of Munipack: http://munipack.physics.muni.cz
  !         Modified, Matthias Cuntz, Jul 2012 - function, k=n/2+1
  !         Modified, Matthias Cuntz, Juliane Mai, Jul 2012 - real median for even n
  INTERFACE qmedian
     MODULE PROCEDURE qmedian_sp, qmedian_dp
  END INTERFACE qmedian

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION median_dp(arrin, mask)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp)                                        :: median_dp

    INTEGER(i4) :: n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: arr
    REAL(dp) :: tmp

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)

       if (n < 2) stop 'median_dp: n < 2'

       if (mod(n,2) == 0) then ! Even
          median_dp = n_element(arr,n/2+1,previous=tmp)
          median_dp = 0.5_dp*(median_dp+tmp)
       else ! Odd
          median_dp = n_element(arr,(n+1)/2)
       endif

       deallocate(arr)
    else
       n = size(arrin)
       if (n < 2) stop 'median_dp: n < 2'

       if (mod(n,2) == 0) then ! Even
          median_dp = n_element(arrin,n/2+1,previous=tmp)
          median_dp = 0.5_dp*(median_dp+tmp)
       else ! Odd
          median_dp = n_element(arrin,(n+1)/2)
       endif
    endif

  END FUNCTION median_dp


  FUNCTION median_sp(arrin, mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp)                                        :: median_sp

    INTEGER(i4) :: n
    REAL(sp), DIMENSION(:), ALLOCATABLE :: arr
    REAL(sp) :: tmp

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)

       if (n < 2) stop 'median_sp: n < 2'

       if (mod(n,2) == 0) then ! Even
          median_sp = n_element(arr,n/2+1,previous=tmp)
          median_sp = 0.5_sp*(median_sp+tmp)
       else ! Odd
          median_sp = n_element(arr,(n+1)/2)
       endif

       deallocate(arr)
    else
       n = size(arrin)
       if (n < 2) stop 'median_sp: n < 2'

       if (mod(n,2) == 0) then ! Even
          median_sp = n_element(arrin,n/2+1,previous=tmp)
          median_sp = 0.5_sp*(median_sp+tmp)
       else ! Odd
          median_sp = n_element(arrin,(n+1)/2)
       endif
    endif

  END FUNCTION median_sp

  ! ------------------------------------------------------------------

  function n_element_dp(idat, n, mask, before, after, previous, next)

    IMPLICIT NONE

    real(dp),    dimension(:), intent(in)            :: idat
    integer(i4),               intent(in)            :: n
    logical,     dimension(:), optional, intent(in)  :: mask
    real(dp),                  optional, intent(out) :: before
    real(dp),                  optional, intent(out) :: after
    real(dp),                  optional, intent(out) :: previous
    real(dp),                  optional, intent(out) :: next
    real(dp)                                         :: n_element_dp

    real(dp), dimension(:), allocatable :: dat
    real(dp) :: w
    integer(i4) :: nn, k
    integer(i4) :: l,r,i,j

    if (present(mask)) then
       nn = count(mask)
       allocate(dat(nn))
       dat = pack(idat,mask)
    else
       nn = size(idat)
       allocate(dat(nn))
       dat = idat
    endif

    if (n < 1)  stop 'n_element_dp: n < 1'
    if (n > nn) stop 'n_element_dp: n > size(pack(dat,mask))'

    dat = idat
    nn = size(dat)
    k = n !nn/2 + 1
    l = 1
    r = nn
    do while( l < r )
       n_element_dp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < n_element_dp )
             i = i + 1
          enddo
          do while( n_element_dp < dat(j) )
             j = j - 1
          enddo
          if ( i <= j ) then
             w      = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1
             j = j - 1
          endif
          if ( i > j ) exit
       enddo
       if ( j < k ) l = i
       if ( k < i ) r = j
    enddo
    ! if (mod(nn,2) == 0) then
    !    n_element_dp = 0.5_dp*(dat(k) + maxval(dat(:k-1)))
    ! else
    !    n_element_dp = dat(k)
    ! endif
    n_element_dp = dat(k)

    if (present(before))   before   = maxval(dat(:k-1))
    if (present(previous)) previous = maxval(dat(:k-1))
    if (present(after))    after    = minval(dat(k+1:))
    if (present(next))     next     = minval(dat(k+1:))

    deallocate(dat)

  end function n_element_dp

  function n_element_sp(idat, n, mask, before, after, previous, next)

    IMPLICIT NONE

    real(sp),    dimension(:), intent(in)            :: idat
    integer(i4),               intent(in)            :: n
    logical,     dimension(:), optional, intent(in)  :: mask
    real(sp),                  optional, intent(out) :: before
    real(sp),                  optional, intent(out) :: after
    real(sp),                  optional, intent(out) :: previous
    real(sp),                  optional, intent(out) :: next
    real(sp)                                         :: n_element_sp

    real(sp), dimension(:), allocatable :: dat
    real(sp) :: w
    integer(i4) :: nn, k
    integer(i4) :: l,r,i,j

    if (present(mask)) then
       nn = count(mask)
       allocate(dat(nn))
       dat = pack(idat,mask)
    else
       nn = size(idat)
       allocate(dat(nn))
       dat = idat
    endif

    if (n < 1)  stop 'n_element_sp: n < 1'
    if (n > nn) stop 'n_element_sp: n > size(pack(dat,mask))'

    dat = idat
    nn = size(dat)
    k = n !nn/2 + 1
    l = 1
    r = nn
    do while( l < r )
       n_element_sp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < n_element_sp )
             i = i + 1
          enddo
          do while( n_element_sp < dat(j) )
             j = j - 1
          enddo
          if ( i <= j ) then
             w      = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1
             j = j - 1
          endif
          if ( i > j ) exit
       enddo
       if ( j < k ) l = i
       if ( k < i ) r = j
    enddo
    ! if (mod(nn,2) == 0) then
    !    n_element_sp = 0.5_sp*(dat(k) + maxval(dat(:k-1)))
    ! else
    !    n_element_sp = dat(k)
    ! endif
    n_element_sp = dat(k)

    if (present(before))   before   = maxval(dat(:k-1))
    if (present(previous)) previous = maxval(dat(:k-1))
    if (present(after))    after    = minval(dat(k+1:))
    if (present(next))     next     = minval(dat(k+1:))

    deallocate(dat)

  end function n_element_sp

  ! ------------------------------------------------------------------

  FUNCTION percentile_0d_dp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in
    REAL(dp)                                        :: percentile_0d_dp

    INTEGER(i4)                         :: n, nn1, nn2
    INTEGER(i4)                         :: mode
    REAL(dp)                            :: kk, ks1, ks2
    REAL(dp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    if (n < 2) stop 'percentile_0d_dp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
    case(1_i4)
       kk = k/100._dp*real(n,dp)
       nn1 = min(n, max(1_i4,ceiling(kk,kind=i4)))
       nn2 = nn1

       ! Linear interpolation (California method)
    case(2_i4)
       kk  = k/100._dp*real(n,dp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Element numbered closest
    case(3_i4)
       kk = 0.5_dp+k/100._dp*real(n,dp)
       nn1 = min(n, max(1_i4,floor(kk,kind=i4)))
       nn2 = nn1

       ! Linear interpolation (hydrologist method)
    case(4_i4)
       kk  = 0.5_dp+k/100._dp*(real(n,dp))
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
    case(5_i4)
       kk  = k/100._dp*(real(n,dp)+1._dp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mode-based estimate
    case(6_i4)
       kk  = 1.0_dp+k/100._dp*(real(n,dp)-1._dp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Median-based estimate
    case(7_i4)
       kk  = 1.0_dp/3.0_dp+k/100._dp*(real(n,dp)+1.0_dp/3.0_dp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Normal distribution estimate
    case(8_i4)
       kk  = 3.0_dp/8.0_dp+k/100._dp*(real(n,dp)+1.0_dp/4.0_dp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! No valid mode
    case default
       stop 'percentile_0d_dp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_dp = n_element(arr,nn1)
       else
          ! interpolation
          ks2 = n_element(arr,nn2,previous=ks1)
          percentile_0d_dp = ks1 + (ks2-ks1)*(kk-real(nn1,dp))
       end if
       deallocate(arr)
    else
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_dp = n_element(arrin,nn1)
       else
          ! interpolation
          ks2 = n_element(arrin,nn2,previous=ks1)
          percentile_0d_dp = ks1 + (ks2-ks1)*(kk-real(nn1,dp))
       end if
    endif

  END FUNCTION percentile_0d_dp


  FUNCTION percentile_0d_sp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in
    REAL(sp)                                        :: percentile_0d_sp

    INTEGER(i4)                         :: n, nn1, nn2
    INTEGER(i4)                         :: mode
    REAL(sp)                            :: kk, ks1, ks2
    REAL(sp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    if (n < 2) stop 'percentile_0d_sp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
    case(1_i4)
       kk = k/100._sp*real(n,sp)
       nn1 = min(n, max(1_i4,ceiling(kk,kind=i4)))
       nn2 = nn1

       ! Linear interpolation (California method)
    case(2_i4)
       kk  = k/100._sp*real(n,sp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Element numbered closest
    case(3_i4)
       kk = 0.5_sp+k/100._sp*real(n,sp)
       nn1 = min(n, max(1_i4,floor(kk,kind=i4)))
       nn2 = nn1

       ! Linear interpolation (hydrologist method)
    case(4_i4)
       kk  = 0.5_sp+k/100._sp*(real(n,sp))
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
    case(5_i4)
       kk  = k/100._sp*(real(n,sp)+1._sp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mode-based estimate
    case(6_i4)
       kk  = 1.0_sp+k/100._sp*(real(n,sp)-1._sp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Median-based estimate
    case(7_i4)
       kk  = 1.0_sp/3.0_sp+k/100._sp*(real(n,sp)+1.0_sp/3.0_sp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Normal distribution estimate
    case(8_i4)
       kk  = 3.0_sp/8.0_sp+k/100._sp*(real(n,sp)+1.0_sp/4.0_sp)
       nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
       nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! No valid mode
    case default
       stop 'percentile_0d_sp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_sp = n_element(arr,nn1)
       else
          ! interpolation
          ks2 = n_element(arr,nn2,previous=ks1)
          percentile_0d_sp = ks1 + (ks2-ks1)*(kk-real(nn1,sp))
       end if
       deallocate(arr)
    else
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_sp = n_element(arrin,nn1)
       else
          ! interpolation
          ks2 = n_element(arrin,nn2,previous=ks1)
          percentile_0d_sp = ks1 + (ks2-ks1)*(kk-real(nn1,sp))
       end if
    endif

  END FUNCTION percentile_0d_sp


  function percentile_1d_dp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),    DIMENSION(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in

    REAL(dp),    DIMENSION(size(k))                 :: percentile_1d_dp

    INTEGER(i4)                            :: i, n
    INTEGER(i4)                            :: mode
    INTEGER(i4), DIMENSION(size(k))        :: nn1, nn2
    REAL(dp),    DIMENSION(size(k))        :: kk
    REAL(dp)                               :: ks1, ks2
    REAL(dp),    DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_dp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_dp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
    case(1_i4)
       kk(:) = k(:)/100._dp*real(n,dp)
       nn1(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))
       nn2 = nn1

       ! Linear interpolation (California method)
    case(2_i4)
       kk(:)  = k(:)/100._dp*real(n,dp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Element numbered closest
    case(3_i4)
       kk(:) = 0.5_dp+k(:)/100._dp*real(n,dp)
       nn1(:) = min(n, max(1_i4,floor(kk(:),kind=i4)))
       nn2 = nn1

       ! Linear interpolation (hydrologist method)
    case(4_i4)
       kk(:)  = 0.5_dp+k(:)/100._dp*(real(n,dp))
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
    case(5_i4)
       kk(:)  = k(:)/100._dp*(real(n,dp)+1._dp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mode-based estimate
    case(6_i4)
       kk(:)  = 1.0_dp+k(:)/100._dp*(real(n,dp)-1._dp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Median-based estimate
    case(7_i4)
       kk(:)  = 1.0_dp/3.0_dp+k(:)/100._dp*(real(n,dp)+1.0_dp/3.0_dp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Normal distribution estimate
    case(8_i4)
       kk(:)  = 3.0_dp/8.0_dp+k(:)/100._dp*(real(n,dp)+1.0_dp/4.0_dp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! No valid mode
    case default
       stop 'percentile_1d_dp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_dp(i) = n_element(arr,nn1(i))
          else
             ! interpolation
             ks2 = n_element(arr,nn2(i),previous=ks1)
             percentile_1d_dp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),dp))
          end if
       end do
       deallocate(arr)
    else
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_dp(i) = n_element(arrin,nn1(i))
          else
             ! interpolation
             ks2 = n_element(arrin,nn2(i),previous=ks1)
             percentile_1d_dp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),dp))
          end if
       end do
    endif

  END function percentile_1d_dp


  function percentile_1d_sp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),    DIMENSION(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in

    REAL(sp),    DIMENSION(size(k))                 :: percentile_1d_sp

    INTEGER(i4)                            :: i, n
    INTEGER(i4)                            :: mode
    INTEGER(i4), DIMENSION(size(k))        :: nn1, nn2
    REAL(sp),    DIMENSION(size(k))        :: kk
    REAL(sp)                               :: ks1, ks2
    REAL(sp),    DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_sp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_sp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
    case(1_i4)
       kk(:) = k(:)/100._sp*real(n,sp)
       nn1(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))
       nn2 = nn1

       ! Linear interpolation (California method)
    case(2_i4)
       kk(:)  = k(:)/100._sp*real(n,sp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Element numbered closest
    case(3_i4)
       kk(:) = 0.5_sp+k(:)/100._sp*real(n,sp)
       nn1(:) = min(n, max(1_i4,floor(kk(:),kind=i4)))
       nn2 = nn1

       ! Linear interpolation (hydrologist method)
    case(4_i4)
       kk(:)  = 0.5_sp+k(:)/100._sp*(real(n,sp))
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
    case(5_i4)
       kk(:)  = k(:)/100._sp*(real(n,sp)+1._sp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mode-based estimate
    case(6_i4)
       kk(:)  = 1.0_sp+k(:)/100._sp*(real(n,sp)-1._sp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Median-based estimate
    case(7_i4)
       kk(:)  = 1.0_sp/3.0_sp+k(:)/100._sp*(real(n,sp)+1.0_sp/3.0_sp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Normal distribution estimate
    case(8_i4)
       kk(:)  = 3.0_sp/8.0_sp+k(:)/100._sp*(real(n,sp)+1.0_sp/4.0_sp)
       nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
       nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! No valid mode
    case default
       stop 'percentile_1d_sp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_sp(i) = n_element(arr,nn1(i))
          else
             ! interpolation
             ks2 = n_element(arr,nn2(i),previous=ks1)
             percentile_1d_sp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),sp))
          end if
       end do
       deallocate(arr)
    else
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_sp(i) = n_element(arrin,nn1(i))
          else
             ! interpolation
             ks2 = n_element(arrin,nn2(i),previous=ks1)
             percentile_1d_sp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),sp))
          end if
       end do
    endif

  END function percentile_1d_sp

  ! ------------------------------------------------------------------

  function qmedian_dp(dat)

    IMPLICIT NONE

    real(dp), dimension(:), intent(inout) :: dat
    real(dp) :: qmedian_dp

    real(dp) :: w
    integer(i4) :: n,k
    integer(i4) :: l,r,i,j

    n = size(dat)
    k = n/2 + 1
    l = 1
    r = n
    do while( l < r )
       qmedian_dp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < qmedian_dp )
             i = i + 1
          enddo
          do while( qmedian_dp < dat(j) )
             j = j - 1
          enddo
          if ( i <= j ) then
             w      = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1
             j = j - 1
          endif
          if ( i > j ) exit
       enddo
       if ( j < k ) l = i
       if ( k < i ) r = j
    enddo
    if (mod(n,2) == 0) then
       qmedian_dp = 0.5_dp*(dat(k) + maxval(dat(:k-1)))
    else
       qmedian_dp = dat(k)
    endif

  end function qmedian_dp


  function qmedian_sp(dat)

    IMPLICIT NONE

    real(sp), dimension(:), intent(inout) :: dat
    real(sp) :: qmedian_sp

    real(sp) :: w
    integer(i4) :: n,k
    integer(i4) :: l,r,i,j

    n = size(dat)
    k = n/2 + 1
    l = 1
    r = n
    do while( l < r )
       qmedian_sp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < qmedian_sp )
             i = i + 1
          enddo
          do while( qmedian_sp < dat(j) )
             j = j - 1
          enddo
          if ( i <= j ) then
             w      = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1
             j = j - 1
          endif
          if ( i > j ) exit
       enddo
       if ( j < k ) l = i
       if ( k < i ) r = j
    enddo
    if (mod(n,2) == 0) then
       qmedian_sp = 0.5_sp*(dat(k) + maxval(dat(:k-1)))
    else
       qmedian_sp = dat(k)
    endif

  end function qmedian_sp

  ! ------------------------------------------------------------------

END MODULE mo_percentile
