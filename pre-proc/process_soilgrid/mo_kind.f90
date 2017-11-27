MODULE mo_kind

  ! This module declares the desired ranges and precisions of the number representations,
  ! such as single precision or double precision, 32-bit or 46-bit integer, etc.
  ! It confirms mostly with the nrtype module of Numerical Recipes in f90.

  ! Number model from which the SELECTED_REAL_KIND are requested:
  !                   4 byte REAL      8 byte REAL
  !          IEEE:    precision =  6   precision =   15
  !                   exponent  = 37   exponent  =  307
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465

  ! Written  Juliane Mai, Nov 2011
  ! Modified Matthias Cuntz, Nov 2011 - private/public
  !                                   - documentation
  !                                   - removed tab characters

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: i1, i2, i4, i8, sp, dp, spc, dpc, lgt, sprs2_sp, sprs2_dp

  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(18)       ! 8 Byte Integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)        ! 4 Byte Integer
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(4)        ! 2 Byte Integer
  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(2)        ! 1 Byte Integer
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)    ! Single Precision Real
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)  ! Double Precision Real
  ! Note: Intrinsic KIND is an extension to the Fortran 90 standard.
  ! Matthias C thinks that spc=sp and dpc=dp would be enough.
  INTEGER, PARAMETER :: spc = KIND((1.0_sp,1.0_sp))      ! Single Precision Complex
  INTEGER, PARAMETER :: dpc = KIND((1.0_dp,1.0_dp))      ! Double Precision Complex
  INTEGER, PARAMETER :: lgt = KIND(.true.)               ! Logical

  ! Numerical recipes types for sparse arrays
  TYPE sprs2_sp
     INTEGER(I4) :: n, len
     REAL(SP),    DIMENSION(:), POINTER :: val
     INTEGER(I4), DIMENSION(:), POINTER :: irow
     INTEGER(I4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4) :: n, len
     REAL(DP),    DIMENSION(:), POINTER :: val
     INTEGER(I4), DIMENSION(:), POINTER :: irow
     INTEGER(I4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

END MODULE mo_kind
