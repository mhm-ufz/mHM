!> \file mo_kind.f90

!> \brief Define number representations

!> \details This module declares the desired ranges and precisions of the number representations,
!> such as single precision or double precision, 32-bit or 46-bit integer, etc.
!> It confirms mostly with the nrtype module of Numerical Recipes in f90.

!> \authors Juliane Mai, Matthias Cuntz, Nov 2011
!> \date Nov 2011
!> \copyright GNU Lesser General Public License http://www.gnu.org/licenses/

!  Number model from which the SELECTED_REAL_KIND are requested:
!                   4 byte REAL      8 byte REAL
!          IEEE:    precision =  6   precision =   15
!                   exponent  = 37   exponent  =  307
!          CRAY:        -            precision =   13
!                                    exponent  = 2465

! Written  Juliane Mai, Nov 2011
! Modified Matthias Cuntz, Nov 2011 - private/public
!                                   - documentation
!                                   - removed tab characters

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

! Copyright 2011 Matthias Cuntz, Juliane Mai

MODULE mo_kind

  IMPLICIT NONE

  !> 8 Byte Integer Kind
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18)
  !> 4 Byte Integer Kind
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)
  !> 2 Byte Integer Kind
  INTEGER, PARAMETER :: i2  = SELECTED_INT_KIND(4)
  !> 1 Byte Integer Kind
  INTEGER, PARAMETER :: i1  = SELECTED_INT_KIND(2)
  !> Single Precision Real Kind
  INTEGER, PARAMETER :: sp  = SELECTED_REAL_KIND(6,37)
  !> Double Precision Real Kind
  INTEGER, PARAMETER :: dp  = SELECTED_REAL_KIND(15,307)
  ! Note: Intrinsic KIND is an extension to the Fortran 90 standard.
  ! Matthias C thinks that spc=sp and dpc=dp would be enough.
  !> Single Precision Complex Kind
  INTEGER, PARAMETER :: spc = KIND((1.0_sp,1.0_sp))
  !> Double Precision Complex Kind
  INTEGER, PARAMETER :: dpc = KIND((1.0_dp,1.0_dp))
  !> Logical Kind
  INTEGER, PARAMETER :: lgt = KIND(.true.)

  ! Types have to be in a public section for doxygen

  !> Single Precision Numerical Recipes types for sparse arrays
  TYPE sprs2_sp
     INTEGER(I4) :: n, len
     REAL(SP),    DIMENSION(:), POINTER :: val
     INTEGER(I4), DIMENSION(:), POINTER :: irow
     INTEGER(I4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  !> Double Precision Numerical Recipes types for sparse arrays
  TYPE sprs2_dp
     INTEGER(I4) :: n, len
     REAL(DP),    DIMENSION(:), POINTER :: val
     INTEGER(I4), DIMENSION(:), POINTER :: irow
     INTEGER(I4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

END MODULE mo_kind
