MODULE mo_boxcox

  ! This module contains routines to calculate the Box-Cox transformation
  ! as well as estimating the best exponent for the Box-Cox transformation

  ! Usage:
  ! USE mo_boxcox, ONLY: boxcox, invboxcox
  ! new_data = boxcox(data, lmbda)
  ! data     = invboxcox(new_data, lmbda)

  ! Written March 2011, Matthias Cuntz
  !   - modified Python code of Travis Oliphant (2002): boxcox, llf_boxcox

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

  ! Copyright 2011-2012 Matthias Cuntz, Juliane Mai

  ! History:
  !   Dec 2019: Robert Schweppe: - removed NR code (get_boxcox)

  USE mo_kind, ONLY : i4, sp, dp
  USE mo_utils, only : eq, le, ge

  IMPLICIT NONE

  PUBLIC :: boxcox     ! Calculate Box-Cox power transformed values given the original values and the exponent lambda
  PUBLIC :: invboxcox  ! Calculate the inverse Box-Cox given the transformed values and the exponent lambda

  ! ------------------------------------------------------------------

  !     NAME
  !         boxcox

  !     PURPOSE
  !         Transform a positive dataset with a Box-Cox power transformation.
  !             boxcox(x) = (x**lambda - 1)/lambda    if lambda <> 0
  !                         ln(x)                     if lambda = 0
  !
  !         If an optinal mask is given, then the Box-Cox transformation is only performed on
  !         those locations that correspond to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = boxcox(x, lmbda, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp) :: x(:)       1D-array with input numbers (>0.)
  !         real(sp/dp) :: lmbda      Exponent power of Box-Cox transform (>= 0.)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: boxcox     Power transformed values (at mask=.True.)

  !     INTENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(x).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be positive, i.e. > 0.

  !     EXAMPLE
  !         lmbda    = get_boxcox(data, mask=(data>0.))
  !         new_data = boxcox(data, lmbda, mask=(data>0.))
  !         idata    = invboxcox(new_data, lmbda, mask=(data>0.))
  !         -> see also example in test directory

  !     HISTORY
  !         Written,  Matthias Cuntz, March 2011
  !            - Modified Python code of Travis Oliphant (2002): boxcox, llf_boxcox, get_boxcox
  !            - Modified numerical recipes: brent, mnbrak, swap, shft
  INTERFACE boxcox
    MODULE PROCEDURE boxcox_sp, boxcox_dp
  END INTERFACE boxcox

  ! ------------------------------------------------------------------

  !     NAME
  !         invboxcox

  !     PURPOSE
  !         Back-transformation of Box-Cox-transformed data.
  !             boxcox(x)    = (x**lambda - 1)/lambda        if lambda <> 0
  !                            ln(x)                         if lambda = 0
  !             invboxcox(y) = (lambda*y + 1)**(1/lambda)    if lambda <> 0
  !                            exp(y)                        if lambda = 0
  !
  !         If an optinal mask is given, then the inverse Box-Cox transformation is only performed on
  !         those locations that correspond to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.
  !         x can be scalar or vector

  !     CALLING SEQUENCE
  !         out = invboxcox(x, lmbda)                   ! scalar x
  !         out = invboxcox(x, lmbda, mask=mask)        ! vector x

  !     INTENT(IN)
  !         real(sp/dp) :: x / x(:)   scalar/1D-array with input numbers (>0.)
  !         real(sp/dp) :: lmbda      Exponent power of Box-Cox transform (>= 0.)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: boxcox     Back transformed values (at mask=.True.)

  !     INTENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(x).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.
  !                                   Only applicable if x is a 1D-array

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         lmbda    = get_boxcox(data, mask=(data>0.))
  !         new_data = boxcox(data, lmbda, mask=(data>0.))
  !         idata    = invboxcox(new_data, lmbda, mask=(data>0.))
  !         -> see also example in test directory

  !     HISTORY
  !         Written,  Matthias Cuntz, March 2011
  !            - Modified MC: Python code of Travis Oliphant (2002): boxcox, llf_boxcox, get_boxcox
  !            - Modified MC: numerical recipes: brent, mnbrak, swap, shft
  !            - Modified JM: scalar version of invboxcox
  INTERFACE invboxcox
     MODULE PROCEDURE invboxcox_0d_sp, invboxcox_0d_dp, invboxcox_1d_sp, invboxcox_1d_dp
  END INTERFACE invboxcox

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION boxcox_sp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(in) :: x
    REAL(sp), INTENT(in) :: lmbda
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(sp), DIMENSION(size(x)) :: boxcox_sp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(sp) :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error boxcox_sp: size(mask) /= size(x)'
      maske = mask
    endif
    if (any((le(x, 0.0_sp)) .and. maske)) stop 'Error boxcox_sp: x <= 0'
    if (abs(lmbda) < tiny(0.0_sp)) then
      where (maske)
        boxcox_sp = log(x)
      elsewhere
        boxcox_sp = x
      end where
    else
      lmbda1 = 1.0_sp / lmbda
      boxcox_sp = merge((exp(lmbda * log(x)) - 1.0_sp) * lmbda1, x, maske)
    endif

  END FUNCTION boxcox_sp


  FUNCTION boxcox_dp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in) :: x
    REAL(dp), INTENT(in) :: lmbda
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(dp), DIMENSION(size(x)) :: boxcox_dp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(dp) :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error boxcox_dp: size(mask) /= size(x)'
      maske = mask
    endif
    if (any((le(x, 0.0_dp)) .and. maske)) then
      stop 'Error boxcox_dp: x <= 0'
    end if
    if (abs(lmbda) < tiny(0.0_dp)) then
      where (maske)
        boxcox_dp = log(x)
      elsewhere
        boxcox_dp = x
      end where
    else
      lmbda1 = 1.0_dp / lmbda
      boxcox_dp = merge((exp(lmbda * log(x)) - 1.0_dp) * lmbda1, x, maske)
    endif

  END FUNCTION boxcox_dp

  ! ------------------------------------------------------------------

  FUNCTION invboxcox_1d_sp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(in) :: x
    REAL(sp), INTENT(in) :: lmbda
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(sp), DIMENSION(size(x)) :: invboxcox_1d_sp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(sp) :: lmbda1
    REAL(sp), DIMENSION(size(x)) :: temp

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error invboxcox_1d_sp: size(mask) /= size(x)'
      maske = mask
    endif
    if (abs(lmbda) < tiny(0.0_sp)) then
      where (maske)
        invboxcox_1d_sp = exp(x)
      elsewhere
        invboxcox_1d_sp = x
      end where
    else
      lmbda1 = 1.0_sp / lmbda
      temp = lmbda * x + 1.0_sp
      where (temp > 0.0_sp)
        temp = exp(lmbda1 * log(temp))
      elsewhere
        temp = 0.0_sp
      end where
      invboxcox_1d_sp = merge(temp, x, maske)
    endif

  END FUNCTION invboxcox_1d_sp

  FUNCTION invboxcox_1d_dp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in) :: x
    REAL(dp), INTENT(in) :: lmbda
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(dp), DIMENSION(size(x)) :: invboxcox_1d_dp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(dp) :: lmbda1
    REAL(dp), DIMENSION(size(x)) :: temp

    maske(:) = .true.
    if (present(mask)) then
      if (size(mask) /= size(x)) stop 'Error invboxcox_1d_dp: size(mask) /= size(x)'
      maske = mask
    endif
    if (abs(lmbda) < tiny(0.0_dp)) then
      where (maske)
        invboxcox_1d_dp = exp(x)
      elsewhere
        invboxcox_1d_dp = x
      end where
    else
      lmbda1 = 1.0_dp / lmbda
      temp = lmbda * x + 1.0_dp
      where (temp > 0.0_dp)
        temp = exp(lmbda1 * log(temp))
      elsewhere
        temp = 0.0_dp
      end where
      invboxcox_1d_dp = merge(temp, x, maske)
    endif

  END FUNCTION invboxcox_1d_dp

  FUNCTION invboxcox_0d_sp(x, lmbda)

    IMPLICIT NONE

    REAL(sp), INTENT(in) :: x
    REAL(sp), INTENT(in) :: lmbda
    REAL(sp) :: invboxcox_0d_sp

    REAL(sp) :: lmbda1
    REAL(sp) :: temp

    if (abs(lmbda) < tiny(0.0_sp)) then
      invboxcox_0d_sp = exp(x)
    else
      lmbda1 = 1.0_sp / lmbda
      temp = lmbda * x + 1.0_sp
      if (temp > 0.0_sp) then
        temp = exp(lmbda1 * log(temp))
      else
        temp = 0.0_sp
      end if
      invboxcox_0d_sp = temp
    endif

  END FUNCTION invboxcox_0d_sp

  FUNCTION invboxcox_0d_dp(x, lmbda)

    IMPLICIT NONE

    REAL(dp), INTENT(in) :: x
    REAL(dp), INTENT(in) :: lmbda
    REAL(dp) :: invboxcox_0d_dp

    REAL(dp) :: lmbda1
    REAL(dp) :: temp

    if (abs(lmbda) < tiny(0.0_dp)) then
      invboxcox_0d_dp = exp(x)
    else
      lmbda1 = 1.0_dp / lmbda
      temp = lmbda * x + 1.0_dp
      if (temp > 0.0_dp) then
        temp = exp(lmbda1 * log(temp))
      else
        temp = 0.0_dp
      end if
      invboxcox_0d_dp = temp
    endif

  END FUNCTION invboxcox_0d_dp

END MODULE mo_boxcox
