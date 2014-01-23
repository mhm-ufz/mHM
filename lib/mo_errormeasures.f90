MODULE mo_errormeasures

  ! This module contains routines for the masked calculation of
  ! error measures like MSE, RMSE, BIAS, SSE, NSE, ...

  ! Note: all except variance and standard deviation are population and not sample moments,
  !       i.e. they are normally divided by n and not (n-1)

  ! Written Aug 2012, Matthias Zink

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

  ! Copyright 2012 Matthias Zink

  USE mo_kind,   ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: BIAS                         ! bias
  PUBLIC :: LNNSE                        ! Logarithmic Nash Sutcliffe efficiency
  PUBLIC :: MAE                          ! Mean of absolute errors
  PUBLIC :: MSE                          ! Mean of squared errors
  PUBLIC :: NSE                          ! Nash Sutcliffe efficiency
  PUBLIC :: SSE                          ! Sum of squared errors
  PUBLIC :: SAE                          ! Sum of absolute errors
  PUBLIC :: RMSE                         ! Root mean squared error

  ! ------------------------------------------------------------------

  !     NAME
  !         BIAS

  !     PURPOSE
  !         Calculates the bias
  !             BIAS = mean(y) - mean(x)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = BIAS(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: BIAS        bias

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = BIAS(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE BIAS
     MODULE PROCEDURE BIAS_sp_1d, BIAS_dp_1d, BIAS_sp_2d, BIAS_dp_2d, BIAS_sp_3d, BIAS_dp_3d
  END INTERFACE BIAS

  ! ------------------------------------------------------------------

  !     NAME
  !         LNNSE

  !     PURPOSE
  !         Calculates the Logarithmic Nash Sutcliffe Efficiency
  !             LNNSE = sum((ln(y) - ln(x))**2) / sum( (ln(x) - ln(mean(x)))**2 )
  !         where x is the observation and y is the modelled data.
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         Note that the mask is intent inout, since values which are less or equal zero will be masked additionally.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = LNNSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: LNNSE         Logarithmic Nash Sutcliffe Efficiency

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.
  !         The mask will be updated if non-masked values are less equal zero.

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = LNNSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Juliane Mai, May 2013
  !         updated,  Rohin Kumar, May 2013  ! for mean of logQ
  INTERFACE LNNSE
     MODULE PROCEDURE LNNSE_sp_1d, LNNSE_dp_1d, LNNSE_dp_2d, LNNSE_sp_2d, LNNSE_sp_3d, LNNSE_dp_3d
  END INTERFACE LNNSE

  ! ------------------------------------------------------------------

  !     NAME
  !         MAE

  !     PURPOSE
  !         Calculates the mean absolute error
  !             MAE = sum(abs(y - x)) / count(mask)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = MAE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: MAE         Mean Absolute Error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = MAE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE MAE
     MODULE PROCEDURE MAE_sp_1d, MAE_dp_1d, MAE_sp_2d, MAE_dp_2d, MAE_sp_3d, MAE_dp_3d
  END INTERFACE MAE

  ! ------------------------------------------------------------------

  !     NAME
  !         MSE

  !     PURPOSE
  !         Calculates the mean squared error
  !             MSE = sum((y - x)**2) / count(mask)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = MSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: MSE         Mean squared error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = MSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE MSE
     MODULE PROCEDURE MSE_sp_1d, MSE_dp_1d, MSE_sp_2d, MSE_dp_2d, MSE_sp_3d, MSE_dp_3d
  END INTERFACE MSE

  ! ------------------------------------------------------------------

  !     NAME
  !         NSE

  !     PURPOSE
  !         Calculates the Nash Sutcliffe Efficiency
  !             NSE = sum((y - x)**2) / sum( (x - mean(x))**2)
  !         where x is the observation and y is the modelled data.
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = NSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: NSE         Nash Sutcliffe Efficiency

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = NSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         NASH, J., & SUTCLIFFE, J. (1970). River flow forecasting through conceptual models part I: A discussion of
  !                  principles. Journal of Hydrology, 10(3), 282-290. doi:10.1016/0022-1694(70)90255-6

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE NSE
     MODULE PROCEDURE NSE_sp_1d, NSE_dp_1d, NSE_dp_2d, NSE_sp_2d, NSE_sp_3d, NSE_dp_3d
  END INTERFACE NSE

  ! ------------------------------------------------------------------

  !     NAME
  !         SAE

  !     PURPOSE
  !         Calculates the sum of absolute errors
  !             SAE = sum(abs(y - x))
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = SAE(x, y, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: SAE         sum of absolute errors

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = SAE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         none

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE SAE
     MODULE PROCEDURE SAE_sp_1d, SAE_dp_1d, SAE_sp_2d, SAE_dp_2d, SAE_sp_3d, SAE_dp_3d
  END INTERFACE SAE

  ! ------------------------------------------------------------------

  !     NAME
  !         SSE

  !     PURPOSE
  !         Calculates the sum of squared errors
  !             SSE = sum((y - x)**2)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = SSE(x, y, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: SSE         sum of squared errors

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = SSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         none

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE SSE
     MODULE PROCEDURE SSE_sp_1d, SSE_dp_1d, SSE_sp_2d, SSE_dp_2d, SSE_sp_3d, SSE_dp_3d
  END INTERFACE SSE

  ! ------------------------------------------------------------------

  !     NAME
  !         RMSE

  !     PURPOSE
  !         Calculates the root-mean-square error
  !             RMSE = sqrt(sum((y - x)**2) / count(mask))
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = RMSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: RMSE        Root-mean-square error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = RMSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE RMSE
     MODULE PROCEDURE RMSE_sp_1d, RMSE_dp_1d, RMSE_sp_2d, RMSE_dp_2d, RMSE_sp_3d, RMSE_dp_3d
  END INTERFACE RMSE

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION BIAS_sp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: BIAS_sp_1d

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_sp_1d: number of arguments must be at least 2'
    !
    BIAS_sp_1d = average(y, mask=maske) - average(x, mask=maske)

  END FUNCTION BIAS_sp_1d

  FUNCTION BIAS_dp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: BIAS_dp_1d

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'BIAS_dp_1d: number of arguments must be at least 2'
    !
    BIAS_dp_1d = average(y, mask=maske) - average(x, mask=maske)

  END FUNCTION BIAS_dp_1d

  FUNCTION BIAS_sp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: BIAS_sp_2d

    INTEGER(i4)                                       :: n

    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_sp_2d: number of arguments must be at least 2'
    !
    BIAS_sp_2d = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2)/)))
    !
  END FUNCTION BIAS_sp_2d

  FUNCTION BIAS_dp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: BIAS_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_dp_2d: number of arguments must be at least 2'
    !
    BIAS_dp_2d = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2)/)))
    !
  END FUNCTION BIAS_dp_2d

  FUNCTION BIAS_sp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: BIAS_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    ! not really sopisticated, it has to be checked if the 3 numbers of x and y are matching in arry position
    if (n .LE. 1_i4) stop 'BIAS_sp_3d: number of arguments must be at least 2'
    !
    BIAS_sp_3d = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)))
    !
  END FUNCTION BIAS_sp_3d

  FUNCTION BIAS_dp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: BIAS_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    ! not really sopisticated, it has to be checked if the 3 numbers of x and y are matching in arry position
    if (n .LE. 1_i4) stop 'BIAS_dp_3d: number of arguments must be at least 2'
    !
    BIAS_dp_3d = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)))
    !
  END FUNCTION BIAS_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_sp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(INOUT)   :: mask
    REAL(sp)                                          :: LNNSE_sp_1d

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(sp)                               :: xmean
    REAL(sp), DIMENSION(size(x))           :: v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_sp_1d: number of arguments must be at least 2'

    ! mean of x
    xmean = average(log(x), mask=maske)
    !
    v1 = merge(log(y) - log(x), 0.0_sp, maske)
    v2 = merge(log(x) - xmean,  0.0_sp, maske)
    !
    LNNSE_sp_1d = 1.0_sp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION LNNSE_sp_1d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_dp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(INOUT)   :: mask
    REAL(dp)                                          :: LNNSE_dp_1d

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(dp)                               :: xmean
    REAL(dp), DIMENSION(size(x))           :: v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_dp_1d: number of arguments must be at least 2'

    ! mean of x
    xmean = average(log(x), mask=maske)
    !
    v1 = merge(log(y) - log(x), 0.0_dp, maske)
    v2 = merge(log(x) -  xmean, 0.0_dp, maske)
    !
    LNNSE_dp_1d = 1.0_dp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION LNNSE_dp_1d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_sp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(sp)                                          :: LNNSE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif
    !
    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_sp_2d: number of arguments must be at least 2'
    !
    ! mean of x
    xmean = average(reshape(log(x(:,:)), (/size(x, dim=1)*size(x, dim=2)/)), &
               mask=reshape(maske(:,:),  (/size(x, dim=1)*size(x, dim=2)/)))
    !
    LNNSE_sp_2d = 1.0_sp - &
         sum((log(y)-log(x))*(log(y)-log(x)), mask=maske) / &
         sum((log(x)- xmean)*(log(x)- xmean), mask=maske)
    !
  END FUNCTION LNNSE_sp_2d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_dp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(dp)                                          :: LNNSE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif
    !
    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_dp_2d: number of arguments must be at least 2'
    !
    ! mean of x
    xmean = average(reshape( log(x(:,:)), (/size(x, dim=1)*size(x, dim=2)/)), &
                 mask=reshape(maske(:,:), (/size(x, dim=1)*size(x, dim=2)/)))
    !
    LNNSE_dp_2d = 1.0_dp - &
         sum((log(y)-log(x))*(log(y)-log(x)), mask=maske) / &
         sum((log(x)-xmean )*(log(x)- xmean), mask=maske)
    !
  END FUNCTION LNNSE_dp_2d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_sp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: LNNSE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif
    !
    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_dp_2d: number of arguments must be at least 2'
    !
    ! mean of x
    xmean = average(reshape(log(x(:,:,:)), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
                mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    LNNSE_sp_3d = 1.0_sp - &
         sum((log(y)-log(x))*(log(y)-log(x)), mask=maske) / &
         sum((log(x)-xmean )*(log(x)-xmean ), mask=maske)
    !
  END FUNCTION LNNSE_sp_3d

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_dp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: LNNSE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif
    !
    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_dp_2d: number of arguments must be at least 2'
    !
    ! mean of x
    xmean = average(reshape(log(x(:,:,:)), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    LNNSE_dp_3d = 1.0_dp - &
         sum((log(y)-log(x))*(log(y)-log(x)), mask=maske) / &
         sum((log(x)-xmean )*(log(x)- xmean), mask=maske)
    !
  END FUNCTION LNNSE_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION MAE_sp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: MAE_sp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MAE_sp_1d: number of arguments must be at least 2'
    !
    MAE_sp_1d = SAE_sp_1d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_sp_1d

  FUNCTION MAE_dp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: MAE_dp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MAE_dp_1d: number of arguments must be at least 2'
    !
    MAE_dp_1d = SAE_dp_1d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_dp_1d

  FUNCTION MAE_sp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: MAE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MAE_sp_2d: number of arguments must be at least 2'
    !
    MAE_sp_2d = SAE_sp_2d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_sp_2d

  FUNCTION MAE_dp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: MAE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MAE_dp_2d: number of arguments must be at least 2'
    !
    MAE_dp_2d = SAE_dp_2d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_dp_2d

  FUNCTION MAE_sp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: MAE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MAE_sp_3d: number of arguments must be at least 2'
    !
    MAE_sp_3d = SAE_sp_3d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_sp_3d

  FUNCTION MAE_dp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: MAE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MAE_dp_3d: number of arguments must be at least 2'
    !
    MAE_dp_3d = SAE_dp_3d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION MSE_sp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: MSE_sp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MSE_sp_1d: number of arguments must be at least 2'
    !
    MSE_sp_1d = SSE_sp_1d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_sp_1d

  FUNCTION MSE_dp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: MSE_dp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MSE_dp_1d: number of arguments must be at least 2'
    !
    MSE_dp_1d = SSE_dp_1d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_dp_1d

  FUNCTION MSE_sp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: MSE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MSE_sp_2d: number of arguments must be at least 2'
    !
    MSE_sp_2d = SSE_sp_2d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_sp_2d

  FUNCTION MSE_dp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: MSE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MSE_dp_2d: number of arguments must be at least 2'
    !
    MSE_dp_2d = SSE_dp_2d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_dp_2d

  FUNCTION MSE_sp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: MSE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MSE_sp_3d: number of arguments must be at least 2'
    !
    MSE_sp_3d = SSE_sp_3d(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_sp_3d

  FUNCTION MSE_dp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: MSE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MSE_dp_3d: number of arguments must be at least 2'
    !
    MSE_dp_3d = SSE_dp_3d(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION NSE_sp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: NSE_sp_1d

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(sp)                               :: xmean
    REAL(sp), DIMENSION(size(x))           :: v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'NSE_sp_1d: number of arguments must be at least 2'
    ! mean of x
    xmean = average(x, mask=maske)
    !
    v1 = merge(y - x    , 0.0_sp, maske)
    v2 = merge(x - xmean, 0.0_sp, maske)
    !
    NSE_sp_1d = 1.0_sp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION NSE_sp_1d

  FUNCTION NSE_dp_1d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: NSE_dp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    REAL(dp), DIMENSION(size(x))                      :: v1, v2
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'NSE_dp_1d: number of arguments must be at least 2'
    ! mean of x
    xmean = average(x, mask=maske)
    !
    v1 = merge(y - x    , 0.0_dp, maske)
    v2 = merge(x - xmean, 0.0_dp, maske)
    !
    NSE_dp_1d = 1.0_dp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION NSE_dp_1d

  FUNCTION NSE_sp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: NSE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_sp_2d: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:), (/size(x, dim=1)*size(x, dim=2)/)))
    !
    NSE_sp_2d = 1.0_sp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_sp_2d

  FUNCTION NSE_dp_2d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: NSE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_dp_2d: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:), (/size(x, dim=1)*size(x, dim=2)/)))
    !
    NSE_dp_2d = 1.0_dp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_dp_2d

  FUNCTION NSE_sp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: NSE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_sp_3d: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    NSE_sp_3d = 1.0_sp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_sp_3d

  FUNCTION NSE_dp_3d(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: NSE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_dp_3d: number of arguments must be at least 2'
    ! Average of x
    xmean = average(reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    NSE_dp_3d = 1.0_dp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_dp_3d


  ! ------------------------------------------------------------------

  FUNCTION SAE_sp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: SAE_sp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SAE_sp_1d: number of arguments must be at least 2'
    !
    SAE_sp_1d = sum(abs(y - x) ,mask = maske)

  END FUNCTION SAE_sp_1d

  FUNCTION SAE_dp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: SAE_dp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SAE_dp_1d: number of arguments must be at least 2'
    !
    SAE_dp_1d = sum(abs(y - x) ,mask = maske)

  END FUNCTION SAE_dp_1d

  FUNCTION SAE_sp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: SAE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'SAE_sp_2d: number of arguments must be at least 2'
    !
    SAE_sp_2d = SAE_sp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SAE_sp_2d

  FUNCTION SAE_dp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: SAE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'SAE_dp_2d: number of arguments must be at least 2'
    !
    SAE_dp_2d = SAE_dp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SAE_dp_2d

  FUNCTION SAE_sp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: SAE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'SAE_sp_3d: number of arguments must be at least 2'
    !
    SAE_sp_3d = SAE_sp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2) &
         * size(maske, dim=3)/)) )

  END FUNCTION SAE_sp_3d

  FUNCTION SAE_dp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: SAE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'SAE_dp_3d: number of arguments must be at least 2'
    !
    SAE_dp_3d = SAE_dp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2) &
         * size(maske, dim=3)/)) )

  END FUNCTION SAE_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION SSE_sp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: SSE_sp_1d

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_sp_1d: number of arguments must be at least 2'
    !
    SSE_sp_1d = sum((y - x)**2_i4 ,mask = maske)

  END FUNCTION SSE_sp_1d

  FUNCTION SSE_dp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: SSE_dp_1d

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_dp_1d: number of arguments must be at least 2'
    !
    SSE_dp_1d = sum((y - x)**2_i4 ,mask = maske)

  END FUNCTION SSE_dp_1d

  FUNCTION SSE_sp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: SSE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1), size(x,dim=2)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_sp_2d: number of arguments must be at least 2'
    !
    SSE_sp_2d = SSE_sp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SSE_sp_2d

  FUNCTION SSE_dp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: SSE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1), size(x,dim=2)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_dp_2d: number of arguments must be at least 2'
    !
    SSE_dp_2d = SSE_dp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SSE_dp_2d

  FUNCTION SSE_sp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: SSE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL, DIMENSION(size(x,dim=1), size(x,dim=2),&
         size(x,dim=3))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_sp_3d: number of arguments must be at least 2'
    !
    SSE_sp_3d = SSE_sp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)   &
         * size(maske, dim=3)/)))

  END FUNCTION SSE_sp_3d

  FUNCTION SSE_dp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: SSE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL, DIMENSION(size(x,dim=1), size(x,dim=2),&
         size(x,dim=3))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_dp_3d: number of arguments must be at least 2'
    !
    SSE_dp_3d = SSE_dp_1d(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)   &
         * size(maske, dim=3)/)))

  END FUNCTION SSE_dp_3d

  ! ------------------------------------------------------------------

  FUNCTION RMSE_sp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: RMSE_sp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_sp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'RMSE_sp_1d: number of arguments must be at least 2'
    !
    RMSE_sp_1d = sqrt(MSE_sp_1d(x,y,mask=maske))

  END FUNCTION RMSE_sp_1d

  FUNCTION RMSE_dp_1d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: RMSE_dp_1d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_dp_1d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'RMSE_dp_1d: number of arguments must be at least 2'
    !
    RMSE_dp_1d = sqrt(MSE_dp_1d(x,y,mask=maske))

  END FUNCTION RMSE_dp_1d

  FUNCTION RMSE_sp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: RMSE_sp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_sp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'RMSE_sp_2d: number of arguments must be at least 2'
    !
    RMSE_sp_2d = sqrt(MSE_sp_2d(x,y,mask=maske))

  END FUNCTION RMSE_sp_2d

  FUNCTION RMSE_dp_2d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: RMSE_dp_2d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_dp_2d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'RMSE_dp_2d: number of arguments must be at least 2'
    !
    RMSE_dp_2d = sqrt(MSE_dp_2d(x,y,mask=maske))

  END FUNCTION RMSE_dp_2d

  FUNCTION RMSE_sp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: RMSE_sp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_sp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'RMSE_sp_3d: number of arguments must be at least 2'
    !
    RMSE_sp_3d = sqrt(MSE_sp_3d(x,y,mask=maske))

  END FUNCTION RMSE_sp_3d

  FUNCTION RMSE_dp_3d(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: RMSE_dp_3d

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_dp_3d: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'RMSE_dp_3d: number of arguments must be at least 2'
    !
    RMSE_dp_3d = sqrt(MSE_dp_3d(x,y,mask=maske))

  END FUNCTION RMSE_dp_3d

END MODULE mo_errormeasures
