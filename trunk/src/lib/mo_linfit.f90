!> \file mo_linfit.f90

!> \brief  Fitting a straight line.

!> \details This module provides a routine to fit a straight line with model I or model II regression.

!> \authors Matthias Cuntz
!> \date Mar 2011

MODULE mo_linfit

  ! This module provides a routine to fit a straight line with model I or model II regression.
  ! Written  Matthias Cuntz, Mar 2011

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

  USE mo_kind, ONLY: sp, dp

  Implicit NONE

  PUBLIC :: linfit ! Fitting straight line (without error bars on input), Model I or Model II

  ! ------------------------------------------------------------------

  !     NAME
  !         linfit

  !     PURPOSE
  !>        \brief Fits a straight line to input data by minimizing chi^2.

  !>        \details Given a set of data points x(1:ndata), y(1:ndata),
  !>         fit them to a straight line \f$ y = a+bx \f$ by minimizing chi2.\n
  !>         Model I minimizes y vs. x while Model II takes the geometric mean of y vs. x and x vs. y.
  !>         Returned is the fitted line at x.
  !>         Optional returns are a, b and their respective probable uncertainties siga and sigb,
  !>         and the chi-square chi2.

  !     CALLING SEQUENCE
  !         out = linfit(x, y, a=a, b=b, siga=siga, sigb=sigb, chi2=chi2, model2=model2)

  !     INTENT(IN)
  !>         \param[in] "real(sp/dp) :: x(:)"                1D-array with input x
  !>         \param[in] "real(sp/dp) :: y(:)"                1D-array with input y

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>         \param[in] "logical, optional :: model2"        If present, use geometric mean regression
  !>                                                         instead of ordinary least square

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(sp/dp), dimension(M) :: a"      intercept
  !>        \param[out] "real(sp/dp), dimension(M) :: b"      slope
  !>        \param[out] "real(sp/dp), dimension(M) :: siga"   error on intercept
  !>        \param[out] "real(sp/dp), dimension(M) :: sigb"   error on slope
  !>        \param[out] "real(sp/dp)               :: chisq"  Minimum chi^2

  !     RETURN
  !>       \return real(sp/dp), dimension(:), allocatable :: out &mdash; fitted values at x(:).

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ytmp = linfit(x, y, a=inter, b=slope, model2=.true.)

  !     LITERATURE
  !     Model I follows closely
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996
  !     For Model II (geometric mean regression)
  !         Robert R. Sokal and F. James Rohlf, Biometry: the principle and practice of statistics
  !             in biological research, Freeman & Co., ISBN 0-7167-2411-1

  !     HISTORY
  !         Written Matthias Cuntz, March 2011 - following closely the routine fit of Numerical Recipes
  !                                            - linfit Model II: geometric mean regression
  !         Modified Matthias Cuntz, Nov 2011  - sp, dp
  !                                            - documentation
  INTERFACE linfit
     MODULE PROCEDURE linfit_sp, linfit_dp
  END INTERFACE linfit

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION linfit_dp(x, y, a, b, siga, sigb, chi2, model2)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: x, y
    REAL(dp), OPTIONAL,     INTENT(OUT) :: a, b, siga, sigb, chi2
    LOGICAL , OPTIONAL,     INTENT(IN)  :: model2
    REAL(dp), DIMENSION(:), allocatable :: linfit_dp

    REAL(dp) :: sigdat, nx, sx, sxoss, sy, st2
    REAL(dp), DIMENSION(size(x)), TARGET :: t
    REAL(dp) :: aa, bb, cchi2
    LOGICAL :: mod2
    REAL(dp) :: mx, my, sx2, sy2, sxy, sxy2, syx2, ssigb
    !REAL(dp) :: r

    if (size(x) /= size(y))   stop 'linfit_dp: size(x) /= size(y)'
    if (.not. allocated(linfit_dp)) allocate(linfit_dp(size(x)))
    if (present(model2)) then
       mod2 = model2
    else
       mod2 = .false.
    endif

    if (mod2) then
       nx = real(size(x),dp)
       mx = sum(x)/nx
       my = sum(y)/nx
       sx2 = sum((x-mx)*(x-mx))
       sy2 = sum((y-my)*(y-my))
       sxy = sum((x-mx)*(y-my))
       !r = sxy / sqrt(sx2) / sqrt(sy2)
       bb = sign(sqrt(sy2/sx2),sxy)
       aa = my - bb*mx
       linfit_dp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2)) then
          t(:) = y(:)-linfit_dp(:)
          chi2 = dot_product(t,t)
       endif
       if (present(siga) .or. present(sigb)) then
          syx2 = (sy2 - sxy*sxy/sx2) / (nx-2.0_dp)
          sxy2 = (sx2 - sxy*sxy/sy2) / (nx-2.0_dp)
          ! syx2 should be >0
          ! ssigb = sqrt(syx2/sx2)
          ssigb = sqrt(abs(syx2)/sx2)
          if (present(sigb)) sigb = ssigb
          if (present(siga)) then
             ! siga = sqrt(syx2*(1.0_dp/nX+mx*mx/sx2))
             siga = sqrt(abs(syx2)*(1.0_dp/nX+mx*mx/sx2))
             ! Add Extra Term for Error in xmean which is not in Sokal & Rohlf.
             ! They take the error estimate of the chi-squared error for a.
             ! siga = sqrt(siga*siga + bb*bb*sxy2/nx)
             siga = sqrt(siga*siga + bb*bb*abs(sxy2)/nx)
          endif
       endif
    else
       nx = real(size(x),dp)
       sx = sum(x)
       sy = sum(y)
       sxoss = sx/nx
       t(:) = x(:)-sxoss
       bb = dot_product(t,y)
       st2 = dot_product(t,t)
       bb = bb/st2
       aa = (sy-sx*bb)/nx
       linfit_dp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2) .or. present(siga) .or. present(sigb)) then
          t(:) = y(:)-linfit_dp(:)
          cchi2 = dot_product(t,t)
          if (present(chi2)) chi2 = cchi2
          if (present(siga) .or. present(sigb)) then
             sigdat = sqrt(cchi2/(size(x)-2))
             if (present(siga)) then
                siga = sqrt((1.0_dp+sx*sx/(nx*st2))/nx)
                siga = siga*sigdat
             endif
             if (present(sigb)) then
                sigb = sqrt(1.0_dp/st2)
                sigb = sigb*sigdat
             endif
          endif
       endif
    endif

  END FUNCTION linfit_dp


  FUNCTION linfit_sp(x, y, a, b, siga, sigb, chi2, model2)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN)  :: x, y
    REAL(sp), OPTIONAL,     INTENT(OUT) :: a, b, siga, sigb, chi2
    LOGICAL , OPTIONAL,     INTENT(IN)  :: model2
    REAL(sp), DIMENSION(:), allocatable :: linfit_sp

    REAL(sp) :: sigdat, nx, sx, sxoss, sy, st2
    REAL(sp), DIMENSION(size(x)), TARGET :: t
    REAL(sp) :: aa, bb, cchi2
    LOGICAL :: mod2
    REAL(sp) :: mx, my, sx2, sy2, sxy, sxy2, syx2, ssigb
    !REAL(sp) :: r

    if (size(x) /= size(y))   stop 'linfit_sp: size(x) /= size(y)'
    if (.not. allocated(linfit_sp)) allocate(linfit_sp(size(x)))
    if (present(model2)) then
       mod2 = model2
    else
       mod2 = .false.
    endif

    if (mod2) then
       nx = real(size(x),sp)
       mx = sum(x)/nx
       my = sum(y)/nx
       sx2 = sum((x-mx)*(x-mx))
       sy2 = sum((y-my)*(y-my))
       sxy = sum((x-mx)*(y-my))
       !r = sxy / sqrt(sx2) / sqrt(sy2)
       bb = sign(sqrt(sy2/sx2),sxy)
       aa = my - bb*mx
       linfit_sp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2)) then
          t(:) = y(:)-linfit_sp(:)
          chi2 = dot_product(t,t)
       endif
       if (present(siga) .or. present(sigb)) then
          syx2 = (sy2 - sxy*sxy/sx2) / (nx-2.0_sp)
          sxy2 = (sx2 - sxy*sxy/sy2) / (nx-2.0_sp)
          ! syx2 should be >0
          ! ssigb = sqrt(syx2/sx2)
          ssigb = sqrt(abs(syx2)/sx2)
          if (present(sigb)) sigb = ssigb
          if (present(siga)) then
             ! siga = sqrt(syx2*(1.0_sp/nX+mx*mx/sx2))
             siga = sqrt(abs(syx2)*(1.0_sp/nX+mx*mx/sx2))
             ! Add Extra Term for Error in xmean which is not in Sokal & Rohlf.
             ! They take the error estimate of the chi-squared error for a.
             ! siga = sqrt(siga*siga + bb*bb*sxy2/nx)
             siga = sqrt(siga*siga + bb*bb*abs(sxy2)/nx)
          endif
       endif
    else
       nx = real(size(x),sp)
       sx = sum(x)
       sy = sum(y)
       sxoss = sx/nx
       t(:) = x(:)-sxoss
       bb = dot_product(t,y)
       st2 = dot_product(t,t)
       bb = bb/st2
       aa = (sy-sx*bb)/nx
       linfit_sp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2) .or. present(siga) .or. present(sigb)) then
          t(:) = y(:)-linfit_sp(:)
          cchi2 = dot_product(t,t)
          if (present(chi2)) chi2 = cchi2
          if (present(siga) .or. present(sigb)) then
             sigdat = sqrt(cchi2/(size(x)-2))
             if (present(siga)) then
                siga = sqrt((1.0_sp+sx*sx/(nx*st2))/nx)
                siga = siga*sigdat
             endif
             if (present(sigb)) then
                sigb = sqrt(1.0_sp/st2)
                sigb = sigb*sigdat
             endif
          endif
       endif
    endif

  END FUNCTION linfit_sp

  ! ------------------------------------------------------------------

END MODULE mo_linfit
