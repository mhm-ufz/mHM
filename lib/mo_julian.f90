!> \file mo_julian.f90

!> \brief Julian date conversion routines

!> \details Julian date to and from day, month, year, and also from day, month, year, hour, minute, and second.
!> Also convience routines for Julian dates of IMSL are provided.

!> \note Julian day definition starts at noon of the 1st January 4713.\n
!> julday and caldat start at midnight of the 1st January 4713. So date2dec and julday as well as dec2date and caldat
!> are shifted by half a day.\n
!> Use date2dec with dec2date together for fractional Julian dates
!> and use julday with caldat together for integer Julian days.

!> \author Matthias Cuntz
!> \date Dec 2011

MODULE mo_julian

  ! This module provides Julian day conversion routines

  ! Written  Matthias Cuntz, Dec 2011
  ! Modified Matthias Cuntz, Jan 2013 - added date2dec and dec2date

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

  ! Copyright 2011-2013 Matthias Cuntz


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

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: caldat   ! Day , month and year from Julian day
  PUBLIC :: date2dec ! Fractional Julian day from day, month, year, hour, minute, and second
  PUBLIC :: dec2date ! Day, month, year, hour, minute, and second from fractional Julian day
  PUBLIC :: julday   ! Julian day from day, month and year
  PUBLIC :: ndays    ! IMSL Julian day from day, month and year
  PUBLIC :: ndyin    ! Day, month and year from IMSL Julian day

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         caldat

  !     PURPOSE
  !>        \brief Day, month and year from Julian day

  !>        \details Inverse of the function julday. Here julian is input as a Julian Day Number,
  !>        and the routine outputs id, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day is 01.01.-4713, i.e. the 1st January 4713 BC.

  !     CALLING SEQUENCE
  !         call caldat(Julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: Julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         call caldat(2415021, dd, mm, yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Written,  Matthias Cuntz - modified caldat from Numerical Recipes
  !>        \date Dec 2011

  SUBROUTINE caldat(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4) :: ja, jalpha, jb, jc, jd, je
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    if (julian >= IGREG) then
       jalpha = int((real(julian-1867216_i4,dp)-0.25_dp)/36524.25_dp, i4)
       ja     = julian + 1 + jalpha - int(0.25_dp*real(jalpha,dp), i4)
    else
       ja     = julian
    end if
    jb = ja + 1524_i4
    jc = int(6680.0_dp+(real(jb-2439870_i4,dp)-122.1_dp)/365.25_dp, i4)
    jd = 365*jc + int(0.25_dp*real(jc,dp), i4)
    je = int(real(jb-jd,dp)/30.6001_dp, i4)
    dd = jb - jd - int(30.6001_dp*real(je,dp), i4)
    mm = je - 1
    if (mm > 12) mm = mm - 12
    yy = jc - 4715_i4
    if (mm > 2)  yy = yy - 1
    if (yy <= 0) yy = yy - 1

  END SUBROUTINE caldat

  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second

  !>        \details In this routine date2dec returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.
  !>        Positive year signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.

  !>        The zeroth Julian Day is 01.01.-4713 at noon, i.e. the 1st January 4713 BC 12:00:00 h.

  !     CALLING SEQUENCE
  !         date2dec = date2dec(dd, mm, yy, hh, ii, ss)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ii"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec &mdash;     Fractional Julian day

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.
  !>        \note Julian day definition starts at noon of the 1st January 4713.\n
  !>        julday starts at midnight of the 1st January 4713. So date2dec and julday are
  !>        shifted by half a day.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         julian = date2dec(01,01,1990)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         julian = date2dec(01,01,1990,12,00)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996
  !         IDL routine julday.pro  Copyright (c) 1988-2011, ITT Visual Information Solutions.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013

  ELEMENTAL PURE FUNCTION date2dec(dd,mm,yy,hh,ii,ss)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(IN), OPTIONAL :: hh, ii, ss
    REAL(dp)                          :: date2dec

    INTEGER(i4), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(i4) :: idd, imm, iyy
    INTEGER(i4) :: ja, jm, jy
    INTEGER(i4) :: idate2dec
    REAL(dp)    :: ihh, iii, iss
    REAL(dp)    :: eps

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    ihh = 0.0_dp
    if (present(hh)) ihh = real(hh,dp)
    iii = 0.0_dp
    if (present(ii)) iii = real(ii,dp)
    iss = 0.0_dp
    if (present(ss)) iss = real(ss,dp)

    ! Julian day as in julday
    jy = iyy
    !if (jy == 0) stop 'date2dec: there is no year zero'
    if (jy < 0) jy = jy+1
    if (imm > 2) then
       jm = imm+1
    else
       jy = jy-1
       jm = imm+13
    end if

    idate2dec = 365*jy + int(0.25_dp*real(jy,dp)+2000._dp, i4) + int(30.6001_dp*real(jm,dp), i4) + idd + 1718995_i4
    if (idd+31*(imm+12*iyy) >= IGREG) then
       ja = int(0.01_dp*jy, i4)
       idate2dec = idate2dec + 2 - ja + int(0.25_dp*ja, i4)
    end if
    date2dec = real(idate2dec,dp)

    ! Fractional Julian day starts at noon

    ! Add a small offset (propor. julian dat) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec), eps)
    ! For Hours, divide by 24, then subtract 0.5 because Julian date starts at noon.
    date2dec = date2dec + ihh/24.0_dp - 0.5_dp + iii/1440.0_dp + iss/86400.0_dp + eps

  END FUNCTION date2dec

  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day

  !>        \details Inverse of the function date2dec. Here dec2date is input as a fractional Julian Day,
  !>        which starts at noon of the 1st January 4713 BC.
  !>        The routine outputs dd, mm, yy, hh, ii, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.

  !>        The zeroth Julian Day is 01.01.-4713 at noon, i.e. the 1st January 4713 BC at noon.

  !     CALLING SEQUENCE
  !         call dec2date(fJulian, dd, mm, yy, hh, ii, ss)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: ii"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.
  !>        \note Julian day definition starts at noon of the 1st January 4713.\n
  !>        caldat starts at midnight of the 1st January 4713. So caldat and dec2date are
  !>        shifted by half a day.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         call caldat(2415020.5, dd, mm, yy, hh, ii)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         call caldat(2415021., dd, mm, yy, hh, ii, ss)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996
  !         IDL routine caldat.pro  Copyright (c) 1988-2011, ITT Visual Information Solutions.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013

  ELEMENTAL PURE SUBROUTINE dec2date(julian,dd,mm,yy,hh,ii,ss)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN)            :: julian
    INTEGER(i4), INTENT(OUT), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(OUT), OPTIONAL :: hh, ii, ss

    INTEGER(i4) :: ja, jalpha, jb, jc, jd, je
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4
    INTEGER(i4) :: day, month, year, hour, minute, second
    INTEGER(i4) :: ijulian
    REAL(dp)    :: fraction
    REAL(dp)    :: eps

    ! Julian of caldat
    ijulian = floor(julian + 0.5_dp)

    ! Normal caldat
    if (ijulian >= IGREG) then
       jalpha = int((real(ijulian-1867216_i4,dp)-0.25_dp)/36524.25_dp, i4)
       ja     = ijulian + 1 + jalpha - int(0.25_dp*real(jalpha,dp), i4)
    else
       ja     = ijulian
    end if
    jb = ja + 1524_i4
    jc = int(6680.0_dp+(real(jb-2439870_i4,dp)-122.1_dp)/365.25_dp, i4)
    jd = 365*jc + int(0.25_dp*real(jc,dp), i4)
    je = int(real(jb-jd,dp)/30.6001_dp, i4)
    day = jb - jd - int(30.6001_dp*real(je,dp), i4)
    month = je - 1
    if (month > 12) month = month - 12
    year = jc - 4715_i4
    if (month > 2)  year = year - 1
    if (year <= 0) year = year - 1

    ! Fractional part
    eps = 1e-12_dp ! ~ 5000*epsilon(1.0_dp)
    eps = max(eps * abs(real(ijulian,dp)), eps)
    fraction = julian + 0.5_dp - ijulian
    hour     = min(max(floor(fraction * 24.0_dp + eps), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp + eps), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(ii)) ii = minute
    if (present(ss)) ss = second

  END SUBROUTINE dec2date

  ! ------------------------------------------------------------------

  !     NAME
  !         julday

  !     PURPOSE
  !>        \brief Julian day from day, month and year

  !>        \details In this routine julday returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables. Positive year
  !>        signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.

  !>        The zeroth Julian Day is 01.01.-4713, i.e. the 1st January 4713 BC.

  !     CALLING SEQUENCE
  !         julian = julday(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian &mdash;     Julian day

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         julian = julday(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Written,  Matthias Cuntz - modified julday from Numerical Recipes
  !>        \date Dec 2011

  FUNCTION julday(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: julday
    INTEGER(i4), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(i4) :: ja,jm,jy

    jy=yy
    if (jy == 0) stop 'julday: there is no year zero'
    if (jy < 0) jy=jy+1
    if (mm > 2) then
       jm=mm+1
    else
       jy=jy-1
       jm=mm+13
    end if

    !Bug julday=int(365.25_dp*jy)+int(30.6001_dp*jm)+dd+1720995_i4
    julday = 365*jy + int(0.25_dp*real(jy,dp)+2000._dp, i4) + int(30.6001_dp*real(jm,dp), i4) + dd + 1718995_i4
    if (dd+31*(mm+12*yy) >= IGREG) then
       ja = int(0.01_dp*jy, i4)
       julday=julday + 2 - ja + int(0.25_dp*ja, i4)
    end if

  END FUNCTION julday

  ! ------------------------------------------------------------------

  !     NAME
  !         ndays

  !     PURPOSE
  !>        \brief IMSL Julian day from day, month and year

  !>        \details In this routine ndays returns the IMSL Julian Day Number. Julian days begin at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables. IMSL treats 01.01.1900
  !>        as a reference and assigns a Julian day 0 to it.
  !>            ndays = julday(dd,mm,yy) - julday(01,01,1900)

  !     CALLING SEQUENCE
  !         julian = ndays(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of IMSL Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of IMSL Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of IMSL Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian &mdash;     IMSL Julian day, i.e. days before or after 01.01.1900

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 0 is 01.01.1900
  !         julian = ndays(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Dec 2011

  FUNCTION ndays(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: ndays

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    ndays = julday(dd, mm, yy) - IMSLday

  END FUNCTION ndays

  ! ------------------------------------------------------------------

  !     NAME
  !         ndyin

  !     PURPOSE
  !>        \brief Day, month and year from IMSL Julian day

  !>        \details Inverse of the function ndys. Here ISML Julian is input as a Julian Day Number
  !>        minus the Julian Day Number of 01.01.1900, and the routine outputs id, mm, and yy
  !>        as the day, month, and year on which the specified Julian Day started at noon.
  !>          ndyin is caldat(IMSLJulian + 2415021, dd, mm, yy)

  !     CALLING SEQUENCE
  !         call ndyin(julian, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julian"     IMSL Julian day, i.e. days before or after 01.01.1900

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of IMSL Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of IMSL Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of IMSL Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 0 is 01.01.1900
  !         call ndyin(0,dd,mm,yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Dec 2011

  SUBROUTINE ndyin(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    call caldat(julian+IMSLday, dd, mm, yy)

  END SUBROUTINE ndyin

  ! ------------------------------------------------------------------

END MODULE mo_julian
