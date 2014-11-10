!> \file mo_julian.f90

!> \brief Julian date conversion routines

!> \details Julian date to and from day, month, year, and also from day, month, year, hour, minute, and second.
!> Also convience routines for Julian dates of IMSL are provided.

!> \note Julian day definition starts at noon of the 1st January 4713 BC.\n
!> Here, the astronomical definition is used,
!> i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
!> This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
!> \n
!> julday and caldat start at midnight of the 1st January 4713 BC.
!> So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
!> Use date2dec with dec2date together for fractional Julian dates
!> and use julday with caldat together for integer Julian days.

!> \author Matthias Cuntz
!> \date Dec 2011

MODULE mo_julian

  ! This module provides Julian day conversion routines

  ! Written  Matthias Cuntz, Dec 2011
  ! Modified Matthias Cuntz, Jan 2013 - added date2dec and dec2date
  ! Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !                                     removed numerical recipes

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

  ! Copyright 2011-2014 Matthias Cuntz

  USE mo_kind, ONLY: i4, i8, dp

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

  !>        The zeroth Julian Day is 01.01.-4712, i.e. the 1st January 4713 BC.

  !>        Julian day definition starts at 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

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
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        julday and caldat start at midnight of the 1st January 4713 BC.
  !>        So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
  !>        Use date2dec with dec2date together for fractional Julian dates
  !>        and use julday with caldat together for integer Julian days.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         call caldat(2415021, dd, mm, yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified julday from Numerical Recipes
  !>        \date Dec 2011
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes

  ELEMENTAL PURE SUBROUTINE caldat(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i8) :: A, B, C, D, E, g
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    if (julian < IGREG) then
       A = int(julian, i8) ! julian
    else
       g = int((real(julian,dp)-1867216.25_dp)/36524.25_dp, i8) ! gregorian
       A = julian + 1_i8 + g - g/4_i8
    endif

    B = A + 1524_i8
    C = int((real(B,dp)-122.1_dp) / 365.25_dp, i8)
    D = int(365.25_dp * real(C,dp), i8)
    E = int(real(B-D,dp) / 30.6001_dp, i8)

    dd = int(B - D - int(30.6001_dp*real(E,dp), i8), i4)

    if (E<14_i8) then
       mm = int(E-1_i8, i4)
    else
       mm = int(E-13_i8, i4)
    endif

    if (mm > 2) then
       yy = int(C - 4716_i8, i4)
    else
       yy = int(C - 4715_i8, i4)
    endif

  END SUBROUTINE caldat

  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second

  !>        \details In this routine date2dec returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC 12:00:00 h.

  !>        Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         date2dec = date2dec(dd, mm, yy, hh, nn, ss)

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
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec &mdash;     Fractional Julian day

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        julday and caldat start at midnight of the 1st January 4713 BC.
  !>        So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
  !>        Use date2dec with dec2date together for fractional Julian dates
  !>        and use julday with caldat together for integer Julian days.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         julian = date2dec(01,01,1990)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         julian = date2dec(01,01,1990,12,00)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.
  !         Numerical regulation of fractions is after
  !             IDL routine julday.pro. Copyright (c) 1988-2011, ITT Visual Information Solutions.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes

  ELEMENTAL PURE FUNCTION date2dec(dd,mm,yy,hh,nn,ss)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(IN), OPTIONAL :: hh, nn, ss
    REAL(dp)                          :: date2dec

    INTEGER(i4), PARAMETER :: IGREG2 = 15 + 31*(10+12*1582)
    INTEGER(i4), PARAMETER :: IGREG1 =  4 + 31*(10+12*1582)
    INTEGER(i4) :: idd, imm, iyy
    REAL(dp)    :: ihh, inn, iss
    INTEGER(i8) :: jm, jy
    REAL(dp)    :: jd, H, eps
    INTEGER(i8) :: A, B

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    ihh = 0.0_dp
    if (present(hh)) ihh = real(hh,dp)
    inn = 0.0_dp
    if (present(nn)) inn = real(nn,dp)
    iss = 0.0_dp
    if (present(ss)) iss = real(ss,dp)

    if (imm > 2) then
       jm = int(imm, i8)
       jy = int(iyy, i8)
    else
       jm = int(imm+12, i8)
       jy = int(iyy-1, i8)
    endif

    jd = real(idd, dp)
    
    H = ihh/24._dp + inn/1440._dp + iss/86400._dp

    if (dd+31*(mm+12*yy) >= IGREG2) then ! gregorian
       A = jy/100_i8
       B = 2_i8 - A + A/4_i8
    else if (dd+31*(mm+12*yy) <= IGREG1) then ! julian
       B = 0_i8
    ! else
    !    stop 'No Gregorian dates between 04.10.1582 and 15.10.1582'
    endif
    
    ! Fractional Julian day starts at noon
    date2dec = floor(365.25_dp*real(jy+4716_i8,dp)) + floor(30.6001_dp*real(jm+1_i8,dp)) + jd + H + real(B,dp) - 1524.5_dp

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec), eps)
    date2dec = date2dec + eps

  END FUNCTION date2dec

  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day

  !>        \details Inverse of the function date2dec. Here dec2date is input as a fractional Julian Day,
  !>        which starts at noon of the 1st January 4713 BC, i.e. 01.01.-4712.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC at noon.

  !>        Julian day definition starts at 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         call dec2date(fJulian, dd, mm, yy, hh, nn, ss)

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
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        julday and caldat start at midnight of the 1st January 4713 BC.
  !>        So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
  !>        Use date2dec with dec2date together for fractional Julian dates
  !>        and use julday with caldat together for integer Julian days.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         call caldat(2415020.5, dd, mm, yy, hh, nn)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         call caldat(2415021., dd, mm, yy, hh, nn, ss)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.
  !         Here the sometimes 60 sec as output are corrected at the end.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes

  ELEMENTAL PURE SUBROUTINE dec2date(julian,dd,mm,yy,hh,nn,ss)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN)            :: julian
    INTEGER(i4), INTENT(OUT), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(OUT), OPTIONAL :: hh, nn, ss

    INTEGER(i4) :: day, month, year, hour, minute, second
    REAL(dp)    :: fraction
    ! REAL(dp)    :: eps

    INTEGER(i8) :: A, B, C, D, E, g, Z
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    Z = int(julian + 0.5, i8)

    if (Z < IGREG) then
       A = Z ! julian
    else
       g = int((real(Z,dp)-1867216.25_dp)/36524.25_dp, i8) ! gregorian
       A = Z + 1_i8 + g - g/4_i8
    endif

    B = A + 1524_i8
    C = int((real(B,dp)-122.1_dp) / 365.25_dp, i8)
    D = int(365.25_dp * real(C,dp), i8)
    E = int(real(B-D,dp) / 30.6001_dp, i8)

    day = int(B - D - int(30.6001_dp*real(E,dp), i8), i4)

    if (E<14_i8) then
       month = int(E-1_i8, i4)
    else
       month = int(E-13_i8, i4)
    endif

    if (month > 2) then
       year = int(C - 4716_i8, i4)
    else
       year = int(C - 4715_i8, i4)
    endif

    ! ! Fractional part
    ! eps = 1e-12_dp ! ~ 5000*epsilon(1.0_dp)
    ! eps = max(eps * abs(real(Z,dp)), eps)
    ! fraction = julian + 0.5_dp - real(Z,dp)
    ! hour     = min(max(floor(fraction * 24.0_dp + eps), 0), 23)
    ! fraction = fraction - real(hour,dp)/24.0_dp
    ! minute   = min(max(floor(fraction*1440.0_dp + eps), 0), 59)
    ! second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! Fractional part
    fraction = julian + 0.5_dp - real(Z,dp)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat(julday(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second

  END SUBROUTINE dec2date

  ! ------------------------------------------------------------------

  !     NAME
  !         julday

  !     PURPOSE
  !>        \brief Julian day from day, month and year

  !>        \details In this routine julday returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC 12:00:00 h.

  !>        Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

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
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        julday and caldat start at midnight of the 1st January 4713 BC.
  !>        So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
  !>        Use date2dec with dec2date together for fractional Julian dates
  !>        and use julday with caldat together for integer Julian days.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         ! 2440588 is 01.01.1970
  !         julian = julday(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified julday from Numerical Recipes
  !>        \date Dec 2011
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes

  ELEMENTAL PURE FUNCTION julday(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: julday

    INTEGER(i4), PARAMETER :: IGREG2 = 15 + 31*(10+12*1582)
    INTEGER(i4), PARAMETER :: IGREG1 =  4 + 31*(10+12*1582)
    INTEGER(i8) :: jd, jm, jy
    INTEGER(i8) :: A, B

    if (mm > 2) then
       jm = int(mm, i8)
       jy = int(yy, i8)
    else
       jm = int(mm+12, i8)
       jy = int(yy-1, i8)
    endif

    jd = int(dd, i8)

    if (dd+31*(mm+12*yy) >= IGREG2) then ! gregorian
       A = jy/100_i8
       B = 2_i8 - A + A/4_i8
    else if (dd+31*(mm+12*yy) <= IGREG1) then ! julian
       B = 0_i8
    ! else
    !    stop 'No Gregorian dates between 04.10.1582 and 15.10.1582'
    endif
    
    ! add 0.5 to Wiki formula because formula was for fractional day
    ! julday = int(365.25_dp*real(jy+4716_i8,dp) + real(int(30.6001*real(jm+1_i8,dp),i8),dp) + real(jd+B,dp) - 1524.5_dp, i4)
    julday = int(365.25_dp*real(jy+4716_i8,dp) + real(int(30.6001*real(jm+1_i8,dp),i8),dp) + real(jd+B,dp) - 1524.5_dp + 0.5_dp, i4)

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
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2011

  ELEMENTAL PURE FUNCTION ndays(dd,mm,yy)

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
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2011

  ELEMENTAL PURE SUBROUTINE ndyin(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    call caldat(julian+IMSLday, dd, mm, yy)

  END SUBROUTINE ndyin

  ! ------------------------------------------------------------------

END MODULE mo_julian
