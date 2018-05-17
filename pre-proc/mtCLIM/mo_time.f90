module mo_time
! * Version:   1.0 2015-12-22
! download from http://fortranwiki.org/fortran/show/m_time
! 2017-06-15, Brenner Johannes
! added get_time mHM subroutine
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------
!EPOCH TIME (UT starts at 0000 on 1 Jan. 1970)
   public d2u            !(dat,UNIXTIME,IERR)                 ! Convert date array to Unix Time
   public u2d            !(DAT,unixtime,IERR)                 ! Convert Unix Time to date array
   public du             !(dat) result (UNIXTIME)             ! Convert date array to Unix Time
   public ud             !(unixtime) result (DAT)             ! Convert Unix Time to date array
!JULIAN
   public j2d            !(DAT,julian,IERR)                   ! Convert Julian date to date array
   public d2j            !(dat,JULIAN,IERR)                   ! Convert date array to Julian date
   public dj             !(dat) result (JULIAN)               ! Convert date array to Julian Day
   public jd             !(julian) result (DAT)               ! Convert Julian Day to date array
!DAY OF WEEK
   public dow            !(dat,[WEEKDAY],[DAY],IERR)          ! Convert date array to day of the week as number(Sun=1) and name
!WEEK OF YEAR
   public woy  !(dat,ISO_YEAR,ISO_WEEK,ISO_WEEKDAY,ISO_NAME)  ! Calculate iso-8601 Week-numbering year date yyyy-Www-d
!ORDINAL DAY
   public d2o            !(dat) result(ORDINAL)               ! given date array return ordinal day of year, Jan 1st=1
!PRINTING DATES
   public fmtdate        !(dat,format) result (TIMESTRING)    ! Convert date array to string using format
   public fmtdate_usage  !()                                  ! display macros recognized by fmtdate(3f)
   public now            !(format) result (NOW)               ! return string representing current time given format
!MONTH NAME
   public v2mo           !(month_number) result (MONTH_NAME)  ! given month number return month name
!C INTERFACE
   public sys_sleep      !(wait_seconds)                      ! Call sleep(3c)
   public get_time
!-----------------------------------------------------------------------------------------------------------------------------------
!INTERNAL
   integer,parameter,private :: dp=kind(0.0d0)
   real(kind=dp)             :: secday=86400.0d0              ! 24:00:00 hours as seconds
!-----------------------------------------------------------------------------------------------------------------------------------
!C INTERFACE
 contains
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine d2j(dat,julian,ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! * Author:    John S. Urban
! * Version:   1.0 2015-12-21
! * Reference: From Wikipedia, the free encyclopedia 2015-12-19
! * There is no year zero
! * Julian Day must be non-negative
! * Julian Day starts at noon; while Civil Calendar date starts at midnight
!-----------------------------------------------------------------------------------------------------------------------------------
character(len=:),parameter :: ident="@(#)d2j(3f): Converts proleptic Gregorian date array to Julian Day -JSU 2015-12-21"
integer,intent(in)         :: dat(8)   ! array like returned by DATE_AND_TIME(3f)
real(kind=dp),intent(out)  :: julian   ! Julian Day (non-negative, but may be non-integer)
integer,intent(out)        :: ierr     ! Error return, 0 for successful execution,-1=invalid year,-2=invalid month,-3=invalid day,
                                       ! -4=invalid date (29th Feb, non leap-year)
   integer                 :: year, month, day, utc, hour, minute
   real(kind=dp)           :: second
   integer                 :: A, Y, M, JDN
!-----------------------------------------------------------------------------------------------------------------------------------
   year   = dat(1)                        ! Year
   month  = dat(2)                        ! Month
   day    = dat(3)                        ! Day
   utc    = dat(4)*60                     ! Delta from UTC, convert from minutes to seconds
   hour   = dat(5)                        ! Hour
   minute = dat(6)                        ! Minute
   second = dat(7)-utc+dat(8)/1000.0d0    ! Second   ! correction for time zone and milliseconds
!-----------------------------------------------------------------------------------------------------------------------------------
   julian = -HUGE(99999)                  ! this is the date if an error occurs and IERR is < 0
!-----------------------------------------------------------------------------------------------------------------------------------
   if(year==0 .or. year .lt. -4713) then
      ierr=-1
      return
   endif
!-----------------------------------------------------------------------------------------------------------------------------------
!  You must compute first the number of years (Y) and months (M) since March 1st -4800 (March 1, 4801 BC)
   A=(14-month)/12 ! A will be 1 for January or Febuary, and 0 for other months, with integer truncation
   Y=year+4800-A
   M=month+12*A-3  ! M will be 0 for March and 11 for Febuary
!  All years in the BC era must be converted to astronomical years, so that 1BC is year 0, 2 BC is year "-1", etc.
!  Convert to a negative number, then increment towards zero
!  Staring from a Gregorian calendar date
   JDN=day + (153*M+2)/5 + 365*Y + Y/4 - Y/100 + Y/400 - 32045  !  with integer truncation
!  Finding the Julian date given the JDN (Julian day number) and time of day
   julian=JDN + dble(hour-12)/24.0d0 + dble(minute)/1440.0d0 + second/86400.0d0
!-----------------------------------------------------------------------------------------------------------------------------------
   if(julian.lt.0.d0) then                  ! Julian Day must be non-negative
      ierr=1
   else
      ierr=0
   endif
!-----------------------------------------------------------------------------------------------------------------------------------
end subroutine d2j
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine j2d(dat,julian,ierr)
character(len=:),parameter :: ident="@(#)j2d(3f): Converts Julian Day to date array"
integer,intent(out)        :: dat(8)
integer                    :: timezone(8), tz
real(kind=dp),intent(in)   :: julian            ! Julian Day (non-negative)
integer,intent(out)        :: ierr              ! 0 for successful execution, otherwise 1
   real(kind=dp)           :: second
   integer                 :: year
   integer                 :: month
   integer                 :: day
   integer                 :: hour
   integer                 :: minute
   integer                 :: jalpha,ja,jb,jc,jd,je,ijul

   if(julian.lt.0.d0) then                      ! Negative Julian Day not allowed
      ierr=1
      return
   else
      ierr=0
   endif
   call date_and_time(values=timezone)
   tz=timezone(4)

   ijul=idint(julian)                           ! Integral Julian Day
   second=sngl((julian-dble(ijul))*secday)      ! Seconds from beginning of Jul. Day
   second=second+(tz*60)

   if(second.ge.(secday/2.0d0)) then            ! In next calendar day
      ijul=ijul+1
      second=second-(secday/2.0d0)              ! Adjust from noon to midnight
   else                                         ! In same calendar day
      second=second+(secday/2.0d0)              ! Adjust from noon to midnight
   endif

   if(second.ge.secday) then                    ! Final check to prevent time 24:00:00
      ijul=ijul+1
      second=second-secday
   endif

   minute=int(second/60.0)                      ! Integral minutes from beginning of day
   second=second-float(minute*60)               ! Seconds from beginning of minute
   hour=minute/60                               ! Integral hours from beginning of day
   minute=minute-hour*60                        ! Integral minutes from beginning of hour

   !---------------------------------------------
   jalpha=idint((dble(ijul-1867216)-0.25d0)/36524.25d0) ! Correction for Gregorian Calendar
   ja=ijul+1+jalpha-idint(0.25d0*dble(jalpha))
   !---------------------------------------------

   jb=ja+1524
   jc=idint(6680.d0+(dble(jb-2439870)-122.1d0)/365.25d0)
   jd=365*jc+idint(0.25d0*dble(jc))
   je=idint(dble(jb-jd)/30.6001d0)
   day=jb-jd-idint(30.6001d0*dble(je))
   month=je-1

   if(month.gt.12)then
      month=month-12
   endif

   year=jc-4715
   if(month.gt.2)then
      year=year-1
   endif

   if(year.le.0)then
      year=year-1
   endif

   dat(1)=year
   dat(2)=month
   dat(3)=day
   dat(4)=tz
   dat(5)=hour
   dat(6)=minute
   dat(7)=int(second)
   dat(8)=int((second-int(second))*1000.0)
   ierr=0

end subroutine j2d
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine d2u(dat,unixtime,ierr)
character(len=:),parameter :: ident="@(#)d2u(3f): Converts date array to Unix Time (UT starts at 0000 on 1 Jan. 1970)"
integer,intent(in)         :: dat(8)                  ! date time array similar to that returned by DATE_AND_TIME
real(kind=dp),intent(out)  :: unixtime                ! Unix time (seconds)
integer,intent(out)        :: ierr                    ! return 0 on successful, otherwise 1
   real(kind=dp)           :: julian
   real(kind=dp),save      :: julian_at_epoch
   logical,save            :: first=.true.
!-----------------------------------------------------------------------------------------------------------------------------------
if(first) then                                        ! Compute zero of Unix Time in Julian days and save
   call d2j([1970,1,1,0,0,0,0,0],julian_at_epoch,ierr)
   if(ierr.ne.0) return                               ! Error
   first=.false.
endif
!-----------------------------------------------------------------------------------------------------------------------------------
   call d2j(dat,julian,ierr)
   if(ierr.ne.0) return                               ! Error
   unixtime=(julian-julian_at_epoch)*secday
end subroutine d2u
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine u2d(dat,unixtime,ierr)
! REF:JRH:1991-05-23
! REF:JSU:2015-12-12
!-----------------------------------------------------------------------------------------------------------------------------------
character(len=:),parameter :: ident="@(#)u2d(3f): Converts Unix Time to date array"
integer,intent(out)        :: dat(8)                           ! date and time array
real(kind=dp),intent(in)   :: unixtime                         ! Unix time (seconds)
integer,intent(out)        :: ierr                             ! 0 for successful execution, otherwise 1
real(kind=dp)              :: julian                           ! Unix time converted to a Julian date
real(kind=dp),save         :: Unix_Origin_as_Julian            ! start of Unix Time as Julian date
logical,save               :: first=.TRUE.
integer                    :: v(8)                             ! date and time array used to get time zone
!-----------------------------------------------------------------------------------------------------------------------------------
if(first)then                                                  ! Initialize calculated constants on first call
   call d2j([1970,1,1,0,0,0,0,0],Unix_Origin_as_Julian,ierr)   ! Compute start of Unix Time in Julian days
   if(ierr.ne.0) return                                        ! Error
   first=.FALSE.
endif
!-----------------------------------------------------------------------------------------------------------------------------------
   call date_and_time(values=v)                                ! need to get time zone
   julian=(unixtime/secday)+Unix_Origin_as_Julian              ! convert seconds from Unix Epoch to Julian days
   call j2d(dat,julian,ierr)                                   ! calculate date array from Julian date
   dat(4)=v(4)
end subroutine u2d
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
FUNCTION d2o(dat) RESULT (ordinal)
CHARACTER(LEN=:),PARAMETER :: ident="@(#)d2o(3f): Converts date-time array to Ordinal day -JSU 2015-12-13"
INTEGER,INTENT(IN)         :: dat(8)                  ! date time array similar to that returned by DATE_AND_TIME
INTEGER                    :: ordinal                 ! the returned number of days
   REAL(KIND=dp)           :: unixtime                ! Unix time (seconds)
   REAL(KIND=dp)           :: unix_first_day
   INTEGER                 :: ierr                    ! return 0 on successful, otherwise 1 from d2u(3f)
   CALL d2u(dat,unixtime,ierr)                        ! convert date to Unix Epoch Time
   IF(ierr.NE.0)THEN
      write(*,*)'*d2o* bad date array'
      ordinal=-1                                      ! initialize to bad value
   ELSE
      CALL d2u([dat(1),1,1,dat(4),0,0,0,0],unix_first_day,ierr)
      ordinal=int((unixtime-unix_first_day)/secday)+1
   ENDIF
END FUNCTION d2o
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
FUNCTION v2mo(imonth) RESULT(month_name)
CHARACTER(LEN=:),PARAMETER   :: ident="@(#)v2mo(3f): returns the month name of a Common month -JSU 2015-12-13"
CHARACTER(LEN=:),ALLOCATABLE :: month_name                                        ! string containing month name or abbreviation.
INTEGER,INTENT(IN)           :: imonth                                            ! the number of the month(1-12)
CHARACTER(LEN=:),PARAMETER   :: names(12)=[                                    &
&'January  ', 'February ', 'March    ', 'April    ', 'May      ', 'June     ', &
&'July     ', 'August   ', 'September', 'October  ', 'November ', 'December ']
   SELECT CASE(imonth)
   CASE (1:12);  month_name=TRIM(names(imonth))
   CASE DEFAULT; month_name='UNKNOWN'
   END SELECT
END FUNCTION v2mo
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
FUNCTION now(format)
CHARACTER(LEN=:),PARAMETER           :: ident="@(#)now(3f): return string representing current time given format - JSU 2015-10-24"
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: format
CHARACTER(LEN=:),ALLOCATABLE         :: now
   INTEGER                           :: values(8)
!-----------------------------------------------------------------------------------------------------------------------------------
   CALL DATE_AND_TIME(VALUES=values)
   IF(PRESENT(format))THEN
      IF(format.NE.' ')THEN
         now=fmtdate(values,format)
      ELSE
         now=fmtdate(values,'%Y-%M-%D %h:%m:%s %Z')
      ENDIF
   ELSE
      NOW=fmtdate(values,'%Y-%M-%D %h:%m:%s %Z Julian date is %J Epoch time is %E ')
   ENDIF
END FUNCTION now
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
FUNCTION fmtdate(values,format) RESULT (timestring)
! Read the FORMAT string and replace the "%" strings per the following rules:
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=:),PARAMETER :: ident="@(#)fmtdate(3f): given date array return date as string using format -JSU 2015-10-24 "
CHARACTER(LEN=*),INTENT(IN)     :: format    ! input format string
INTEGER,DIMENSION(8),INTENT(IN) :: values    ! numeric time values as DATE_AND_TIME(3f) intrinsic returns
CHARACTER(LEN=:),ALLOCATABLE    :: timestring
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   INTEGER              :: i10
   LOGICAL              :: keyword   ! flag that previous character was a % character
   CHARACTER(LEN=9)     :: day       ! day of week
   CHARACTER(LEN=1)     :: chara     ! character being looked at in format string
   CHARACTER(LEN=4096)  :: text      ! character array
   INTEGER              :: iout
   INTEGER              :: weekday
   INTEGER              :: ierr
   INTEGER,SAVE         :: called=0
   LOGICAL,SAVE         :: since=.FALSE.
   REAL(KIND=dp)        :: julian
   REAL(KIND=dp)        :: cputime
   INTEGER              :: ii
   REAL(KIND=dp)        :: unixtime
   REAL(KIND=dp),save   :: unixtime_last
   INTEGER              :: systemclock, countrate
   INTEGER              :: iso_year, iso_week, iso_weekday
   CHARACTER(LEN=10)    :: iso_name
   CHARACTER(LEN=2)     :: dayend

   text=' '
!  write string, when encounter a percent character do a substitution
   keyword=.FALSE.
   iout=1
   DO i10=1,LEN(format)
      chara=format(i10:i10)
      IF(chara.eq.'%'.and..not.keyword)THEN
            keyword=.TRUE.
            CYCLE
      ENDIF
      IF(keyword)THEN
         keyword=.FALSE.
         SELECT CASE(chara)
         !=====================================================================================
         CASE('%'); WRITE(text(iout:),'(A1)')chara                        ! literal percent character
         !=====================================================================================
         CASE('b'); WRITE(text(iout:),'(A1)')' '                          ! space character
         !=====================================================================================
         CASE('c'); CALL cpu_time(cputime)                                ! CPU_TIME()
                    WRITE(text(iout:),'(G0)')cputime
         !=====================================================================================
         CASE('C'); called = called + 1                                   ! number of times this routine called
                    WRITE(text(iout:),'(I0)')called
         !=====================================================================================
         CASE('d');                                                       ! the day of the month 1st..31st
                    dayend='  '
                    select case(values(3))
                    case(1,21,31); dayend='st'
                    case(2,22); dayend='nd'
                    case(3,23); dayend='rd'
                    case(4:20,24:30); dayend='th'
                    case default
                    end select
                    WRITE(text(iout:),'(I2,a)')values(3),dayend
         !=====================================================================================
         CASE('D'); WRITE(text(iout:),'(I2.2)')values(3)                  ! the day of the month 1..31
         !=====================================================================================
         CASE('e'); CALL d2u(values,unixtime,ierr)                        ! integer Unix Epoch time in seconds
                    WRITE(text(iout:),'(G0)')int(unixtime)
         !=====================================================================================
         CASE('E'); CALL d2u(values,unixtime,ierr)                        ! Unix Epoch time in seconds
                    WRITE(text(iout:),'(G0)')unixtime
         !=====================================================================================
         CASE('h'); WRITE(text(iout:),'(I2.2)')values(5)                  ! the hour of the day, in the range of 0 to 23
         !=====================================================================================
         CASE('H'); ii=mod(values(5),12)                                  ! hour of day in range 1..12
                    if(ii.eq.0)then
                       ii=12
                    endif
                    WRITE(text(iout:),'(I2.2)')ii
         !=====================================================================================
         CASE('i'); CALL woy(values,iso_year,iso_week,iso_weekday,iso_name) ! ISO week of year
                    WRITE(text(iout:),'(I0)')iso_week
         !=====================================================================================
         CASE('I'); CALL woy(values,iso_year,iso_week,iso_weekday,iso_name) ! iso-8601 Week-numbering year date
                    WRITE(text(iout:),'(a)')iso_name
         !=====================================================================================
         CASE('j'); CALL d2j(values,julian,ierr)                          ! integer Julian date (truncated to integer)
                    WRITE(text(iout:),'(I0)')int(julian)
         !=====================================================================================
         CASE('J'); CALL d2j(values,julian,ierr)                          ! Julian date to milliseconds
                    WRITE(text(iout:),'(I0,".",i3.3)')int(julian),int((julian-int(julian))*1000.0)
         !=====================================================================================
         CASE('k'); call system_clock(count=systemclock,count_rate=countrate)  ! systemclock/countrate
                    WRITE(text(iout:),'(G0)')real(systemclock)/countrate
         !=====================================================================================
         CASE('l'); WRITE(text(iout:),'(A3)')v2mo(values(2))              ! three characters of the name of the month of the year
         !=====================================================================================
         CASE('L'); WRITE(text(iout:),'(A)')v2mo(values(2))               ! name of the month of the year
         !=====================================================================================
         CASE('m'); WRITE(text(iout:),'(I2.2)')values(6)                  ! the minutes of the hour, in the range 0 to 59
         !=====================================================================================
         CASE('M'); WRITE(text(iout:),'(I2.2)')values(2)                  ! month of year (1..12)
         !=====================================================================================
         CASE('N'); if( values(5).ge.12)then                              ! AM||PM
                       WRITE(text(iout:),'("PM")')
                    else
                       WRITE(text(iout:),'("AM")')
                    endif
         !=====================================================================================
         CASE('O'); WRITE(text(iout:),'(I3.3)')d2o(values)                ! Ordinal day of year
         !=====================================================================================
         CASE('s'); WRITE(text(iout:),'(I2.2)')values(7)                  ! the seconds of the minute, in the range 0 to 60
         !=====================================================================================
         CASE('S'); IF(.NOT.since)THEN                                    ! seconds since last called
                       since=.TRUE.
                       CALL d2u(values,unixtime_last,ierr)
                    ENDIF
                    CALL d2u(values,unixtime,ierr)
                    WRITE(text(iout:),'(G0)')unixtime-unixtime_last
                    unixtime_last=unixtime
         !=====================================================================================
         CASE('t'); WRITE(text(iout:),'(A1)')CHAR(9)                      ! tab character
         !=====================================================================================
         CASE('U'); CALL dow(values,weekday,day,ierr)                     ! Return the day of the week, 1..7 Sunday=1
                    WRITE(text(iout:),'(I1)')weekday
         !=====================================================================================
         CASE('u'); CALL dow(values,weekday,day,ierr)                     ! Return the day of the week, 1..7 Monday=1
                    WRITE(text(iout:),'(I1)')mod(weekday+5,7)+1
         !=====================================================================================
         CASE('W'); CALL dow(values,weekday,day,ierr)                     ! Return the name of the day of the week
                    WRITE(text(iout:),'(a)')day
         !=====================================================================================
         CASE('w'); CALL dow(values,weekday,day,ierr)                     ! Return the first three characters of the day of the week
                    WRITE(text(iout:),'(A3)')day(1:3)
         !=====================================================================================
         CASE('x'); WRITE(text(iout:),'(I3.3)')values(8)                  ! the milliseconds of the second, in the range 0 to 999
         !=====================================================================================
         CASE('Y'); WRITE(text(iout:),'(I4.4)')values(1)                  ! the year, including the century (for example, 1990)
         !=====================================================================================
         CASE('Z'); WRITE(text(iout:),'(SP,I5.4)')values(4)               ! time difference with respect to UTC in minutes
         !=====================================================================================
         CASE('z'); WRITE(text(iout:),'(I3.2,":",I2.2)')int(values(4)/60),abs(mod(values(4),60)) ! time from UTC as +-hh:mm
         !=====================================================================================
         CASE DEFAULT
            WRITE(text(iout:),'(A1)')chara
         !=====================================================================================
         END SELECT
         !=====================================================================================
         iout=len_trim(text)+1
      ELSE
         WRITE(text(iout:),'(A1)')chara;iout=iout+1
      ENDIF
   ENDDO
   timestring=trim(text)
END FUNCTION fmtdate
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine fmtdate_usage(ii)
character(len=:),parameter :: ident="@(#)fmtdate_usage(3f): display macros recognized by fmtdate(3f) - JSU 2015-10-24"
character(len=51),allocatable :: usage(:)
integer                       :: i,ii
character(len=ii)             :: blanks
usage=[ &                                               !date(1) COMMAND
&' Base time array:                                  ',&
&' (1) %Y -- year, yyyy                              ',&
&' (2) %M -- month of year, 01 to 12                 ',&
&' (3) %D -- day of month, 01 to 31                  ',&
&'     %d -- day of month, with suffix (1st, 2nd,...)',&
&' (4) %Z -- minutes from UTC                        ',&
&'     %z -- -+hh:mm from UTC                        ',&
&' (5) %h -- hours, 00 to 23                         ',&
&'     %H -- hour (1 to 12, or twelve-hour clock)    ',&
&'     %N -- AM (before noon) PM (>=after noon)      ',&
&' (6) %m -- minutes, 00 to 59                       ',&
&' (7) %s -- sec, 00 to 60                           ',&
&' (8) %x -- milliseconds 000 to 999                 ',&
&'Conversions                                        ',&
&'     %E -- Unix Epoch time                         ',&
&'     %e -- integer value of Unix Epoch time        ',&
&'     %J -- Julian  date                            ',&
&'     %j -- integer value of Julian date            ',&
&'     %O -- Ordinal day (day of year)               ',&
&'     %U -- day of week, 1..7 Sunday=1              ',&
&'     %u -- day of week, 1..7 Monday=1              ',&
&'     %i -- ISO week of year 1..53                  ',&
&'     %I -- iso-8601 week-numbering date(yyyy-Www-d)',&
&' Names                                             ',&
&'     %l -- abbreviated month name                  ',&
&'     %L -- full month name                         ',&
&'     %w -- first three characters of weekday       ',&
&'     %W -- weekday name                            ',&
&' Literals                                          ',&
&'     %% -- a literal %%                            ',&
&'     %t -- tab character                           ',&
&'     %b -- blank character                         ',&
&' Program timing:                                   ',&
&'     %c -- CPU_TIME(3f) output                     ',&
&'     %C -- number of times this routine is used    ',&
&'     %k -- time in seconds from SYSTEM_CLOCK(3f)   ',&
&'                                                   ']
   blanks=' '
   WRITE(*,'(a,a)')(blanks(:ii),usage(i),i=1,size(usage))
end subroutine fmtdate_usage
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine dow(values, weekday, day, ierr)
character(len=:),parameter :: ident="@(#)dow(3f): Return the day of the week"
real(kind=dp)                      :: julian    ! the julian day for which the weekday is required,
integer,intent(in)                 :: values(8) ! date and time array used to get time zone
integer,intent(out),optional       :: weekday   ! The day of the week, 1 = Sunday
character*(*),intent(out),optional :: day       ! The name of the day of the week, e.g. 'Sunday'. Minimum length = 9
integer,intent(out)                :: ierr      ! Error return,0=correct,-1=invalid Julian day,-2=neither day nor weekday specified
   integer                         :: iweekday

   call d2j(values,julian,ierr)                 ! need julian date to calculate day of week for first day of month
   ierr = 0

   if(julian < 0) then
      ierr = -1
      return
   endif

   if(.not.present(day).and. .not.present(weekday)) then
      ierr=-2
      return
   endif

   ! julian day starts at noon so add 1/2 day
   ! add time zone
   iweekday = mod(int((julian+dble(values(4)/60.0d0/24.0d0)+0.5d0)+1.0d0), 7)
   iweekday = iweekday +1

   if(present(day)) then
      select case(iweekday)
      case(1)     ;day = 'Sunday'
      case(2)     ;day = 'Monday'
      case(3)     ;day = 'Tuesday'
      case(4)     ;day = 'Wednesday'
      case(5)     ;day = 'Thursday'
      case(6)     ;day = 'Friday'
      case(7)     ;day = 'Saturday'
      case default;day = 'E-R-R-O-R'
      end select
   endif

   if(present(weekday))then
      weekday=iweekday
   endif

end subroutine dow
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
subroutine woy(dat,iso_year,iso_week,iso_weekday,iso_name)
!-----------------------------------------------------------------------------------------------------------------------------------
!  The ISO-8601 date and time standard was issued by the International Organization for Standardization (ISO).
!  It is used (mainly) in government and business for fiscal years, as well as in timekeeping.
!  The system specifies a week year atop the Gregorian calendar by defining a notation for ordinal weeks of the year.
!
!  An ISO week-numbering year (also called ISO year informally) has 52 or 53 full weeks.
!  That is 364 or 371 days instead of the usual 365 or 366 days.
!  The extra week is referred to here as a leap week, although ISO 8601 does not use this term.
!  Weeks start with Monday.
!  The first week of a year is the week that contains the first Thursday of the year (and, hence, always contains 4 January).
!  ISO week year numbering therefore slightly deviates from the Gregorian for some days close to 1 January.
!-----------------------------------------------------------------------------------------------------------------------------------
!CALCULATION:
!  The ISO-8601 week number of any date can be calculated, given its ordinal date (i.e. position within the year)
!  and its day of the week.

!METHOD:
!   Using ISO weekday numbers (running from 1 for Monday to 7 for Sunday),
!   subtract the weekday from the ordinal date, then add 10. Divide the result
!   by 7. Ignore the remainder; the quotient equals the week number. If
!   the week number thus obtained equals 0, it means that the given date
!   belongs to the preceding (week-based) year. If a week number of 53 is
!   obtained, one must check that the date is not actually in week 1 of the
!   following year.
! These two statements are assumed true when correcting the dates around January 1st ...
!   o  The number of weeks in a given year is equal to the corresponding week number of 28 December.
!   o  January 4th is always in the first week.
!
!ISO_NAME:
!  Week date representations are in the format YYYYWww-D.
!  o [YYYY] indicates the ISO week-numbering year which is slightly different from the traditional Gregorian calendar year.
!  o [Www] is the week number prefixed by the letter W, from W01 through W53.
!  o [D] is the weekday number, from 1 through 7, beginning with Monday and ending with Sunday.
!
!  For example, the Gregorian date 31 December 2006 corresponds to the Sunday of the 52nd week of 2006, and is written
!     2006-W52-7 (extended form)
!  or 2006W527 (compact form).
!
!REFERENCE:
!  From Wikipedia, the free encyclopedia 2015-12-19
!AUTHOR:
!  John S. Urban, 2015-12-19
!-----------------------------------------------------------------------------------------------------------------------------------
character(len=:),parameter      :: ident="@(#)woy(3f): Calculate iso-8601 Week-numbering year date yyyy-Www-d"
integer,parameter               :: dp=kind(0.0d0)
integer,intent(in)              :: dat(8)     ! input date array
integer,intent(out)             :: iso_year, iso_week, iso_weekday
character(len=10),intent(out)   :: iso_name
integer                         :: shared_weekday
integer                         :: last_week_this_year
integer                         :: dec28_lastyear(8)   ! December 28th is always in last week
integer                         :: dec28_thisyear(8)   ! December 28th is always in last week
character(len=9)                :: day
integer                         :: ierr
   iso_year=dat(1)                                               ! initially assume the iso_year is the same as the data array year
   iso_week=uncorrected_week_of_year(dat)                        ! this is the week number unless around January 1st
   iso_weekday=shared_weekday                                    ! this is the number of the day of the week assuming Monday=1
   dec28_thisyear=[dat(1),12,28,dat(4),0,0,0,0]                  ! Dec 28th is always in last week; use this to get number of weeks
   last_week_this_year=uncorrected_week_of_year(dec28_thisyear)  ! get the number of the last week of the year (52 or 53)
   ! correct dates around January 1st
   if(iso_week  < 1)then                                         ! if week < 1 then week = lastWeek(year -1)
      dec28_lastyear=[dat(1)-1,12,28,dat(4),0,0,0,0]             ! Dec 28th is always in last week, we want its week number
      iso_week=uncorrected_week_of_year(dec28_lastyear)          ! got the week number for the last week of last year (52 or 53)
      iso_year=dat(1)-1                                          ! our date belongs to last year
   elseif(iso_week >last_week_this_year)then                     ! if week > lastweek(year) then week = 1
      iso_week=iso_week-last_week_this_year                      ! our date belongs to next year
      iso_year=dat(1)+1
   endif

   write(iso_name,'(i4.4,"-W",i2.2,"-",i1)')iso_year,iso_week,iso_weekday ! create ISO string designation for our date

contains
   function uncorrected_week_of_year(datin)
   implicit none
   integer            :: uncorrected_week_of_year
   integer,intent(in) :: datin(8)
      integer         :: ordinal
      call dow(datin,shared_weekday,day,ierr)                 ! formula needs day of week 1..7 where Monday=1
      shared_weekday=mod(shared_weekday+5,7)+1                ! change from Sunday=1 to Monday=1
      ordinal=d2o(datin)                                      ! formula needs ordinal day of year where Jan 1st=1
      uncorrected_week_of_year=(ordinal-shared_weekday+10)/7
   end function uncorrected_week_of_year

end subroutine woy
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
function dj(dat) result (julian)
character(len=:),parameter :: ident="@(#)dj(3f): Given date array returns Julian Day"
real(kind=dp)              :: julian
integer,intent(in)         :: dat(8)
   integer                 :: ierr
call d2j(dat,julian,ierr)
end function dj

function jd(julian) result (dat)
character(len=:),parameter :: ident="@(#)jd(3f): Given Julian Day returns date array"
real(kind=dp),intent(in)   :: julian
integer                    :: dat(8)
   integer                 :: ierr
call j2d(dat,julian,ierr)
end function jd

function du(dat) result (unixtime)
character(len=:),parameter :: ident="@(#)du(3f): Given date array returns Unix Epoch time "
real(kind=dp)              :: unixtime
integer,intent(in)         :: dat(8)
   integer                 :: ierr
call d2u(dat,unixtime,ierr)
end function du

function ud(unixtime) result (dat)
character(len=:),parameter :: ident="@(#)ud(3f): Given Unix Epoch Time returns date array"
real(kind=dp),intent(in)   :: unixtime
integer                    :: dat(8)
   integer                 :: ierr
call u2d(dat,unixtime,ierr)
end function ud
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================
!
!   XXXX
!  X    X
! X
! X
! X
! X
! X
!  X    X
!   XXXX
!
subroutine sys_sleep(wait_seconds)
use, intrinsic  :: iso_c_binding, only: c_int
character(len=:),parameter :: ident="@(#)sys_sleep(3f): call sleep(3c)"
integer (c_int) :: wait_seconds, how_long
interface
      function c_sleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: c_sleep !  should be unsigned int (not available in Fortran).  OK until highest bit gets set.
          integer (c_int), intent (in), VALUE :: seconds
      end function c_sleep
end interface
   if(wait_seconds.gt.0)then
      how_long = c_sleep(wait_seconds)
   endif
end subroutine sys_sleep
!===================================================================================================================================
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!===================================================================================================================================

!     PORPOSE
!         Determine data time interval & check timesteps

!     RESTRICTIONS
!         Input values must be floating points.

!     LITERATURE
!         None

!     HISTORY
!         Written,  Matthias Zink, Oct 2012
!         Modified  Matthias Cuntz & Juliane Mai Nov 2014 - time int or double
!                   Stephan Thober Sep 2015 - added read for hourly data

subroutine get_time(fName, vName, julStart, julEnd, nctimestep)
  !
  use mo_kind,         only: i4, dp
  use mo_julian,       only: date2dec
  use mo_message,      only: message
  use mo_NcRead,       only: Get_NcVar, Get_NcDim, Get_NcVarAtt
  use mo_utils,        only: ne
  use mo_string_utils, only: DIVIDE_STRING
  use netcdf,          only: NF90_INT, NF90_DOUBLE, NF90_NOWRITE
  use netcdf,          only: nf90_open, nf90_close, nf90_inq_varid, nf90_inquire_variable
  !
  implicit none
  !
  character(len=*)            , intent(in)  :: fName               ! name of NetCDF file
  character(len=*)            , intent(in)  :: vName               ! name of variable
  integer(i4)                 , intent(out) :: julStart
  integer(i4)                 , intent(out) :: julEnd
  integer(i4),        optional, intent(in)  :: nctimestep          ! -1: daily (default); -2:monthly; -3:yearly
  !
  integer(i4)                               :: i
  integer(i4)                               :: yRef, dRef, mRef    ! reference time of NetCDF (unit attribute of
  integer(i4)                               :: datatype            ! datatype of attribute
  integer(i4),    dimension(5)              :: dimen
  !
  integer(i4),   dimension(:), allocatable  :: timesteps           ! time variable of NetCDF
  real(dp),      dimension(:), allocatable  :: itimesteps          ! time variable of NetCDF
  !
  character(256)                            :: AttValues           ! netcdf attribute values
  character(256), dimension(:), allocatable :: strArr              ! dummy for netcdf attribute handling
  character(256), dimension(:), allocatable :: date                ! dummy for netcdf attribute handling date
  real(dp)                                  :: jday_frac           ! julian day from dec2date
  integer(i4)                               :: inctimestep         ! local nctimestep
  integer(i4) :: ncid    ! id of input stream
  integer(i4) :: varid   ! id of variable to be read
  integer(i4) :: status  ! netcdf inquire return
  real(dp) :: deltaT ! diff between single time steps in NetCDF in days
  real(dp) :: nTStepDay ! number of timesteps per day
  !
  ! check optional nctimestep
  inctimestep = -1 !ST: to do inctimestep should be determined from file
  nTStepDay = -1.
  if (present(nctimestep)) inctimestep = nctimestep

  dimen = Get_NCDim(fName, trim(vName))
  ! get unit attribute of variable 'time'
  call Get_NcVarAtt(fName, 'time', 'units', AttValues, dtype=datatype)
  ! AttValues looks like "<unit> since YYYY-MM-DD HH:MM:SS"
  call DIVIDE_STRING(trim(AttValues), ' ', strArr)
  !
  ! determine reference time and convert to integer
  call DIVIDE_STRING(trim(strArr(3)), '-', date)
  read(date(1),*) yRef
  read(date(2),*) mRef
  read(date(3),*) dRef
  jday_frac = date2dec(dd=dRef, mm=mRef, yy=yRef)
  !
  status = nf90_open(fName, NF90_NOWRITE, ncid)
  status = nf90_inq_varid(ncid, 'time', varid)
  status = nf90_inquire_variable(ncid, varid, xtype=datatype)
  status = nf90_close(ncid)
  status = status + 0  ! this is only to make the variable used
  if ((datatype .eq. NF90_INT) .or. (datatype .eq. NF90_DOUBLE)) then
     if (datatype .eq. NF90_INT) then
        allocate(timesteps(dimen(3)))
        call Get_NcVar(fName, 'time', timesteps)
     else
        allocate(itimesteps(dimen(3)))
        allocate(timesteps(dimen(3)))
        call Get_NcVar(fName, 'time', itimesteps)
        timesteps = nint(itimesteps, i4)
     endif
  else
     call message('***ERROR: data type of time must be NF90_INT or NF90_DOUBLE in netcdf file.')
     stop
  end if
  !
  ! strArr(1) is <unit>
  if (strArr(1) .EQ. 'days') then
     nTStepDay = 1
  else if (strArr(1) .eq. 'hours') then
     nTStepDay = 24
  else
     call message('***ERROR: Please provide the input data in (days or hours) since YYYY-MM-DD HH:MM:SS in ', trim(vName))
     stop
  end if
  !
  ! check consistency of timesteps
  select case(inctimestep)
  case(-1) ! daily
     do i = 2, dimen(3)
        deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
        ! deltaT has to be one but with integer conversion
        if ( ne(deltaT, 1._dp) ) then
           call message('***ERROR: ', trim(vName),' must have daily time steps.')
           stop
        end if
     end do
  case(-2) ! monthly
     do i = 2, dimen(3)
        deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
        ! deltaT has to be one but with integer conversion
        if (( deltaT .lt. 28._dp) .or. ( deltaT .gt. 31._dp)) then
           call message('***ERROR: ', trim(vName),' must have monthly time steps.')
           stop
        end if
     end do
  case(-3) ! yearly
     do i = 2, dimen(3)
        deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
        ! deltaT has to be one but with integer conversion
        if ( deltaT .lt. 360._dp) then
           call message('***ERROR: ', trim(vName),' must have yearly time steps.')
           stop
        end if
     end do
  case(-4) ! hourly
     do i = 2, dimen(3)
        deltaT = real(timesteps(i) - timesteps(i-1), dp)
        ! deltaT has to be one but with integer conversion
        if ( ne(deltaT, 1._dp) ) then
           call message('***ERROR: ', trim(vName),' must have hourly time steps.')
           stop
        end if
     end do
  case default
     call message('***ERROR: get_time: unknown nctimestep switch.')
     stop
  end select
  !
  ! determine starting and ending julian day of the dataset
  i  = nint(jday_frac, i4 )
  ! calculate days
  julStart = i + timesteps(1) / int(nTStepDay, i4)
  julEnd   = i + timesteps(dimen(3)) / int(nTStepDay, i4)
  !
end subroutine get_time

end module mo_time
