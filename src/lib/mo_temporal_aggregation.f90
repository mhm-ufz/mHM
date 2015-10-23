!> \file mo_temporal_aggregation.f90

!> \brief Temporal aggregation for time series (averaging)

!> \details This module does temporal aggregation (averaging) of time series

!> \authors Oldrich Rakovec, Rohini Kumar
!> \date October 2015

MODULE mo_temporal_aggregation

  ! This module calculates does temporal aggregation of ot time series

  ! Written  Oldrich Rakovec, October 2015

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

  ! Copyright 2015 Oldrich Rakovec, Rohini Kumar

  use mo_kind,   ONLY: i4, dp
  use mo_julian, ONLY: julday, dec2date
  use mo_constants,    ONLY: eps_dp


  IMPLICIT NONE

  PUBLIC :: day2mon_average    ! converts daily time series to monthly
  PUBLIC :: hour2day_average   ! converts hourly time series to daily 

  ! ------------------------------------------------------------------
  !     NAME
  !         day2mon_average

  !     PURPOSE
  !         Calculates monthly average values from daily values
  !
  !>        \brief Day-to-month average (day2mon_average)
  !
  !>        \details converts daily time series to monthly
  !
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: daily_data(:)"   array of daily time series
  !>        \param[in] "integer(i4) :: year"            year of the starting time
  !>        \param[in] "integer(i4) :: month"           month of the starting time
  !>        \param[in] "integer(i4) :: day"             day of the starting time
  !
  !     INTENT(INOUT)
  !>        \param[in] "real(sp/dp) :: mon_average(:)"  array of monthly averaged values
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: misval"          missing value definition
  !>        \param[in] "logical     :: rm_misval"       switch to exclude missing values 
  !
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note
  !
  !     EXAMPLE
  !         -> see example in test directory

  !     LITERATURE
  !
  
  !     HISTORY
  !>        \author Oldrich Rakovec, Rohini Kumar
  !>        \date Oct 2015
  !         Modified,

   INTERFACE day2mon_average
      MODULE PROCEDURE day2mon_average_dp
   END INTERFACE day2mon_average

  ! ------------------------------------------------------------------
  !     NAME
  !         hour2day_average

  !     PURPOSE
  !         Calculates daily average values from hourly values
  !
  !>        \brief Hour-to-day average (hour2day_average)
  !
  !>        \details converts hourly time series to daily
  !
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: hourly_data(:)"  array of hourly time series
  !>        \param[in] "integer(i4) :: year"            year of the starting time
  !>        \param[in] "integer(i4) :: month"           month of the starting time
  !>        \param[in] "integer(i4) :: day"             day of the starting time
  !>        \param[in] "integer(i4) :: hour"            hour of the starting time
  !
  !     INTENT(INOUT)
  !>        \param[in] "real(sp/dp) :: day_average(:)"  array of daily averaged values
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: misval"          missing value definition
  !>        \param[in] "logical     :: rm_misval"       switch to exclude missing values 
  !
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note Hours values should be from 0 to 23 (NOT from 1 to 24!)
  !
  !     EXAMPLE
  !         -> see example in test directory

  !     LITERATURE
  !
  
  !     HISTORY
  !>        \author Oldrich Rakovec, Rohini Kumar
  !>        \date Oct 2015
  !         Modified,

   INTERFACE hour2day_average
      MODULE PROCEDURE hour2day_average_dp
   END INTERFACE hour2day_average

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  SUBROUTINE day2mon_average_dp(daily_data,yearS,monthS,dayS,mon_avg, misval, rm_misval )

    IMPLICIT NONE

    REAL(dp), dimension(:),                         INTENT(IN)  :: daily_data      ! array of daily data
    INTEGER(i4),                                    INTENT(IN)  :: yearS           ! year of the initial time step
    INTEGER(i4),                                    INTENT(IN)  :: monthS          ! month of the initial time step
    INTEGER(i4),                                    INTENT(IN)  :: dayS            ! day of the initial time step

    REAL(dp), dimension(:), allocatable,          INTENT(INOUT) :: mon_avg         ! array of the monthly averages

    REAL(dp),                             optional, INTENT(IN)  :: misval          ! missing value definition
    logical,                              optional, INTENT(IN)  :: rm_misval       ! switch to remove missing values

    ! local variables 
    INTEGER(i4)                                                 :: ndays, tt, kk      ! number of days, indices
    INTEGER(i4)                                                 :: start_day, end_day ! size of input array, size of days  
    INTEGER(i4)                                                 :: y, m                
    INTEGER(i4)                                                 :: year, month, day    ! variables for date
    INTEGER(i4)                                                 :: yearE, monthE, dayE ! vatiables for End date    
    REAL(dp)                                                    :: newTime 

    REAL(dp), dimension(:,:), allocatable                       :: nCounter_m       ! counter number of days in months (w/ data)
    REAL(dp), dimension(:,:), allocatable                       :: nCounter_m_full  ! counter number of days in months (complete) 
    REAL(dp), dimension(:,:), allocatable                       :: mon_sum          ! monthly sum

    INTEGER(i4)                                                 :: nmonths     ! number of days, number of months
    LOGICAL                                                     :: remove      ! switch for considering missing data
    REAL(dp)                                                    :: missing  ! switch for reading missing value or default -9999.

    if (present(misval)) then
       missing = misval
    else
       missing = -9999._dp
    endif   
   
    if (present(rm_misval)) then
       remove = rm_misval
    else
       remove = .FALSE.
    endif   
    
    ! get total number of days 
    ndays   = SIZE(daily_data)
  
    ! assign initial julian day
    start_day = julday(dayS,monthS,yearS) 

    ! calculate last julian day
    end_day = start_day + ndays - 1_i4

    ! get year, month and day of the end date:
    call dec2date( real(end_day, dp), yy=yearE, mm=monthE, dd=dayE)
 
    ! get number of days with data for each month
    allocate( nCounter_m(yearS:yearE,12) )
    allocate( nCounter_m_full(yearS:yearE,12) )
    allocate(mon_sum(yearS:yearE,12))
    nCounter_m(:,:)      = 0
    nCounter_m_full(:,:) = 0
    mon_sum(:,:)         = 0.0_dp

    newTime = real(start_day, dp)
    ! calculate monthly sums
    do tt = 1, (end_day - start_day + 1)          
       call dec2date( (newTime + tt - 1), yy=year, mm=month, dd=day)
       nCounter_m_full(year,month) = nCounter_m_full(year,month) + 1.0_dp
       if ( abs( daily_data(tt) - missing ) .lt. eps_dp ) cycle
       mon_sum(year,month) = mon_sum(year,month) + daily_data(tt)
       nCounter_m(year,month) = nCounter_m(year,month) + 1.0_dp
    end do

    ! calculate number of months
    nmonths = 0
    do y = yearS, yearE
       do m = 1, 12
          if ( (y .EQ. yearS) .AND. ( m .LT. monthS ) ) cycle
          if ( (y .EQ. yearE) .AND. ( m .GT. monthE ) ) cycle         
          nmonths = nmonths + 1
       end do
    end do
   
    ! calculate monthly averages
    allocate(mon_avg(nmonths))
    mon_avg(:) = missing
    kk = 0
    do y = yearS, yearE
       do m = 1, 12
          if ( (y .EQ. yearS) .AND. ( m .LT. monthS ) ) cycle
          if ( (y .EQ. yearE) .AND. ( m .GT. monthE ) ) cycle         
          kk = kk + 1
          if ( ( nCounter_m(y,m) .GT. 0 ) .AND. &
               ( abs( nCounter_m_full(y,m) - nCounter_m(y,m)) .LT. eps_dp ) ) then
             mon_avg(kk) = mon_sum(y,m) / real( nCounter_m(y,m), dp)
          else if ( ( nCounter_m(y,m) .GT. 0 ) .AND. remove) then
             mon_avg(kk) = mon_sum(y,m) / real( nCounter_m(y,m), dp)
          endif    
       end do
    end do
    
    deallocate(nCounter_m_full)
    deallocate(nCounter_m)
    deallocate(mon_sum)
    
  END SUBROUTINE day2mon_average_dp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE hour2day_average_dp(hourly_data,yearS,monthS,dayS,hourS,day_avg, misval, rm_misval )

    IMPLICIT NONE

    REAL(dp), dimension(:),                         INTENT(IN)  :: hourly_data     ! array of hourly data
    INTEGER(i4),                                    INTENT(IN)  :: yearS           ! year of the initial time step
    INTEGER(i4),                                    INTENT(IN)  :: monthS          ! month of the initial time step
    INTEGER(i4),                                    INTENT(IN)  :: dayS            ! day of the initial time step
    INTEGER(i4),                                    INTENT(IN)  :: hourS           ! hour of the initial time step   

    REAL(dp), dimension(:), allocatable,          INTENT(INOUT) :: day_avg         ! array of the daily averages

    REAL(dp),                             optional, INTENT(IN)  :: misval          ! missing value definition
    logical,                              optional, INTENT(IN)  :: rm_misval       ! switch to remove missing values

    ! local variables 
    INTEGER(i4)                                                 :: nhours, ndays_dummy, tt, dd, kk 
    REAL(dp)                                                    :: start_day, end_day   ! assign julian values
    INTEGER(i4)                                                 :: yearE, monthE, dayE, hourE, hourEd ! End dates, incl. Dummy

    REAL(dp), dimension(:), allocatable                         :: nCounter_h       ! counter number of hours in day (w/ data)
    REAL(dp), dimension(:), allocatable                         :: nCounter_h_full  ! counter number of hours in day (complete) 
    REAL(dp), dimension(:), allocatable                         :: day_sum          ! daily sum
    
    LOGICAL                                                     :: remove   ! switch for considering missing data
    REAL(dp)                                                    :: missing  ! switch for reading missing value or default -9999.

    if (present(misval)) then
       missing = misval
    else
       missing = -9999._dp
    endif   

    if (present(rm_misval)) then
       remove = rm_misval
    else
       remove = .FALSE.
    endif   
   
    ! get total number of hours 
    nhours   = SIZE(hourly_data)
    ! assign initial julian day
    start_day = julday(dayS,monthS,yearS) - 0.5_dp + real(hourS,dp)/24._dp
    
    ! calculate last julian day
    end_day = start_day + real( nhours - 1._dp ,dp )/ 24._dp 

    ! get year, month and day of the end date 
    call dec2date( end_day , yy=yearE, mm=monthE, dd=dayE, hh=hourE)

    ! get largerst possible number of calendar days
    ndays_dummy = ceiling( real(nhours,dp) / 24._dp + 2._dp ) 
   
    allocate( day_sum(ndays_dummy))
    allocate( nCounter_h(ndays_dummy) )
    allocate( nCounter_h_full(ndays_dummy) )
    day_sum(:)         = 0.0_dp
    nCounter_h(:)      = 0
    nCounter_h_full(:) = 0

    ! calculate daily sums
    dd = 1
    do tt = 1, nhours
       call dec2date( start_day + real(tt-1,dp)/24._dp , hh=hourEd)
       nCounter_h_full(dd) = nCounter_h_full(dd) + 1
       if ( abs( hourly_data(tt) - missing ) .lt. eps_dp ) then
          day_sum(dd) = day_sum(dd)          
       else
          day_sum(dd) = day_sum(dd) + hourly_data(tt)
          nCounter_h(dd) = nCounter_h(dd) + 1
       endif  
       if ( (hourEd .EQ. 23) .AND. (tt .LT. nhours)) dd = dd + 1    
    end do

    ! dd is the total number of calendar days, between hourS and hourE
    allocate(day_avg(dd))
    day_avg(:) = missing

    ! calculate daily average
    do kk = 1,dd      
       if ( ( nCounter_h(kk) .GT. 0 ) .AND. &
            ( abs( nCounter_h_full(kk) - nCounter_h(kk)) .LT. eps_dp ) ) then
          day_avg(kk) = day_sum(kk) / real( nCounter_h(kk), dp)
       else if ( ( nCounter_h(kk) .GT. 0 ) .AND. remove) then
          day_avg(kk) = day_sum(kk) / real( nCounter_h(kk), dp)
       endif
    end do   
    
    deallocate(nCounter_h_full)
    deallocate(nCounter_h)
    deallocate(day_sum)
    
  END SUBROUTINE hour2day_average_dp
   
END MODULE mo_temporal_aggregation
