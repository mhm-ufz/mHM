!> \file mo_signatures.f90

!> \brief Module with calculations for several hydrological signatures.

!> \details This module contains calculations for hydrological signatures.
!> It contains:\n
!> * Mean annual runoff
!> * Autocorrelation
!> * Autocorrelation for low flows
!> * Autocorrelation for high flows
!> * Rising and declining lomb densities
!> * Flow duration curves
!> * Flow duration curves low flows
!> * Flow duration curves high flows
!> * Peak distribution parameter
!> * Peaks for low flows
!> * Peaks for high flows

!> \authors Remko Nijzink,
!> \date March 2014
!  modified, Stephan Thober, Jan 2015 - uncommented unused subroutines that are incompatible
!                                       with multiple time periods

MODULE mo_signatures

  ! Module with calculations for several hydrological signatures.

  ! Written  Remko Nijzink, March 2014
  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: Autocorrelation         ! Autocorrelation function
  !PUBLIC :: Autocorrelation_low     ! Autocorrelation for low flows
  !PUBLIC :: Autocorrelation_high    ! Autocorrelation for high flows
  PUBLIC :: FlowDurationCurves      ! Flow duration curves
  !PUBLIC :: FlowDurationCurves_high ! Flow duration curves high flows
  !PUBLIC :: FlowDurationCurves_low  ! Flow duration curves low flows
  PUBLIC :: Limb_densities          ! Rising and declining limb densities
  PUBLIC :: PeakDistribution        ! Peak distribution parameter
  !PUBLIC :: PeakDistribution_high   ! Peaks for high flows
  !PUBLIC :: PeakDistribution_low    ! Peaks for low flows
  !PUBLIC :: Q_MeanAnnual            ! Autocorrelation function


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  !     NAME
  !         Q_MeanAnnual

  !     PURPOSE
  !>        \brief The mean annual runoff

  !>        \details Calculates the mean annual runoff

  !     CALLING SEQUENCE
  !         call Q_Meannual(Q_serie,Q_MA)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"        ! Discharge series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: Q_MA"        Mean annual runoff

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE Q_MeanAnnual(Q_serie,Q_MA)

  !   USE mo_julian,              only : caldat
  !   USE mo_global_variables,    only : evalPer

  !   IMPLICIT NONE

  !   real(dp), dimension(:), intent(in)     :: Q_serie    ! simulated  daily river flow, 
  !   real(dp),               intent(out)    :: Q_MA       ! Mean annual discharge

  !   ! local variables
  !   integer(i4)                            :: dd         ! day
  !   integer(i4)                            :: mm         ! month
  !   integer(i4)                            :: numYears   ! number of years
  !   real(dp)                               :: Qsum       ! sum of discharge 
  !   real(dp), allocatable, dimension(:)    :: Q_Ann      ! Annual discharge
  !   integer(i4)                            :: ii,jj      ! Counters
  !   integer(i4)                            :: year_new   ! new year
  !   integer(i4)                            :: year_old   ! old year
  !   integer(i4)                            :: yy         ! year

  !   !initialize
  !   year_old    = evalPer%yStart
  !   year_new    = year_old
  !   Qsum        = 0.0_dp
  !   jj          = 0_i4
  !   numYears    = evalPer%yEnd-evalPer%yStart+1_i4 
  !   allocate(Q_Ann(numYears))
  !   ! evalPerDate = evalPer%julStart        

  !   !loop over time series and sum discharge per year
  !   do ii=1, size(Q_serie,1)

  !      if(year_new .eq. year_old) then
  !         Qsum=Qsum+Q_serie(ii)
  !         if(ii .eq. size(Q_serie,1)) then !exception for end of timeseries
  !            jj=jj+1_i4
  !            Q_Ann(jj)=Qsum
  !         end if
  !      else
  !         !store value when next year is reached and initialize for new year
  !         jj=jj+1_i4
  !         Q_Ann(jj)=Qsum
  !         Qsum=Q_serie(ii)
  !      end if

  !      !update counters
  !      year_old=year_new
  !      call caldat(evalPer%julStart+ii, dd, mm, yy)
  !      year_new=yy
  !   end do

  !   !determine the mean of the annual discharge sums
  !   Q_MA=sum(Q_Ann)/real(size(Q_Ann,1),dp)

  !   deallocate(Q_Ann)

  ! END SUBROUTINE Q_MeanAnnual

  !-------------------------------------------------------------------------------
  !     NAME
  !         Autocorrelation

  !     PURPOSE
  !>        \brief Calculates autocorrelation

  !>        \details Calculates autocorrelation serie with lag times to 200 timesteps

  !     CALLING SEQUENCE
  !         call Autocorrelation(Q_serie, maxlag, AC)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"        ! Discharge series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), allocatable, dimension(:) :: AC"        ! Autocorrelation function

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4) ,optional                :: maxlag"    ! maximum lag of series

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  SUBROUTINE Autocorrelation(Q_serie, AC, maxlag_in)

    USE mo_template

    IMPLICIT NONE

    real(dp),              dimension(:),           intent(in)  :: Q_serie    ! Simulated  daily river flow
    real(dp), allocatable, dimension(:),           intent(out) :: AC         ! Autocorrelation serie
    integer(i4),                         optional, intent(in)  :: maxlag_in  ! maximum lag in autocorrelation

    ! local variables
    real(dp)    :: AC_up          ! Nominator 
    real(dp)    :: AC_down        ! Denominator
    integer(i4) :: tlag           ! lag time
    integer(i4) :: ii             ! Counter
    integer(i4) :: maxlag         ! maximum lag in autocorrelation
    real(dp)    :: Qmean          ! Mean discharge

    ! check optional
    if( present(maxlag_in) ) then
       maxlag = maxlag_in
    else
       maxlag = 200_i4
    end if

    !initialize
    AC_up=0.0_dp
    if(.not. allocated(AC)) allocate( AC(  maxlag) )

    !determine mean
    Qmean=mean(Q_serie(:))

    do tlag=1, maxlag
       AC_down=sum((Q_serie-mean(Q_serie(1:(size(Q_serie,1)-tlag))))**2)

       do ii=1, (size(Q_serie,1)-tlag)
          AC_up=AC_up+(Q_serie(ii)-Qmean)*(Q_serie(ii+tlag)-Qmean)
       end do

       AC(tlag)=AC_up/AC_down
       AC_up=0.0_dp
    end do

  END SUBROUTINE Autocorrelation

  !-------------------------------------------------------------------------------
  !     NAME
  !         Autocorrelation_low

  !     PURPOSE
  !>        \brief Calculates autocorrelation for low flows

  !>        \details Calculates autocorrelation for low flows with delay of 1

  !     CALLING SEQUENCE
  !         call Autocorrelation_low(Q_serie, mm_low, mm_high, AC_low)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"       ! Discharge series
  !>        \param[in] "real(dp), dimension(:) :: mm_low"        ! month where low flow starts
  !>        \param[in] "real(dp), dimension(:) :: mm_high"       ! month where high flow starts

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:) :: AC_low"        Autocorrelation function

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         NONE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE Autocorrelation_low(Q_serie, mm_low, mm_high, AC_low)

  !   USE mo_global_variables,            only : evalPer
  !   USE mo_julian,                      only : caldat

  !   IMPLICIT NONE

  !   real(dp), dimension(:), intent(in)  :: Q_serie ! River flow
  !   integer(i4),            intent(in)  :: mm_low  ! month where low flow starts
  !   integer(i4),            intent(in)  :: mm_high ! month where high flow starts
  !   real(dp),               intent(out) :: AC_low  ! Autocorrelation 1 day low flows

  !   ! local variables
  !   real(dp)    :: AC_up   ! Nominator 
  !   real(dp)    :: AC_down ! Denominator
  !   integer(i4) :: date    ! date in julian day
  !   integer(i4) :: dd      ! day
  !   integer(i4) :: dd1     ! another day
  !   integer(i4) :: mm      ! month
  !   integer(i4) :: mm1     ! another month
  !   integer(i4) :: jj      ! Counters
  !   integer(i4) :: n_low   ! Number of low flow days
  !   real(dp)    :: Qmean   ! Mean discharge
  !   real(dp)    :: Qsum    ! Sum of discharge
  !   integer(i4) :: yy      ! year
  !   integer(i4) :: yy1     ! another year

  !   !initialize
  !   AC_down = 0.0_dp
  !   AC_up   = 0.0_dp
  !   Qsum    = 0.0_dp
  !   n_low   = 0_i4

  !   !calculate the total mean of the low flows

  !   if(mm_low .lt. mm_high) then

  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            n_low=n_low+1_i4
  !            Qsum=Qsum+Q_serie(jj)
  !         end if
  !         date=date+1
  !      end do

  !      Qmean=Qsum/real(n_low,dp)

  !      !calculate AClow
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            call caldat(date+1, dd=dd1, mm=mm1, yy=yy1)
  !            if( (mm1 .ge. mm_low) .and. (mm1 .lt. mm_high) ) then
  !               AC_down=AC_down+(Q_serie(jj)-Qmean)**2
  !               AC_up=AC_up+(Q_serie(jj)-Qmean)*(Q_serie(jj+1)-Qmean)
  !            end if
  !         end if
  !         date=date+1
  !      end do

  !   else

  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            n_low=n_low+1_i4
  !            Qsum=Qsum+Q_serie(jj)
  !         end if
  !         date=date+1
  !      end do

  !      Qmean=Qsum/real(n_low,dp)

  !      !calculate AClow
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            call caldat(date+1, dd=dd1, mm=mm1, yy=yy1)
  !            if( (mm1 .ge. mm_low) .or. (mm1 .lt. mm_high) ) then
  !               AC_down=AC_down+(Q_serie(jj)-Qmean)**2
  !               AC_up=AC_up+(Q_serie(jj)-Qmean)*(Q_serie(jj+1)-Qmean)
  !            end if
  !         end if
  !         date=date+1
  !      end do

  !   end if

  !   AC_low=AC_up/AC_down

  ! END SUBROUTINE Autocorrelation_low

  !-------------------------------------------------------------------------------
  !     NAME
  !         Autocorrelation_high

  !     PURPOSE
  !>        \brief Calculates autocorrelation for high flows

  !>        \details Calculates autocorrelation for high flows with delay of 1

  !     CALLING SEQUENCE
  !         call Autocorrelation_high(Q_serie, mm_low, mm_high, AC_high)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"       ! Discharge series
  !>        \param[in] "real(dp), dimension(:) :: mm_low"        ! month where low flow starts
  !>        \param[in] "real(dp), dimension(:) :: mm_high"       ! month where high flow starts

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:) :: AC_high"        Autocorrelation function

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE Autocorrelation_high(Q_serie, mm_low, mm_high, AC_high)

  !   USE mo_julian,                      only : caldat
  !   USE mo_global_variables,            only : evalPer

  !   IMPLICIT NONE

  !   real(dp), dimension(:), intent(in)  :: Q_serie ! River flow
  !   integer(i4),            intent(in)  :: mm_low  ! month where low flow starts
  !   integer(i4),            intent(in)  :: mm_high ! month where high flow starts   
  !   real(dp),               intent(out) :: AC_high ! Autocorrelation 1 day low flows

  !   ! local variables
  !   real(dp)    :: AC_up   ! Nominator 
  !   real(dp)    :: AC_down ! Denominator
  !   integer(i4) :: date    ! date in julian day
  !   integer(i4) :: dd      ! day
  !   integer(i4) :: dd1     ! another day
  !   integer(i4) :: mm      ! month
  !   integer(i4) :: mm1     ! another month

  !   integer(i4) :: jj      ! Counters
  !   integer(i4) :: n_high  ! Number of low flow days
  !   real(dp)    :: Qmean   ! Mean discharge
  !   real(dp)    :: Qsum    ! Sum of discharge
  !   integer(i4) :: yy      ! year
  !   integer(i4) :: yy1     ! another year

  !   !initialize
  !   AC_down = 0.0_dp
  !   AC_up   = 0.0_dp
  !   Qsum    = 0.0_dp
  !   n_high  = 0_i4

  !   !calculate the total mean of the high flows

  !   if(mm_low .lt. mm_high) then  !middle of year is low flow period

  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .or. (mm .ge. mm_high) ) then
  !            n_high=n_high+1_i4
  !            Qsum=Qsum+Q_serie(jj)
  !         end if
  !         date=date+1
  !      end do

  !      Qmean=Qsum/real(n_high,dp)

  !      !calculate AChigh
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)-1_i4
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .or. (mm .ge. mm_high) ) then
  !            call caldat(date+1, dd=dd1, mm=mm1, yy=yy1)
  !            if( (mm .lt. mm_low .and. mm1 .lt. mm_low) .or. (mm .ge. mm_high .and. mm1 .ge. mm_high) ) then
  !               AC_down=AC_down+(Q_serie(jj)-Qmean)**2
  !               AC_up=AC_up+(Q_serie(jj)-Qmean)*(Q_serie(jj+1)-Qmean)
  !            end if
  !         end if
  !         date=date+1
  !      end do

  !   else  !low flows in begin and end of year

  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .and. (mm .ge. mm_high) ) then
  !            n_high=n_high+1_i4
  !            Qsum=Qsum+Q_serie(jj)
  !         end if
  !         date=date+1
  !      end do

  !      Qmean=Qsum/real(n_high,dp)

  !      !calculate AChigh
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)-1_i4
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .and. (mm .ge. mm_high) ) then
  !            call caldat(date+1, dd=dd1, mm=mm1, yy=yy1)
  !            if( (mm .lt. mm_low .and. mm1 .lt. mm_low) .and. (mm .ge. mm_high .and. mm1 .ge. mm_high) ) then
  !               AC_down=AC_down+(Q_serie(jj)-Qmean)**2
  !               AC_up=AC_up+(Q_serie(jj)-Qmean)*(Q_serie(jj+1)-Qmean)
  !            end if
  !         end if
  !         date=date+1
  !      end do

  !   end if

  !   AC_high=AC_up/AC_down

  ! END SUBROUTINE Autocorrelation_high

  !-------------------------------------------------------------------------------
  !     NAME
  !         Limb_densities

  !     PURPOSE
  !>        \brief Calculates limb densities

  !>        \details Calculates rising and declinging limb densities
  !>                     \f[ RLD=t_rise/n_peak \f]
  !>                 Duration the hydrograph increases divided by the number of peaks, 
  !>                 i.e. timestep before and after have lower values
  !>                     \f[ DLD=t_fall/n_peak \f]
  !>                 Duration the hydrograph decreases divided by the number of peaks,
  !>                 i.e. timestep before and after have lower values

  !     CALLING SEQUENCE
  !         call Limb_densities(Q_serie,RLD,DLD)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"        ! Discharge series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: RLD"        Rising    limb density
  !>        \param[out] "real(dp) :: DLD"        Declining limb density

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  SUBROUTINE Limb_densities(Q_serie, RLD, DLD)

    IMPLICIT NONE

    real(dp), dimension(:), intent(in)  :: Q_serie ! River flow,    
    real(dp),               intent(out) :: RLD     ! Rising limb density
    real(dp),               intent(out) :: DLD     ! Declining limb density

                                                   ! local variables
    integer(i4)                         :: jj      ! Counter
    integer(i4)                         :: n_peak  ! Number of peaks
    real(dp)                            :: t_fall  ! Duration of the hydrograph declining
    real(dp)                            :: t_rise  ! Duration of the hydrograph rising 

    !initialize
    t_fall=0.0_dp
    t_rise=0.0_dp
    n_peak=0_i4

    !calculate the total mean of the high flows

    do jj=2, size(Q_serie,1)-1

       !check for peak
       if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) ) then
          n_peak=n_peak+1_i4
       end if

       !check if hydrograph rises
       if(Q_serie(jj-1) .lt. Q_serie(jj)) then
          t_rise=t_rise+1.0_dp
       end if

       !check if hydrograph falls
       if(Q_serie(jj-1) .gt. Q_serie(jj)) then
          t_fall=t_fall+1.0_dp
       end if

    end do

    RLD = t_rise/real(n_peak,dp)
    DLD = t_fall/real(n_peak,dp)

  END SUBROUTINE Limb_densities

  !-------------------------------------------------------------------------------
  !     NAME
  !         FlowDurationCurves

  !     PURPOSE
  !>        \brief Calculates the flow duration curves

  !>        \details Calculates the flow duration curves for discharge time series

  !     CALLING SEQUENCE
  !         call FlowDurationCurves(Q_serie, Quantiles, FDC_serie, Q_quantiles)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"     ! Discharge series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)  
  !>        \param[out] "real(dp), allocatable, dimension(:,:) :: FDC_serie" !serie (ordened discharge)


  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), dimension(:) :: Quantiles"    ! Percentages of occurrence

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), allocatable, dimension(:)   :: Q_quantiles"  !Q with p=10% and 50%   

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  SUBROUTINE FlowDurationCurves(Q_serie, FDC_serie, Quantiles, Q_quantiles)

    use mo_orderpack, only: sort_index

    IMPLICIT NONE

    real(dp), dimension(:),                        intent(in)  :: Q_serie      ! River flow,
    real(dp), allocatable, dimension(:,:),         intent(out) :: FDC_serie    ! ordened Q
    real(dp),              dimension(:), optional, intent(in)  :: Quantiles    ! Percentages of occurrence,
    real(dp), allocatable, dimension(:), optional, intent(out) :: Q_quantiles  ! Q with p=10% and 50%

    !local variables
    integer(i4)                             :: ii         ! Counter
    integer(i4), dimension(size(Q_serie,1)) :: ind_sort   ! indexes
    integer(i4)                             :: n          ! Size of array
    real(dp),    dimension(size(Q_serie,1)) :: per        ! percentages and sorted Q
    real(dp),    dimension(size(Q_serie,1)) :: Q_sort     ! percentages and sorted Q

    ind_sort = sort_index(Q_serie)
    Q_sort   = Q_serie(ind_sort)

    do ii=1, size(Q_serie,1)
       per(ii)=real(ii,dp)/size(Q_serie,1)*100_dp
    end do

    n=size(Q_sort,1)
    Q_sort=Q_sort(n:1:-1)

    if ( present( Quantiles ) .and. present(Q_quantiles) )  call Quantile(Q_sort, Quantiles, Q_quantiles)

    allocate( FDC_serie(n, 2) )

    FDC_serie(:,1)=per(:)
    FDC_serie(:,2)=Q_sort

  END SUBROUTINE FlowDurationCurves

  !-------------------------------------------------------------------------------
  !     NAME
  !         FlowDurationCurves_low

  !     PURPOSE
  !>        \brief Calculates the flow duration curves for low flow periods

  !>        \details Calculates the flow duration curves for discharge time series during low flow periods

  !     CALLING SEQUENCE
  !         call FlowDurationCurves_low(Q_serie, Quantiles, mm_low, mm_high, FDClow_serie, Q_quantiles)

  !     INTENT(IN)
  !>         \param[in]  real(dp), dimension(:)          :: Q_serie         ! River flow,
  !>         \param[in]  integer(i4)                     :: mm_low          ! low flow month
  !>         \param[in]  integer(i4)                     :: mm_high         ! high flow month 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>         \param[out]  real(dp),allocatable, dimension(:,:) :: FDClow_serie    ! ordened Q

  !     INTENT(IN), OPTIONAL
  !>         \param[in]  real(dp), dimension(:)          :: Quantiles       ! Percentages of occurrence,

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>         \param[out]  real(dp),allocatable, dimension(:)   :: Q_quantiles     ! Q with different 

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE FlowDurationCurves_low(Q_serie, mm_low, mm_high, FDClow_serie, Quantiles, Q_quantiles)

  !   USE mo_julian,           only : caldat
  !   USE mo_global_variables, only : evalPer
  !   USE mo_orderpack,        only : sort_index

  !   IMPLICIT NONE

  !   real(dp), dimension(:),                       intent(in)  :: Q_serie         ! River flow,
  !   integer(i4),                                  intent(in)  :: mm_low          ! low flow month
  !   integer(i4),                                  intent(in)  :: mm_high         ! high flow month 
  !   real(dp), allocatable,dimension(:,:),         intent(out) :: FDClow_serie    ! ordened Q
  !   real(dp), dimension(:),             optional, intent(in)  :: Quantiles       ! Percentages of occurrence,
  !   real(dp), allocatable,dimension(:), optional, intent(out) :: Q_quantiles     ! Q with different occurences

  !   !local variables
  !   integer(i4)                            :: date         ! date
  !   integer(i4)                            :: dd           ! day
  !   integer(i4)                            :: ii           ! Counters
  !   integer(i4), allocatable, dimension(:) :: ind_sort     ! indexes
  !   integer(i4)                            :: jj           ! Counters
  !   integer(i4)                            :: kk           ! Counters
  !   integer(i4)                            :: mm           ! month
  !   integer(i4)                            :: m_low        ! Array sizes
  !   integer(i4)                            :: n_low        ! Array sizes
  !   real(dp),    allocatable, dimension(:) :: per          ! percentages
  !   real(dp),    allocatable, dimension(:) :: Qlow         ! discharge in low flow period
  !   real(dp),    allocatable, dimension(:) :: Q_sort       ! percentages and sorted Q
  !   integer(i4)                            :: yy           ! year

  !   n_low=0_i4
  !   m_low=0_i4

  !   if(mm_low .lt. mm_high) then

  !      !determine number of low flow days
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            n_low=n_low+1_i4
  !         end if
  !         date=date+1_i4
  !      end do

  !      !Create low flow array
  !      allocate(Qlow(n_low))

  !      date=evalPer%julStart
  !      do ii=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            m_low=m_low+1
  !            Qlow(m_low)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   else

  !      !determine number of low flow days
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            n_low=n_low+1_i4
  !         end if
  !         date=date+1_i4
  !      end do

  !      !Create low flow array
  !      allocate(Qlow(n_low))

  !      date=evalPer%julStart
  !      do ii=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            m_low=m_low+1
  !            Qlow(m_low)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   end if

  !   allocate( ind_sort(n_low) ) 
  !   allocate( per(n_low) ) 
  !   allocate( Q_sort(n_low) ) 
  !   allocate( FDClow_serie(n_low,2)) 

  !   !Create low flow duration curve
  !   ind_sort=sort_index(Qlow)
  !   Q_sort=Qlow(ind_sort)

  !   do kk=1, n_low
  !      per(kk)=real(kk,dp)/real(n_low,dp)*100_dp
  !   end do

  !   Q_sort=Q_sort(n_low:1:-1)

  !   FDClow_serie(:,1)=per
  !   FDClow_serie(:,2)=Q_sort

  !   if(present( Quantiles )) call Quantile(Q_sort, Quantiles, Q_quantiles)

  !   deallocate( ind_sort ) 
  !   deallocate( per ) 
  !   deallocate( Q_sort ) 

  ! END SUBROUTINE FlowDurationCurves_low

  !-------------------------------------------------------------------------------
  !     NAME
  !         FlowDurationCurves_high

  !     PURPOSE
  !>        \brief Calculates the flow duration curves for high flow periods

  !>        \details Calculates the flow duration curves for discharge time series during high flow periods

  !     CALLING SEQUENCE
  !         call FlowDurationCurves_high(Q_serie, Quantiles, mm_low, mm_high, FDClow_serie, Q_quantiles)

  !     INTENT(IN)
  !>         \param[in]  real(dp), dimension(:)          :: Q_serie         ! River flow,
  !>         \param[in]  real(dp), dimension(:)          :: Quantiles       ! Percentages of occurrence,
  !>         \param[in]  integer(i4)                     :: mm_low          ! low flow month
  !>         \param[in]  integer(i4)                     :: mm_high         ! high flow month 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>         \param[out]  real(dp), allocatable, dimension(:,:) :: FDChigh_serie ! ordened Q
  !>         \param[out]  real(dp), allocatable, dimension(:)   :: Q_quantiles   ! Q with different probabilities   

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE FlowDurationCurves_high(Q_serie, mm_low, mm_high, FDChigh_serie, Quantiles, Q_quantiles)

  !   USE mo_julian,                      only : caldat
  !   USE mo_global_variables,            only : evalPer
  !   USE mo_orderpack,                   only: sort_index

  !   IMPLICIT NONE

  !   real(dp),              dimension(:),           intent(in)  :: Q_serie         ! River flow,
  !   integer(i4),                                   intent(in)  :: mm_low          ! low flow month
  !   integer(i4),                                   intent(in)  :: mm_high         ! high flow month 
  !   real(dp), allocatable, dimension(:,:),         intent(out) :: FDChigh_serie   ! ordened Q
  !   real(dp),              dimension(:), optional, intent(in)  :: Quantiles       ! Percentages of occurrence,
  !   real(dp), allocatable, dimension(:), optional, intent(out) :: Q_quantiles     ! Q with different occurences

  !   !local variables
  !   integer(i4)                            :: date         ! date
  !   integer(i4)                            :: dd           ! day
  !   integer(i4)                            :: ii           ! Counters
  !   integer(i4), allocatable, dimension(:) :: ind_sort     ! indexes
  !   integer(i4)                            :: jj           ! Counters
  !   integer(i4)                            :: kk           ! Counters
  !   integer(i4)                            :: mm           ! month
  !   integer(i4)                            :: m_high       ! Array sizes
  !   integer(i4)                            :: n_high       ! Array sizes
  !   real(dp),   allocatable, dimension(:)  :: per          ! percentages
  !   real(dp),   allocatable, dimension(:)  :: Qhigh        ! discharge in high flow period
  !   real(dp),   allocatable, dimension(:)  :: Q_sort       ! percentages and sorted Q
  !   integer(i4)                            :: yy           ! year

  !   n_high=0_i4
  !   m_high=0_i4

  !   if(mm_low .lt. mm_high) then

  !      !determine number of high flow days
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .or. (mm .ge. mm_high) ) then
  !            n_high=n_high+1_i4
  !         end if
  !         date=date+1_4
  !      end do

  !      !Create high flow array
  !      allocate(Qhigh(n_high))

  !      date=evalPer%julStart
  !      do ii=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .or. (mm .ge. mm_high) ) then
  !            m_high=m_high+1_i4
  !            Qhigh(m_high)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   else

  !      !determine number of high flow days
  !      date=evalPer%julStart
  !      do jj=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .and. (mm .ge. mm_high) ) then
  !            n_high=n_high+1_i4
  !         end if
  !         date=date+1_4
  !      end do

  !      !Create high flow array
  !      allocate(Qhigh(n_high))

  !      date=evalPer%julStart
  !      do ii=1, size(Q_serie,1)
  !         call caldat(date, dd, mm, yy)
  !         if( (mm .lt. mm_low) .and. (mm .ge. mm_high) ) then
  !            m_high=m_high+1_i4
  !            Qhigh(m_high)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   end if

  !   allocate( ind_sort(n_high) ) 
  !   allocate( per(n_high) ) 
  !   allocate( Q_sort(n_high) ) 
  !   allocate( FDChigh_serie(n_high,2)) 

  !   !Create high flow duration curve
  !   ind_sort=sort_index(Qhigh)
  !   Q_sort=Qhigh(ind_sort)

  !   do kk=1, n_high
  !      per(kk)=real(kk,dp)/real(n_high,dp)*100_dp
  !   end do

  !   Q_sort=Q_sort(n_high:1:-1)

  !   if(present( Quantiles ))     call Quantile(Q_sort, Quantiles, Q_quantiles)

  !   FDChigh_serie(:,1) = per
  !   FDChigh_serie(:,2) = Q_sort

  !   deallocate( ind_sort ) 
  !   deallocate( per ) 
  !   deallocate( Q_sort ) 

  ! END SUBROUTINE FlowDurationCurves_high

  !-------------------------------------------------------------------------------
  !     NAME
  !         PeakDistribution

  !     PURPOSE
  !>        \brief Calculates the peak distributions

  !>        \details Calculates the peak distributions according to 
  !>                    peaks= (Qpeaks_10%-Qpeaks_50%)/(0.9-0.5)
  !>                 It therefore represents the slope between the peak flow 
  !>                 duration curve between the 10th and 50th percentile.

  !     CALLING SEQUENCE
  !         call PeakDistribution(Q_serie, Quantiles, peaks, Q_quantiles)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"        ! Discharge series
  !>        \param[in] "real(dp), dimension(:) :: Quantiles"      ! Probabilities

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         \param[out] real(dp)                          :: peaks       ! slope between 50th and 10t peal flow percentile
  !>        \param[out] real(dp),allocatable,dimension(:) :: Q_quantiles ! 50% peak river flow 

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  SUBROUTINE PeakDistribution(Q_serie, Quantiles, peaks, Q_quantiles)

    USE mo_orderpack, only: sort_index

    IMPLICIT NONE

    real(dp),              dimension(:), intent(in)  :: Q_serie     ! River flow, 
    real(dp),              dimension(:), intent(in)  :: Quantiles   ! Percentages of occurence,       
    real(dp),                            intent(out) :: peaks       ! slope of peak flow duration curve
    real(dp), allocatable, dimension(:), intent(out) :: Q_quantiles ! Q with p=10% and 50%

    ! local variables
    integer(i4)                              :: ii          ! counters
    integer(i4), allocatable, dimension(:)   :: ind_sort    ! indexes
    integer(i4)                              :: jj          ! counters
    integer(i4)                              :: kk          ! counters
    integer(i4)                              :: nn          ! counters
    integer(i4)                              :: n_peak      ! Number of peaks
    real(dp),                 dimension(2)   :: p           ! prob. of p=10% and 50%
    real(dp),    allocatable, dimension(:)   :: per         ! percentages
    real(dp),    allocatable, dimension(:)   :: Q_p         ! Q with p=10% and 50%
    real(dp) ,   allocatable, dimension(:)   :: Qpeak       ! peak distribution discharge
    real(dp),    allocatable, dimension(:)   :: Q_sort      ! sorted Q

    n_peak=0_i4

    !calculate the total mean of the high flows and count peaks
    do jj=2, size(Q_serie,1)-1
       if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) ) then
          n_peak=n_peak+1_i4
       end if
    end do

    allocate(Qpeak(n_peak))

    ! find peaks
    kk=0
    do ii=2, size(Q_serie,1)-1
       if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) ) then
          kk=kk+1_i4       
          Qpeak(kk)=Q_serie(ii)
       end if
    end do

    allocate( ind_sort(n_peak) ) 
    allocate( per(n_peak) ) 
    allocate( Q_sort(n_peak) ) 

    !Create peak flow duration curve
    ind_sort=sort_index(Qpeak)
    Q_sort=Qpeak(ind_sort)

    do nn=1, size(Q_sort,1)
       per(nn)=real(nn,dp)/real(n_peak,dp)*100_dp
    end do

    Q_sort=Q_sort(n_peak:1:-1)

    !calculate slope between 10 and 50 percent quantiles, per definition
    p=(/0.1_dp,0.5_dp/)
    call Quantile(Q_sort,p,Q_p)

    peaks=(Q_p(1)-Q_p(2))/(0.9_dp-0.5_dp)

    call Quantile(Q_sort,Quantiles,Q_quantiles)

    deallocate( ind_sort ) 
    deallocate( per ) 
    deallocate( Q_sort ) 
    deallocate( Qpeak)

  END SUBROUTINE PeakDistribution

  !-------------------------------------------------------------------------------
  !     NAME
  !         PeakDistribution_low

  !     PURPOSE
  !>        \brief Calculates the peak distributions

  !>        \details Calculates the peak distributions according to 
  !>                    peaks= (Qpeaks_10%-Qpeaks_50%)/(0.9-0.5)
  !>                 It therefore represents the slope between the peak flow
  !>                 duration curve between the 10th and 50th percentile.

  !     CALLING SEQUENCE
  !         call PeakDistribution_low(Q_serie, Quantiles, mm_low, mm_high, peaks_low, Q_quantiles)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"    !River flow, 
  !>        \param[in] "real(dp),allocatable, dimension(:) :: Quantiles"  !Percentages of occurence, 
  !>        \param[in] "integer(i4)            :: mm_low      !start month of low flow period  
  !>        \param[in] "integer(i4)            :: mm_high     !start month of high flow

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] real(dp) :: peaks_low        ! slope between 50th and 10th low flow percentile, 
  !>        \param[out] real(dp) :: Q_quantiles      ! Q for prob. p 

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE PeakDistribution_low(Q_serie, Quantiles, mm_low, mm_high, peaks_low, Q_quantiles)

  !   USE mo_julian,                      only: caldat
  !   USE mo_global_variables,            only: evalPer
  !   USE mo_orderpack,                   only: sort_index

  !   IMPLICIT NONE

  !   real(dp),              dimension(:), intent(in)  :: Q_serie     ! River flow, 
  !   real(dp),              dimension(:), intent(in)  :: Quantiles   ! Percentages of occurence, 
  !   integer(i4),                         intent(in)  :: mm_low      ! start month of low flow period  
  !   integer(i4),                         intent(in)  :: mm_high     ! start month of high flow period     
  !   real(dp),                            intent(out) :: peaks_low   ! slope of peak flow duration curve
  !   real(dp), allocatable, dimension(:), intent(out) :: Q_quantiles ! Q with p=10% and 50%

  !   !local variables
  !   integer(i4)                              :: date        ! date
  !   integer(i4)                              :: dd          ! day
  !   integer(i4)                              :: ii          ! counters
  !   integer(i4), allocatable, dimension(:)   :: ind_sort    ! indexes
  !   integer(i4)                              :: jj          ! counters
  !   integer(i4)                              :: kk          ! counters
  !   integer(i4)                              :: nn          ! counters
  !   integer(i4)                              :: n_peak      ! Number of peaks
  !   integer(i4)                              :: mm          ! month
  !   real(dp),                 dimension(2)   :: p           ! prob. of p=10% and 50%
  !   real(dp),    allocatable, dimension(:)   :: per         ! percentages
  !   real(dp),    allocatable, dimension(:)   :: Q_p         ! Q with p=10% and 50%
  !   real(dp),    allocatable, dimension(:)   :: Qpeak       ! peak distribution discharge
  !   real(dp),    allocatable, dimension(:)   :: Q_sort      ! sorted Q
  !   integer(i4)                              :: yy          ! day, month, year

  !   n_peak=0_i4
  !   date=evalPer%julStart+1_i4

  !   if(mm_low .lt. mm_high) then

  !      do jj=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
  !              .and. (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            n_peak=n_peak+1_i4
  !         end if
  !         date=date+1_i4
  !      end do

  !      allocate(Qpeak(n_peak))

  !      kk=0_i4
  !      date=evalPer%julStart+1_i4

  !      ! now fill Qpeak
  !      do ii=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
  !              .and. (mm .ge. mm_low) .and. (mm .lt. mm_high) ) then
  !            kk=kk+1_i4       
  !            Qpeak(kk)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   else

  !      do jj=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
  !              .and. (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            n_peak=n_peak+1_i4
  !         end if
  !         date=date+1_i4
  !      end do

  !      allocate(Qpeak(n_peak))

  !      kk=0_i4
  !      date=evalPer%julStart+1_i4

  !      ! now fill Qpeak
  !      do ii=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
  !              .and. (mm .ge. mm_low) .or. (mm .lt. mm_high) ) then
  !            kk=kk+1_i4       
  !            Qpeak(kk)=Q_serie(ii)
  !         end if
  !         date=date+1_i4
  !      end do

  !   end if

  !   allocate( ind_sort(n_peak) ) 
  !   allocate( per(n_peak) ) 
  !   allocate( Q_sort(n_peak) ) 

  !   !Create peak flow duration curve

  !   ind_sort=sort_index(Qpeak)
  !   Q_sort=Qpeak(ind_sort)

  !   do nn=1, n_peak
  !      per(nn)=real(nn,dp)/real(n_peak,dp)*100_dp
  !   end do

  !   Q_sort=Q_sort(n_peak:1:-1)

  !   !calculate slope between 10 and 50 percent quantiles, per definition
  !   p=(/0.1_dp,0.5_dp/)
  !   call Quantile(Q_sort, p, Q_p)

  !   peaks_low=(Q_p(1)-Q_p(2))/(0.9_dp-0.5_dp)

  !   call Quantile(Q_sort, Quantiles, Q_quantiles)

  !   deallocate( ind_sort ) 
  !   deallocate( per ) 
  !   deallocate( Q_sort ) 
  !   deallocate( Qpeak )

  ! END SUBROUTINE PeakDistribution_low

  !-------------------------------------------------------------------------------
  !     NAME
  !         PeakDistribution_high

  !     PURPOSE
  !>        \brief   Calculates the peak distributions

  !>        \details Calculates the peak distributions according to 
  !>                    peaks= (Qpeaks_10%-Qpeaks_50%)/(0.9-0.5)
  !>                 It therefore represents the slope between the peak flow 
  !>                 duration curve between the 10th and 50th percentile.

  !     CALLING SEQUENCE
  !         call PeakDistribution_high(Q_serie, Quantiles, mm_low, mm_high, peaks_high, Q_quantiles)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: Q_serie"    !River flow, 
  !>        \param[in] "real(dp), allocatable, dimension(:) :: Quantiles"  !Percentages of occurence, 
  !>        \param[in] "integer(i4)            :: mm_low      !start month of low flow period  
  !>        \param[in] "integer(i4)            :: mm_high     !start month of high flow

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] real(dp) :: peaks_high       ! slope between 50th and 10th low flow percentile, 
  !>        \param[out] real(dp) :: Q_quantiles      ! Q for prob. p 

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  ! SUBROUTINE PeakDistribution_high(Q_serie, Quantiles, mm_low, mm_high, peaks_high, Q_quantiles)

  !   USE mo_julian,                      only : caldat
  !   USE mo_global_variables,            only : evalPer
  !   USE mo_orderpack,                   only: sort_index

  !   IMPLICIT NONE       

  !   real(dp),              dimension(:), intent(in)  :: Q_serie     ! River flow, 
  !   real(dp),              dimension(:), intent(in)  :: Quantiles   ! Percentages of occurence, 
  !   integer(i4),                         intent(in)  :: mm_low      ! start month of low flow period  
  !   integer(i4),                         intent(in)  :: mm_high     ! start month of high flow period     
  !   real(dp),                            intent(out) :: peaks_high  ! slope of peak flow duration curve
  !   real(dp), allocatable, dimension(:), intent(out) :: Q_quantiles ! Q with p=10% and 50%

  !   !local variables
  !   integer(i4)                              :: date        ! date
  !   integer(i4)                              :: dd          ! day
  !   integer(i4)                              :: ii          ! counters
  !   integer(i4), allocatable, dimension(:)   :: ind_sort    ! indexes
  !   integer(i4)                              :: jj          ! counters
  !   integer(i4)                              :: kk          ! counters
  !   integer(i4)                              :: nn          ! counters
  !   integer(i4)                              :: n_peak      ! Number of peaks
  !   integer(i4)                              :: mm          ! month
  !   real(dp),                 dimension(2)   :: p           ! prob. of p=10% and 50%
  !   real(dp),    allocatable, dimension(:)   :: per         ! percentages
  !   real(dp),    allocatable, dimension(:)   :: Q_p         ! Q with p=10% and 50%
  !   real(dp),    allocatable, dimension(:)   :: Qpeak       ! peak distribution discharge
  !   real(dp),    allocatable, dimension(:)   :: Q_sort      ! sorted Q
  !   integer(i4)                              :: yy          ! day, month, year

  !   n_peak=0_i4
  !   date=evalPer%julStart+1_i4

  !   if(mm_low .lt. mm_high) then
  !      ! count peaks
  !      do jj=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
  !              .and. ((mm .lt. mm_low) .or. (mm .ge. mm_high)) ) then
  !            n_peak=n_peak+1_i4
  !         end if
  !         date=date+1_4
  !      end do

  !      allocate(Qpeak(n_peak))

  !      kk=0_i4
  !      date=evalPer%julStart+1

  !      ! now fill Qpeak
  !      do ii=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
  !              .and. ((mm .lt. mm_low) .or. (mm .ge. mm_high)) ) then
  !            kk=kk+1_i4       
  !            Qpeak(kk)=Q_serie(ii)
  !         end if
  !         date=date+1
  !      end do

  !   else
  !      ! count peaks
  !      do jj=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
  !              .and. ((mm .lt. mm_low) .and. (mm .ge. mm_high)) ) then
  !            n_peak=n_peak+1_i4
  !         end if
  !         date=date+1_4
  !      end do

  !      allocate(Qpeak(n_peak))

  !      kk=0_i4
  !      date=evalPer%julStart+1

  !      ! now fill Qpeak
  !      do ii=2, size(Q_serie,1)-1
  !         call caldat(date, dd, mm, yy)
  !         if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
  !              .and. ((mm .lt. mm_low) .and. (mm .ge. mm_high)) ) then
  !            kk=kk+1_i4       
  !            Qpeak(kk)=Q_serie(ii)
  !         end if
  !         date=date+1
  !      end do

  !   end if

  !   allocate( ind_sort(n_peak) ) 
  !   allocate( per(n_peak) ) 
  !   allocate( Q_sort(n_peak) ) 

  !   !Create peak flow duration curve

  !   ind_sort=sort_index(Qpeak)
  !   Q_sort=Qpeak(ind_sort)

  !   do nn=1, n_peak
  !      per(nn)=real(nn,dp)/real(n_peak,dp)*100_dp
  !   end do

  !   Q_sort=Q_sort(n_peak:1:-1)

  !   !calculate slope between 10 and 50 percent quantiles
  !   p=(/0.1_dp,0.5_dp/)
  !   call Quantile(Q_sort,p,Q_p)

  !   peaks_high=(Q_p(1)-Q_p(2))/(0.9_dp-0.5_dp)

  !   call Quantile(Q_sort,Quantiles,Q_quantiles)

  !   deallocate( ind_sort ) 
  !   deallocate( per ) 
  !   deallocate( Q_sort ) 
  !   deallocate( Qpeak )

  ! END SUBROUTINE PeakDistribution_high

  !-------------------------------------------------------------------------------
  !     NAME
  !         Quantile

  !     PURPOSE
  !>        \brief   Calculates quantiles of a series

  !>        \details Calculates the quantiles of a series
  !>         

  !     CALLING SEQUENCE
  !         call PeakDistribution_high(Q_serie, peaks_high, Qhighpeak10, Qhighpeak50)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: x        ! Input series
  !>        \param[in] "real(dp), dimension(:) :: p        ! Array of probabilities

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] real(dp),allocatable, dimension(:)  :: Q     ! Quantiles 


  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014

  SUBROUTINE Quantile(x,p,Q)

    IMPLICIT NONE

    real(dp),    dimension(:),             intent(in)  ::x   ! input serie
    real(dp),    dimension(:),             intent(in)  ::p   ! probability array
    real(dp),    dimension(:),allocatable, intent(out) ::Q   ! quantiles

    ! local variables
    real(dp),    dimension(size(p,1)) :: g   ! temp variables
    integer(i4), dimension(size(p,1)) :: j   ! indexes
    real(dp),    dimension(size(p,1)) :: m   ! temp variables
    integer(i4)                       :: n   ! length of x

    allocate( Q(size(p,1)))

    n=size(x,1)

    m=1-p
    j=floor(n*p+m)
    g=n*p+m-real(j,dp)
    Q=(1-g)*x(j)+g*x(j+1) 

  END SUBROUTINE Quantile

END MODULE mo_signatures
