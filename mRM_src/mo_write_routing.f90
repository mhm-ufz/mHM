module mo_write_routing

  use mo_kind, only: i4, dp
  
  implicit none

  public :: write_routing

contains

  subroutine write_routing(runoff)
    
    use mo_global_variables_routing, only: gauge, nGaugesTotal, basin_mrm, nBasins, evalPer, warmingDays, simPer, &
         ntstepday, write_restart, dirRestartOut, &
         nMeasPerDay, optimize
    use mo_restart_routing,  only: write_restart_routing

    implicit none

    ! input variables
    real(dp), dimension(:,:), allocatable, optional, intent(in) :: runoff       ! dim1=time dim2=gauge
    
    ! local variables
    integer(i4) :: iBasin
    integer(i4) :: iDay, iS, iE
    integer(i4) :: ii
    integer(i4) :: tt
    integer(i4) :: gg
    integer(i4) :: nTimeSteps
    real(dp), dimension(:,:), allocatable :: d_Qmod

    ! --------------------------------------------------------------------------
    ! CHECK CONDITIONS FOR WRITING ROUTING
    ! --------------------------------------------------------------------------
    if ((optimize) .or. (nMeasPerDay .ne. 1)) return
    
    ! --------------------------------------------------------------------------
    ! WRITE RESTART
    ! --------------------------------------------------------------------------
    if (write_restart) then
       do iBasin = 1, nBasins
          call write_restart_routing(iBasin, dirRestartOut)
       end do
    end if
    
    ! --------------------------------------------------------------------------
    ! STORE DAILY DISCHARGE TIMESERIES OF EACH GAUGING STATION 
    ! FOR SIMULATIONS DURING THE EVALUATION PERIOD
    !
    !  **** AT DAILY TIME STEPS ****
    ! Note:: Observed Q are stored only for the evaluation period and not for
    !        the warming days
    ! --------------------------------------------------------------------------
    ii = maxval( evalPer(1:nBasins)%julEnd - evalPer(1:nBasins)%julStart + 1 )
    allocate( d_Qmod(ii, nGaugesTotal) ) 
    d_Qmod = 0.0_dp

    ! loop over basins
    do ii = 1, nBasins
       nTimeSteps = ( simPer(ii)%julEnd - simPer(ii)%julStart + 1 ) * NTSTEPDAY
       iDay = 0
       ! loop over timesteps
       do tt = warmingDays(ii)*NTSTEPDAY+1, nTimeSteps, NTSTEPDAY
          iS = tt
          iE = tt + NTSTEPDAY - 1
          iDay = iDay + 1
          ! over gauges
          do gg = 1, basin_mrm%nGauges(ii)
             d_Qmod(iDay, basin_mrm%gaugeIndexList(ii,gg) ) = &
                  sum( runoff(iS:iE, basin_mrm%gaugeIndexList(ii,gg)) )/ real(NTSTEPDAY,dp)
          end do
          !
       end do
    end do
    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_discharge( gauge%Q(:,:), d_Qmod(:,:) )
    ! free space
    deallocate(d_Qmod)        
  end subroutine write_routing
  ! ------------------------------------------------------------------

  !     NAME
  !         write_daily_obs_sim_discharge

  !     PURPOSE
  !>        \brief Write a file for the daily observed and simulated discharge timeseries 
  !>                during the evaluation period for each gauging station

  !>        \details Write a file for the daily observed and simulated discharge timeseries 
  !>                during the evaluation period for each gauging station

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in]  "real(dp), dimension(:,:)    :: Qobs"    daily time series of observed discharge 
  !>                                                             dims = (nModeling_days , nGauges_total)
  !>        \param[in]  "real(dp), dimension(:,:)    :: Qsim"    daily time series of modeled discharge
  !>                                                             dims = (nModeling_days , nGauges_total)

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
  !         None 

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date August 2013

  subroutine write_daily_obs_sim_discharge(Qobs, Qsim)


    use mo_errormeasures,       only: kge, nse
    use mo_julian,              only: dec2date
    use mo_message,             only: message
    use mo_string_utils,        only: num2str
    use mo_utils,               only: ge
    use mo_mrm_file,            only: file_daily_discharge, udaily_discharge
    use mo_global_variables_routing, only: &
         nBasins, &
         basin_mrm, &
          dirOut, evalPer, &
         gauge

    implicit none

    ! input arguments
    real(dp), dimension(:,:), intent(in) :: Qobs      ! observed time series  [nModeling_days X nGauges_total]
    real(dp), dimension(:,:), intent(in) :: Qsim      ! simulated time series [nModeling_days X nGauges_total]

    ! local vars
    character(256) :: fName, formHeader, formData, dummy
    integer(i4) :: bb, gg, tt, err
    integer(i4) :: igauge_start, igauge_end
    integer(i4) :: day, month, year
    real(dp) :: newTime


    ! initalize igauge_start
    igauge_start = 1

    ! basin loop
    do bb = 1, nBasins
       if( basin_mrm%nGauges(bb) .lt. 1 ) cycle

       ! estimate igauge_end
       igauge_end = igauge_start + basin_mrm%nGauges(bb) - 1

       ! check the existance of file
       fName = trim(adjustl(dirOut(bb))) // trim(adjustl(file_daily_discharge))
       open(udaily_discharge, file=trim(fName), status='unknown', action='write', iostat=err)
       if( err .ne. 0 ) then
          call message ('  IOError while openening ',trim(fName))
          call message ('  Error-Code ', num2str(err))
          stop
       end if

       ! header
       write(formHeader, *) '( 4a8, ' , basin_mrm%nGauges(bb),'(2X, a5, i10.10, 2X, a5, i10.10) )' 
       write(udaily_discharge, formHeader) 'No', 'Day', 'Mon', 'Year', &
            ( 'Qobs_', gauge%gaugeId(gg), &
            'Qsim_', gauge%gaugeId(gg), gg=igauge_start, igauge_end )

       ! form data
       write(formData, *) '( 4I8, ' , basin_mrm%nGauges(bb),'(2X,   f15.7 , 2X,  f15.7  ) )' 

       ! write data
       newTime  = real(evalPer(bb)%julStart,dp) - 0.5_dp

       do tt = 1, (evalPer(bb)%julEnd - evalPer(bb)%julStart + 1)          
          call dec2date(newTime, yy=year, mm=month, dd=day)
          write(udaily_discharge, formData) tt, day, month, year, ( Qobs(tt,gg), Qsim(tt,gg) , gg=igauge_start, igauge_end )
          newTime = newTime + 1.0_dp
       end do

       ! close file
       close(udaily_discharge)

       ! screen output
       call message()
       write(dummy,'(I3)') bb
       call message('  OUTPUT: saved daily discharge file for basin ', trim(adjustl(dummy)))
       call message('    to ',trim(fname))
       do gg=igauge_start, igauge_end
          call message('    KGE of daily discharge (gauge #',trim(adjustl(num2str(gg))),'): ', &
               trim(adjustl(num2str(kge(Qobs(:,gg), Qsim(:,gg), mask=(ge(Qobs(:,gg), 0.0_dp)))))) )
          call message('    NSE of daily discharge (gauge #',trim(adjustl(num2str(gg))),'): ', &
               trim(adjustl(num2str(nse(Qobs(:,gg), Qsim(:,gg), mask=(ge(Qobs(:,gg), 0.0_dp)))))) )
       end do

       ! update igauge_start
       igauge_start = igauge_end + 1
       !
    end do
    !
  end subroutine write_daily_obs_sim_discharge

end module mo_write_routing
