module mo_write_routing

  use mo_kind, only: i4, dp
  
  implicit none

  public write_routing

contains

  subroutine write_routing(runoff)

    use mo_write_ascii,         only : write_daily_obs_sim_discharge
    !ST: the following dependency has to be removed
    use mo_global_variables, only: basin, nBasins, gauge, nGaugesTotal, &
         evalPer, warmingDays, ntstepday, simPer
    
    implicit none

    ! input variables
    real(dp), dimension(:,:), allocatable, optional, intent(in) :: runoff       ! dim1=time dim2=gauge
    
    ! local variables
    integer(i4) :: iDay, iS, iE
    integer(i4) :: ii
    integer(i4) :: tt
    integer(i4) :: gg
    integer(i4) :: nTimeSteps
    real(dp), dimension(:,:), allocatable :: d_Qmod
    
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
          do gg = 1, basin%nGauges(ii)
             d_Qmod(iDay, basin%gaugeIndexList(ii,gg) ) = &
                  sum( runoff(iS:iE, basin%gaugeIndexList(ii,gg)) )/ real(NTSTEPDAY,dp)
          end do
          !
       end do
    end do
    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_discharge( gauge%Q(:,:), d_Qmod(:,:) )
    ! free space
    deallocate(d_Qmod)        
  end subroutine write_routing

end module mo_write_routing
