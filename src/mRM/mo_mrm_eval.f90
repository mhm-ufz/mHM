!>       \file mo_mrm_eval.f90

!>       \brief Runs mrm with a specific parameter set and returns required variables, e.g. runoff.

!>       \details Runs mrm with a specific parameter set and returns required variables, e.g. runoff.

!>       \authors Stephan Thober

!>       \date Sep 2015

! Modifications:

module mo_mrm_eval

  implicit none

  public :: mrm_eval

contains

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_eval

  !    PURPOSE
  !>       \brief Runs mrm with a specific parameter set and returns required variables, e.g. runoff.

  !>       \details Runs mrm with a specific parameter set and returns required variables, e.g. runoff.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset" a set of global parameter (gamma) to run mHM, DIMENSION
  !>       [no. of global_Parameters]

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), dimension(:, :), optional :: runoff"        returns runoff time series, DIMENSION
  !>       [nTimeSteps, nGaugesTotal]
  !>       \param[out] "real(dp), dimension(:, :), optional :: sm_opti"       dim1=ncells, dim2=time
  !>       \param[out] "real(dp), dimension(:, :), optional :: basin_avg_tws" dim1=time dim2=nDomains
  !>       \param[out] "real(dp), dimension(:, :), optional :: neutrons_opti" dim1=ncells, dim2=time
  !>       \param[out] "real(dp), dimension(:, :), optional :: et_opti"       dim1=ncells, dim2=time

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Sep 2015

  ! Modifications:
  ! Stephan Thober Nov 2016 - implemented second routing process i.e. adaptive timestep
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_eval(parameterset, runoff, sm_opti, basin_avg_tws, neutrons_opti, et_opti)

    use mo_common_constants, only : HourSecs
    use mo_common_mHM_mRM_variables, only : LCYearId, dirRestartIn, nTStepDay, optimize, read_restart, resolutionRouting, simPer, &
                                            timestep, warmingDays
    use mo_common_variables, only : level1, domainMeta, processMatrix, resolutionHydrology
    use mo_julian, only : caldat, julday
    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_mrm_global_variables, only : InflowGauge, &
                                        L11_C1, L11_C2, L11_L1_ID, L11_TSrout, L11_fromN, L11_length, L11_nLinkFracFPimp, &
                                        L11_nOutlets, L11_netPerm, L11_qMod, L11_qOUT, L11_qTIN, L11_qTR, L11_slope, &
                                        L11_toN, L1_L11_ID, L1_total_runoff_in, basin_mrm, level11, mRM_runoff, &
                                        outputFlxState_mrm, timeStep_model_outputs_mrm, gw_coupling, L0_river_head_mon_sum
    use mo_mrm_init, only : variables_default_init_routing
    use mo_mrm_mpr, only : mrm_update_param
    use mo_mrm_restart, only : mrm_read_restart_states
    use mo_mrm_routing, only : mrm_routing
    use mo_mrm_river_head, only : calc_river_head, avg_and_write_timestep
    use mo_mrm_write, only : mrm_write_output_fluxes
    use mo_utils, only : ge

    implicit none

    ! a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), intent(in) :: parameterset

    ! returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: runoff

    ! dim1=ncells, dim2=time
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: sm_opti

    ! dim1=time dim2=nDomains
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: basin_avg_tws

    ! dim1=ncells, dim2=time
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: neutrons_opti

    ! dim1=ncells, dim2=time
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: et_opti

    integer(i4) :: domainID, iDomain

    integer(i4) :: jj

    integer(i4) :: tt

    integer(i4) :: day

    integer(i4) :: month, prev_month

    integer(i4) :: year

    integer(i4) :: hour

    integer(i4) :: nTimeSteps

    ! Land cover year ID
    integer(i4) :: Lcover_yID

    ! start and end index at level 1 for current domain
    integer(i4) :: s1, e1

    ! start and end index at L11
    integer(i4) :: s11, e11

    ! discharge timestep
    integer(i4) :: iDischargeTS

    ! factor between routing and hydrological modelling resolution
    real(dp) :: tsRoutFactor

    ! factor between routing and hydrological modelling resolution (temporary variable)
    real(dp) :: tsRoutFactorIn

    ! timestep of runoff to rout [h]
    ! - identical to timestep of input if tsRoutFactor is less than 1
    ! - tsRoutFactor * timestep if tsRoutFactor is greater than 1
    integer(i4) :: timestep_rout

    real(dp) :: newTime

    ! Runoff that is input for routing
    real(dp), allocatable, dimension(:) :: RunToRout

    ! inflowing discharge
    real(dp), allocatable, dimension(:) :: InflowDischarge

    logical, pointer, dimension(:, :) :: mask11

    ! flag for performing routing
    logical :: do_rout

    ! flag for monthly mean of river head
    logical :: is_new_month = .false.

    ! initialize variables
    month = 0_i4
    
    if (present(sm_opti) .or. present(basin_avg_tws) .or. present(neutrons_opti) .or. present(et_opti)) then
      call message("Error during initialization of mrm_eval, incorrect call from optimization routine.")
      stop 1
    end if
    if (.not. read_restart) then
      !-------------------------------------------
      ! L11 ROUTING STATE VARIABLES, FLUXES AND
      !             PARAMETERS
      !-------------------------------------------
      call variables_default_init_routing()
    end if

    ! ----------------------------------------
    ! initialize factor between routing resolution and hydrologic model resolution
    ! ----------------------------------------
    tsRoutFactor = 1_i4
    allocate(InflowDischarge(size(InflowGauge%Q, dim = 2)))
    InflowDischarge = 0._dp
    ! ----------------------------------------
    ! loop over domains
    ! ----------------------------------------
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      ! read states from restart
      if (read_restart) call mrm_read_restart_states(iDomain, domainID, dirRestartIn(iDomain))
      !
      ! get domain information at L11 and L1 if routing is activated
      s1 = level1(iDomain)%iStart
      e1 = level1(iDomain)%iEnd
      s11 = level11(iDomain)%iStart
      e11 = level11(iDomain)%iEnd
      mask11 => level11(iDomain)%mask
      !
      ! initialize routing parameters (has to be called only for Routing option 2)
      if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) &
          call mrm_update_param(iDomain, parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)))
      ! calculate NtimeSteps for this domain
      nTimeSteps = (simPer(iDomain)%julEnd - simPer(iDomain)%julStart + 1) * NTSTEPDAY
      ! initialize timestep
      newTime = real(simPer(iDomain)%julStart, dp)
      ! initialize variable for runoff for routing
      allocate(RunToRout(e1 - s1 + 1))
      RunToRout = 0._dp
      ! ----------------------------------------
      ! loop over time
      ! ----------------------------------------
      hour = -timestep
      do tt = 1, nTimeSteps
        ! set discharge timestep
        iDischargeTS = ceiling(real(tt, dp) / real(NTSTEPDAY, dp))
        ! calculate current timestep
        is_new_month = .false.
        prev_month = month
        call caldat(int(newTime), yy = year, mm = month, dd = day)
        if (prev_month .ne. month) is_new_month = .true.
        ! -------------------------------------------------------------------
        ! PERFORM ROUTING
        ! -------------------------------------------------------------------
        !
        ! set input variables for routing
        if (processMatrix(8, 1) .eq. 1) then
          ! >>>
          ! >>> original Muskingum routing, executed every time
          ! >>>
          !
          ! initialize land cover year id
          Lcover_yID = LCyearId(year, iDomain)
          !
          do_rout = .True.
          L11_tsRout(iDomain) = (timestep * HourSecs)
          tsRoutFactorIn = 1._dp
          timestep_rout = timestep
          RunToRout = L1_total_runoff_in(s1 : e1, tt) ! runoff [mm TST-1] mm per timestep
          InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
          !
        else if ((processMatrix(8, 1) .eq. 2) .or. &
                 (processMatrix(8, 1) .eq. 3)) then
          ! >>>
          ! >>> adaptive timestep
          ! >>>
          !
          ! set dummy lcover_yID
          Lcover_yID = 1_i4
          !
          do_rout = .False.
          ! calculate factor
          tsRoutFactor = L11_tsRout(iDomain) / (timestep * HourSecs)
          ! print *, 'routing factor: ', tsRoutFactor
          ! prepare routing call
          if (tsRoutFactor .lt. 1._dp) then
            ! ----------------------------------------------------------------
            ! routing timesteps are shorter than hydrologic time steps
            ! ----------------------------------------------------------------
            ! set all input variables
            tsRoutFactorIn = tsRoutFactor
            RunToRout = L1_total_runoff_in(s1 : e1, tt) ! runoff [mm TST-1] mm per timestep
            InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
            timestep_rout = timestep
            do_rout = .True.
          else
            ! ----------------------------------------------------------------
            ! routing timesteps are longer than hydrologic time steps
            ! ----------------------------------------------------------------
            ! set all input variables
            tsRoutFactorIn = tsRoutFactor
            RunToRout = RunToRout + L1_total_runoff_in(s1 : e1, tt)
            InflowDischarge = InflowDischarge + InflowGauge%Q(iDischargeTS, :)
            ! reset tsRoutFactorIn if last period did not cover full period
            if ((tt .eq. nTimeSteps) .and. (mod(tt, nint(tsRoutFactorIn)) .ne. 0_i4)) &
                    tsRoutFactorIn = mod(tt, nint(tsRoutFactorIn))
            if ((mod(tt, nint(tsRoutFactorIn)) .eq. 0_i4) .or. (tt .eq. nTimeSteps)) then
              InflowDischarge = InflowDischarge / tsRoutFactorIn
              timestep_rout = nint(real(timestep, dp) * tsRoutFactorIn)
              do_rout = .True.
            end if
          end if
        end if
        ! -------------------------------------------------------------------
        ! execute routing
        ! -------------------------------------------------------------------
        if (do_rout) call mRM_routing(&
                ! general INPUT variables
                read_restart, &
                processMatrix(8, 1), & ! parse process Case to be used
                parameterset, & ! routing par.
                RunToRout, & ! runoff [mm TST-1] mm per timestep old: L1_total_runoff_in(s1:e1, tt), &
                level1(iDomain)%CellArea * 1.E-6_dp, &
                L1_L11_Id(s1 : e1), &
                level11(iDomain)%CellArea * 1.E-6_dp, &
                L11_L1_Id(s11 : e11), &
                L11_netPerm(s11 : e11), & ! routing order at L11
                L11_fromN(s11 : e11), & ! link source at L11
                L11_toN(s11 : e11), & ! link target at L11
                L11_nOutlets(iDomain), & ! number of outlets
                timestep_rout, & ! timestep of runoff to rout [h]
                tsRoutFactorIn, & ! Factor between routing and hydrologic resolution
                level11(iDomain)%nCells, & ! number of Nodes
                basin_mrm(iDomain)%nInflowGauges, &
                basin_mrm(iDomain)%InflowGaugeIndexList(:), &
                basin_mrm(iDomain)%InflowGaugeHeadwater(:), &
                basin_mrm(iDomain)%InflowGaugeNodeList(:), &
                InflowDischarge, &
                basin_mrm(iDomain)%nGauges, &
                basin_mrm(iDomain)%gaugeIndexList(:), &
                basin_mrm(iDomain)%gaugeNodeList(:), &
                ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)), &
                ! original routing specific input variables
                L11_length(s11 : e11 - 1), & ! link length
                L11_slope(s11 : e11 - 1), &
                L11_nLinkFracFPimp(s11 : e11, Lcover_yID), & ! fraction of impervious layer at L11 scale
                ! general INPUT/OUTPUT variables
                L11_C1(s11 : e11), & ! first muskingum parameter
                L11_C2(s11 : e11), & ! second muskigum parameter
                L11_qOUT(s11 : e11), & ! routed runoff flowing out of L11 cell
                L11_qTIN(s11 : e11, :), & ! inflow water into the reach at L11
                L11_qTR(s11 : e11, :), & !
                L11_qMod(s11 : e11), &
                mRM_runoff(tt, :) &
                )

        ! -------------------------------------------------------------------
        ! groundwater coupling
        ! -------------------------------------------------------------------
            if (gw_coupling) then
                call calc_river_head(iDomain, L11_Qmod, L0_river_head_mon_sum)
                if (is_new_month) then
                    call avg_and_write_timestep(iDomain, tt, L0_river_head_mon_sum)
                end if
            end if
        ! -------------------------------------------------------------------
        ! reset variables
        ! -------------------------------------------------------------------
        if (processMatrix(8, 1) .eq. 1) then
          ! reset Input variables
          InflowDischarge = 0._dp
          RunToRout = 0._dp
        else if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) then
          if ((.not. (tsRoutFactorIn .lt. 1._dp)) .and. do_rout) then
            do jj = 1, nint(tsRoutFactorIn)
              mRM_runoff(tt - jj + 1, :) = mRM_runoff(tt, :)
            end do
            ! reset Input variables
            InflowDischarge = 0._dp
            RunToRout = 0._dp
          end if
        end if
        ! -------------------------------------------------------------------
        ! INCREMENT TIME
        ! -------------------------------------------------------------------
        hour = mod(hour + timestep, 24)
        newTime = julday(day, month, year) + real(hour + timestep, dp) / 24._dp
        ! -------------------------------------------------------------------
        ! WRITE OUTPUT
        ! -------------------------------------------------------------------
        if (.not. optimize .and. any(outputFlxState_mrm)) then
          call mrm_write_output_fluxes(&
                  ! domain id
                  iDomain, &
                  ! nCells in domain
                  level11(iDomain)%nCells, &
                  ! output specification
                  timeStep_model_outputs_mrm, &
                  ! time specification
                  warmingDays(iDomain), newTime, nTimeSteps, nTStepDay, &
                  tt, day, month, year, timestep, &
                  ! mask specification
                  mask11, &
                  ! output variables
                  L11_qmod(s11 : e11))
        end if
      end do
      ! clean runoff variable
      deallocate(RunToRout)
    end do

    ! =========================================================================
    ! SET RUNOFF OUTPUT VARIABLE IF REQUIRED
    ! =========================================================================
    if (present(runoff)) runoff = mRM_runoff

    ! free memory
    deallocate(InflowDischarge)
  end subroutine mrm_eval

end module mo_mrm_eval
