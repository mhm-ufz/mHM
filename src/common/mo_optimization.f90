!> \file mo_optimization.f90

!> \brief Wrapper subroutine for optimization against runoff and sm.

!> \details This module provides a wrapper subroutine for optimization of mRM/mHM
!>          against runoff or soil moisture.\n

!> \authors Stephan Thober
!> \date Oct 2015
module mo_optimization
  implicit none
  public :: optimization
  private
contains
  ! ------------------------------------------------------------------

  !      NAME
  !          optimization

  !>        \brief Wrapper for optimization.

  !>        \details This subroutine selects the optimization defined in a namelist, 
  !>        i.e. the global variable \e opti\_method.\n
  !>        It return the objective function value for a specific parameter set.

  !     INTENT(IN)
  !>        \param[in] "real(dp)         :: objective"    - objective function used in the optimization
  !>        \param[in] "character(len=*) :: dirConfigOut" - directory where to write ascii output

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)             :: funcBest"     - best objective function value obtained during optimization
  !>        \param[out] "logical, allocatable :: maskpara(:) " - mask of optimized parameters

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !        None
  
  !     RESTRICTIONS
  !        None
  
  !     EXAMPLE
  !         None
  
  !     LITERATURE

  !     HISTORY
  !>        \author Matthias Cuntz, Luis Samaniego, Juliane Mai, Matthias Zink and Stephan Thober
  !>        \date Oct 2015

  subroutine optimization(objective, dirConfigOut, funcBest, maskpara)
    use mo_kind,                          only: i4, i8, dp
    use mo_anneal,                        only: anneal                             ! Optimize with Simulated Annealing SA
    use mo_dds,                           only: dds                                ! Optimize with Dynam. Dimens. Search DDS
    use mo_string_utils,                  only: num2str
    use mo_common_variables,              only: &
         opti_method,                           &                                  ! Optimization algorithm used
         opti_function,                         &                                  ! Objective function used
         optimize_restart,                      &                                  ! Optimization will be restarted from
         !                                                                         ! mo_<opti_method>.restart file (.true.)
         global_parameters,                     &                                  ! Matrix of global parameters (former: gamma)
         !                                                                         !     col1: min,  col2: max, col3: initial, 
         !                                                                         !     col4: flag, col5: scaling
         nIterations,                           &                                  ! number of iterations for optimization
         seed,                                  &                                  ! seed used for optimization
         mcmc_opti,                             &                                  ! MCMC: Optimization (.true. ) or
         !                                                                         !       Only parameter uncertainty (.false.)
         mcmc_error_params,                     &                                  !       Parameters of error model if mcmc_opti=.false.
         !                                                                         !       e.g. for opti_function=8: 0.01, 0.3
         sa_temp,                               &                                  ! SA:  initial temperature
         dds_r,                                 &                                  ! DDS: perturbation rate
         sce_ngs,                               &                                  ! SCE: # of complexes
         sce_npg,                               &                                  ! SCE: # of points per complex
         sce_nps                                                                   ! SCE: # of points per subcomplex
    USE mo_mcmc,                          only: mcmc, mcmc_stddev                  ! Monte Carlo Markov Chain method
    use mo_sce,                           only: sce                                ! Optimize with Shuffled Complex evolution
    use mo_message,                       only: message
    use mo_finish,                        only: finish
    use mo_timer,                         only: timer_start, timer_stop, timer_get ! Timing of processes
    use mo_xor4096,                       only: get_timeseed                       ! generating a seed from clock
    ! objective functions and likelihood for runoff only
    use mo_mrm_objective_function_runoff, only: loglikelihood, loglikelihood_stddev
    
    implicit none
    
    ! -------------------------------------------------------------------------
    ! INPUT VARIABLES
    ! -------------------------------------------------------------------------
    interface 
       function objective (pp)
         use mo_kind, only: dp
         implicit none
         real(dp), intent (in) :: pp(:)
         real(dp) :: objective
      end function objective
    end interface 
    character(len=*), intent(in) :: dirConfigOut
    
    ! -------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    ! -------------------------------------------------------------------------
    real(dp), intent(out)              :: funcbest    ! best objective function achivied during optimization
    logical,  intent(out), allocatable :: maskpara(:) ! true  = parameter will be optimized     = parameter(i,4) = 1
    !                                                 ! false = parameter will not be optimized = parameter(i,4) = 0
    
    ! -------------------------------------------------------------------------
    ! LOCAL VARIABLES
    ! -------------------------------------------------------------------------
    integer(i4)                             :: ii
    integer(i4)                             :: iTimer                ! current timer number
                                                                     ! mcmc
    real(dp), allocatable                   :: burnin_paras(:,:)     ! parameter sets sampled during burnin
    real(dp), allocatable                   :: mcmc_paras(:,:)       ! parameter sets sampled during proper mcmc
                                                                     ! setting step sizes manually
    ! real(dp), dimension(:),   allocatable :: step                  ! pre-determined stepsize 
    integer(i4)                             :: npara
    real(dp), allocatable                   :: local_parameters(:,:) ! global_parameters but includes a and b for likelihood
    logical,  allocatable                   :: local_maskpara(:)     ! maskpara but includes a and b for likelihood
    integer(i8)                             :: iseed                 ! local seed used for optimization
    character(256)                          :: tFile                 ! file for temporal optimization outputs
    character(256)                          :: pFile                 ! file for temporal SCE optimization outputs

    ! -------------------------------------------------------------------------
    ! START
    ! -------------------------------------------------------------------------
    call message('  Start optimization')
    iTimer = 1
    call timer_start(iTimer)

    ! mask parameter which have a FLAG=0 in mhm_parameter.nml
    ! maskpara = true : parameter will be optimized
    ! maskpara = false : parameter is discarded during optimization
    npara = size(global_parameters,1)
    allocate(maskpara(npara))
    maskpara = .true.
    do ii=1, npara
       if ( nint(global_parameters(ii,4),i4) .eq. 0_i4 ) then
          maskpara(ii) = .false.
       end if
    end do

    ! add two extra parameter for optimisation of likelihood
    if (opti_function == 8) then
       allocate(local_parameters(npara+2,size(global_parameters,2)))
       allocate(local_maskpara(npara+2))
       ! setting step sizes manually
       ! allocate(step(npara+2))
       local_parameters(1:npara,:) = global_parameters(:,:)
       local_maskpara(1:npara)     = maskpara(:)
       local_parameters(npara+1,1) = 0.001_dp
       local_parameters(npara+1,2) = 100._dp
       local_parameters(npara+1,3) = 1._dp
       local_parameters(npara+1,4) = 1._dp
       local_parameters(npara+1,5) = 0._dp
       local_parameters(npara+2,1) = 0.001_dp
       local_parameters(npara+2,2) = 10._dp
       local_parameters(npara+2,3) = 0.1_dp
       local_parameters(npara+2,4) = 1._dp
       local_parameters(npara+2,5) = 0._dp
       local_maskpara(npara+1:)    = .true.
       if ((opti_method == 0) .and. (.not. mcmc_opti)) then ! MCMC but only for parameter uncertainties
          local_parameters(npara+1,3) = mcmc_error_params(1)
          local_parameters(npara+2,3) = mcmc_error_params(2)
          local_maskpara(npara+1:)    = .false.
       endif
    else
       allocate(local_parameters(npara,size(global_parameters,2)))
       allocate(local_maskpara(npara))
       local_parameters = global_parameters
       local_maskpara   = maskpara
    endif

    ! Seed for random numbers in optimisation
    if (seed .gt. 0_i8) then    ! fixed user-defined seed
       iseed = seed
    else                        ! flexible clock-time seed
       call get_timeseed(iseed)
    endif

    select case (opti_method)
    case (0)
       call message('    Use MCMC')

       tFile = trim(adjustl(dirConfigOut)) // 'mcmc_tmp_parasets.nc'

       ! setting step sizes manually
       ! step=(/ 1.00000000000000,  ... , /)

       if (opti_function == 8) then
          call message('    Use MCMC')
          call mcmc(loglikelihood, local_parameters(:,3), local_parameters(:,1:2), mcmc_paras, burnin_paras,               &
               ParaSelectMode_in=2_i4, tmp_file=tFile,                                                                     &
               maskpara_in=local_maskpara,                                                                                 &
               restart=optimize_restart, restart_file='mo_mcmc.restart',                                                   &
               ! stepsize_in=step,                                                                                         &
               seed_in=iseed, loglike_in=.true., printflag_in=.true.)
       else
          if (optimize_restart) then
             call message('ERROR: A restart of this optimization method is not implemented yet!')
             stop
          end if
          call message('    Use MCMC_STDDEV')
          call mcmc_stddev(loglikelihood_stddev, local_parameters(:,3), local_parameters(:,1:2), mcmc_paras, burnin_paras, &
               ParaSelectMode_in=2_i4, tmp_file=tFile,                                                                     &
               maskpara_in=local_maskpara,                                                                                 &
               seed_in=iseed, loglike_in=.true., printflag_in=.true.)
       endif

    case (1)
       call message('    Use DDS')

       tFile = trim(adjustl(dirConfigOut)) // 'dds_results.out'

       if (optimize_restart) then
          call message('ERROR: A restart of this optimization method is not implemented yet!')
          stop
       end if
       
       ! use fixed user-defined seed
       local_parameters(:,3) = dds(objective, local_parameters(:,3), local_parameters(:,1:2),                              &
            maxiter=int(nIterations,i8), r=dds_r, seed=iseed,                                                              &
            tmp_file=tFile, mask=local_maskpara,                                                                           &
            funcbest=funcbest)
    case (2)
       call message('    Use Simulated Annealing')

       tFile = trim(adjustl(dirConfigOut)) // 'anneal_results.out'

       if (optimize_restart) then
          call message('ERROR: A restart of this optimization method is not implemented yet!')
          stop
       end if

       if (sa_temp .gt. 0.0_dp) then
          ! use fixed user-defined seed and user-defined initial temperature
          local_parameters(:,3) = anneal(objective, local_parameters(:,3), local_parameters(:,1:2),                        &
               temp=sa_temp, seeds=(/iseed, iseed+1000_i8, iseed+2000_i8/), nITERmax=nIterations,                          &
               tmp_file=tFile, maskpara=local_maskpara,                                                                    &
               funcbest=funcbest)
       else
          ! use fixed user-defined seed and adaptive initial temperature
          local_parameters(:,3) = anneal(objective, local_parameters(:,3), local_parameters(:,1:2),                        &
               seeds=(/iseed, iseed+1000_i8, iseed+2000_i8/), nITERmax=nIterations,                                        &
               tmp_file=tFile, maskpara=local_maskpara,                                                                    &
               funcbest=funcbest)
       end if
    case (3)
       call message('    Use SCE')

       tFile = trim(adjustl(dirConfigOut)) // 'sce_results.out'
       pFile =  trim(adjustl(dirConfigOut)) // 'sce_population.out'

       ! use fixed user-defined seed
       local_parameters(:,3) = sce(objective, local_parameters(:,3), local_parameters(:,1:2),                              &
            mymaxn=int(nIterations,i8), myseed=iseed, myngs=sce_ngs, mynpg=sce_npg, mynps=sce_nps,                         &
            parallel=.false., mymask=local_maskpara,                                                                       &
            restart=optimize_restart, restart_file='mo_sce.restart',                                                       &
            tmp_file=tFile, popul_file=pFile,                                                                              &
            bestf=funcbest)
    case default
       call finish('mRM','This optimization method is not implemented.')
    end select
    call timer_stop(iTimer)
    call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

    global_parameters(:,:) = local_parameters(1:npara,:)
    maskpara(:)    = local_maskpara(1:npara)

    deallocate(local_parameters)
    deallocate(local_maskpara)

  end subroutine optimization
  
end module mo_optimization
