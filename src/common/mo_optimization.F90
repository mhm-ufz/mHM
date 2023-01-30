!> \file mo_optimization.f90
!> \brief \copybrief mo_optimization
!> \details \copydetails mo_optimization

!> \brief Wrapper subroutine for optimization against runoff and sm.
!> \details This module provides a wrapper subroutine for optimization of mRM/mHM against runoff or soil moisture.
!> \authors Stephan Thober
!> \date Oct 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_optimization
  use mo_kind, only : i4, i8, dp
  use mo_optimization_utils, only : eval_interface, objective_interface
  implicit none
  private
  public :: optimization

contains
  ! ------------------------------------------------------------------

  !    NAME
  !        optimization

  !    PURPOSE
  !>       \brief Wrapper for optimization.

  !>       \details This subroutine selects the optimization defined in a namelist,
  !>       i.e. the global variable \e opti\_method.
  !>       It return the objective function value for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "procedure(eval_interface) :: eval"
  !>       \param[in] "procedure(objective_interface) :: objective" - objective function used in the optimization
  !>       \param[in] "character(len = *) :: dirConfigOut"          - directory where to write ascii output

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: funcbest"              - best objective function value obtained during optimization
  !>       \param[out] "logical, dimension(:) :: maskpara" true  = parameter will be optimized     = parameter(i,4) = 1
  !>       false = parameter will not be optimized = parameter(i,4) = 0

  !    HISTORY
  !>       \authors Matthias Cuntz, Luis Samaniego, Juliane Mai, Matthias Zink and Stephan Thober

  !>       \date Oct 2015

  ! Modifications:

  subroutine optimization(eval, objective, dirConfigOut, funcBest, maskpara)

    use mo_anneal, only : anneal
    use mo_common_mHM_mRM_variables, only : dds_r, mcmc_error_params, mcmc_opti, nIterations, opti_function, opti_method, &
                                            optimize_restart, sa_temp, sce_ngs, sce_npg, sce_nps, seed
    use mo_common_variables, only : global_parameters
#ifdef MPI
    use mo_common_variables, only : domainMeta
#endif
    use mo_dds, only : dds
    use mo_mcmc, only : mcmc, mcmc_stddev
    use mo_message, only : message, error_message
    use mo_sce, only : sce
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, &
                         timer_stop
    use mo_xor4096, only : get_timeseed

    implicit none

    procedure(eval_interface), intent(in), pointer :: eval

    ! - objective function used in the optimization
    procedure(objective_interface), intent(in), pointer :: objective

    ! - directory where to write ascii output
    character(len = *), intent(in) :: dirConfigOut

    ! - best objective function value obtained during optimization
    real(dp), intent(out) :: funcbest

    ! true  = parameter will be optimized     = parameter(i,4) = 1
    ! false = parameter will not be optimized = parameter(i,4) = 0
    logical, intent(out), allocatable, dimension(:) :: maskpara

    integer(i4) :: ii

    ! current timer number
    integer(i4) :: iTimer

    ! parameter sets sampled during burnin
    real(dp), allocatable, dimension(:, :) :: burnin_paras

    ! parameter sets sampled during proper mcmc
    real(dp), allocatable, dimension(:, :) :: mcmc_paras

    ! pre-determined stepsize
    ! real(dp), dimension(:), allocatable :: step

    integer(i4) :: npara

    ! global_parameters but includes a and b for likelihood
    real(dp), allocatable, dimension(:, :) :: local_parameters

    ! maskpara but includes a and b for likelihood
    logical, allocatable, dimension(:) :: local_maskpara

    ! local seed used for optimization
    integer(i8) :: iseed

    ! file for temporal optimization outputs
    character(256) :: tFile

    ! file for temporal SCE optimization outputs
    character(256) :: pFile


    ! -------------------------------------------------------------------------
    ! START
    ! -------------------------------------------------------------------------
    call message('  Start optimization')
    iTimer = 1
    call timer_start(iTimer)

    ! mask parameter which have a FLAG=0 in mhm_parameter.nml
    ! maskpara = true : parameter will be optimized
    ! maskpara = false : parameter is discarded during optimization
    npara = size(global_parameters, 1)
    allocate(maskpara(npara))
    maskpara = .true.
    do ii = 1, npara
      if (nint(global_parameters(ii, 4), i4) .eq. 0_i4) then
        maskpara(ii) = .false.
      end if
    end do

    ! add two extra parameter for optimisation of likelihood
    if (opti_function == 8) then
      allocate(local_parameters(npara + 2, size(global_parameters, 2)))
      allocate(local_maskpara(npara + 2))
      ! setting step sizes manually
      ! allocate(step(npara+2))
      local_parameters(1 : npara, :) = global_parameters(:, :)
      local_maskpara(1 : npara) = maskpara(:)
      local_parameters(npara + 1, 1) = 0.001_dp
      local_parameters(npara + 1, 2) = 100._dp
      local_parameters(npara + 1, 3) = 1._dp
      local_parameters(npara + 1, 4) = 1._dp
      local_parameters(npara + 1, 5) = 0._dp
      local_parameters(npara + 2, 1) = 0.001_dp
      local_parameters(npara + 2, 2) = 10._dp
      local_parameters(npara + 2, 3) = 0.1_dp
      local_parameters(npara + 2, 4) = 1._dp
      local_parameters(npara + 2, 5) = 0._dp
      local_maskpara(npara + 1 :) = .true.
      if ((opti_method == 0) .and. (.not. mcmc_opti)) then ! MCMC but only for parameter uncertainties
        local_parameters(npara + 1, 3) = mcmc_error_params(1)
        local_parameters(npara + 2, 3) = mcmc_error_params(2)
        local_maskpara(npara + 1 :) = .false.
      end if
    else
      allocate(local_parameters(npara, size(global_parameters, 2)))
      allocate(local_maskpara(npara))
      local_parameters = global_parameters
      local_maskpara = maskpara
    end if

    ! Seed for random numbers in optimisation
    if (seed .gt. 0_i8) then    ! fixed user-defined seed
      iseed = seed
    else                        ! flexible clock-time seed
      call get_timeseed(iseed)
    end if

    select case (opti_method)
    case (0)
      call message('    Use MCMC')

      tFile = trim(adjustl(dirConfigOut)) // 'mcmc_tmp_parasets.nc'

      ! setting step sizes manually
      ! step=(/ 1.00000000000000,  ... , /)

      select case (opti_function)
      case (8)
        call message('    Use MCMC')
        call mcmc(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), mcmc_paras, burnin_paras, &
                ParaSelectMode_in = 2_i4, tmp_file = tFile, &
                maskpara_in = local_maskpara, &
                restart = optimize_restart, restart_file = 'mo_mcmc.restart', &
                ! stepsize_in=step,                                                                                         &
                seed_in = iseed, loglike_in = .true., printflag_in = .true.)
      case (4)
        if (optimize_restart) then
          call error_message('ERROR: A restart of this optimization method is not implemented yet!')
        end if
        call message('    Use MCMC_STDDEV')
        call mcmc_stddev(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), mcmc_paras, burnin_paras, &
                ParaSelectMode_in = 2_i4, tmp_file = tFile, &
                maskpara_in = local_maskpara, &
                seed_in = iseed, loglike_in = .true., printflag_in = .true.)
      case default
        call error_message("Error objective: This opti_function is either not implemented yet.")
      end select

    case (1)
      call message('    Use DDS')

      tFile = trim(adjustl(dirConfigOut)) // 'dds_results.out'

      if (optimize_restart) then
        call error_message('ERROR: A restart of this optimization method is not implemented yet!')
      end if
      ! use fixed user-defined seed
#ifdef MPI
      local_parameters(:, 3) = dds(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), &
              maxiter = int(nIterations, i8), r = dds_r, seed = iseed, &
              tmp_file = tFile, comm = domainMeta%comMaster, mask = local_maskpara, &
              funcbest = funcbest)
#else
      local_parameters(:, 3) = dds(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), &
              maxiter = int(nIterations, i8), r = dds_r, seed = iseed, &
              tmp_file = tFile, mask = local_maskpara, &
              funcbest = funcbest)
#endif

    case (2)
      call message('    Use Simulated Annealing')

      tFile = trim(adjustl(dirConfigOut)) // 'anneal_results.out'

      if (optimize_restart) then
        call error_message('ERROR: A restart of this optimization method is not implemented yet!')
      end if

      if (sa_temp .gt. 0.0_dp) then
        ! use fixed user-defined seed and user-defined initial temperature
        local_parameters(:, 3) = anneal(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), &
                temp = sa_temp, seeds = (/iseed, iseed + 1000_i8, iseed + 2000_i8/), nITERmax = nIterations, &
                tmp_file = tFile, maskpara = local_maskpara, &
                funcbest = funcbest)
      else
        ! use fixed user-defined seed and adaptive initial temperature
        local_parameters(:, 3) = anneal(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), &
                seeds = (/iseed, iseed + 1000_i8, iseed + 2000_i8/), nITERmax = nIterations, &
                tmp_file = tFile, maskpara = local_maskpara, &
                funcbest = funcbest)
      end if
    case (3)
      call message('    Use SCE')

      tFile = trim(adjustl(dirConfigOut)) // 'sce_results.out'
      pFile = trim(adjustl(dirConfigOut)) // 'sce_population.out'

      ! use fixed user-defined seed
      local_parameters(:, 3) = sce(eval, objective, local_parameters(:, 3), local_parameters(:, 1 : 2), &
              mymaxn = int(nIterations, i8), myseed = iseed, myngs = sce_ngs, mynpg = sce_npg, mynps = sce_nps, &
              parallel = .false., mymask = local_maskpara, &
              restart = optimize_restart, restart_file = 'mo_sce.restart', &
#ifdef MPI
              comm = domainMeta%comMaster, &
#endif
              tmp_file = tFile, popul_file = pFile, &
              bestf = funcbest)
    case default
      call error_message('mRM', 'This optimization method is not implemented.')
    end select
    call timer_stop(iTimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

    global_parameters(:, :) = local_parameters(1 : npara, :)
    maskpara(:) = local_maskpara(1 : npara)

    deallocate(local_parameters)
    deallocate(local_maskpara)

  end subroutine optimization

end module mo_optimization
