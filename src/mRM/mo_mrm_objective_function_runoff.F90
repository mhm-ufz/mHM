!> \file mo_mrm_objective_function_runoff.f90
!> \brief \copybrief mo_mrm_objective_function_runoff
!> \details \copydetails mo_mrm_objective_function_runoff

!> \brief Objective Functions for Optimization of mHM/mRM against runoff.
!> \details This module provides a wrapper for several objective functions used to optimize mRM/mHM against
!> runoff.
!!
!> If the objective contains besides runoff another variable like TWS move it to mHM/mo_objective_function.f90.
!> If it is only regarding runoff implement it here.
!!
!! All the objective functions are supposed to be minimized!
!! 1.  SO: Q:        1.0 - NSE
!! 2.  SO: Q:        1.0 - lnNSE
!! 3.  SO: Q:        1.0 - 0.5*(NSE+lnNSE)
!! 4.  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
!! 5.  SO: Q:        ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
!! 6.  SO: Q:        SSE
!! 7.  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors
!! 8.  SO: Q:       -1.0 * loglikelihood with trend removed from the relative errors and then lag(1)-autocorrelation removed
!! 9.  SO: Q:        1.0 - KGE (Kling-Gupta efficiency measure)
!! 14. SO: Q:        sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) > combination of KGE of every gauging station based on a power-6 norm
!! 16. MO: Q:        1st objective: 1.0 - NSE
!!         Q:        2nd objective: 1.0 - lnNSE
!! 18. MO: Q:        1st objective: 1.0 - lnNSE(Q_highflow)  (95% percentile)
!!         Q:        2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)
!! 19. MO: Q:        1st objective: 1.0 - lnNSE(Q_highflow)  (non-low flow)
!!         Q:        2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)eshold for Q
!! 20. MO: Q:        1st objective: absolute difference in FDC's low-segment volume
!!         Q:        2nd objective: 1.0 - NSE of discharge of months DJF
!! 31. SO: Q:        1.0 - wNSE - weighted NSE
!! 32. SO: Q:        SSE of boxcox-transformed streamflow
!> \changelog
!! - Stephan Thober             Oct 2015
!!   - adapted for mRM
!! - Juliane Mai                Nov 2015
!!   - introducing multi
!!   - and single-objective
!!   - first multi-objective function (16), but not used yet
!! - Juliane Mai                Feb 2016
!!   - multi-objective function (18) using lnNSE(highflows) and lnNSE(lowflows)
!!   - multi-objective function (19) using lnNSE(highflows) and lnNSE(lowflows)
!!   - multi-objective function (20) using FDC and discharge of months DJF
!! - Stephan Thober,Bjoern Guse May 2018
!!   - single objective function (21) using weighted NSE following (Hundecha and Bardossy, 2004)
!! - Robert Schweppe            Jun 2018
!!   - refactoring and reformatting
!! - Stephan Thober             Aug 2019
!!   - added OF 32: SSE of boxcox-transformed streamflow
!> \authors Juliane Mai
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
MODULE mo_mrm_objective_function_runoff

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 - removed all none runoff objectives,
  !                                     these can be found mo_objective_functions_sm

  USE mo_kind, ONLY : i4, dp
  use mo_optimization_utils, only : eval_interface
  use mo_message, only: message, error_message
  use mo_string_utils, only : num2str

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: single_objective_runoff ! single-objective function wrapper
  PUBLIC :: multi_objective_runoff  ! multi-objective function wrapper
  PUBLIC :: extract_runoff          ! extract runoff period specified in mhm.nml from available runoff time series
#ifdef MPI
  PUBLIC :: single_objective_runoff_master, single_objective_runoff_subprocess ! objective function wrapper for soil moisture only
#endif

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        single_objective_runoff

  !    PURPOSE
  !>       \brief Wrapper for objective functions optimizing agains runoff.

  !>       \details The functions selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It return the objective function value for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp), DIMENSION(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: arg1"

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: arg2"
  !>       \param[out] "real(dp), optional :: arg3"

  !    RETURN
  !>       \return real(dp) :: objective &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober              Oct 2015 - only runoff objective functions
  ! Stephan Thober, Bjoern Guse May 2018 - added weighted objective function
  ! Robert Schweppe             Jun 2018 - refactoring and reformatting


  FUNCTION single_objective_runoff(parameterset, eval, arg1, arg2, arg3)

    use mo_common_mHM_mRM_variables, only : opti_function, opti_method

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: single_objective_runoff


    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
    case (1)
      ! 1.0-nse
      single_objective_runoff = objective_nse(parameterset, eval)
    case (2)
      ! 1.0-lnnse
      single_objective_runoff = objective_lnnse(parameterset, eval)
    case (3)
      ! 1.0-0.5*(nse+lnnse)
      single_objective_runoff = objective_equal_nse_lnnse(parameterset, eval)
    case (4)
      if (opti_method .eq. 0_i4) then
        ! MCMC
        single_objective_runoff = loglikelihood_stddev(parameterset, eval, arg1, arg2, arg3)
      else
        ! -loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
        single_objective_runoff = - loglikelihood_stddev(parameterset, eval, 1.0_dp)
      end if
    case (5)
      ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
      single_objective_runoff = objective_power6_nse_lnnse(parameterset, eval)
    case (6)
      ! SSE
      single_objective_runoff = objective_sse(parameterset, eval)
    case (7)
      ! -loglikelihood with trend removed from absolute errors
      single_objective_runoff = -loglikelihood_trend_no_autocorr(parameterset, eval, 1.0_dp)
    case (8)
      if (opti_method .eq. 0_i4) then
        ! MCMC
        single_objective_runoff = loglikelihood_evin2013_2(parameterset, eval, regularize = .true.)
      else
        ! -loglikelihood of approach 2 of Evin et al. (2013),
        !  i.e. linear error model with lag(1)-autocorrelation on relative errors
        single_objective_runoff = -loglikelihood_evin2013_2(parameterset, eval)
      end if

    case (9)
      ! KGE
      single_objective_runoff = objective_kge(parameterset, eval)
    case (14)
      ! combination of KGE of every gauging station based on a power-6 norm \n
      ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6)
      single_objective_runoff = objective_multiple_gauges_kge_power6(parameterset, eval)
    case (31)
      ! weighted NSE with observed streamflow
      single_objective_runoff = objective_weighted_nse(parameterset, eval)
    case (32)
      ! sum of squared errors (SSE) of boxcox_transformed streamflow
      single_objective_runoff = objective_sse_boxcox(parameterset, eval)
    case default
      call error_message("Error objective: This opti_function is either not implemented yet or is not a single-objective one.")
    end select

  END FUNCTION single_objective_runoff


  ! ------------------------------------------------------------------

  !    NAME
  !        single_objective_runoff_master

  !    PURPOSE
  !>       \brief Wrapper for objective functions optimizing agains runoff.

  !>       \details The functions selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It return the objective function value for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp), DIMENSION(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: arg1"

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: arg2"
  !>       \param[out] "real(dp), optional :: arg3"

  !    RETURN
  !>       \return real(dp) :: objective &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober              Oct 2015 - only runoff objective functions
  ! Stephan Thober, Bjoern Guse May 2018 - added weighted objective function
  ! Robert Schweppe             Jun 2018 - refactoring and reformatting

#ifdef MPI
  FUNCTION single_objective_runoff_master(parameterset, eval, arg1, arg2, arg3)

    use mo_common_mHM_mRM_variables, only : opti_function, opti_method
    use mo_common_mpi_tools, only : distribute_parameterset
    use mo_mrm_global_variables, only: nGaugesTotal
    use mo_common_variables, only : domainMeta
    use mpi_f08

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: single_objective_runoff_master

    REAL(dp) :: partial_objective

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    integer(i4) :: iproc, nproc

    integer(i4) :: ierror

    type(MPI_Status) :: status

    !write(*,*) 'parameterset: ',parameterset(:)
    call distribute_parameterset(parameterset)
    select case (opti_function)
    case(1 : 3, 5, 6, 9, 31)
      call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
      single_objective_runoff_master = 0.0_dp
      do iproc = 1, nproc - 1
        call MPI_Recv(partial_objective, 1, MPI_DOUBLE_PRECISION, iproc, 0, domainMeta%comMaster, status, ierror)
        single_objective_runoff_master = single_objective_runoff_master + partial_objective
      end do
      single_objective_runoff_master = 1.0_dp - single_objective_runoff_master / real(nGaugesTotal, dp)
    case(14)
      call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
      single_objective_runoff_master = 0.0_dp
      do iproc = 1, nproc - 1
        call MPI_Recv(partial_objective, 1, MPI_DOUBLE_PRECISION, iproc, 0, domainMeta%comMaster, status, ierror)
        single_objective_runoff_master = single_objective_runoff_master + partial_objective
      end do
      single_objective_runoff_master = single_objective_runoff_master**onesixth
    case(4, 7, 8)
      call message("case 4, 7, 8 are not implemented in parallel yet")
    case default
      call error_message("Error single_objective_runoff_master:", raise=.false.)
      call error_message("This opti_function is either not implemented yet or is not a single-objective one.")
    end select

    select case (opti_function)
    case (1)
      ! 1.0-nse
      call message('objective_nse (i.e., 1 - NSE) = ', num2str(single_objective_runoff_master))
    case (2)
      ! 1.0-lnnse
      call message('objective_lnnse = ', num2str(single_objective_runoff_master))
    case (3)
      ! 1.0-0.5*(nse+lnnse)
      call message('objective_equal_nse_lnnse = ', num2str(single_objective_runoff_master))
    case (4)
      call error_message("case 4, loglikelihood_stddev not implemented in parallel yet")
    case (5)
      ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
      call message('objective_power6_nse_lnnse = ', num2str(single_objective_runoff_master))
    case (6)
      ! SSE
      call message('objective_sse = ', num2str(single_objective_runoff_master))
    case (7)
      ! -loglikelihood with trend removed from absolute errors
      call error_message("case 7, single_objective_runoff_master not implemented in parallel yet")
    case (8)
      call error_message("case 8, loglikelihood_evin2013_2 not implemented in parallel yet")
    case (9)
      ! KGE
      call message('objective_kge (i.e., 1 - KGE) = ', num2str(single_objective_runoff_master))
    case (14)
      ! combination of KGE of every gauging station based on a power-6 norm \n
      ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6)
      call message('objective_multiple_gauges_kge_power6 = ', num2str(single_objective_runoff_master))
    case (31)
      ! weighted NSE with observed streamflow
      call message('objective_weighted_nse (i.e., 1 - wNSE) = ', num2str(single_objective_runoff_master))
    case (32)
      ! SSE of boxcox-transformed streamflow
      call message('sse_boxcox_streamflow = ', num2str(single_objective_runoff_master))
    case default
      call error_message("Error single_objective_runoff_master:", raise=.false.)
      call error_message("This opti_function is either not implemented yet or is not a single-objective one.", raise=.false.)
      call error_message("This part of the code should never be executed.")
    end select

  END FUNCTION single_objective_runoff_master


  ! ------------------------------------------------------------------

  !    NAME
  !        single_objective_runoff_subprocess

  !    PURPOSE
  !>       \brief Wrapper for objective functions optimizing agains runoff.

  !>       \details The functions selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It return the objective function value for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp), DIMENSION(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: arg1"

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: arg2"
  !>       \param[out] "real(dp), optional :: arg3"

  !    RETURN
  !>       \return real(dp) :: objective &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober              Oct 2015 - only runoff objective functions
  ! Stephan Thober, Bjoern Guse May 2018 - added weighted objective function
  ! Robert Schweppe             Jun 2018 - refactoring and reformatting


  subroutine single_objective_runoff_subprocess(eval, arg1, arg2, arg3)

    use mo_common_mHM_mRM_variables, only : opti_function, opti_method
    use mo_common_mpi_tools, only : get_parameterset
    use mo_common_variables, only : domainMeta
    use mpi_f08

    implicit none

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: partial_single_objective_runoff

    REAL(dp), DIMENSION(:), allocatable :: parameterset

    integer(i4) :: ierror

    type(MPI_Status) :: status

    logical :: do_obj_loop


    do ! a do loop without condition runs until exit
      call MPI_Recv(do_obj_loop, 1, MPI_LOGICAL, 0, 0, domainMeta%comMaster, status, ierror)

      if (.not. do_obj_loop) exit
      !write(*,*) 'parameterset: ',parameterset(:)
      call get_parameterset(parameterset)
      select case (opti_function)
      case (1)
        ! 1.0-nse
        partial_single_objective_runoff = objective_nse(parameterset, eval)
      case (2)
        ! 1.0-lnnse
        partial_single_objective_runoff = objective_lnnse(parameterset, eval)
      case (3)
        ! 1.0-0.5*(nse+lnnse)
        partial_single_objective_runoff = objective_equal_nse_lnnse(parameterset, eval)
      case (4)
        if (opti_method .eq. 0_i4) then
          ! MCMC
          ! partial_single_objective_runoff = loglikelihood_stddev(parameterset, eval, arg1, arg2, arg3)
          call error_message("Error single_objective_runoff_subprocess: case 4 with optimethod 0 not supported")
        else
          ! -loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
          partial_single_objective_runoff = - loglikelihood_stddev(parameterset, eval, 1.0_dp)
        end if
      case (5)
        ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
        partial_single_objective_runoff = objective_power6_nse_lnnse(parameterset, eval)
      case (6)
        ! SSE
        partial_single_objective_runoff = objective_sse(parameterset, eval)
      case (7)
        ! -loglikelihood with trend removed from absolute errors
        ! partial_single_objective_runoff = -loglikelihood_trend_no_autocorr(parameterset, eval, 1.0_dp)
        call error_message("Error single_objective_runoff_subprocess: case 7 not supported")
      case (8)
        if (opti_method .eq. 0_i4) then
          ! MCMC
          ! partial_single_objective_runoff = loglikelihood_evin2013_2(parameterset, eval, regularize = .true.)
          call error_message("Error single_objective_runoff_subprocess: case 8 with optimethod 0 not supported")
        else
          ! -loglikelihood of approach 2 of Evin et al. (2013),
          !  i.e. linear error model with lag(1)-autocorrelation on relative errors
          partial_single_objective_runoff = -loglikelihood_evin2013_2(parameterset, eval)
        end if

      case (9)
        ! KGE
        partial_single_objective_runoff = objective_kge(parameterset, eval)
      case (14)
        ! combination of KGE of every gauging station based on a power-6 norm \n
        ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6)
        partial_single_objective_runoff = objective_multiple_gauges_kge_power6(parameterset, eval)
      case (31)
         ! weighted NSE with observed streamflow
         partial_single_objective_runoff = objective_weighted_nse(parameterset, eval)
      case (32)
         ! SSE of transformed streamflow
         partial_single_objective_runoff = objective_sse_boxcox(parameterset, eval)
      case default
        call error_message("Error single_objective_runoff_subprocess:", raise=.false.)
        call error_message("This opti_function is either not implemented yet or is not a single-objective one.")
      end select

      select case (opti_function)
      case (1 : 3, 5, 6, 9, 14, 31)
        call MPI_Send(partial_single_objective_runoff,1, MPI_DOUBLE_PRECISION,0,0,domainMeta%comMaster,ierror)
      case default
        call error_message("Error objective_subprocess: this part should not be executed -> error in the code.")
      end select
      deallocate(parameterset)
    end do

  END subroutine single_objective_runoff_subprocess

#endif

  ! ------------------------------------------------------------------

  !    NAME
  !        multi_objective_runoff

  !    PURPOSE
  !>       \brief Wrapper for multi-objective functions where at least one is regarding runoff.

  !>       \details The functions selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It return the multiple objective function values for a specific parameter set.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp), DIMENSION(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp), DIMENSION(:) :: multi_objectives"

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Oct 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE multi_objective_runoff(parameterset, eval, multi_objectives)

    use mo_common_mHM_mRM_variables, only : opti_function

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: multi_objectives


    select case (opti_function)
    case (16)
      ! 1st objective: 1.0-nse
      ! 2nd objective: 1.0-lnnse
      multi_objectives = multi_objective_nse_lnnse(parameterset, eval)
    case (18)
      ! 1st objective: 1.0 - lnNSE(Q_highflow)  (95% percentile)
      ! 2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)
      multi_objectives = multi_objective_lnnse_highflow_lnnse_lowflow(parameterset, eval)
    case (19)
      ! 1st objective: 1.0 - lnNSE(Q_highflow)  (non-low flow)
      ! 2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)
      multi_objectives = multi_objective_lnnse_highflow_lnnse_lowflow_2(parameterset, eval)
    case (20)
      ! 1st objective: 1.0 - difference in FDC's low-segment volume
      ! 2nd objective: 1.0 - NSE of discharge of months DJF
      multi_objectives = multi_objective_ae_fdc_lsv_nse_djf(parameterset, eval)
    case default
      call error_message("Error objective: Either this opti_function is not implemented yet or it is not a multi-objective one.")
    end select

  END SUBROUTINE multi_objective_runoff

  ! ------------------------------------------------------------------

  !    NAME
  !        loglikelihood_stddev

  !    PURPOSE
  !>       \brief Logarithmic likelihood function with removed linear trend and Lag(1)-autocorrelation.

  !>       \details The logarithmis likelihood function is used when mHM runs in MCMC mode.
  !>       It can also be used for optimization when selecting the likelihood in the
  !>       namelist as \e opti\_function.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"
  !>       \param[in] "real(dp) :: stddev"                     standard deviation of data

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: stddev_new" standard deviation of errors with removed trend and
  !>       correlationbetween model run using parameter set and observation
  !>       \param[out] "real(dp), optional :: likeli_new" logarithmic likelihood determined with stddev_new instead of
  !>       stddev

  !    RETURN
  !>       \return real(dp) :: loglikelihood_stddev &mdash; logarithmic likelihood using given stddev
  !>       but remove optimal trend and lag(1)-autocorrelation in errors
  !>       (absolute between running model with parameterset and observation)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION loglikelihood_stddev(parameterset, eval, stddev, stddev_new, likeli_new)

    use mo_append, only : append
    use mo_linfit, only : linfit
    use mo_moment, only : correlation, mean

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    ! standard deviation of data
    real(dp), intent(in) :: stddev

    ! standard deviation of errors with removed trend and correlationbetween model run using parameter set and
    ! observation
    real(dp), intent(out), optional :: stddev_new

    ! logarithmic likelihood determined with stddev_new instead of stddev
    real(dp), intent(out), optional :: likeli_new

    real(dp) :: loglikelihood_stddev

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), dimension(:, :), allocatable :: runoff

    ! gauges counter
    integer(i4) :: gg

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    integer(i4) :: nmeas

    real(dp), dimension(:), allocatable :: errors

    real(dp), dimension(:), allocatable :: obs, calc, out

    real(dp) :: a, b, c

    real(dp) :: stddev_tmp


    call eval(parameterset, runoff = runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim = 2) ! second dimension equals nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! append it to variables
      call append(obs, runoff_obs)
      call append(calc, runoff_agg)

    end do
    ! ----------------------------------------

    nmeas = size(obs, dim = 1)

    allocate(out(nmeas), errors(nmeas))
    errors(:) = abs(calc(:) - obs(:))

    ! remove linear trend of errors - must be model NOT obs
    ! out = linfit(obs, errors, a=a, b=b, model2=.False.)
    ! errors(:) = errors(:) - (a + obs(:)*b)
    out = linfit(calc, errors, a = a, b = b, model2 = .False.)
    errors(:) = errors(:) - (a + calc(:) * b)

    ! remove lag(1)-autocorrelation of errors
    c = correlation(errors(2 : nmeas), errors(1 : nmeas - 1))
    errors(1 : nmeas - 1) = errors(2 : nmeas) - c * errors(1 : nmeas - 1)
    errors(nmeas) = 0.0_dp

    ! you have to take stddev=const because otherwise loglikelihood is always N
    ! in MCMC stddev gets updated only when a better likelihood is found.
    loglikelihood_stddev = sum(errors(:) * errors(:) / stddev**2)
    loglikelihood_stddev = -0.5_dp * loglikelihood_stddev

    call message('-loglikelihood_stddev = ', num2str(-loglikelihood_stddev))

    stddev_tmp = sqrt(sum((errors(:) - mean(errors)) * (errors(:) - mean(errors))) / real(nmeas - 1, dp))
    if (present(stddev_new)) then
      stddev_new = stddev_tmp
    end if
    if (present(likeli_new)) then
      likeli_new = sum(errors(:) * errors(:) / stddev_tmp**2)
      likeli_new = -0.5_dp * likeli_new
    end if

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs)
    deallocate(obs, calc, out, errors)

  END FUNCTION loglikelihood_stddev

  ! ------------------------------------------------------------------

  !    NAME
  !        loglikelihood_evin2013_2

  !    PURPOSE
  !>       \brief Logarithmised likelihood with linear error model and lag(1)-autocorrelation
  !>       of the relative errors.

  !>       \details This loglikelihood uses a linear error model and a lag(1)-autocorrelation
  !>       on the relative errors. This is approach 2 of the paper Evin et al. (WRR, 2013).

  !>       This is opti_function = 8.

  !>       mHM then adds two extra (local) parameters for the error model in mhm_driver,
  !>       which get optimised together with the other, global parameters.
  !>       ADDITIONAL INFORMATION
  !>       Evin et al., WRR 49, 4518-4524, 2013

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: regularize"

  !    RETURN
  !>       \return real(dp) :: loglikelihood_evin2013_2 &mdash; logarithmic likelihood using given stddev
  !>       but remove optimal trend and lag(1)-autocorrelation in errors
  !>       (absolute between running model with parameterset and observation)

  !    HISTORY
  !>       \authors Juliane Mai and Matthias Cuntz

  !>       \date Mar 2014

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION loglikelihood_evin2013_2(parameterset, eval, regularize)

    use mo_append, only : append
    use mo_common_variables, only : global_parameters
    use mo_constants, only : pi_dp
    use mo_moment, only : correlation
    use mo_utils, only : eq

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    logical, optional, intent(in) :: regularize

    real(dp) :: loglikelihood_evin2013_2

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), dimension(:, :), allocatable :: runoff

    ! penalty term due to a parmeter set out of bound
    real(dp) :: penalty

    ! gauges counter
    integer(i4) :: gg

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    integer(i4) :: nmeas

    real(dp), dimension(:), allocatable :: errors, sigma, eta, y

    real(dp), dimension(:), allocatable :: obs, calc, out

    real(dp) :: a, b, c, vary, vary1, ln2pi, tmp

    integer(i4) :: npara

    logical :: iregularize


    iregularize = .false.
    if (present(regularize)) iregularize = regularize

    npara = size(parameterset)
    call eval(parameterset(1 : npara - 2), runoff = runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim = 2) ! second dimension equals nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! append it to variables
      call append(obs, runoff_obs)
      call append(calc, runoff_agg)

    end do

    ! ----------------------------------------

    nmeas = size(obs, dim = 1)

    allocate(out(nmeas), errors(nmeas), sigma(nmeas), eta(nmeas), y(nmeas))
    ! residual errors
    errors(:) = calc(:) - obs(:)

    ! linear error model
    a = parameterset(npara - 1)
    b = parameterset(npara)
    sigma(:) = a + b * calc(:)
    ! standardized residual errors (SRE)
    eta(:) = errors(:) / sigma(:)

    ! remove lag(1)-autocorrelation of SRE
    c = correlation(eta(2 : nmeas), eta(1 : nmeas - 1))
    y(1) = 0.0_dp ! only for completeness
    y(2 : nmeas) = eta(2 : nmeas) - c * eta(1 : nmeas - 1)

    ! likelihood of residual errors
    ln2pi = log(sqrt(2.0_dp * pi_dp))
    vary = 1.0_dp - c * c
    vary1 = 1.0_dp / vary
    loglikelihood_evin2013_2 = -ln2pi - 0.5_dp * eta(1) * eta(1) - log(sigma(1)) & ! li(eta(1))/sigma(1)
            - real(nmeas - 1, dp) * log(sqrt(2.0_dp * pi_dp * vary)) &
            - sum(0.5_dp * y(2 : nmeas) * y(2 : nmeas) * vary1) - sum(log(sigma(2 : nmeas)))

    if (iregularize) then
      ! Regularistion term as deviation from initial parameter value
      penalty = parameter_regularization(&
              parameterset(1 : npara - 2), &        ! current parameter set
              global_parameters(1 : npara - 2, 3), &        ! prior/initial parameter set
              global_parameters(1 : npara - 2, 1 : 2), &        ! bounds
              eq(global_parameters(1 : npara - 2, 4), 1.0_dp)) ! used/unused

      tmp = loglikelihood_evin2013_2 + penalty
      call message( &
        '-loglikelihood_evin2013_2, + penalty, chi^2: ', &
        num2str(-loglikelihood_evin2013_2), &
        num2str(-tmp), &
        num2str(-tmp / real(nmeas, dp)))
      loglikelihood_evin2013_2 = tmp
    else
      call message( &
        '-loglikelihood_evin2013_2, chi^2: ', &
        num2str(-loglikelihood_evin2013_2), &
        num2str(-loglikelihood_evin2013_2 / real(nmeas, dp)))
    end if

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs)
    deallocate(obs, calc, out, errors, sigma, eta, y)

  END FUNCTION loglikelihood_evin2013_2

  !    NAME
  !        parameter_regularization

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: paraset"
  !>       \param[in] "real(dp), dimension(size(paraset)) :: prior"
  !>       \param[in] "real(dp), dimension(size(paraset), 2) :: bounds" (min, max)
  !>       \param[in] "logical, dimension(size(paraset)) :: mask"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  FUNCTION parameter_regularization(paraset, prior, bounds, mask)

    use mo_constants, only : pi_dp

    implicit none

    real(dp), dimension(:), intent(in) :: paraset

    real(dp), dimension(size(paraset)), intent(in) :: prior

    ! (min, max)
    real(dp), dimension(size(paraset), 2), intent(in) :: bounds

    logical, dimension(size(paraset)), intent(in) :: mask

    real(dp) :: parameter_regularization

    integer(i4) :: ipara

    integer(i4) :: npara

    real(dp), parameter :: onetwelveth = 1._dp / 12._dp

    real(dp), dimension(size(paraset)) :: sigma


    npara = size(paraset, 1)

    sigma = sqrt(onetwelveth * (bounds(:, 2) - bounds(:, 1))**2) ! standard deviation of uniform distribution
    parameter_regularization = -sum(log(sqrt(2.0_dp * pi_dp) * sigma), mask = mask)

    do ipara = 1, npara
      if (mask(ipara)) then
        ! if ((paraset(ipara) .lt. bounds(ipara,1)) .or. (paraset(ipara) .gt. bounds(ipara,2))) then
        !    ! outside bounds
        parameter_regularization = parameter_regularization - &
                0.5_dp * ((paraset(ipara) - prior(ipara)) / sigma(ipara))**2
        ! else
        !    ! in bound
        !    parameter_regularization = 0.0_dp
        ! end if
      end if
    end do

  END FUNCTION parameter_regularization

  ! ------------------------------------------------------------------

  !    NAME
  !        loglikelihood_trend_no_autocorr

  !    PURPOSE
  !>       \brief Logarithmic likelihood function with linear trend removed.

  !>       \details The logarithmis likelihood function is used when mHM runs in MCMC mode.
  !>       It can also be used for optimization when selecting the likelihood in the
  !>       namelist as \e opti\_function.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"
  !>       \param[in] "real(dp) :: stddev_old"                 standard deviation of data

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), optional :: stddev_new" standard deviation of errors with removed trendbetween model
  !>       run using parameter set and observation
  !>       \param[out] "real(dp), optional :: likeli_new" logarithmic likelihood determined with stddev_new instead of
  !>       stddev

  !    RETURN
  !>       \return real(dp) :: loglikelihood_trend_no_autocorr &mdash; logarithmic likelihood using given stddev
  !>       but remove optimal trend in errors
  !>       (absolute between running model with parameterset and observation)

  !    HISTORY
  !>       \authors Juliane Mai and Matthias Cuntz

  !>       \date Mar 2014

  ! Modifications:
  ! Stephan Thober  Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                            to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION loglikelihood_trend_no_autocorr(parameterset, eval, stddev_old, stddev_new, likeli_new)

    use mo_append, only : append
    use mo_linfit, only : linfit
    use mo_moment, only : stddev

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    ! standard deviation of data
    real(dp), intent(in) :: stddev_old

    ! standard deviation of errors with removed trendbetween model run using parameter set and observation
    real(dp), intent(out), optional :: stddev_new

    ! logarithmic likelihood determined with stddev_new instead of stddev
    real(dp), intent(out), optional :: likeli_new

    real(dp) :: loglikelihood_trend_no_autocorr

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), dimension(:, :), allocatable :: runoff

    ! gauges counter
    integer(i4) :: gg

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    integer(i4) :: nmeas

    real(dp), dimension(:), allocatable :: errors

    real(dp), dimension(:), allocatable :: obs, calc, out

    real(dp) :: a, b

    real(dp) :: stddev_tmp


    call eval(parameterset, runoff = runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim = 2) ! second dimension equals nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! append it to variables
      call append(obs, runoff_obs)
      call append(calc, runoff_agg)

    end do

    ! ----------------------------------------
    nmeas = size(obs, dim = 1)

    ! allocate output variables
    allocate(out(nmeas), errors(nmeas))
    errors(:) = abs(calc(:) - obs(:))

    ! remove linear trend of errors - must be model NOT obs
    out = linfit(calc, errors, a = a, b = b, model2 = .False.)
    errors(:) = errors(:) - (a + calc(:) * b)

    ! you have to take stddev_old=const because otherwise loglikelihood_trend_no_autocorr is always N
    ! in MCMC stddev_old gets updated only when a better likelihood is found.
    loglikelihood_trend_no_autocorr = sum(errors(:) * errors(:) / stddev_old**2)
    loglikelihood_trend_no_autocorr = -0.5_dp * loglikelihood_trend_no_autocorr

    call message('-loglikelihood_trend_no_autocorr = ', num2str(-loglikelihood_trend_no_autocorr))

    stddev_tmp = 1.0_dp  ! initialization
    if (present(stddev_new) .or. present(likeli_new)) then
      stddev_tmp = stddev(errors(:))
    end if
    if (present(stddev_new)) then
      stddev_new = stddev_tmp
    end if
    if (present(likeli_new)) then
      likeli_new = sum(errors(:) * errors(:) / stddev_tmp**2)
      likeli_new = -0.5_dp * likeli_new
    end if

    deallocate(runoff, runoff_agg, runoff_obs, runoff_obs_mask)
    deallocate(obs, calc, out, errors)

  END FUNCTION loglikelihood_trend_no_autocorr

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_lnnse

  !    PURPOSE
  !>       \brief Objective function of logarithmic NSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>       \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>       is calculated.
  !>       \f[ obj\_value = lnNSE \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_lnnse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date May 2013

  ! Modifications:
  ! Stephan Thober  Jan 2015 - introduced extract_runoff
  ! Stephan Thober  Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_lnnse(parameterset, eval)

    use mo_errormeasures, only : lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_lnnse

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_lnnse = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! lnNSE
      objective_lnnse = objective_lnnse + &
              lnnse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
#ifndef MPI
    ! objective function value which will be minimized
    objective_lnnse = 1.0_dp - objective_lnnse / real(nGaugesTotal, dp)

    call message('objective_lnnse = ', num2str(objective_lnnse))
    ! pause
#endif

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_lnnse

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sse

  !    PURPOSE
  !>       \brief Objective function of SSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the sum squared errors
  !>       \f[ SSE = \sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2 \f]
  !>       is calculated and the objective function is
  !>       \f[ obj\_value = SSE \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai and Matthias Cuntz

  !>       \date March 2014

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sse(parameterset, eval)

    use mo_errormeasures, only : sse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_sse

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_sse = 0.0_dp
    do gg = 1, nGaugesTotal
      !
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      !
      objective_sse = objective_sse + &
              sse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
#ifndef MPI
    ! objective_sse = objective_sse + sse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_sse = objective_sse / real(nGaugesTotal, dp)

    call message('objective_sse = ', num2str(objective_sse))
    ! pause
#endif

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_sse

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_nse

  !    PURPOSE
  !>       \brief Objective function of NSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>       \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>       is calculated and the objective function is
  !>       \f[ obj\_value = 1-NSE \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_nse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date May 2013

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_nse(parameterset, eval)

    use mo_errormeasures, only : nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_nse

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_nse = 0.0_dp
    do gg = 1, nGaugesTotal
      !
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      !
      objective_nse = objective_nse + &
              nse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
#ifndef MPI
    ! objective_nse = objective_nse + nse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_nse = 1.0_dp - objective_nse / real(nGaugesTotal, dp)

    call message('objective_nse (i.e., 1 - NSE) = ', num2str(objective_nse))
    ! pause
#endif

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_nse

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_equal_nse_lnnse

  !    PURPOSE
  !>       \brief Objective function equally weighting NSE and lnNSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>       \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>       and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>       \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>       are calculated and added up equally weighted:
  !>       \f[ obj\_value = \frac{1}{2} (NSE + lnNSE) \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_equal_nse_lnnse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date May 2013

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_equal_nse_lnnse(parameterset, eval)

    use mo_errormeasures, only : lnnse, nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_equal_nse_lnnse

    ! modelled runoff for a given parameter set
    ! dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_equal_nse_lnnse = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      !
      ! NSE
      objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
              0.5_dp * nse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      ! lnNSE
      objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
              0.5_dp * lnnse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
#ifndef MPI
    ! objective function value which will be minimized
    objective_equal_nse_lnnse = 1.0_dp - objective_equal_nse_lnnse / real(nGaugesTotal, dp)

    call message('objective_equal_nse_lnnse = ', num2str(objective_equal_nse_lnnse))
#endif

    ! clean up
    deallocate(runoff_agg, runoff_obs)
    deallocate(runoff_obs_mask)

  END FUNCTION objective_equal_nse_lnnse


  ! ------------------------------------------------------------------

  !    NAME
  !        multi_objective_nse_lnnse

  !    PURPOSE
  !>       \brief Multi-objective function with NSE and lnNSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>       \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>       and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>       \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>       are calculated and both returned.
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.
  !>       To calibrate this objective you need a multi-objective optimizer like PA-DDS.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp), dimension(2) :: multi_objective_nse_lnnse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like PA-DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Oct 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION multi_objective_nse_lnnse(parameterset, eval)

    use mo_errormeasures, only : lnnse, nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp), dimension(2) :: multi_objective_nse_lnnse

    ! modelled runoff for a given parameter set
    ! dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    ! call mhm_eval(parameterset, runoff=runoff)
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    multi_objective_nse_lnnse = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      !
      ! NSE
      multi_objective_nse_lnnse(1) = multi_objective_nse_lnnse(1) + &
              nse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      ! lnNSE
      multi_objective_nse_lnnse(2) = multi_objective_nse_lnnse(2) + &
              lnnse(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
    ! objective function value which will be minimized
    multi_objective_nse_lnnse(:) = 1.0_dp - multi_objective_nse_lnnse(:) / real(nGaugesTotal, dp)

    ! write(*,*) 'multi_objective_nse_lnnse = ',multi_objective_nse_lnnse

    ! clean up
    deallocate(runoff_agg, runoff_obs)
    deallocate(runoff_obs_mask)

  END FUNCTION multi_objective_nse_lnnse

  ! ------------------------------------------------------------------

  !    NAME
  !        multi_objective_lnnse_highflow_lnnse_lowflow

  !    PURPOSE
  !>       \brief Multi-objective function with NSE and lnNSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       A timepoint \f$t\f$ of the observed data is marked as a lowflow timepoint \f$t_{low}\f$ if
  !>       \f[ Q_{obs}(t) < min(Q_{obs}) + 0.05 * ( max(Q_{obs}) - min(Q_{obs}) )\f]
  !>       and a timepoint \f$t\f$ of the observed data is marked as a highflow timepoint \f$t_{high}\f$ if
  !>       \f[ t_{high} if Q_{obs}(i) > percentile(Q_{obs},95.)\f]
  !>       This timepoint identification is only performed for the observed data.

  !>       The first objective is the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE_{high} \f$
  !>       of discharge values at high-flow timepoints
  !>       \f[ lnNSE_{high} = 1 - \frac{\sum_{i=1}^{N_{high}} (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f] .
  !>       The second objective is the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE_{low} \f$
  !>       of discharge values at low-flow timepoints
  !>       \f[ lnNSE_{low} = 1 - \frac{\sum_{i=1}^{N_{low}} (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f] .
  !>       Both objectives are returned.
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.
  !>       To calibrate this objective you need a multi-objective optimizer like PA-DDS.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp), dimension(2) :: multi_objective_lnnse_highflow_lnnse_lowflow &mdash; objective function
  !>       value
  !>       (which will be e.g. minimized by an optimization routine like PA-DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Oct 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION multi_objective_lnnse_highflow_lnnse_lowflow(parameterset, eval)

    use mo_errormeasures, only : lnnse
    use mo_percentile, only : percentile

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp), dimension(2) :: multi_objective_lnnse_highflow_lnnse_lowflow

    ! modelled runoff for a given parameter set
    real(dp), dimension(:, :), allocatable :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! upper discharge value to determine lowflow timepoints
    real(dp) :: q_low

    ! lower discharge value to determine highflow timepoints
    real(dp) :: q_high

    ! total number of discharge values
    integer(i4) :: nrunoff

    ! timepoint counter
    integer(i4) :: tt

    ! mask to get lowflow values
    logical, dimension(:), allocatable :: lowflow_mask

    ! mask to get highflow values
    logical, dimension(:), allocatable :: highflow_mask


    ! call mhm_eval(parameterset, runoff=runoff)
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    multi_objective_lnnse_highflow_lnnse_lowflow = 0.0_dp
    do gg = 1, nGaugesTotal
      !
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      nrunoff = size(runoff_obs, dim = 1)
      !print*, 'nrunoff = ',nrunoff
      !
      ! mask for highflow timepoints
      if (allocated(highflow_mask)) deallocate(highflow_mask)
      allocate(highflow_mask(nrunoff))
      highflow_mask = .false.
      q_high = percentile(runoff_obs, 95._dp, mask = runoff_obs_mask)
      forall(tt = 1 : nrunoff) highflow_mask(tt) = (runoff_obs(tt) > q_high .and. runoff_obs_mask(tt))
      !print*, 'nhigh = ',count(highflow_mask)
      !
      ! mask for lowflow timepoints
      if (allocated(lowflow_mask)) deallocate(lowflow_mask)
      allocate(lowflow_mask(nrunoff))
      lowflow_mask = .false.
      q_low = minval(runoff_obs, mask = runoff_obs_mask) &
              + 0.05_dp * (maxval(runoff_obs, mask = runoff_obs_mask) - minval(runoff_obs, mask = runoff_obs_mask))
      forall(tt = 1 : nrunoff) lowflow_mask(tt) = (runoff_obs(tt) < q_low .and. runoff_obs_mask(tt))
      !print*, 'nlow  = ',count(lowflow_mask)
      !
      ! lnNSE highflows
      multi_objective_lnnse_highflow_lnnse_lowflow(1) = multi_objective_lnnse_highflow_lnnse_lowflow(1) + &
              lnnse(runoff_obs, runoff_agg, mask = highflow_mask)
      !
      ! lnNSE lowflows
      multi_objective_lnnse_highflow_lnnse_lowflow(2) = multi_objective_lnnse_highflow_lnnse_lowflow(2) + &
              lnnse(runoff_obs, runoff_agg, mask = lowflow_mask)
    end do
    ! objective function value which will be minimized
    multi_objective_lnnse_highflow_lnnse_lowflow(:) = 1.0_dp &
            - multi_objective_lnnse_highflow_lnnse_lowflow(:) / real(nGaugesTotal, dp)

    ! write(*,*) 'multi_objective_lnnse_highflow_lnnse_lowflow = ',multi_objective_lnnse_highflow_lnnse_lowflow

    ! clean up
    deallocate(runoff_agg, runoff_obs)
    deallocate(runoff_obs_mask)

  END FUNCTION multi_objective_lnnse_highflow_lnnse_lowflow

  ! ------------------------------------------------------------------

  !    NAME
  !        multi_objective_lnnse_highflow_lnnse_lowflow_2

  !    PURPOSE
  !>       \brief Multi-objective function with NSE and lnNSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       A timepoint \f$t\f$ of the observed data is marked as a lowflow timepoint \f$t_{low}\f$ if
  !>       \f[ Q_{obs}(t) < min(Q_{obs}) + 0.05 * ( max(Q_{obs}) - min(Q_{obs}) )\f]
  !>       and all other timepoints are marked as a highflow timepoints \f$t_{high}\f$.
  !>       This timepoint identification is only performed for the observed data.

  !>       The first objective is the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE_{high} \f$
  !>       of discharge values at high-flow timepoints
  !>       \f[ lnNSE_{high} = 1 - \frac{\sum_{i=1}^{N_{high}} (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f] .
  !>       The second objective is the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE_{low} \f$
  !>       of discharge values at low-flow timepoints
  !>       \f[ lnNSE_{low} = 1 - \frac{\sum_{i=1}^{N_{low}} (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f] .
  !>       Both objectives are returned.
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.
  !>       To calibrate this objective you need a multi-objective optimizer like PA-DDS.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp), dimension(2) :: multi_objective_lnnse_highflow_lnnse_lowflow_2 &mdash; objective function
  !>       value
  !>       (which will be e.g. minimized by an optimization routine like PA-DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Oct 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION multi_objective_lnnse_highflow_lnnse_lowflow_2(parameterset, eval)

    use mo_errormeasures, only : lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp), dimension(2) :: multi_objective_lnnse_highflow_lnnse_lowflow_2

    ! modelled runoff for a given parameter set
    real(dp), dimension(:, :), allocatable :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! upper discharge value to determine lowflow timepoints
    real(dp) :: q_low

    ! total number of discharge values
    integer(i4) :: nrunoff

    ! timepoint counter
    integer(i4) :: tt

    ! mask to get lowflow values
    logical, dimension(:), allocatable :: lowflow_mask

    ! mask to get highflow values
    logical, dimension(:), allocatable :: highflow_mask


    ! call mhm_eval(parameterset, runoff=runoff)
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    multi_objective_lnnse_highflow_lnnse_lowflow_2 = 0.0_dp
    do gg = 1, nGaugesTotal
      !
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      nrunoff = size(runoff_obs, dim = 1)
      !print*, 'nrunoff = ',nrunoff
      !
      ! mask for lowflow timepoints
      if (allocated(lowflow_mask)) deallocate(lowflow_mask)
      allocate(lowflow_mask(nrunoff))
      lowflow_mask = .false.
      q_low = minval(runoff_obs, mask = runoff_obs_mask) &
              + 0.05_dp * (maxval(runoff_obs, mask = runoff_obs_mask) - minval(runoff_obs, mask = runoff_obs_mask))
      forall(tt = 1 : nrunoff) lowflow_mask(tt) = (runoff_obs(tt) < q_low .and. runoff_obs_mask(tt))
      !print*, 'nlow  = ',count(lowflow_mask)
      !
      ! mask for highflow timepoints
      if (allocated(highflow_mask)) deallocate(highflow_mask)
      allocate(highflow_mask(nrunoff))
      highflow_mask = (.not. lowflow_mask) .and. runoff_obs_mask
      !print*, 'nhigh = ',count(highflow_mask)
      !
      ! lnNSE highflows
      multi_objective_lnnse_highflow_lnnse_lowflow_2(1) = multi_objective_lnnse_highflow_lnnse_lowflow_2(1) + &
              lnnse(runoff_obs, runoff_agg, mask = highflow_mask)
      !
      ! lnNSE lowflows
      multi_objective_lnnse_highflow_lnnse_lowflow_2(2) = multi_objective_lnnse_highflow_lnnse_lowflow_2(2) + &
              lnnse(runoff_obs, runoff_agg, mask = lowflow_mask)
    end do
    ! objective function value which will be minimized
    multi_objective_lnnse_highflow_lnnse_lowflow_2(:) = 1.0_dp &
            - multi_objective_lnnse_highflow_lnnse_lowflow_2(:) / real(nGaugesTotal, dp)

    ! write(*,*) 'multi_objective_lnnse_highflow_lnnse_lowflow_2 = ',multi_objective_lnnse_highflow_lnnse_lowflow_2

    ! clean up
    deallocate(runoff_agg, runoff_obs)
    deallocate(runoff_obs_mask)

  END FUNCTION multi_objective_lnnse_highflow_lnnse_lowflow_2

  ! ------------------------------------------------------------------

  !    NAME
  !        multi_objective_ae_fdc_lsv_nse_djf

  !    PURPOSE
  !>       \brief Multi-objective function with absolute error of Flow Duration Curves
  !>       low-segment volume and nse of DJF's discharge.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       The first objective is using the routine "FlowDurationCurves" from "mo_signatures" to determine the
  !>       low-segment volume of the FDC. The objective is the absolute difference between the observed volume
  !>       and the simulated volume.

  !>       For the second objective the discharge of the winter months December, January and February are extracted
  !>       from the time series. The objective is then the Nash-Sutcliffe efficiency NSE of the observed winter
  !>       discharge against the simulated winter discharge.

  !>       The observed data \f$ Q_{obs} \f$ are global in this module.
  !>       To calibrate this objective you need a multi-objective optimizer like PA-DDS.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp), dimension(2) :: multi_objective_ae_fdc_lsv_nse_djf &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like PA-DDS)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Feb 2016

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION multi_objective_ae_fdc_lsv_nse_djf(parameterset, eval)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_errormeasures, only : nse
    use mo_julian, only : dec2date
    use mo_mrm_global_variables, only : gauge, nMeasPerDay
    use mo_mrm_signatures, only : FlowDurationCurve

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp), dimension(2) :: multi_objective_ae_fdc_lsv_nse_djf

    ! modelled runoff for a given parameter set
    real(dp), dimension(:, :), allocatable :: runoff

    ! gauges counter
    integer(i4) :: gg

    ! total number of gauges
    integer(i4) :: nGaugesTotal

    ! domain ID of gauge
    integer(i4) :: iDomain

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! total number of discharge values
    integer(i4) :: nrunoff

    ! timepoint counter
    integer(i4) :: tt

    ! month of current time step
    integer(i4) :: month

    ! Fractional Julian day of current time step
    real(dp) :: current_time

    ! mask to get lowflow values
    logical, dimension(:), allocatable :: djf_mask

    ! quantiles for FDC
    real(dp), dimension(10) :: quantiles

    ! number of quantiles
    integer(i4) :: nquantiles

    ! FDC of simulated or observed discharge
    real(dp), dimension(size(quantiles)) :: fdc

    ! low-segment volume of FDC of simulated discharge
    real(dp) :: lsv_mod

    ! low-segment volume of FDC of observed  discharge
    real(dp) :: lsv_obs


    ! call mhm_eval(parameterset, runoff=runoff)
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)
    nquantiles = size(quantiles)
    forall(tt = 1 : nquantiles) quantiles(tt) = real(tt, dp) / real(nquantiles, dp)

    multi_objective_ae_fdc_lsv_nse_djf = 0.0_dp
    do gg = 1, nGaugesTotal
      !
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      nrunoff = size(runoff_obs, dim = 1)
      !
      ! mask DJF timepoints
      if (allocated(djf_mask)) deallocate(djf_mask)
      allocate(djf_mask(nrunoff))
      djf_mask = .false.

      iDomain = gauge%domainId(gg)
      do tt = 1, nrunoff
        current_time = evalPer(iDomain)%julStart + (tt - 1) * 1.0_dp / real(nMeasPerDay, dp)
        call dec2date(current_time, mm = month)
        if ((month == 1 .or. month == 2 .or. month == 12) .and. runoff_obs_mask(tt)) djf_mask(tt) = .True.
      end do
      !
      ! Absolute error of low-segment volume of FDC
      fdc = FlowDurationCurve(runoff_obs, quantiles, mask = runoff_obs_mask, low_segment_volume = lsv_obs)
      fdc = FlowDurationCurve(runoff_agg, quantiles, mask = runoff_obs_mask, low_segment_volume = lsv_mod)
      fdc = fdc * 1.0_dp ! only to avoid warning of unused variable
      !
      ! Absolute distance between low-segment volumes
      multi_objective_ae_fdc_lsv_nse_djf(1) = multi_objective_ae_fdc_lsv_nse_djf(1) + &
              abs(lsv_obs - lsv_mod)
      !
      ! NSE of DJF discharge
      multi_objective_ae_fdc_lsv_nse_djf(2) = multi_objective_ae_fdc_lsv_nse_djf(2) + &
              nse(runoff_obs, runoff_agg, mask = djf_mask)
    end do
    ! objective function value which will be minimized
    multi_objective_ae_fdc_lsv_nse_djf(1) = &
            multi_objective_ae_fdc_lsv_nse_djf(1) / real(nGaugesTotal, dp)
    multi_objective_ae_fdc_lsv_nse_djf(2) = 1.0_dp &
            - multi_objective_ae_fdc_lsv_nse_djf(2) / real(nGaugesTotal, dp)

    call message('multi_objective_ae_fdc_lsv_nse_djf = [', &
      num2str(multi_objective_ae_fdc_lsv_nse_djf(1)), ', ', &
      num2str(multi_objective_ae_fdc_lsv_nse_djf(2)), ']')

    ! clean up
    deallocate(runoff_agg, runoff_obs)
    deallocate(runoff_obs_mask)

  END FUNCTION multi_objective_ae_fdc_lsv_nse_djf

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_power6_nse_lnnse

  !    PURPOSE
  !>       \brief Objective function of combined NSE and lnNSE with power of 5
  !>       i.e. the p-norm with p=5.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>       \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>       and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>       \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>       {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>       are calculated and added up equally weighted:
  !>       \f[ obj\_value = \sqrt[6]{(1-NSE)^6 + (1-lnNSE)^6} \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_power6_nse_lnnse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Juliane Mai and Matthias Cuntz

  !>       \date March 2014

  ! Modifications:
  ! Stephan Thober Jan 2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_power6_nse_lnnse(parameterset, eval)

    use mo_errormeasures, only : lnnse, nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_power6_nse_lnnse

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp


    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_power6_nse_lnnse = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! NSE + lnNSE
      objective_power6_nse_lnnse = objective_power6_nse_lnnse + &
              ((1.0_dp - nse(runoff_obs, runoff_agg, mask = runoff_obs_mask))**6 + &
                      (1.0_dp - lnnse(runoff_obs, runoff_agg, mask = runoff_obs_mask))**6)**onesixth
    end do
#ifndef MPI
    ! objective function value which will be minimized
    objective_power6_nse_lnnse = objective_power6_nse_lnnse / real(nGaugesTotal, dp)

    call message('objective_power6_nse_lnnse = ', num2str(objective_power6_nse_lnnse))
    ! pause
#endif

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_power6_nse_lnnse

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_kge

  !    PURPOSE
  !>       \brief Objective function of KGE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency coefficient \f$ KGE \f$
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of similated mean to observed mean
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function is
  !>       \f[ obj\_value = 1.0 - KGE \f]
  !>       \f$(1-KGE)\f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$(1-KGE)\f$ is 0 for the optimal KGE of 1.0.

  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date August 2014

  ! Modifications:
  ! Stephan Thober Jan  2015 - introduced extract_runoff
  ! Stephan Thober Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2) to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge(parameterset, eval)

    use mo_errormeasures, only : kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_kge

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    !
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_kge = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! KGE
      objective_kge = objective_kge + &
              kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)
    end do
#ifndef MPI
    ! objective_kge = objective_kge + kge(gauge%Q, runoff_model_agg, runoff_model_agg_mask)
    objective_kge = 1.0_dp - objective_kge / real(nGaugesTotal, dp)

    call message('objective_kge (i.e., 1 - KGE) = ', num2str(objective_kge))
#endif
    ! pause

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_kge

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_multiple_gauges_kge_power6

  !    PURPOSE
  !>       \brief combined objective function based on KGE raised to the power 6

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency coefficient \f$ KGE \f$ for a given gauging station
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of similated mean to observed mean
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given gauging station (\f$ i \f$) is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall \f$ OF \f$ is estimated based on the power-6 norm to
  !>       combine the \f$ \phi_{i} \f$ from all gauging stations (\f$ N \f$).
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }  \f].

  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_multiple_gauges_kge_power6 &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date March 2015

  ! Modifications:
  ! Stephan Thober  Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                            to not interfere with mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_multiple_gauges_kge_power6(parameterset, eval)

    use mo_errormeasures, only : kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_multiple_gauges_kge_power6

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    !
    call eval(parameterset, runoff = runoff)
    nGaugesTotal = size(runoff, dim = 2)

    objective_multiple_gauges_kge_power6 = 0.0_dp
    do gg = 1, nGaugesTotal
      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)
      ! KGE
      objective_multiple_gauges_kge_power6 = objective_multiple_gauges_kge_power6 + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6
    end do
#ifndef MPI
    objective_multiple_gauges_kge_power6 = objective_multiple_gauges_kge_power6**onesixth
    call message('objective_multiple_gauges_kge_power6 = ', num2str(objective_multiple_gauges_kge_power6))
#endif

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

  END FUNCTION objective_multiple_gauges_kge_power6

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_weighted_nse

  !    PURPOSE
  !>       \brief Objective function of weighted NSE.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the weighted Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>       \f[ wNSE = 1 - \frac{\sum_{i=1}^N Q_{obs}(i) * (Q_{obs}(i) - Q_{model}(i))^2}
  !>       {\sum_{i=1}^N Q_{obs}(i) * (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>       is calculated and the objective function is
  !>       \f[ obj\_value = 1- wNSE \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_weighted_nse &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Stephan Thober, Bjoern Guse

  !>       \date May 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_weighted_nse(parameterset, eval)

    use mo_errormeasures, only : wnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_weighted_nse

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:,:) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask


    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_weighted_nse = 0.0_dp
    do gg=1, nGaugesTotal
       !
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       !
       objective_weighted_nse = objective_weighted_nse + &
            wnse( runoff_obs, runoff_agg, mask=runoff_obs_mask)
    end do
#ifndef MPI
    ! objective_nse = objective_nse + nse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_weighted_nse = 1.0_dp - objective_weighted_nse / real(nGaugesTotal,dp)

    call message('objective_weighted_nse (i.e., 1 - wNSE) = ', num2str(objective_weighted_nse))
    ! pause
#endif

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_weighted_nse

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sse_boxcox

  !    PURPOSE
  !>       \brief Objective function of sum of squared errors of transformed streamflow.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.
  !>       Therefore, the sum of squared error \f$ tSSE \f$
  !>       \f[ tSSE = \sum_{i=1}^N (z(Q_{obs}(i), \lambda) - z(Q_{model}(i), \lambda))^2 \f]
  !>       is calculated where \f$ z \f$ is the transform and given by
  !>       \f[ z(x, \lambda) = \frac{x^\lambda -1}{\lambda} \f] for \f$ \lambda \f$ unequal to zero and
  !>       \f[ z(x, \lambda) = log x \f] for \f$ \lambda \f$ equal to zero.
  !>       The objective function is
  !>       \f[ obj\_value = tSSE \f]
  !>       The observed data \f$ Q_{obs} \f$ are global in this module.
  !>
  !>       The boxcox transformation uses a parameter of 0.2, suggested by
  !>       Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278.
  !>       "Evaluating post-processing approaches for monthly and seasonal streamflow forecasts."
  !>       https://www.hydrol-earth-syst-sci.net/22/6257/2018/


  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sse_boxcox &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Stephan Thober, Dmitri Kavetski

  !>       \date Aug 2019

  FUNCTION objective_sse_boxcox(parameterset, eval)

    use mo_errormeasures, only: SSE
    use mo_boxcox, only: boxcox

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), pointer :: eval

    real(dp) :: objective_sse_boxcox

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:,:) :: runoff

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for aggregated measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! boxcox parameter
    real(dp) :: lambda

    ! set box cox transformation to 0.2
    ! suggested by:
    ! Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278.
    ! "Evaluating post-processing approaches for monthly and seasonal streamflow forecasts."
    ! https://www.hydrol-earth-syst-sci.net/22/6257/2018/
    lambda = 0.2_dp

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_sse_boxcox = 0.0_dp
    do gg=1, nGaugesTotal
       !
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       !
       objective_sse_boxcox = objective_sse_boxcox + &
            SSE( boxcox(runoff_obs, lambda, mask=runoff_obs_mask), boxcox(runoff_agg, lambda), mask=runoff_obs_mask)
    end do
#ifndef MPI
    ! objective_nse = objective_nse + nse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_sse_boxcox = objective_sse_boxcox / real(nGaugesTotal,dp)

    call message('objective_sse_boxcox = ', num2str(objective_sse_boxcox))
    ! pause
#endif

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_sse_boxcox


  ! ------------------------------------------------------------------

  !    NAME
  !        extract_runoff

  !    PURPOSE
  !>       \brief extracts runoff data from global variables

  !>       \details extracts simulated and measured runoff from global variables,
  !>       such that they overlay exactly. For measured runoff, only the runoff
  !>       during the evaluation period are cut, not succeeding nodata values.
  !>       For simulated runoff, warming days as well as succeeding nodata values
  !>       are neglected and the simulated runoff is aggregated to the resolution
  !>       of the observed runoff.
  !>       see use in this module above

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: gaugeId"              current gauge Id
  !>       \param[in] "real(dp), dimension(:, :) :: runoff" simulated runoff

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: runoff_agg"     aggregated simulated
  !>       \param[out] "real(dp), dimension(:) :: runoff_obs"     extracted measured
  !>       \param[out] "logical, dimension(:) :: runoff_obs_mask" mask of no data values

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Jan 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine extract_runoff(gaugeId, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

    use mo_common_mhm_mrm_variables, only : evalPer, nTstepDay, warmingDays
    use mo_mrm_global_variables, only : gauge, nMeasPerDay
    use mo_utils, only : ge

    implicit none

    ! current gauge Id
    integer(i4), intent(in) :: gaugeId

    ! simulated runoff
    real(dp), dimension(:, :), intent(in) :: runoff

    ! aggregated simulated
    real(dp), dimension(:), allocatable, intent(out) :: runoff_agg

    ! extracted measured
    real(dp), dimension(:), allocatable, intent(out) :: runoff_obs

    ! mask of no data values
    logical, dimension(:), allocatable, intent(out) :: runoff_obs_mask

    ! domain id
    integer(i4) :: iDomain

    ! timestep counter
    integer(i4) :: tt

    ! length of extracted time series
    integer(i4) :: length

    ! between simulated and measured time scale
    integer(i4) :: factor

    ! simulated Timesteps per Day
    integer(i4) :: TPD_sim

    ! observed Timesteps per Day
    integer(i4) :: TPD_obs

    real(dp), dimension(:), allocatable :: dummy


    ! copy time resolution to local variables
    TPD_sim = nTstepDay
    TPD_obs = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if (modulo(TPD_sim, TPD_obs) .eq. 0) then
      factor = TPD_sim / TPD_obs
    else
      call error_message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
    end if

    ! extract domain Id from gauge Id
    iDomain = gauge%domainId(gaugeId)

    ! get length of evaluation period times TPD_obs
    length = (evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1) * TPD_obs

    ! extract measurements
    if (allocated(runoff_obs)) deallocate(runoff_obs)
    allocate(runoff_obs(length))
    runoff_obs = gauge%Q(1 : length, gaugeId)

    ! create mask of observed runoff
    if (allocated(runoff_obs_mask)) deallocate(runoff_obs_mask)
    allocate(runoff_obs_mask(length))
    runoff_obs_mask = .false.
    forall(tt = 1 : length) runoff_obs_mask(tt) = ge(runoff_obs(tt), 0.0_dp)

    ! extract and aggregate simulated runoff
    if (allocated(runoff_agg)) deallocate(runoff_agg)
    allocate(runoff_agg(length))
    ! remove warming days
    length = (evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1) * TPD_sim
    allocate(dummy(length))
    dummy = runoff(warmingDays(iDomain) * TPD_sim + 1 : warmingDays(iDomain) * TPD_sim + length, gaugeId)
    ! aggregate runoff
    length = (evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1) * TPD_obs
    forall(tt = 1 : length) runoff_agg(tt) = sum(dummy((tt - 1) * factor + 1 : tt * factor)) / &
            real(factor, dp)
    ! clean up
    deallocate(dummy)

  end subroutine extract_runoff


END MODULE mo_mrm_objective_function_runoff
