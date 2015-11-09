!> \file mo_objective_function.f90

!> \brief Objective Functions for Optimization of mHM/mRM against runoff.\n

!> \details This module provides a wrapper for several objective functions used to optimize mRM/mHM against runoff.\n
!>          If the objective contains besides runoff another variable like TWS move it to mHM/mo_objective_function.f90.
!>          If it is only regarding runoff implement it here.\n
!>
!>          All the objective functions are supposed to be minimized! \n
!>               (1)  SO: Q:        1.0 - NSE  \n
!>               (2)  SO: Q:        1.0 - lnNSE  \n
!>               (3)  SO: Q:        1.0 - 0.5*(NSE+lnNSE)  \n
!>               (4)  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors and
!>                                        then lag(1)-autocorrelation removed  \n
!>               (5)  SO: Q:        ((1-NSE)**6+(1-lnNSE)**6)**(1/6)  \n
!>               (6)  SO: Q:        SSE  \n
!>               (7)  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors  \n
!>               (8)  SO: Q:       -1.0 * loglikelihood with trend removed from the relative errors and
!>                                        then lag(1)-autocorrelation removed \n
!>               (9)  SO: Q:        1.0 - KGE (Kling-Gupta efficiency measure)  \n
!>               (14) SO: Q:        sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) > combination of KGE of every gauging
!>                                        station based on a power-6 norm \n
!>               (16) (reserved) not used yet
!>                    MO: Q:        1st objective: (1) = 1.0 - NSE  
!>                                  2nd objective: (2) = 1.0 - lnNSE

!> \authors Juliane Mai
!> \date Dec 2012
!  Modified, Oct 2015, Stephan Thober - adapted for mRM
!            Nov 2015, Juliane Mai    - introducing multi- and single-objective
!                                     - first multi-objective function (16), but not used yet

MODULE mo_mrm_objective_function_runoff

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 - removed all none runoff objectives,
  !                                     these can be found mo_objective_functions_sm

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: loglikelihood           ! loglikelihood with errormodel including linear trend and lag(1)-correlation
  PUBLIC :: loglikelihood_stddev    ! loglikelihood where error is computed from difference of obs vs model
  PUBLIC :: single_objective_runoff ! single-objective function wrapper
  PUBLIC :: multi_objective_runoff  ! multi-objective function wrapper
  PUBLIC :: extract_runoff          ! extract runoff period specified in mhm.nml from available runoff time series 

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_runoff

  !>        \brief Wrapper for objective functions optimizing agains runoff.

  !>        \details The functions selects the objective function case defined in a namelist, 
  !>        i.e. the global variable \e opti\_function.\n
  !>        It return the objective function value for a specific parameter set.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_runoff(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified,
  !               Oct 2015, Stephan Thober - only runoff objective functions

  FUNCTION single_objective_runoff(parameterset)

    USE mo_common_variables, ONLY: opti_function

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: parameterset
    REAL(dp)                            :: single_objective_runoff

    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
    case (1)
       ! 1.0-nse
       single_objective_runoff = objective_nse(parameterset)
    case (2)
       ! 1.0-lnnse
       single_objective_runoff = objective_lnnse(parameterset)
    case (3)
       ! 1.0-0.5*(nse+lnnse)
       single_objective_runoff = objective_equal_nse_lnnse(parameterset)
    case (4)
       ! -loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
       single_objective_runoff = - loglikelihood_stddev(parameterset, 1.0_dp)
    case (5)
       ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
       single_objective_runoff = objective_power6_nse_lnnse(parameterset)
    case (6)
       ! SSE
       single_objective_runoff = objective_sse(parameterset)
    case (7)
       ! -loglikelihood with trend removed from absolute errors
       single_objective_runoff = -loglikelihood_trend_no_autocorr(parameterset, 1.0_dp)
    case (8)
       ! -loglikelihood of approach 2 of Evin et al. (2013),
       !  i.e. linear error model with lag(1)-autocorrelation on relative errors
       single_objective_runoff = -loglikelihood_evin2013_2(parameterset)
    case (9)
       ! KGE
       single_objective_runoff = objective_kge(parameterset)
    case (14)
       ! combination of KGE of every gauging station based on a power-6 norm \n
       ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) 
       single_objective_runoff = objective_multiple_gauges_kge_power6(parameterset)
    case default
       stop "Error objective: This opti_function is either not implemented yet or is not a single-objective one."
    end select

  END FUNCTION single_objective_runoff

  ! ------------------------------------------------------------------ 

  !      NAME 
  !          multi_objective_runoff 

  !>        \brief Wrapper for multi-objective functions where at least one is regarding runoff. 

  !>        \details The functions selects the objective function case defined in a namelist,  
  !>        i.e. the global variable \e opti\_function.\n 
  !>        It return the multiple objective function values for a specific parameter set. 

  !     INTENT(IN) 
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with 

  !     INTENT(INOUT) 
  !         None 

  !     INTENT(OUT) 
  !>        \param[out] "real(dp) :: multi_objectives(:)"   1D-array with multiple objective function values 

  !     INTENT(IN), OPTIONAL 
  !         None 

  !     INTENT(INOUT), OPTIONAL 
  !         None 

  !     INTENT(OUT), OPTIONAL 
  !         None 

  !     RETURN 
  !        None 

  !     RESTRICTIONS 
  !>       \note Input values must be floating points. 

  !     EXAMPLE 
  !         para = (/ 1., 2, 3., -999., 5., 6. /) 
  !         obj_value = objective(para) 

  !     LITERATURE 

  !     HISTORY 
  !>        \author Juliane Mai 
  !>        \date Oct 2015 
  !         Modified,  

  SUBROUTINE multi_objective_runoff(parameterset, multi_objectives) 

    USE mo_common_variables, ONLY: opti_function 

    IMPLICIT NONE 

    REAL(dp), DIMENSION(:),              INTENT(IN)  :: parameterset 
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: multi_objectives 

    select case (opti_function) 
    case (16) 
       ! 1st objective: 1.0-nse 
       ! 2nd objective: 1.0-lnnse 
       multi_objectives = multi_objective_nse_lnnse(parameterset) 
    case default 
       stop "Error objective: Either this opti_function is not implemented yet or it is not a multi-objective one." 
    end select

  END SUBROUTINE multi_objective_runoff

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood

  !>        \brief Wrapper for loglikelihood functions.

  !>        \details This wrapper picks the loglikelihood function selected from the namelist parameter
  !>                 \e opti\_function. \n
  !>                 It returns the return value of the selected loglikelihood function.
  !>
  !>                 This routine assumes that the wrapped function is a real likelihood
  !>                 where the errors are known or modelled.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: loglikelihood &mdash; loglikelihood function value

  !     RESTRICTIONS
  !>       \note The wrapped functions must return real loglikelihoods, i.e. 
  !>        errors are either known from observations or modelled.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = loglikelihood(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, 
  FUNCTION loglikelihood(parameterset)

    USE mo_common_variables, ONLY: opti_function, opti_method

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset
    REAL(dp)                           :: loglikelihood

    select case (opti_function)
    case (8)
       if (opti_method == 0) then
          ! Approach 2 of Evin et al. (2013), i.e. linear error model with lag(1)-autocorrelation on relative errors
          loglikelihood = loglikelihood_evin2013_2(parameterset, regularize=.true.)
       else
          loglikelihood = loglikelihood_evin2013_2(parameterset)
       endif
    case default
       stop "Error loglikelihood: chosen opti_function is no loglikelihood."
    end select

  END FUNCTION loglikelihood

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood_stddev

  !>        \brief Logarithmic likelihood function with removed linear trend and Lag(1)-autocorrelation.

  !>        \details The logarithmis likelihood function is used when mHM runs in MCMC mode.\n
  !>                 It can also be used for optimization when selecting the likelihood in the
  !>                 namelist as \e opti\_function.\n

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with
  !>        \param[in] "real(dp) :: stddev"                 standard deviation of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional :: stddev_new"   standard deviation of errors with removed trend and correlation 
  !>                                                         between model run using parameter set and observation
  !>        \param[out] "real(dp), optional :: likeli_new"   logarithmic likelihood determined with stddev_new instead of stddev

  !     RETURN
  !>       \return     real(dp) :: loglikelihood_stddev &mdash; logarithmic likelihood using given stddev 
  !>                                                     but remove optimal trend and lag(1)-autocorrelation in errors 
  !>                                                     (absolute between running model with parameterset and observation) 

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !     EXAMPLE
  !         para = (/ 1._dp, 2._dp, 3._dp, -999._dp, 5._dp, 6._dp /)
  !         stddev = 0.5_dp
  !         log_likeli = loglikelihood_stddev(para, stddev, stddev_new=stddev_new, likeli_new=likeli_new)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION loglikelihood_stddev(parameterset, stddev, stddev_new, likeli_new)
    use mo_moment,           only: mean, correlation
    use mo_linfit,           only: linfit
    use mo_append,           only: append

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev           ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood_stddev

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for aggregated measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b, c
    real(dp)                              :: stddev_tmp

    call eval(parameterset, runoff=runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim=2) ! second dimension equals nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! append it to variables
       call append( obs,  runoff_obs )
       call append( calc, runoff_agg )

    end do
    ! ----------------------------------------

    nmeas     = size(obs, dim = 1)

    allocate(out(nmeas), errors(nmeas))
    errors(:) = abs( calc(:) - obs(:) )

    ! remove linear trend of errors - must be model NOT obs
    ! out = linfit(obs, errors, a=a, b=b, model2=.False.)
    ! errors(:) = errors(:) - (a + obs(:)*b)
    out = linfit(calc, errors, a=a, b=b, model2=.False.)
    errors(:) = errors(:) - (a + calc(:)*b)

    ! remove lag(1)-autocorrelation of errors
    c = correlation(errors(2:nmeas),errors(1:nmeas-1))
    errors(1:nmeas-1) = errors(2:nmeas) - c*errors(1:nmeas-1)
    errors(nmeas)     = 0.0_dp

    ! you have to take stddev=const because otherwise loglikelihood is always N
    ! in MCMC stddev gets updated only when a better likelihood is found.
    loglikelihood_stddev = sum( errors(:) * errors(:) / stddev**2 )
    loglikelihood_stddev = -0.5_dp * loglikelihood_stddev

    write(*,*) '-loglikelihood_stddev = ', -loglikelihood_stddev

    stddev_tmp = sqrt(sum( (errors(:) - mean(errors)) * (errors(:) - mean(errors))) / real(nmeas-1,dp))
    if (present(stddev_new)) then
       stddev_new = stddev_tmp
    end if
    if (present(likeli_new)) then
       likeli_new = sum( errors(:) * errors(:) / stddev_tmp**2 )
       likeli_new = -0.5_dp * likeli_new
    end if

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs )
    deallocate(obs, calc, out, errors)

  END FUNCTION loglikelihood_stddev

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood_evin2013_2

  !>        \brief Logarithmised likelihood with linear error model and lag(1)-autocorrelation
  !>               of the relative errors.

  !>        \details This loglikelihood uses a linear error model and a lag(1)-autocorrelation
  !>                 on the relative errors. This is approach 2 of the paper Evin et al. (WRR, 2013).
  !>
  !>                 This is opti_function = 8.
  !>
  !>                 mHM then adds two extra (local) parameters for the error model in mhm_driver,
  !>                 which get optimised together with the other, global parameters.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with
  !>        \param[in] "real(dp) :: stddev_old"                 standard deviation of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional :: stddev_new"   standard deviation of errors with removed trend and correlation 
  !>                                                         between model run using parameter set and observation
  !>        \param[out] "real(dp), optional :: likeli_new"   logarithmic likelihood determined with stddev_new instead of stddev

  !     RETURN
  !>       \return     real(dp) :: loglikelihood_evin2013_2 &mdash; logarithmic likelihood using given stddev 
  !>                                                     but remove optimal trend and lag(1)-autocorrelation in errors 
  !>                                                     (absolute between running model with parameterset and observation) 

  !     RESTRICTIONS
  !>       \note Does not work with MCMC yet.

  !     EXAMPLE
  !         para = (/ 1._dp, 2._dp, 3._dp, -999._dp, 5._dp, 6._dp /)
  !         log_likeli = loglikelihood_evin2013_2(para, 1.0_dp)

  !     LITERATURE
  !         Evin et al., WRR 49, 4518-4524, 2013

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date Mar 2014
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION loglikelihood_evin2013_2(parameterset, regularize)

    use mo_constants,        only: pi_dp
    use mo_moment,           only: correlation
    use mo_common_variables, only: global_parameters ! for parameter ranges --> col1=min, col2=max
    use mo_utils,            only: eq
    use mo_append,           only: append

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    logical,  optional,     intent(in)            :: regularize
    real(dp)                                      :: loglikelihood_evin2013_2

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    real(dp)                              :: penalty                  ! penalty term due to a parmeter set out of bound
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors, sigma, eta, y
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b, c, vary, vary1, ln2pi, tmp
    integer(i4)                           :: npara
    logical                               :: iregularize

    iregularize = .false.
    if (present(regularize)) iregularize = regularize

    npara = size(parameterset)
    call eval(parameterset(1:npara-2), runoff=runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim=2) ! second dimension equals nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! append it to variables
       call append( obs,  runoff_obs )
       call append( calc, runoff_agg )

    end do

    ! ----------------------------------------

    nmeas     = size(obs, dim = 1 )

    allocate( out(nmeas), errors(nmeas), sigma(nmeas), eta(nmeas), y(nmeas))
    ! residual errors
    errors(:) = calc(:) - obs(:)

    ! linear error model
    a = parameterset(npara-1)
    b = parameterset(npara)
    sigma(:) = a + b * calc(:)
    ! standardized residual errors (SRE)
    eta(:)   = errors(:) / sigma(:)

    ! remove lag(1)-autocorrelation of SRE
    c = correlation(eta(2:nmeas),eta(1:nmeas-1))
    y(1) = 0.0_dp ! only for completeness
    y(2:nmeas) = eta(2:nmeas) - c*eta(1:nmeas-1)

    ! likelihood of residual errors
    ln2pi = log(sqrt(2.0_dp*pi_dp))
    vary  = 1.0_dp - c*c
    vary1 = 1.0_dp / vary
    loglikelihood_evin2013_2 = -ln2pi - 0.5_dp*eta(1)*eta(1) - log(sigma(1)) & ! li(eta(1))/sigma(1)
         - real(nmeas-1,dp)*log(sqrt(2.0_dp*pi_dp*vary)) &
         - sum(0.5_dp*y(2:nmeas)*y(2:nmeas)*vary1) - sum(log(sigma(2:nmeas)))

    if (iregularize) then
       ! Regularistion term as deviation from initial parameter value
       penalty = parameter_regularization(        &
            parameterset(1:npara-2),          &        ! current parameter set 
            global_parameters(1:npara-2,3),   &        ! prior/initial parameter set
            global_parameters(1:npara-2,1:2), &        ! bounds
            eq(global_parameters(1:npara-2,4),1.0_dp)) ! used/unused

       tmp = loglikelihood_evin2013_2 + penalty
       write(*,*) '-loglikelihood_evin2013_2, + penalty, chi^2: ', -loglikelihood_evin2013_2, -tmp, -tmp/real(nmeas,dp)
       loglikelihood_evin2013_2 = tmp
    else
       write(*,*) '-loglikelihood_evin2013_2, chi^2: ', -loglikelihood_evin2013_2, -loglikelihood_evin2013_2/real(nmeas,dp)
    endif

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs )
    deallocate(obs, calc, out, errors, sigma, eta, y)

  END FUNCTION loglikelihood_evin2013_2

  ! Regularisation function sum(((para-ini)/sigma)**2)
  FUNCTION parameter_regularization(paraset, prior, bounds, mask)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:),               intent(in) :: paraset
    real(dp), dimension(size(paraset)),   intent(in) :: prior
    real(dp), dimension(size(paraset),2), intent(in) :: bounds                      ! (min, max)
    logical,  dimension(size(paraset)),   intent(in) :: mask
    real(dp)                                         :: parameter_regularization

    ! local variables
    integer(i4) :: ipara
    integer(i4) :: npara
    real(dp), parameter :: onetwelveth = 1._dp/12._dp
    real(dp), dimension(size(paraset)) :: sigma

    npara = size(paraset,1)

    sigma = sqrt(onetwelveth*(bounds(:,2)-bounds(:,1))**2) ! standard deviation of uniform distribution 
    parameter_regularization = -sum(log(sqrt(2.0_dp*pi_dp)*sigma), mask=mask)

    do ipara=1,npara
       if (mask(ipara)) then
          ! if ((paraset(ipara) .lt. bounds(ipara,1)) .or. (paraset(ipara) .gt. bounds(ipara,2))) then
          !    ! outside bounds
          parameter_regularization = parameter_regularization - &
               0.5_dp*((paraset(ipara)-prior(ipara))/sigma(ipara))**2
          ! else
          !    ! in bound
          !    parameter_regularization = 0.0_dp
          ! end if
       endif
    end do

  END FUNCTION parameter_regularization

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood_trend_no_autocorr

  !>        \brief Logarithmic likelihood function with linear trend removed.

  !>        \details The logarithmis likelihood function is used when mHM runs in MCMC mode.\n
  !>                 It can also be used for optimization when selecting the likelihood in the
  !>                 namelist as \e opti\_function.\n

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with
  !>        \param[in] "real(dp) :: stddev"                 standard deviation of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional :: stddev_new"   standard deviation of errors with removed trend 
  !>                                                         between model run using parameter set and observation
  !>        \param[out] "real(dp), optional :: likeli_new"   logarithmic likelihood determined with stddev_new instead of stddev

  !     RETURN
  !>       \return     real(dp) :: loglikelihood_trend_no_autocorr &mdash; logarithmic likelihood using given stddev 
  !>                                                     but remove optimal trend in errors 
  !>                                                     (absolute between running model with parameterset and observation) 

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !     EXAMPLE
  !         para = (/ 1._dp, 2._dp, 3._dp, -999._dp, 5._dp, 6._dp /)
  !         stddev = 0.5_dp
  !         log_likeli = loglikelihood_trend_no_autocorr(para, stddev, stddev_new=stddev_new, likeli_new=likeli_new)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date Mar 2014
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION loglikelihood_trend_no_autocorr(parameterset, stddev_old, stddev_new, likeli_new)
    use mo_moment,           only: stddev
    use mo_linfit,           only: linfit
    use mo_append,           only: append

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev_old       ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood_trend_no_autocorr

    ! local
    real(dp), dimension(:,:), allocatable :: runoff          ! modelled runoff for a given parameter set
    !                                                        ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg              ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg      ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs      ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask ! mask for aggregated measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b
    real(dp)                              :: stddev_tmp

    call eval(parameterset, runoff=runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, size(runoff, dim=2) ! second dimension equals nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! append it to variables
       call append( obs,  runoff_obs )
       call append( calc, runoff_agg )

    end do

    ! ----------------------------------------
    nmeas     = size(obs, dim = 1)

    ! allocate output variables
    allocate(out(nmeas), errors(nmeas))
    errors(:) = abs( calc(:) - obs(:) )

    ! remove linear trend of errors - must be model NOT obs
    out = linfit(calc, errors, a=a, b=b, model2=.False.)
    errors(:) = errors(:) - (a + calc(:)*b)

    ! you have to take stddev_old=const because otherwise loglikelihood_trend_no_autocorr is always N
    ! in MCMC stddev_old gets updated only when a better likelihood is found.
    loglikelihood_trend_no_autocorr = sum( errors(:) * errors(:) / stddev_old**2 )
    loglikelihood_trend_no_autocorr = -0.5_dp * loglikelihood_trend_no_autocorr

    write(*,*) '-loglikelihood_trend_no_autocorr = ', -loglikelihood_trend_no_autocorr

    stddev_tmp = 1.0_dp  ! initialization
    if (present(stddev_new) .or. present(likeli_new)) then
       stddev_tmp = stddev(errors(:))
    end if
    if (present(stddev_new)) then
       stddev_new = stddev_tmp
    end if
    if (present(likeli_new)) then
       likeli_new = sum( errors(:) * errors(:) / stddev_tmp**2 )
       likeli_new = -0.5_dp * likeli_new
    end if

    deallocate(runoff, runoff_agg, runoff_obs, runoff_obs_mask)
    deallocate(obs, calc, out, errors)

  END FUNCTION loglikelihood_trend_no_autocorr

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_lnnse

  !>        \brief Objective function of logarithmic NSE.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>        \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>                             {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>        is calculated.
  !>        \f[ obj\_value = lnNSE \f]
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_lnnse &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ 1-lnnse \f$ will be returned such that it can be minimized.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_lnnse(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date May 2013
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_lnnse(parameterset)

    use mo_errormeasures,    only: lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_lnnse = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! lnNSE
       objective_lnnse = objective_lnnse + &
            lnnse( runoff_obs, runoff_agg, mask=runoff_obs_mask)
    end do
    ! objective function value which will be minimized
    objective_lnnse = 1.0_dp - objective_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_lnnse = ',objective_lnnse
    ! pause

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_lnnse

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_sse

  !>        \brief Objective function of SSE.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the sum squared errors
  !>        \f[ SSE = \sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2 \f]
  !>        is calculated and the objective function is
  !>        \f[ obj\_value = SSE \f]
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_sse &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_sse(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date March 2014
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_sse(parameterset)

    use mo_errormeasures,    only: sse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_sse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_sse = 0.0_dp
    do gg=1, nGaugesTotal
       !
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       !
       objective_sse = objective_sse + &
            sse( runoff_obs, runoff_agg, mask=runoff_obs_mask)
    end do
    ! objective_sse = objective_sse + sse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_sse = objective_sse / real(nGaugesTotal,dp)

    write(*,*) 'objective_sse = ', objective_sse
    ! pause

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_sse

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_nse

  !>        \brief Objective function of NSE.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>        \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>                           {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>        is calculated and the objective function is
  !>        \f[ obj\_value = 1-NSE \f]
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_nse &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ 1-NSE \f$ will be returned such that it can be minimized.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_nse(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date May 2013
  !         Modified, Stephan Thober, Jan 2015 - introduced extract runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_nse(parameterset)

    use mo_errormeasures,    only: nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_nse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for aggregated measured runoff

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_nse = 0.0_dp
    do gg=1, nGaugesTotal
       !
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       !
       objective_nse = objective_nse + &
            nse( runoff_obs, runoff_agg, mask=runoff_obs_mask)
    end do
    ! objective_nse = objective_nse + nse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_nse = 1.0_dp - objective_nse / real(nGaugesTotal,dp)

    write(*,*) 'objective_nse = ',objective_nse
    ! pause

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_nse

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_equal_nse_lnnse

  !>        \brief Objective function equally weighting NSE and lnNSE.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>        \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>                           {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>        and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>        \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>                             {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>        are calculated and added up equally weighted:
  !>        \f[ obj\_value = \frac{1}{2} (NSE + lnNSE) \f]
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_equal_nse_lnnse &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ 1-0.5*(nse + lnnse) \f$ will be returned such that it can be minimized.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_equal_nse_lnnse(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date May 2013
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_equal_nse_lnnse(parameterset)

    use mo_errormeasures,    only: nse, lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_equal_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff             ! modelled runoff for a given parameter set
    !                                                           ! dim2=nGauges
    integer(i4)                           :: gg                 ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_obs         ! measured runoff
    real(dp), dimension(:),   allocatable :: runoff_agg         ! aggregated simulated runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask    ! mask for aggregated measured runoff

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_equal_nse_lnnse = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       !
       ! NSE
       objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
            nse(   runoff_obs, runoff_agg, mask=runoff_obs_mask )
       ! lnNSE
       objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
            lnnse( runoff_obs, runoff_agg, mask=runoff_obs_mask )
    end do
    ! objective function value which will be minimized
    objective_equal_nse_lnnse = 1.0_dp - 0.5_dp * objective_equal_nse_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_equal_nse_lnnse = ',objective_equal_nse_lnnse

    ! clean up
    deallocate( runoff_agg, runoff_obs )
    deallocate( runoff_obs_mask )

  END FUNCTION objective_equal_nse_lnnse


  ! ------------------------------------------------------------------ 

  !      NAME 
  !          multi_objective_nse_lnnse 

  !>        \brief Multi-objective function with NSE and lnNSE. 

  !>        \details The objective function only depends on a parameter vector.  
  !>        The model will be called with that parameter vector and  
  !>        the model output is subsequently compared to observed data. 
  !>        Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$ 
  !>        \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2} 
  !>                           {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f] 
  !>        and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$ 
  !>        \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2} 
  !>                             {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f] 
  !>        are calculated and both returned. 
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 
  !>        To calibrate this objective you need a multi-objective optimizer like PA-DDS. 

  !     INTENT(IN) 
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with 

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
  !>       \return     real(dp), dimension(2) :: multi_objective_nse_lnnse &mdash; objective function value  
  !>                                             (which will be e.g. minimized by an optimization routine like PA-DDS) 

  !     RESTRICTIONS 
  !>       \note Input values must be floating points. \n 
  !>             Actually, \f$ 1-nse \f$ and \f$1-lnnse\f$ will be returned such that it can be minimized. 

  !     EXAMPLE 
  !         para = (/ 1., 2, 3., -999., 5., 6. /) 
  !         obj_value = objective_equal_nse_lnnse(para) 

  !     LITERATURE 

  !     HISTORY 
  !>        \author Juliane Mai 
  !>        \date Oct 2015 
  !         Modified,  

  FUNCTION multi_objective_nse_lnnse(parameterset) 

    ! use mo_mhm_eval,         only: mhm_eval 
    use mo_errormeasures,    only: nse, lnnse 

    implicit none 

    real(dp), dimension(:), intent(in) :: parameterset 
    real(dp), dimension(2)             :: multi_objective_nse_lnnse 

    ! local 
    real(dp), allocatable, dimension(:,:) :: runoff             ! modelled runoff for a given parameter set 
    !                                                           ! dim2=nGauges 
    integer(i4)                           :: gg                 ! gauges counter 
    integer(i4)                           :: nGaugesTotal 
    real(dp), dimension(:),   allocatable :: runoff_obs         ! measured runoff 
    real(dp), dimension(:),   allocatable :: runoff_agg         ! aggregated simulated runoff 
    logical,  dimension(:),   allocatable :: runoff_obs_mask    ! mask for aggregated measured runoff 

    ! call mhm_eval(parameterset, runoff=runoff) 
    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2) 

    multi_objective_nse_lnnse = 0.0_dp 
    do gg=1, nGaugesTotal 
       ! extract runoff 
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask ) 
       ! 
       ! NSE 
       multi_objective_nse_lnnse(1) = multi_objective_nse_lnnse(1) + & 
            nse(   runoff_obs, runoff_agg, mask=runoff_obs_mask ) 
       ! lnNSE 
       multi_objective_nse_lnnse(2) = multi_objective_nse_lnnse(2) + & 
            lnnse( runoff_obs, runoff_agg, mask=runoff_obs_mask ) 
    end do
    ! objective function value which will be minimized 
    multi_objective_nse_lnnse(:) = 1.0_dp - multi_objective_nse_lnnse(:) / real(nGaugesTotal,dp) 

    ! write(*,*) 'multi_objective_nse_lnnse = ',multi_objective_nse_lnnse 

    ! clean up 
    deallocate( runoff_agg, runoff_obs ) 
    deallocate( runoff_obs_mask ) 

  END FUNCTION multi_objective_nse_lnnse

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_power6_nse_lnnse

  !>        \brief Objective function of combined NSE and lnNSE with power of 5
  !>               i.e. the p-norm with p=5.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>        \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>                           {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>        and the logarithmic Nash-Sutcliffe model efficiency coefficient \f$ lnNSE \f$
  !>        \f[ lnNSE = 1 - \frac{\sum_{i=1}^N (\ln Q_{obs}(i) - \ln Q_{model}(i))^2}
  !>                             {\sum_{i=1}^N (\ln Q_{obs}(i) - \bar{\ln Q_{obs}})^2} \f]
  !>        are calculated and added up equally weighted:
  !>        \f[ obj\_value = \sqrt[6]{(1-NSE)^6 + (1-lnNSE)^6} \f]
  !>        The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_power6_nse_lnnse &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ \sqrt[6]{(1-NSE)^6 + (1-lnNSE)^6} \f$ will be returned such that
  !>             it can be minimized and converges to 0.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_power6_nse_lnnse(para)

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date March 2014
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_power6_nse_lnnse(parameterset)

    use mo_errormeasures,    only: nse, lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_power6_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    real(dp), parameter :: onesixth = 1.0_dp/6.0_dp

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_power6_nse_lnnse = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! NSE + lnNSE
       objective_power6_nse_lnnse = objective_power6_nse_lnnse + &
            ( (1.0_dp-nse(  runoff_obs, runoff_agg, mask=runoff_obs_mask) )**6 + &
            (1.0_dp-lnnse(runoff_obs, runoff_agg, mask=runoff_obs_mask) )**6 )**onesixth
    end do
    ! objective function value which will be minimized
    objective_power6_nse_lnnse = objective_power6_nse_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_power6_nse_lnnse = ', objective_power6_nse_lnnse
    ! pause

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_power6_nse_lnnse

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_kge

  !>        \brief Objective function of KGE.

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore, the Kling-Gupta model efficiency coefficient \f$ KGE \f$
  !>                       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>                 where \n
  !>                       \f$ r \f$ = Pearson product-moment correlation coefficient \n
  !>                       \f$ \alpha \f$ = ratio of similated mean to observed mean \n
  !>                       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation \n
  !>                 is calculated and the objective function is
  !>                       \f[ obj\_value = 1.0 - KGE \f]
  !>                 \f$(1-KGE)\f$ is the objective since we always apply minimization methods. 
  !>                 The minimal value of \f$(1-KGE)\f$ is 0 for the optimal KGE of 1.0. \n
  !>
  !>                 The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_kge &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ KGE \f$ will be returned such that it can be minimized.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_kge(para)

  !     LITERATURE
  !>        Gupta, Hoshin V., et al. "Decomposition of the mean squared error and NSE performance criteria: 
  !>        Implications for improving hydrological modelling." Journal of Hydrology 377.1 (2009): 80-91.


  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date August 2014
  !         Modified, R. Kumar & O. Rakovec, Sep. 2014
  !                   Stephan Thober,        Jan  2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_kge(parameterset)

    use mo_errormeasures,    only: kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_kge

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    !

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_kge = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! KGE
       objective_kge = objective_kge + &
            kge( runoff_obs, runoff_agg, mask=runoff_obs_mask)
    end do
    ! objective_kge = objective_kge + kge(gauge%Q, runoff_model_agg, runoff_model_agg_mask)
    objective_kge = 1.0_dp - objective_kge / real(nGaugesTotal,dp)

    write(*,*) 'objective_kge = ', objective_kge
    ! pause

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_kge

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_multiple_gauges_kge_power6

  !>        \brief combined objective function based on KGE raised to the power 6

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore, the Kling-Gupta model efficiency coefficient \f$ KGE \f$ for a given gauging station
  !>                       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>                 where \n
  !>                       \f$ r \f$ = Pearson product-moment correlation coefficient \n
  !>                       \f$ \alpha \f$ = ratio of similated mean to observed mean \n
  !>                       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation \n
  !>                 is calculated and the objective function for a given gauging station (\f$ i \f$) is
  !>                       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>                 \f$ \phi_{i} \f$ is the objective since we always apply minimization methods. 
  !>                 The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.\n
  !>
  !>                 Finally, the overall \f$ OF \f$ is estimated based on the power-6 norm to 
  !>                 combine the \f$ \phi_{i} \f$ from all gauging stations (\f$ N \f$). 
  !>                 \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }  \f]. \n
  !>                 
  !>                 The observed data \f$ Q_{obs} \f$ are global in this module. 

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !>       \return     real(dp) :: objective_multiple_gauges_kge_power6 &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ OF \f$ will be returned such that it can be minimized.

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_multiple_gauges_kge_power6(para)

  !     LITERATURE
  !>    Gupta, Hoshin V., et al. "Decomposition of the mean squared error and NSE performance criteria: 
  !>    Implications for improving hydrological modelling." Journal of Hydrology 377.1 (2009): 80-91.


  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date March 2015
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM

  FUNCTION objective_multiple_gauges_kge_power6(parameterset)

    use mo_errormeasures,    only: kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_multiple_gauges_kge_power6
    real(dp), parameter                :: onesixth = 1.0_dp/6.0_dp

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    !

    call eval(parameterset, runoff=runoff)
    nGaugesTotal = size(runoff, dim=2)

    objective_multiple_gauges_kge_power6 = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
       ! KGE
       objective_multiple_gauges_kge_power6 = objective_multiple_gauges_kge_power6 + &
            ( (1.0_dp - kge(runoff_obs, runoff_agg, mask=runoff_obs_mask) )/ real(nGaugesTotal,dp) )**6 
    end do
    objective_multiple_gauges_kge_power6 = objective_multiple_gauges_kge_power6**onesixth 
    write(*,*) 'objective_multiple_gauges_kge_power6 = ', objective_multiple_gauges_kge_power6

    deallocate( runoff_agg, runoff_obs, runoff_obs_mask )

  END FUNCTION objective_multiple_gauges_kge_power6


  ! ------------------------------------------------------------------

  ! NAME
  !         extract_runoff

  !>        \brief extracts runoff data from global variables

  !>        \details extracts simulated and measured runoff from global variables,
  !>                 such that they overlay exactly. For measured runoff, only the runoff
  !>                 during the evaluation period are cut, not succeeding nodata values.
  !>                 For simulated runoff, warming days as well as succeeding nodata values
  !>                 are neglected and the simulated runoff is aggregated to the resolution
  !>                 of the observed runoff.\n

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: gaugeID"   - ID of the current gauge to process
  !>        \param[in] "real(dp)    :: runoff(:)" - simulated runoff at this gauge

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)   :: runoff_agg(:)"      - aggregated simulated runoff at this gauge\n
  !>        \param[out] "real(dp)   :: runoff_obs(:)"      - extracted observed runoff\n
  !>        \param[out] "logical    :: runoff_obs_mask(:)" - masking non-negative values in runoff_obs\n

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
  !         see use in this module above

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Jan 2015

  ! ------------------------------------------------------------------
  subroutine extract_runoff( gaugeId, runoff, runoff_agg, runoff_obs, runoff_obs_mask )

    use mo_mrm_global_variables, only: gauge, nMeasPerDay, evalPer, warmingDays_mrm, nTstepDay
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    ! input variables
    integer(i4),               intent(in) :: gaugeId      ! current gauge Id
    real(dp),  dimension(:,:), intent(in) :: runoff       ! simulated runoff

    ! output variables
    real(dp), dimension(:), allocatable, intent(out) :: runoff_agg      ! aggregated simulated
    ! runoff to the resolution
    ! of the measurement
    real(dp), dimension(:), allocatable, intent(out) :: runoff_obs      ! extracted measured 
    ! runoff to exactly the
    ! evaluation period
    logical,  dimension(:), allocatable, intent(out) :: runoff_obs_mask ! mask of no data values
    ! in runoff_obs

    ! local variables
    integer(i4)                         :: iBasin  ! basin id
    integer(i4)                         :: tt      ! timestep counter
    integer(i4)                         :: length  ! length of extracted time series
    integer(i4)                         :: factor  ! between simulated and measured time scale
    integer(i4)                         :: TPD_sim ! simulated Timesteps per Day
    integer(i4)                         :: TPD_obs ! observed Timesteps per Day
    real(dp), dimension(:), allocatable :: dummy

    ! copy time resolution to local variables
    TPD_sim = nTstepDay
    TPD_obs = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo( TPD_sim, TPD_obs) .eq. 0 ) then
       factor = TPD_sim / TPD_obs
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! extract basin Id from gauge Id
    iBasin = gauge%basinId( gaugeId )

    ! get length of evaluation period times TPD_obs
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs

    ! extract measurements
    if ( allocated( runoff_obs ) ) deallocate( runoff_obs )
    allocate( runoff_obs( length ) )
    runoff_obs = gauge%Q( 1 : length, gaugeId )

    ! create mask of observed runoff
    if ( allocated( runoff_obs_mask ) ) deallocate( runoff_obs_mask )
    allocate( runoff_obs_mask( length ) )
    runoff_obs_mask = .false.
    forall(tt=1:length) runoff_obs_mask(tt) = ge( runoff_obs(tt), 0.0_dp)    

    ! extract and aggregate simulated runoff
    if ( allocated( runoff_agg ) ) deallocate( runoff_agg )
    allocate( runoff_agg( length ) )
    ! remove warming days
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_sim
    allocate( dummy( length ) )
    dummy = runoff( warmingDays_mrm(iBasin)*TPD_sim + 1:warmingDays_mrm(iBasin)*TPD_sim + length, gaugeId )
    ! aggregate runoff
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs
    forall(tt=1:length) runoff_agg(tt) = sum( dummy( (tt-1)*factor+1: tt*factor ) ) / &
         real(factor,dp)
    ! clean up
    deallocate( dummy )

  end subroutine extract_runoff

  ! ==================================================================
  ! PRIVATE ROUTINES =================================================
  ! ==================================================================

  ! ------------------------------------------------------------------

  ! NAME
  !         eval

  !>        \brief returns mHM_eval or mRM_eval given preprocessor flag

  !>        \details call mHM_eval if mrm2mhm preprocessor flag is used while
  !>                 compilation or mRM_eval otherwise

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: parameterset" - mHM or mRM parameter set

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !>        \param[out] "real(dp), optional :: runoff(:)" - simulated runoff

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         see use in this module above

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Oct 2015

  ! ------------------------------------------------------------------
  subroutine eval(parameterset, runoff, basin_avg_tws)
#ifdef mrm2mhm
    use mo_mhm_eval, only: mHM_eval
#else
    use mo_mrm_eval, only: mRM_eval
    use mo_message, only: message
#endif
    implicit none
    ! input variables
    real(dp), intent(in) :: parameterset(:)
    ! output variables
    real(dp), allocatable, optional, intent(out) :: runoff(:,:)
    real(dp), allocatable, optional, intent(out) :: basin_avg_tws(:,:)
#ifdef mrm2mhm
    call mHM_eval(parameterset, runoff=runoff, basin_avg_tws=basin_avg_tws)
#else
    call mRM_eval(parameterset, runoff=runoff)
#endif
  end subroutine eval

END MODULE mo_mrm_objective_function_runoff
