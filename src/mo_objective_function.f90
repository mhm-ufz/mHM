!> \file mo_objective_function.f90

!> \brief Objective Functions for Optimization of mHM.

!> \details This module provides a wrapper for several objective functions used to optimize mHM.\n
!>          All the objective functions are supposed to be minimized! \n
!>          (1) 1.0 - NSE  \n
!>          (2) 1.0 - lnNSE  \n
!>          (3) 1.0 - 0.5*(NSE+lnNSE)  \n
!>          (4) -1.0 * loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed  \n
!>          (5) ((1-NSE)**6+(1-lnNSE)**6)**(1/6)  \n
!>          (6) SSE  \n
!>          (7) -1.0 * loglikelihood with trend removed from absolute errors  \n
!>          (8) -1.0 * loglikelihood with linear error model and lag(1)-autocorrelation of the relative errors \n
!>          (9) 1.0 - KGE \n

!> \authors Juliane Mai
!> \date Dec 2012

MODULE mo_objective_function

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified 

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: loglikelihood    ! loglikelihood with errormodel including linear trend and lag(1)-correlation
  PUBLIC :: objective        ! objective function wrapper

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood

  !>        \brief Logarithmic likelihood function with removed linear trend and Lag(1)-autocorrelation.

  !>        \details The logarithmis likelihood function is used when mHM runs in MCMC mode.\n
  !>        It can also be used for optimization when selecting the likelihood in the namelist as \e opti\_function.\n

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
  !>       \return     real(dp) :: loglikelihood &mdash; logarithmic likelihood using given stddev 
  !>                                                     but remove optimal trend and lag(1)-autocorrelation in errors 
  !>                                                     (absolute between running model with parameterset and observation) 

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !     EXAMPLE
  !         para = (/ 1._dp, 2._dp, 3._dp, -999._dp, 5._dp, 6._dp /)
  !         stddev = 0.5_dp
  !         log_likeli = loglikelihood(para, stddev, stddev_new=stddev_new, likeli_new=likeli_new)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, 

  FUNCTION loglikelihood( parameterset, stddev, stddev_new, likeli_new)
    use mo_moment,           only: mean, correlation
    use mo_linfit,           only: linfit
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev           ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b, c
    real(dp)                              :: stddev_tmp

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg),0.0_dp))    

    ! ----------------------------------------

    nmeas     = size(runoff_model_agg,1) * nGaugesTotal
    
    allocate(obs(nmeas), calc(nmeas), out(nmeas), errors(nmeas))

    obs       = reshape(gauge%Q, (/nmeas/))
    calc      = reshape(runoff_model_agg, (/nmeas/))
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
    loglikelihood = sum( errors(:) * errors(:) / stddev**2 )
    loglikelihood = -0.5_dp * loglikelihood

    write(*,*) '-loglikelihood = ', -loglikelihood

    stddev_tmp = sqrt(sum( (errors(:) - mean(errors)) * (errors(:) - mean(errors))) / real(nmeas-1,dp))
    if (present(stddev_new)) then
       stddev_new = stddev_tmp
    end if
    if (present(likeli_new)) then
       likeli_new = sum( errors(:) * errors(:) / stddev_tmp**2 )
       likeli_new = -0.5_dp * likeli_new
    end if

    deallocate(runoff, runoff_model_agg, runoff_model_agg_mask)
    deallocate(obs, calc, out, errors)
    
  END FUNCTION loglikelihood

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood_kavetski

  !>        \brief Logarithmised likelihood with linear error model and lag(1)-autocorrelation
  !>               of the relative errors.

  !>        \details This loglikelihood uses a linear error model and a lag(1)-autocorrelation
  !>                 on the relative errors. This is approach 2 of the paper Evin et al. (WRR, 2013).
  !>
  !>                 This is opti_method = 8.
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
  !>       \return     real(dp) :: loglikelihood_kavetski &mdash; logarithmic likelihood using given stddev 
  !>                                                     but remove optimal trend and lag(1)-autocorrelation in errors 
  !>                                                     (absolute between running model with parameterset and observation) 

  !     RESTRICTIONS
  !>       \note Does not work with MCMC yet.

  !     EXAMPLE
  !         para = (/ 1._dp, 2._dp, 3._dp, -999._dp, 5._dp, 6._dp /)
  !         log_likeli = loglikelihood_kavetski(para, 1.0_dp)

  !     LITERATURE
  !         Evin et al., WRR 49, 4518-4524, 2013

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date Mar 2014
  !         Modified, 

  FUNCTION loglikelihood_kavetski( parameterset, stddev_old, stddev_new, likeli_new)
    use mo_constants,        only: pi_dp
    use mo_moment,           only: stddev, correlation
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev_old       ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood_kavetski

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors, sigma, eta, y
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b, c, vary, vary1
    real(dp)                              :: stddev_tmp
    integer(i4)                           :: npara

    npara = size(parameterset)
    call mhm_eval(parameterset(1:npara-2), runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg),0.0_dp))    

    ! ----------------------------------------

    nmeas     = size(runoff_model_agg,1) * nGaugesTotal
    
    allocate(obs(nmeas), calc(nmeas), out(nmeas), errors(nmeas), sigma(nmeas), eta(nmeas), y(nmeas))

    obs       = reshape(gauge%Q, (/nmeas/))
    calc      = reshape(runoff_model_agg, (/nmeas/))
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

    ! likelihood of residual errors (leave out ln(1/sqrt(2*pi)))
    vary  = 1.0_dp - c*c
    vary1 = 1.0_dp / vary
    loglikelihood_kavetski = real(nmeas-1,dp)*log(1.0_dp/sqrt(2.0_dp*pi_dp*vary)) &
         - log(sigma(1)) - 0.5_dp*eta(1)*eta(1) &
         - sum(0.5_dp*y(2:nmeas)*y(2:nmeas)*vary1 + log(sigma(2:nmeas)))

    write(*,*) '-loglikelihood_kavetski = ', -loglikelihood_kavetski

    ! This is for the interface of MCMC
    stddev_tmp = stddev_old   ! this is for the compiler so that stddev_old gets used
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

    deallocate(runoff, runoff_model_agg, runoff_model_agg_mask)
    deallocate(obs, calc, out, errors, sigma, eta, y)
    
  END FUNCTION loglikelihood_kavetski

  ! ------------------------------------------------------------------

  !      NAME
  !          loglikelihood_trend_no_autocorr

  !>        \brief Logarithmic likelihood function with linear trend removed.

  !>        \details The logarithmis likelihood function is used when mHM runs in MCMC mode.\n
  !>        It can also be used for optimization when selecting the likelihood in the namelist as \e opti\_function.\n

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

  FUNCTION loglikelihood_trend_no_autocorr(parameterset, stddev_old, stddev_new, likeli_new)
    use mo_moment,           only: stddev
    use mo_linfit,           only: linfit
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev_old       ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood_trend_no_autocorr

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b
    real(dp)                              :: stddev_tmp

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg),0.0_dp))    

    ! ----------------------------------------

    nmeas     = size(runoff_model_agg,1) * nGaugesTotal
    
    allocate(obs(nmeas), calc(nmeas), out(nmeas), errors(nmeas))

    obs       = reshape(gauge%Q, (/nmeas/))
    calc      = reshape(runoff_model_agg, (/nmeas/))
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

    deallocate(runoff, runoff_model_agg, runoff_model_agg_mask)
    deallocate(obs, calc, out, errors)
    
  END FUNCTION loglikelihood_trend_no_autocorr

  ! ------------------------------------------------------------------

  !      NAME
  !          objective

  !>        \brief Wrapper for objective functions.

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
  !         obj_value = objective(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, 

  FUNCTION objective(parameterset)

    USE mo_global_variables, ONLY: opti_function

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN)  :: parameterset
    REAL(DP)                            :: objective

    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
    case (1)
       ! 1.0-nse
       objective = objective_nse(parameterset)
    case (2)
       ! 1.0-lnnse
       objective = objective_lnnse(parameterset)
    case (3)
       ! 1.0-0.5*(nse+lnnse)
       objective = objective_equal_nse_lnnse(parameterset)
    case (4)
       ! -loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
       objective = - loglikelihood(parameterset, 1.0_dp)
    case (5)
       ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
       objective = objective_power6_nse_lnnse(parameterset)
    case (6)
       ! SSE
       objective = objective_sse(parameterset)
    case (7)
       ! -loglikelihood with trend removed from absolute errors
       objective = -loglikelihood_trend_no_autocorr(parameterset, 1.0_dp)
    case (8)
       ! -loglikelihood with trend removed from relative errors and then lag(1)-autocorrelation removed
       objective = -loglikelihood_kavetski(parameterset, 1.0_dp)
    case (9)
       ! KGE
       objective = objective_kge(parameterset)
    case default
       stop "Error objective: opti_function not implemented yet."
    end select
    
  END FUNCTION objective

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
  !         Modified, 

  FUNCTION objective_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: lnnse
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_lnnse = 0.0_dp
    do gg=1,nGaugesTotal
       ! lnNSE
       objective_lnnse = objective_lnnse + &
            lnnse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
    end do
    ! objective function value which will be minimized
    objective_lnnse = 1.0_dp - objective_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_lnnse = ',objective_lnnse
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
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

  FUNCTION objective_sse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: sse
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_sse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_sse = 0.0_dp
    do gg=1,nGaugesTotal
       objective_sse = objective_sse + &
            sse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
    end do
    ! objective_sse = objective_sse + sse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_sse = objective_sse / real(nGaugesTotal,dp)

    write(*,*) 'objective_sse = ', objective_sse
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
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
  !         Modified, 

  FUNCTION objective_nse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: nse
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_nse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_nse = 0.0_dp
    do gg=1,nGaugesTotal
       objective_nse = objective_nse + &
            nse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
    end do
    ! objective_nse = objective_nse + nse(gauge%Q, runoff_model_agg) !, runoff_model_agg_mask)
    objective_nse = 1.0_dp - objective_nse / real(nGaugesTotal,dp)

    write(*,*) 'objective_nse = ',objective_nse
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
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
  !         Modified, 

  FUNCTION objective_equal_nse_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: nse, lnnse
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_equal_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_equal_nse_lnnse = 0.0_dp
    do gg=1,nGaugesTotal
       ! NSE
       objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
            nse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
       ! lnNSE
       objective_equal_nse_lnnse = objective_equal_nse_lnnse + &
            lnnse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
    end do
    ! objective function value which will be minimized
    objective_equal_nse_lnnse = 1.0_dp - 0.5_dp * objective_equal_nse_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_equal_nse_lnnse = ',objective_equal_nse_lnnse
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
  END FUNCTION objective_equal_nse_lnnse


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

  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date March 2014

  FUNCTION objective_power6_nse_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: nse, lnnse
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_power6_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff
    real(dp), parameter :: onesixth = 1.0_dp/6.0_dp

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_power6_nse_lnnse = 0.0_dp
    do gg=1, nGaugesTotal
       ! NSE + lnNSE
       objective_power6_nse_lnnse = objective_power6_nse_lnnse + &
            ((1.0_dp-nse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg)))**6 + &
            (1.0_dp-lnnse(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg)))**6)**onesixth
    end do
    ! objective function value which will be minimized
    objective_power6_nse_lnnse = objective_power6_nse_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_power6_nse_lnnse = ', objective_power6_nse_lnnse
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
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
  !>                 where
  !>                       \f[ r \f] = Pearson product-moment correlation coefficient
  !>                       \f[ \alpha \f] = ratio of similated mean to observed mean 
  !>                       \f[ \beta  \f] = ratio of similated standard deviation to observed standard deviation
  !>                 is calculated and the objective function is
  !>                       \f[ obj\_value = 1.0 - KGE \f]
  !>                 (1-KGE) is the objective since we always apply minimization methods. 
  !>                 The minimal value of (1-KGE) is 0 for the optimal KGE of 1.0.\n
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
  !>    Gupta, Hoshin V., et al. "Decomposition of the mean squared error and NSE performance criteria: 
  !>    Implications for improving hydrological modelling." Journal of Hydrology 377.1 (2009): 80-91.


  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date August 2014
  !         Modified, R. Kumar & O. Rakovec, Sep. 2014

  FUNCTION objective_kge(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nTstepDay, nMeasPerDay, nGaugesTotal, warmingDays
    use mo_global_variables, only: gauge
    use mo_errormeasures,    only: kge
    use mo_mhm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_utils,            only: ge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_kge

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: nTimeSteps               ! number of modelled timesteps in total
    integer(i4)                           :: timestepsPerDay_modelled !
    integer(i4)                           :: timestepsPerDay_measured !
    integer(i4)                           :: multiple                 ! timestepsPerDay_modelled = 
    !                                                                 !     multiple * timestepsPerDay_measured
    integer(i4)                           :: tt                       ! timestep counter
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:,:), allocatable :: runoff_model_agg         ! aggregated measured runoff
    logical,  dimension(:,:), allocatable :: runoff_model_agg_mask    ! mask for aggregated measured runoff
    !

    call mhm_eval(parameterset, runoff=runoff)

    ! simulated timesteps per day 
    timestepsPerDay_modelled = nTstepDay
    ! measured timesteps per day 
    timestepsPerDay_measured = nMeasPerDay

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo(timestepsPerDay_modelled,timestepsPerDay_measured) .eq. 0 ) then
       multiple = timestepsPerDay_modelled/timestepsPerDay_measured
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! total number of simulated timesteps
    ntimeSteps = size(runoff,1)
    
    ! allocation and initialization
    allocate( runoff_model_agg((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    allocate( runoff_model_agg_mask((nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple, nGaugesTotal) )
    runoff_model_agg      = nodata_dp
    runoff_model_agg_mask = .false. ! take mask of observation

    ! average <multiple> datapoints
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg(tt,gg) = &
         sum(runoff((tt-1)*multiple+timestepsPerDay_modelled*warmingDays+1: &
         tt*multiple+timestepsPerDay_modelled*warmingDays,gg))/real(multiple,dp)
    ! set mask
    forall(tt=1:(nTimeSteps-timestepsPerDay_modelled*warmingDays)/multiple,gg=1:nGaugesTotal) &
         runoff_model_agg_mask(tt,gg) = (ge(gauge%Q(tt,gg), 0.0_dp))    

    objective_kge = 0.0_dp
    do gg=1,nGaugesTotal
       ! KGE
       objective_kge = objective_kge + &
            kge(gauge%Q(:,gg), runoff_model_agg(:,gg), mask=runoff_model_agg_mask(:,gg))
    end do
    ! objective_kge = objective_kge + kge(gauge%Q, runoff_model_agg, runoff_model_agg_mask)
    objective_kge = 1.0_dp - objective_kge / real(nGaugesTotal,dp)

    write(*,*) 'objective_kge = ', objective_kge
    ! pause

    deallocate( runoff_model_agg )
    deallocate( runoff_model_agg_mask )
    
  END FUNCTION objective_kge
  
END MODULE mo_objective_function
