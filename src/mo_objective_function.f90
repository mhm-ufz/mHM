!> \file mo_objective_function.f90

!> \brief Objective Functions for Optimization of mHM.

!> \details This module provides a wrapper for several objective functions used to optimize mHM.\n
!>          All the objective functions are supposed to be minimized! \n
!> (1)  Q:   1.0 - NSE  \n
!> (2)  Q:   1.0 - lnNSE  \n
!> (3)  Q:   1.0 - 0.5*(NSE+lnNSE)  \n
!> (4)  Q:  -1.0 * loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed  \n
!> (5)  Q:   ((1-NSE)**6+(1-lnNSE)**6)**(1/6)  \n
!> (6)  Q:   SSE  \n
!> (7)  Q:  -1.0 * loglikelihood with trend removed from absolute errors  \n
!> (8)  Q:  -1.0 * loglikelihood with trend removed from the relative errors and then lag(1)-autocorrelation removed \n
!> (9)  Q:  1.0 - KGE (Kling-Gupta efficiency measure)  \n
!> (10) SM: 1.0 - KGE of catchment average soilmoisture \n
!> (11) SM: 1.0 - Pattern dissimilarity (PD) of spatially distributed soil moisture \n
!> (12) SM: Sum of squared errors (SSE) of spatially distributed standard score (normalization) of soil moisture \n
!> (13) SM: 1.0 - average temporal correlation of spatially distributed soil moisture \n
!> (14) Q:  sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) > combination of KGE of every gauging station based on a power-6 norm\n

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
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff

  FUNCTION loglikelihood( parameterset, stddev, stddev_new, likeli_new)
    use mo_moment,           only: mean, correlation
    use mo_linfit,           only: linfit
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_append,           only: append

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev           ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood

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

    call mhm_eval(parameterset, runoff=runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, nGaugesTotal
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

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs )
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
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff

  FUNCTION loglikelihood_kavetski( parameterset, stddev_old, stddev_new, likeli_new)
    use mo_constants,        only: pi_dp
    use mo_moment,           only: stddev, correlation
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_append,           only: append

    implicit none

    real(dp), dimension(:), intent(in)            :: parameterset
    real(dp),               intent(in)            :: stddev_old       ! standard deviation of data
    real(dp),               intent(out), optional :: stddev_new       ! standard deviation of errors using paraset
    real(dp),               intent(out), optional :: likeli_new       ! likelihood using stddev_new, i.e. using new parameter set
    real(dp)                                      :: loglikelihood_kavetski

    ! local
    real(dp), dimension(:,:), allocatable :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    integer(i4)                           :: nmeas
    real(dp), dimension(:),   allocatable :: errors, sigma, eta, y
    real(dp), dimension(:),   allocatable :: obs, calc, out
    real(dp)                              :: a, b, c, vary, vary1
    real(dp)                              :: stddev_tmp
    integer(i4)                           :: npara

    npara = size(parameterset)
    call mhm_eval(parameterset(1:npara-2), runoff=runoff)

  
    ! extract runoff and append it to obs and calc
    do gg = 1, nGaugesTotal
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

    deallocate(runoff, runoff_agg, runoff_obs_mask, runoff_obs )
    deallocate(obs, calc, out, errors, sigma, eta, y)
    
  END FUNCTION loglikelihood_kavetski

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

  FUNCTION loglikelihood_trend_no_autocorr(parameterset, stddev_old, stddev_new, likeli_new)
    use mo_moment,           only: stddev
    use mo_linfit,           only: linfit
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
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

    call mhm_eval(parameterset, runoff=runoff)

    ! extract runoff and append it to obs and calc
    do gg = 1, nGaugesTotal
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
    case (10)
       ! KGE of catchment average SM
       objective = objective_sm_kge_catchment_avg(parameterset)
    case (11)
       ! pattern dissimilarity (PD) of SM fields
       objective = objective_sm_pd(parameterset)
    case (12)
       ! sum of squared errors of standard_score SM
       objective = objective_sm_sse_standard_score(parameterset)
    case (13)
       ! soil moisture correlation - temporal
       objective = objective_sm_corr(parameterset)
    case (14)
       ! combination of KGE of every gauging station based on a power-6 norm \n
       ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) 
       objective = objective_multiple_gauges_kge_power6(parameterset)
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
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff

  FUNCTION objective_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    objective_lnnse = 0.0_dp
    do gg=1,nGaugesTotal
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

  FUNCTION objective_sse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: sse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_sse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    objective_sse = 0.0_dp
    do gg=1,nGaugesTotal
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

  FUNCTION objective_nse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: nse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_nse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    objective_nse = 0.0_dp
    do gg=1,nGaugesTotal
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

  FUNCTION objective_equal_nse_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: nse, lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_equal_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff             ! modelled runoff for a given parameter set
    !                                                           ! dim2=nGauges
    integer(i4)                           :: gg                 ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_obs         ! measured runoff
    real(dp), dimension(:),   allocatable :: runoff_agg         ! aggregated simulated runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask    ! mask for aggregated measured runoff

    call mhm_eval(parameterset, runoff=runoff)

    objective_equal_nse_lnnse = 0.0_dp
    do gg=1,nGaugesTotal
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

  FUNCTION objective_power6_nse_lnnse(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: nse, lnnse

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_power6_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    real(dp), parameter :: onesixth = 1.0_dp/6.0_dp

    call mhm_eval(parameterset, runoff=runoff)

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

  FUNCTION objective_kge(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_kge

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    !

    call mhm_eval(parameterset, runoff=runoff)

    objective_kge = 0.0_dp
    do gg=1,nGaugesTotal
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

  FUNCTION objective_multiple_gauges_kge_power6(parameterset)
    
    use mo_mhm_eval,         only: mhm_eval
    use mo_global_variables, only: nGaugesTotal
    use mo_errormeasures,    only: kge

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_multiple_gauges_kge_power6
    real(dp), parameter                :: onesixth = 1.0_dp/6.0_dp

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    integer(i4)                           :: gg                       ! gauges counter
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    !

    call mhm_eval(parameterset, runoff=runoff)

    objective_multiple_gauges_kge_power6 = 0.0_dp
    do gg=1,nGaugesTotal
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

  !      NAME
  !          objective_sm_kge_catchment_avg

  !>        \brief Objective function for soil moisture.

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>                       soil mloisture (SM) is calculated
  !>                       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>                 where \n
  !>                       \f$ r \f$ = Pearson product-moment correlation coefficient \n
  !>                       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM \n
  !>                       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation \n
  !>                 is calculated and the objective function for a given basin \f$ i \f$ is
  !>                       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>                 \f$ \phi_{i} \f$ is the objective since we always apply minimization methods. 
  !>                 The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.\n
  !>
  !>                 Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6 
  !>                 norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$. 
  !>                 \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f] \n              
  !>                 The observed data L1_sm, L1_sm_mask are global in this module. 

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
  !>       \return     real(dp) :: objective_sm_kge_catchment_avg &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_sm_corr(para)

  !     LITERATURE
  !         none

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    May 2015

  FUNCTION objective_sm_kge_catchment_avg(parameterset)
    
    use mo_mhm_eval,         only : mhm_eval
    use mo_init_states,      only : get_basin_info
    use mo_message,          only : message
    use mo_moment,           only : average
    use mo_errormeasures,    only : KGE
    use mo_string_utils,     only : num2str
    !
    use mo_global_variables, only: nBasins,             & ! number of basins
                                   L1_sm, L1_sm_mask      ! packed measured sm, sm-mask (dim1=ncells, dim2=time)
    use mo_mhm_constants,    only: nodata_dp              ! global nodata value
    
    implicit none

    real(dp), dimension(:), intent(in)      :: parameterset
    real(dp)                                :: objective_sm_kge_catchment_avg

    ! local
    integer(i4)                             :: iBasin                   ! basin loop counter
    integer(i4)                             :: iTime                    ! time loop counter
    integer(i4)                             :: nrows1, ncols1           ! level 1 number of culomns and rows
    integer(i4)                             :: s1, e1                   ! start and end index for the current basin
    integer(i4)                             :: ncells1                  ! ncells1 of level 1
    real(dp), parameter                     :: onesixth = 1.0_dp/6.0_dp ! for sixth root
    real(dp), dimension(:),   allocatable   :: sm_catch_avg_basin       ! spatial average of observed soil moisture
    real(dp), dimension(:),   allocatable   :: sm_opti_catch_avg_basin  ! spatial avergae of modeled  soil moisture
    real(dp), dimension(:,:), allocatable   :: sm_opti                  ! simulated soil moisture
    !                                                                   ! (dim1=ncells, dim2=time)
    logical,  dimension(:),   allocatable   :: mask_times               ! mask for valid sm catchment avg time steps

    call mhm_eval(parameterset, sm_opti=sm_opti)

    ! initialize some variables
    objective_sm_kge_catchment_avg = nodata_dp

    ! loop over basin - for applying power law later on
    do iBasin=1, nBasins

       ! get basin information
       call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, iStart=s1,  iEnd=e1 ) 

       ! allocate
       allocate(mask_times             (size(sm_opti, dim=2)))
       allocate(sm_catch_avg_basin     (size(sm_opti, dim=2)))
       allocate(sm_opti_catch_avg_basin(size(sm_opti, dim=2)))

       ! initalize
       mask_times              = .TRUE.
       sm_catch_avg_basin      = nodata_dp
       sm_opti_catch_avg_basin = nodata_dp

       ! calculate catchment average soil moisture
       do iTime = 1, size(sm_opti, dim=2)

          ! check for enough data points in time for correlation
          if ( all(.NOT. L1_sm_mask(:,iTime)) .OR. (count(L1_sm_mask(:,iTime)) .LE. 10) ) then
             call message('WARNING: objective_sm_kge_catchment_avg: ignored currrent time step since less than')
             call message('         10 valid cells available in soil moisture observation')
             mask_times(iTime) = .FALSE. 
             cycle
          end if
          sm_catch_avg_basin(iTime)      = average(  L1_sm(s1:e1,iTime), mask=L1_sm_mask(s1:e1,iTime))
          sm_opti_catch_avg_basin(iTime) = average(sm_opti(s1:e1,iTime), mask=L1_sm_mask(s1:e1,iTime))
       end do

       ! calculate average soil moisture KGE over all basins with power law
       ! basins are weighted equally ( 1 / real(nBasin,dp))**6
       objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg + &
            ( (1.0_dp-KGE(sm_catch_avg_basin, sm_opti_catch_avg_basin, mask=mask_times)) / real(nBasins,dp) )**6
    end do

    objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg**onesixth
    
    call message('    objective_sm_kge_catchment_avg = ', num2str(objective_sm_kge_catchment_avg,'(F9.5)'))
    
  END FUNCTION objective_sm_kge_catchment_avg

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_sm_corr

  !>        \brief Objective function for soil moisture.

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore the Pearson correlation between observed and modeled soil 
  !>                 moisture on each grid cell \f$ j \f$ is compared
  !>                       \f[ r_j = r^2(SM_{obs}^j, SM_{sim}^j) \f]
  !>                 where \n
  !>                       \f$ r^2\f$        = Pearson correlation coefficient,\n
  !>                       \f$ SM_{obs} \f$  = observed soil moisture,\n
  !>                       \f$ SM_{sim}  \f$ = simulated soil moisture.\n
  !>                 The observed data \f$ SM_{obs} \f$ are global in this module.\n
  !>
  !>                 The the correlation is spatially averaged as 
  !>                 \f[ \phi_{i} = \frac{1}{K} \cdot \sum_{j=1}^K r_j \f] 
  !
  !>                 where \f$ K \f$ denotes the number of valid cells in the study domain.\n
  !
  !>                 Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6 
  !>                 norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$. 
  !>                 \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 }. \f] \n              
  !>                 The observed data L1_sm, L1_sm_mask are global in this module. 

  
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
  !>       \return     real(dp) :: objective_sm_corr &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_sm_corr(para)

  !     LITERATURE
  !         none

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    March 2015

  FUNCTION objective_sm_corr(parameterset)
    
    use mo_mhm_eval,         only : mhm_eval
    use mo_init_states,      only : get_basin_info
    use mo_message,          only : message
    use mo_moment,           only : correlation
    use mo_string_utils,     only : num2str
    !
    use mo_global_variables, only: nBasins,             & ! number of basins
                                   L1_sm, L1_sm_mask      ! packed measured sm, sm-mask (dim1=ncells, dim2=time)


    implicit none

    real(dp), dimension(:), intent(in)      :: parameterset
    real(dp)                                :: objective_sm_corr

    ! local
    integer(i4)                             :: iBasin             ! basin loop counter
    integer(i4)                             :: iCell              ! cell loop counter
    integer(i4)                             :: nrows1, ncols1     ! level 1 number of culomns and rows
    integer(i4)                             :: s1, e1             ! start and end index for the current basin
    integer(i4)                             :: ncells1                 ! ncells1 of level 1
    real(dp)                                :: objective_sm_corr_basin ! basins wise objectives
    real(dp), parameter                     :: onesixth = 1.0_dp/6.0_dp
    real(dp), dimension(:,:), allocatable   :: sm_opti                 ! simulated soil moisture
    !                                                                  ! (dim1=ncells, dim2=time)
 
    call mhm_eval(parameterset, sm_opti=sm_opti)

    ! initialize some variables
    objective_sm_corr          = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin=1, nBasins

       ! init 
       objective_sm_corr_basin = 0.0_dp
       ! get basin information
       call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, iStart=s1,  iEnd=e1 ) 

       ! temporal correlation is calculated on individual gridd cells
       do iCell = s1, e1

          ! check for enough data points in time for correlation
          if ( all(.NOT. L1_sm_mask(iCell,:)) .OR. (count(L1_sm_mask(iCell,:)) .LE. 10) ) then
             call message('WARNING: objective_sm_corr: ignored currrent cell since less than 10 time steps')
             call message('         available in soil moisture observation')
             cycle
          end if
          objective_sm_corr_basin = objective_sm_corr_basin + &
               correlation( L1_sm(iCell,:), sm_opti(iCell,:), mask=L1_sm_mask(iCell,:))
       end do

       ! calculate average soil moisture correlation over all basins with power law
       ! basins are weighted equally ( 1 / real(nBasin,dp))**6
       objective_sm_corr = objective_sm_corr + &
            ( (1.0_dp-objective_sm_corr_basin/ real(nCells1,dp)) / real(nBasins,dp) )**6
    end do

    objective_sm_corr = objective_sm_corr**onesixth
    
    call message('    objective_sm_corr = ', num2str(objective_sm_corr,'(F9.5)'))
    
  END FUNCTION objective_sm_corr

  ! ------------------------------------------------------------------

  !      NAME
  !          objecive_sm_pd

  !>        \brief Objective function for soil moisture.

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore the Pattern Dissimilarity (PD) of observed and modeled soil 
  !>                 moisture fields is calculated - aim: matching spatial patters
  !>                  \f[ E(t) = PD\left( SM_{obs}(t), SM_{sim}(t) \right) \f]
  !>                 where \n
  !>                       \f$ PD \f$        = pattern dissimilarity function,\n
  !>                       \f$ SM_{obs} \f$  = observed soil moisture,\n
  !>                       \f$ SM_{sim}  \f$ = simulated soil moisture.\n
  !>                       \f$ E(t)  \f$     = pattern dissimilarity at timestep \f$ t \f$.\n
  !
  !>                 The the pattern dissimilaity (E) is spatially averaged as 
  !>                 \f[ \phi_{i} = \frac{1}{T} \cdot \sum_{t=1}^T E_t \f] 
  !
  !>                 where \f$ T \f$ denotes the number of time steps.\n
  !
  !>                 Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6 
  !>                 norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$. 
  !>                 \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 } . \f] \n              
  !>                 The observed data L1_sm, L1_sm_mask are global in this module.
  
  !
  !>                  The observed data L1_sm, L1_sm_mask are global in this module.\n

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
  !>       \return     real(dp) :: objecive_sm_pd &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_sm_corr(para)

  !     LITERATURE
  !         none

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    May 2015

  FUNCTION objective_sm_pd(parameterset)
    
    use mo_mhm_eval,          only : mhm_eval
    use mo_init_states,       only : get_basin_info
    use mo_message,           only : message
    use mo_spatialsimilarity, only : PD
    use mo_string_utils,      only : num2str
    !
    use mo_global_variables, only: nBasins,             & ! number of basins
                                   L1_sm, L1_sm_mask      ! packed measured sm, sm-mask (dim1=ncells, dim2=time)
    use mo_mhm_constants,    only: nodata_dp              ! global nodata value
    
    implicit none

    real(dp), dimension(:), intent(in)      :: parameterset
    real(dp)                                :: objective_sm_pd          ! objective function value
 
    ! local
    integer(i4)                             :: iBasin                   ! basin loop counter
    integer(i4)                             :: iTime                    ! time loop counter
    integer(i4)                             :: nrows1, ncols1           ! level 1 number of culomns and rows
    integer(i4)                             :: s1, e1                   ! start and end index for the current basin
    integer(i4)                             :: ncells1                  ! ncells1 of level 1
    real(dp), parameter                     :: onesixth = 1.0_dp/6.0_dp ! for sixth root
    real(dp), dimension(:,:), allocatable   :: mat1, mat2               ! matrices of SM from vectorized arrays
    real(dp), dimension(:),   allocatable   :: pd_time_series           ! pattern dissimilarity (pd) at every time step
    real(dp), dimension(:,:), allocatable   :: sm_opti                  ! simulated soil moisture
    !                                                                   ! (dim1=ncells, dim2=time)
    logical,  dimension(:,:), allocatable   :: mask1                    ! mask of valid cells at level1
    logical,  dimension(:,:), allocatable   :: mask_sm                  ! mask of valid sm cells
    logical,  dimension(:),   allocatable   :: mask_times               ! mask for valid sm catchment avg time steps

    call mhm_eval(parameterset, sm_opti=sm_opti)

    ! initialize some variables
    objective_sm_pd = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin=1, nBasins

       ! get basin information
       call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, iStart=s1,  iEnd=e1, mask=mask1) 

       ! allocate
       allocate(mask_times    (size(sm_opti, dim=2)))
       allocate(pd_time_series(size(sm_opti, dim=2)))
       allocate(mat1   (nrows1, ncols1))
       allocate(mat2   (nrows1, ncols1))
       allocate(mask_sm(nrows1, ncols1))

       ! initalize
       mask_times              = .FALSE.
       pd_time_series          = 0.0_dp

       ! calculate catchment average soil moisture
       do iTime = 1, size(sm_opti, dim=2)
          mat1    = unpack(     L1_sm(s1:e1,iTime), mask1, nodata_dp)
          mat2    = unpack(   sm_opti(s1:e1,iTime), mask1, nodata_dp)
          mask_sm = unpack(L1_sm_mask(s1:e1,iTime), mask1, .FALSE.) 
          pd_time_series = PD(mat1, mat2, mask=mask_sm, valid=mask_times(itime))
       end do

       if (count(mask_times) > 0_i4) then
          ! calculate avergae PD over all basins with power law -basins are weighted equally ( 1 / real(nBasin,dp))**6
          ! print*, 'PD(SM)', sum(pd_time_series, mask=mask_times) / real(count(mask_times), dp) ! MZMZMZMZ
          objective_sm_pd = objective_sm_pd + &
               ((1.0_dp - sum(pd_time_series, mask=mask_times) / real(count(mask_times), dp)) / real(nBasins,dp) )**6
       else
          call message('***ERROR: mo_objective_funtion: objective_sm_pd: No soil moisture observations available!')
          stop
       end if
    end do

    objective_sm_pd = objective_sm_pd**onesixth
    
    call message('    objective_sm_pd = ', num2str(objective_sm_pd,'(F9.5)'))
    
  END FUNCTION objective_sm_pd

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_sm_sse_standard_score

  !>        \brief Objective function for soil moisture.

  !>        \details The objective function only depends on a parameter vector. 
  !>                 The model will be called with that parameter vector and 
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore the sum of squared errors (SSE) of the standard score of observed and 
  !>                 modeled soil moisture is calculated. The standard score or normalization (anomaly)
  !>                 make the objctive function bias insensitive and basically the dynamics of the soil moisture
  !>                 is tried to capture by this obejective function.
  !>                  \f[ phi_i = \sum_{j=1}^K \{ standard\_score( SM_{obs}(j) )- standard\_score(SM_{sim}(j)) \}^2 \f]
  !>                 where \n
  !>                       \f$  standard\_score \f$ = standard score function,\n
  !>                       \f$ SM_{obs} \f$  = observed soil moisture,\n
  !>                       \f$ SM_{sim}  \f$ = simulated soil moisture.\n
  !>                       \f$ K  \f$ = valid elements in study domain.\n
  !
  !>                 Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6 
  !>                 norm to combine the \f$ \phi_{i} \f$ from all basins \f$ N \f$. 
  !>                 \f[ OF = \sqrt[6]{\sum(\phi_{i}/N)^6 }.  \f] \n              
  !>                 The observed data L1_sm, L1_sm_mask are global in this module.

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
  !>       \return     real(dp) :: objective_sm_sse_standard_score &mdash; objective function value 
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_sm_sse_standard_score(para)

  !     LITERATURE
  !         none

  !     HISTORY
  !>        \author  Matthias Zink
  !>        \date    March 2015

  FUNCTION objective_sm_sse_standard_score(parameterset)
    
    use mo_mhm_eval,         only : mhm_eval
    use mo_init_states,      only : get_basin_info
    use mo_message,          only : message
    use mo_errormeasures,    only : SSE
    use mo_standard_score,   only : standard_score
    use mo_string_utils,     only : num2str
    !
    use mo_global_variables, only: nBasins,             & ! number of basins
                                   L1_sm, L1_sm_mask      ! packed measured sm, sm-mask (dim1=ncells, dim2=time)

    implicit none

    real(dp), dimension(:), intent(in)      :: parameterset
    real(dp)                                :: objective_sm_sse_standard_score

    ! local
    integer(i4)                              :: iBasin             ! basin loop counter
    integer(i4)                              :: iCell              ! cell loop counter
    integer(i4)                              :: nrows1, ncols1     ! level 1 number of culomns and rows
    integer(i4)                              :: s1, e1             ! start and end index for the current basin
    integer(i4)                              :: ncells1            ! ncells1 of level 1
    real(dp)                                 :: objective_sm_sse_standard_score_basin ! basins wise objectives
    real(dp),    parameter                   :: onesixth = 1.0_dp/6.0_dp
    real(dp),    dimension(:,:), allocatable :: sm_opti                 ! simulated soil moisture
    !                                                                  ! (dim1=ncells, dim2=time)
 
    call mhm_eval(parameterset, sm_opti=sm_opti)

    ! initialize some variables
    objective_sm_sse_standard_score          = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin=1, nBasins

       ! init 
       objective_sm_sse_standard_score_basin = 0.0_dp
       ! get basin information
       call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, iStart=s1,  iEnd=e1 ) 
      
       ! standard_score signal is calculated on individual grid cells
       do iCell = s1, e1

          ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
          if ( all(.NOT. L1_sm_mask(iCell,:)) .OR. (count(L1_sm_mask(iCell,:)) .LE. 10) ) then
             call message('WARNING: objective_sm_sse_standard_score: ignored currrent cell since less than 10 time steps')
             call message('         available in soil moisture observation')
             cycle
          end if
          objective_sm_sse_standard_score_basin = objective_sm_sse_standard_score_basin + &
               SSE( standard_score(L1_sm(iCell,:), mask=L1_sm_mask(iCell,:)), &
               standard_score(sm_opti(iCell,:), mask=L1_sm_mask(iCell,:)), mask=L1_sm_mask(iCell,:))
         
       end do
       print*, iBasin,  objective_sm_sse_standard_score_basin
       ! calculate average soil moisture correlation over all basins with power law
       ! basins are weighted equally ( 1 / real(nBasin,dp))**6
       objective_sm_sse_standard_score = objective_sm_sse_standard_score + &
            ( objective_sm_sse_standard_score_basin / real(nBasins,dp) )**6
    end do

    objective_sm_sse_standard_score = objective_sm_sse_standard_score**onesixth
    
    call message('    objective_sm_sse_standard_score = ', num2str(objective_sm_sse_standard_score,'(E12.5)'))
    
  END FUNCTION objective_sm_sse_standard_score
  
! private routine

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
  
  use mo_global_variables, only: gauge, evalPer, warmingDays, nTstepDay, nMeasPerDay
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
  dummy = runoff( warmingDays(iBasin)*TPD_sim + 1:warmingDays(iBasin)*TPD_sim + length, gaugeId )
  ! aggregate runoff
  length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs
  forall(tt=1:length) runoff_agg(tt) = sum( dummy( (tt-1)*factor+1: tt*factor ) ) / &
       real(factor,dp)
  ! clean up
  deallocate( dummy )

end subroutine extract_runoff

END MODULE mo_objective_function
