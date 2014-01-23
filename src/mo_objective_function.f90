!> \file mo_objective_function.f90

!> \brief Objective Functions for Optimization of mHM.

!> \details This module provides a wrapper for several objective functions used to optimize mHM.\n
!>          All the objective functions are supposed to be minimized! \n
!>          (1) 1.0 - NSE
!>          (2) 1.0 - lnNSE
!>          (3) 1.0 - 0.5*(NSE+lnNSE)
!>          (4) -1.0 * loglikelihood with removed trend and lag(1)-autocorrelation

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

    print*,'eval likelihood'
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
         runoff_model_agg_mask(tt,gg) = (gauge%Q(tt,gg) .gt. 0.0_dp)    

    ! ----------------------------------------

    nmeas     = size(runoff_model_agg,1) * nGaugesTotal
    
    allocate(obs(nmeas), calc(nmeas), out(nmeas), errors(nmeas))

    obs       = reshape(gauge%Q, (/nmeas/))
    calc      = reshape(runoff_model_agg, (/nmeas/))
    errors(:) = abs( calc(:) - obs(:) )

    ! remove linear trend of errors
    ! out = linfit(obs, errors, a=a, b=b, model2=.False.)
    ! errors(:) = errors(:) - (a + obs(:)*b)
    out = linfit(calc, errors, a=a, b=b, model2=.False.)
    errors(:) = errors(:) - (a + calc(:)*b)

    ! remove lag(1)-autocorrelation of errors
    c = correlation(errors(2:nmeas),errors(1:nmeas-1))
    errors(1:nmeas-1) = errors(2:nmeas) - c*errors(1:nmeas-1)
    errors(nmeas)     = 0.0_dp

    loglikelihood = sum( errors(:) * errors(:) / stddev**2 )
    loglikelihood = -0.5_dp * loglikelihood

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
       ! -loglikelihood with removed trend and lag(1)-autocorrelation
       objective = - loglikelihood( parameterset, 1.0_dp )
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
         runoff_model_agg_mask(tt,gg) = (gauge%Q(tt,gg) .gt. 0.0_dp)    

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
  !          objective_nse

  !>        \brief Objective function of NSE.

  !>        \details The objective function only depends on a parameter vector. 
  !>        The model will be called with that parameter vector and 
  !>        the model output is subsequently compared to observed data.
  !>        Therefore, the Nash-Sutcliffe model efficiency coefficient \f$ NSE \f$
  !>        \f[ NSE = 1 - \frac{\sum_{i=1}^N (Q_{obs}(i) - Q_{model}(i))^2}
  !>                           {\sum_{i=1}^N (Q_{obs}(i) - \bar{Q_{obs}})^2} \f]
  !>        is calculated and added up equally weighted:
  !>        \f[ obj\_value = NSE \f]
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
  !>             Actually, \f$ 1-nse \f$ will be returned such that it can be minimized.

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
         runoff_model_agg_mask(tt,gg) = (gauge%Q(tt,gg) .gt. 0.0_dp)    

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
         runoff_model_agg_mask(tt,gg) = (gauge%Q(tt,gg) .gt. 0.0_dp)    

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


END MODULE mo_objective_function
