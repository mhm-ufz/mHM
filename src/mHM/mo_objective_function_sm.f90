!> \file mo_objective_function_sm.f90

!> \brief Objective Functions for Optimization of mHM against soil moisture.

!> \details This module provides a wrapper for several objective functions used to optimize mHM against soil moisture.\n
!>          Objective functions for optimzation agains runoff can be found in mo_mrm_objective_function_runoff.f90.\n
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

MODULE mo_objective_function_sm

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 moved all runoff only related objectives to mRM

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: objective_sm ! objective function wrapper for soil moisture only

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_sm

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
  !>       \note Input values must be floating points. Numbering of objective functions
  !>             has to be consistent with that of mo_mrm_objective_function_runoff\n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective(para)

  !     LITERATURE

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified,
  !               Oct 2015, Stephan Thober - moved all runoff related objective functions to mRM

  FUNCTION objective_sm(parameterset)

    USE mo_common_variables, ONLY: opti_function

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: parameterset
    REAL(dp)                            :: objective_sm

    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
    case (10)
       ! KGE of catchment average SM
       objective_sm = objective_sm_kge_catchment_avg(parameterset)
    case (11)
       ! pattern dissimilarity (PD) of SM fields
       objective_sm = objective_sm_pd(parameterset)
    case (12)
       ! sum of squared errors of standard_score SM
       objective_sm = objective_sm_sse_standard_score(parameterset)
    case (13)
       ! soil moisture correlation - temporal
       objective_sm = objective_sm_corr(parameterset)
    case default
       stop "Error objective: opti_function not implemented yet."
    end select
    
  END FUNCTION objective_sm

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
    
    use mo_init_states,      only : get_basin_info
    use mo_mhm_eval,         only : mhm_eval
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
    objective_sm_kge_catchment_avg = 0.0_dp

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
    
    use mo_init_states,      only : get_basin_info
    use mo_mhm_eval,         only : mhm_eval
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
    
    use mo_init_states,       only : get_basin_info
    use mo_mhm_eval,         only : mhm_eval
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
    
    use mo_init_states,      only : get_basin_info
    use mo_mhm_eval,         only : mhm_eval
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

END MODULE mo_objective_function_sm
