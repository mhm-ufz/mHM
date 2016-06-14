!> \file mo_objective_function.f90

!> \brief Objective Functions for Optimization of mHM.

!> \details This module provides a wrapper for several objective functions used to optimize mHM against various variables.\n
!>          If the objective is only regarding runoff move it to mRM/mo_mrm_objective_function_runoff.f90.
!>          If it contains besides runoff another variable like TWS implement it here.\n
!>
!>          All the objective functions are supposed to be minimized! \n
!>               (10) SO: SM:       1.0 - KGE of catchment average soilmoisture \n
!>               (11) SO: SM:       1.0 - Pattern dissimilarity (PD) of spatially distributed soil moisture \n
!>               (12) SO: SM:       Sum of squared errors (SSE) of spatially distributed standard score (normalization)
!>                                        of soil moisture \n
!>               (13) SO: SM:       1.0 - average temporal correlation of spatially distributed soil moisture \n
!>               (15) SO: Q + TWS:  [1.0-KGE(Q)]*RMSE(basin_avg_TWS) - objective function using Q and basin average
!>                                        (standard score) TWS\n
!>               (17) SO: N:        1.0 - KGE of spatio-temporal neutron data, catchment-average \n

!> \authors Juliane Mai
!> \date Dec 2012
!  Modified, Oct 2015, Oldrich Rakovec - added obj. func. 15 (objective_kge_q_rmse_tws) and extract_basin_avg_tws routine

MODULE mo_objective_function

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 moved all runoff only related objectives to mRM

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: objective ! objective function wrapper for soil moisture only

  ! ------------------------------------------------------------------

CONTAINS

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

  FUNCTION objective(parameterset)

    USE mo_common_variables, ONLY: opti_function

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: parameterset
    REAL(dp)                            :: objective

    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
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
    case (15)
       ! KGE for Q * RMSE for basin_avg TWS (standarized scored)
       objective = objective_kge_q_rmse_tws(parameterset)
   case (17)
      ! KGE of catchment average SM
      objective = objective_neutrons_kge_catchment_avg(parameterset)
    case default
       stop "Error objective: opti_function not implemented yet."
    end select

  END FUNCTION objective

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
       ! print*, iBasin,  objective_sm_sse_standard_score_basin
       ! calculate average soil moisture correlation over all basins with power law
       ! basins are weighted equally ( 1 / real(nBasin,dp))**6
       objective_sm_sse_standard_score = objective_sm_sse_standard_score + &
            ( objective_sm_sse_standard_score_basin / real(nBasins,dp) )**6
    end do

    objective_sm_sse_standard_score = objective_sm_sse_standard_score**onesixth

    call message('    objective_sm_sse_standard_score = ', num2str(objective_sm_sse_standard_score,'(E12.5)'))

  END FUNCTION objective_sm_sse_standard_score


  ! -----------------------------------------------------------------

  !      NAME
  !          objective_kge_q_rmse_tws

  !>        \brief Objective function of KGE for runoff and RMSE for basin_avg TWS (standarized scores)

  !>        \details Objective function of KGE for runoff and RMSE for basin_avg TWS (standarized scores)

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
  !>       \return     real(dp) :: objective_kge_q_rmse_tws &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n
  !>             Actually, \f$ KGE \f$ will be returned such that it can be minimized. \n
  !>             For TWS required at least 5 years of data!

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_kge(para)

  !     LITERATURE
  !>        Gupta, Hoshin V., et al. "Decomposition of the mean squared error and NSE performance criteria:
  !>        Implications for improving hydrological modelling." Journal of Hydrology 377.1 (2009): 80-91.


  !     HISTORY
  !>        \author Oldrich Rakovec, Rohini Kumar
  !>        \date Oct. 2015
  !         Modified, Oct 2015, Stephan Thober - moved tws optimization from mo_mrm_objective_function_runoff here

  FUNCTION objective_kge_q_rmse_tws(parameterset)

    use mo_mhm_eval,             only: mhm_eval
    use mo_constants,            only: eps_dp
    use mo_errormeasures,        only: kge, rmse
    use mo_global_variables,     only: nBasins, evalPer
    use mo_julian,               only: caldat
    use mo_message,              only: message
    use mo_mrm_constants,        only: nodata_dp
    use mo_moment,               only: mean
    use mo_temporal_aggregation, only: day2mon_average
    use mo_standard_score,       only: classified_standard_score
    use mo_string_utils,         only: num2str
#ifdef mrm2mhm
    use mo_mrm_objective_function_runoff, only: extract_runoff
#endif

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_kge_q_rmse_tws

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:,:) :: tws                      ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges

    integer(i4)                           :: gg, ii, pp, mmm          ! gauges counter, basin counter, month counters
    integer(i4)                           :: year, month, day
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(nBasins)          :: initTime

    real(dp), dimension(:),   allocatable :: runoff_agg            ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs            ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask       ! mask for measured runoff

    real(dp), dimension(:),   allocatable :: tws_sim               ! simulated tws
    real(dp), dimension(:),   allocatable :: tws_obs               ! measured tws
    logical,  dimension(:),   allocatable :: tws_obs_mask          ! mask for measured tws

    integer(i4)                           :: nMonths               ! total number of months
    integer(i4), dimension(:),allocatable :: month_classes         ! vector with months' classes

    !
    real(dp), DIMENSION(:), allocatable   :: tws_sim_m, tws_obs_m           ! monthly values original time series
    real(dp), DIMENSION(:), allocatable   :: tws_sim_m_anom, tws_obs_m_anom ! monthly values anomaly time series
    logical,  DIMENSION(:), allocatable   :: tws_obs_m_mask   !
    real(dp), dimension(:), allocatable   :: rmse_tws, kge_q  ! rmse_tws(nBasins), kge_q(nGaugesTotal)
    real(dp)                              :: rmse_tws_avg, kge_q_avg ! obj. functions

    ! obtain hourly values of runoff and tws:
    call mHM_eval(parameterset, runoff=runoff, basin_avg_tws=tws)

    !--------------------------------------------
    !! TWS
    !--------------------------------------------

    ! allocate
    allocate( rmse_tws(nBasins) )
    rmse_tws(:) = nodata_dp

    do ii=1, nBasins

       ! extract tws the same way as runoff using mrm
       call extract_basin_avg_tws( ii, tws, tws_sim, tws_obs, tws_obs_mask )

       ! check for whether to proceed with this basin or not
       ! potentially 5 years of data
       if (count(tws_obs_mask) .lt.  365*5 ) then
          deallocate (tws_sim, tws_obs, tws_obs_mask)
          call message('objective_kge_q_rmse_tws: Length of TWS data of Basin ', trim(adjustl(num2str(ii)))  , &
                       ' less than 5 years: these data are not considered')
          cycle
       else
          ! get initial time of the evaluation period
          initTime(ii) = real(evalPer(ii)%julStart,dp)

          ! get calendar days, months, year
          call caldat(int(initTime(ii)), yy=year, mm=month, dd=day)

          ! calculate monthly averages from daily values of the model
          call day2mon_average(tws_sim,year,month,day,tws_sim_m, misval=nodata_dp)

          ! remove mean from modelled time series
          tws_sim_m(:) = tws_sim_m(:) - mean(tws_sim_m(:))

          ! calculate monthly averages from daily values of the observations, which already have removed the long-term mean
          call day2mon_average(tws_obs,year,month,day,tws_obs_m, misval=nodata_dp)

          ! get number of months for given basin
          nMonths = size(tws_obs_m)

          allocate ( month_classes(  nMonths ) )
          allocate ( tws_obs_m_mask( nMonths ) )
          allocate ( tws_obs_m_anom( nMonths ) )
          allocate ( tws_sim_m_anom( nMonths ) )

          month_classes(:) = 0
          tws_obs_m_mask(:) = .TRUE.
          tws_obs_m_anom(:) = nodata_dp
          tws_sim_m_anom(:) = nodata_dp

          ! define months' classes
          mmm = month
          do pp = 1, nMonths
             month_classes(pp) = mmm
             if (mmm .LT. 12) then
                mmm = mmm + 1
             else
                mmm = 1
             end if
          end do

          ! define mask for missing data in observations (there are always data for simulations)
          where( abs( tws_obs_m - nodata_dp) .lt. eps_dp ) tws_obs_m_mask = .FALSE.

          ! calculate standard score
          tws_obs_m_anom = classified_standard_score(tws_obs_m, month_classes, mask = tws_obs_m_mask)
          tws_sim_m_anom = classified_standard_score(tws_sim_m, month_classes, mask = tws_obs_m_mask)
          rmse_tws(ii) = rmse(tws_sim_m_anom,tws_obs_m_anom,  mask = tws_obs_m_mask)

          deallocate ( month_classes )
          deallocate ( tws_obs_m )
          deallocate ( tws_sim_m )
          deallocate ( tws_obs_m_mask )
          deallocate ( tws_sim_m_anom )
          deallocate ( tws_obs_m_anom )
          deallocate (tws_sim, tws_obs, tws_obs_mask)
       endif

    end do

    rmse_tws_avg = sum( rmse_tws(:), abs( rmse_tws - nodata_dp) .gt. eps_dp ) / &
                                        real(count( abs( rmse_tws - nodata_dp) .gt. eps_dp),dp)
    deallocate(rmse_tws)

    !--------------------------------------------
    !! RUNOFF
    !--------------------------------------------
#ifdef mrm2mhm
    nGaugesTotal = size(runoff, dim=2)
    allocate( kge_q(nGaugesTotal))
    kge_q(:) = nodata_dp

    do gg=1, nGaugesTotal

       ! extract runoff
       call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )

       ! check for whether to proceed with this basin or not
       ! potentially 5 years of data
       pp = count(runoff_agg .ge. 0.0_dp )
       if (pp .lt.  365*5 ) then
          deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
          cycle
       else
          ! calculate KGE for each basin:
          kge_q(gg) = kge( runoff_obs, runoff_agg, mask=runoff_obs_mask)
          deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
       end if

    end do

    ! calculate average KGE value for runoff
    kge_q_avg = sum( kge_q(:), abs( kge_q - nodata_dp) .gt. eps_dp ) / &
                     real(count( abs( kge_q - nodata_dp) .gt. eps_dp ),dp)
    deallocate(kge_q)
#else
    call message('***ERROR: objective_kge_q_rmse_tws: missing routing module for optimization')
    stop
#endif

    !
    objective_kge_q_rmse_tws = rmse_tws_avg * ( 1._dp - kge_q_avg )

    call message('    objective_kge_q_rmse_tws = ', num2str(objective_kge_q_rmse_tws,'(F9.5)'))

  END FUNCTION objective_kge_q_rmse_tws

  ! ==================================================================
  ! PRIVATE ROUTINES =================================================
  ! ==================================================================


  ! ------------------------------------------------------------------

  ! NAME
  !         extract_basin_avg_tws

  !>        \brief extracts basin average tws data from global variables

  !>        \details extracts simulated and measured basin average tws from global variables,
  !>                 such that they overlay exactly. For measured tws, only the tws
  !>                 during the evaluation period are cut, not succeeding nodata values.
  !>                 For simulated tws, warming days as well as succeeding nodata values
  !>                 are neglected.\n

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: basinID"   - ID of the current basin to process
  !>        \param[in] "real(dp)    :: tws(:)"    - simulated basin average tws

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)   :: tws_sim(:)"      - simulated basin average tws at this basin\n
  !>        \param[out] "real(dp)   :: tws_obs(:)"      - extracted observed basin average tws\n
  !>        \param[out] "logical    :: tws_obs_mask(:)" - masking non-negative values in tws_obs\n

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
  !         Modified from the module extract_runoff by S. Thober
  !>        \author O. Rakovec
  !>        \date Oct 2015
  !         Modified, Oct 2015, Stephan Thober - moved subroutine to objective_function_sm

  ! ------------------------------------------------------------------
  subroutine extract_basin_avg_tws( basinId, tws, tws_sim, tws_obs, tws_obs_mask )

    use mo_global_variables, only: basin_avg_TWS_obs, nMeasPerDay_TWS, evalPer, warmingDays, nTstepDay
    use mo_message,          only: message
    use mo_constants,        only: eps_dp
    use mo_mrm_constants,    only: nodata_dp

    implicit none

    ! input variables
    integer(i4),               intent(in) :: basinId      ! current basin Id
    real(dp),  dimension(:,:), intent(in) :: tws          ! simulated basin average tws

    ! output variables
    real(dp), dimension(:), allocatable, intent(out) :: tws_sim         ! aggregated simulated
    ! tws to the resolution
    ! of the measurement
    real(dp), dimension(:), allocatable, intent(out) :: tws_obs         ! extracted measured
    ! tws to exactly the
    ! evaluation period
    logical,  dimension(:), allocatable, intent(out) :: tws_obs_mask    ! mask of no data values
    ! in tws_obs

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
    TPD_obs = nMeasPerDay_TWS

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo( TPD_sim, TPD_obs) .eq. 0 ) then
       factor = TPD_sim / TPD_obs
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if

    ! extract basin Id
    iBasin = basin_avg_TWS_obs%basinId( basinId )

    ! get length of evaluation period times TPD_obs
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs

    ! extract measurements
    if ( allocated( tws_obs ) ) deallocate( tws_obs )
    allocate( tws_obs( length ) )
    tws_obs = basin_avg_TWS_obs%TWS( 1 : length, basinId )

    ! create mask of observed tws
    if ( allocated( tws_obs_mask ) ) deallocate( tws_obs_mask )
    allocate( tws_obs_mask( length ) )
    tws_obs_mask = .TRUE.
    where( abs( tws_obs - nodata_dp) .lt. eps_dp ) tws_obs_mask = .FALSE.

    ! extract and aggregate simulated tws
    if ( allocated( tws_sim ) ) deallocate( tws_sim )
    allocate( tws_sim( length ) )
    ! remove warming days
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_sim
    allocate( dummy( length ) )
    dummy = tws( warmingDays(iBasin)*TPD_sim + 1:warmingDays(iBasin)*TPD_sim + length, basinId )
    ! aggregate tws
    length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs
    forall(tt=1:length) tws_sim(tt) = sum( dummy( (tt-1)*factor+1: tt*factor ) ) / &
         real(factor,dp)
    ! clean up
    deallocate( dummy )

  end subroutine extract_basin_avg_tws

  ! ------------------------------------------------------------------

  !      NAME
  !          objective_neutrons_kge_catchment_avg

  !>        \brief Objective function for neutrons.

  !>        \details The objective function only depends on a parameter vector.
  !>                 The model will be called with that parameter vector and
  !>                 the model output is subsequently compared to observed data.\n
  !>
  !>                 Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>                       neutrons (N) is calculated
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
  !>                 The observed data L1_neutronsdata, L1_neutronsdata_mask are global in this module.

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
  !>       \return     real(dp) :: objective_neutrons_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         obj_value = objective_neutrons_kge_catchment_avg(para)

  !     LITERATURE
  !         none

  !     HISTORY
  !>        \author  Martin Schroen
  !>        \date    Jun 2015

  FUNCTION objective_neutrons_kge_catchment_avg(parameterset)

    use mo_init_states,      only : get_basin_info
    use mo_mhm_eval,         only : mhm_eval
    use mo_message,          only : message
    use mo_moment,           only : average
    use mo_errormeasures,    only : KGE
    use mo_string_utils,     only : num2str
    !
    use mo_global_variables, only: nBasins,             & ! number of basins
         L1_neutronsdata, L1_neutronsdata_mask      ! packed measured neutrons, neutrons-mask (dim1=ncells, dim2=time)
    use mo_mhm_constants,    only: nodata_dp              ! global nodata value

    implicit none

    real(dp), dimension(:), intent(in)      :: parameterset
    real(dp)                                :: objective_neutrons_kge_catchment_avg

    ! local
    integer(i4)                             :: iBasin                   ! basin loop counter
    integer(i4)                             :: iTime                    ! time loop counter
    integer(i4)                             :: nrows1, ncols1           ! level 1 number of culomns and rows
    integer(i4)                             :: s1, e1                   ! start and end index for the current basin
    integer(i4)                             :: ncells1                  ! ncells1 of level 1
    real(dp), parameter                     :: onesixth = 1.0_dp/6.0_dp ! for sixth root
    real(dp), dimension(:),   allocatable   :: neutrons_catch_avg_basin       ! spatial average of observed neutrons
    real(dp), dimension(:),   allocatable   :: neutrons_opti_catch_avg_basin  ! spatial avergae of modeled  neutrons
    real(dp), dimension(:,:), allocatable   :: neutrons_opti                  ! simulated neutrons
    !                                                                   ! (dim1=ncells, dim2=time)
    logical,  dimension(:),   allocatable   :: mask_times               ! mask for valid neutrons catchment avg time steps

    call mhm_eval(parameterset, neutrons_opti=neutrons_opti)

    ! initialize some variables
    objective_neutrons_kge_catchment_avg = 0.0_dp

    ! loop over basin - for applying power law later on
    do iBasin=1, nBasins

       ! get basin information
       call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, iStart=s1,  iEnd=e1 )

       ! allocate
       allocate(mask_times             (size(neutrons_opti, dim=2)))
       allocate(neutrons_catch_avg_basin     (size(neutrons_opti, dim=2)))
       allocate(neutrons_opti_catch_avg_basin(size(neutrons_opti, dim=2)))

       ! initalize
       mask_times              = .TRUE.
       neutrons_catch_avg_basin      = nodata_dp
       neutrons_opti_catch_avg_basin = nodata_dp

       ! calculate catchment average soil moisture
       do iTime = 1, size(neutrons_opti, dim=2)

          ! check for enough data points in time for correlation
          if ( all(.NOT. L1_neutronsdata_mask(:,iTime))) then
              write (*,*) 'WARNING: neutrons data at time ', iTime, ' is empty.'
             !call message('WARNING: objective_neutrons_kge_catchment_avg: ignored currrent time step since less than')
             !call message('         10 valid cells available in soil moisture observation')
             mask_times(iTime) = .FALSE.
             cycle
          end if
          neutrons_catch_avg_basin(iTime)      = average(  L1_neutronsdata(s1:e1,iTime), mask=L1_neutronsdata_mask(s1:e1,iTime))
          neutrons_opti_catch_avg_basin(iTime) = average(neutrons_opti(s1:e1,iTime), mask=L1_neutronsdata_mask(s1:e1,iTime))
       end do

       ! calculate average neutrons KGE over all basins with power law
       ! basins are weighted equally ( 1 / real(nBasin,dp))**6
       objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg + &
            ( (1.0_dp-KGE(neutrons_catch_avg_basin, neutrons_opti_catch_avg_basin, mask=mask_times)) / real(nBasins,dp) )**6
    end do

    objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg**onesixth

    call message('    objective_neutrons_kge_catchment_avg = ', num2str(objective_neutrons_kge_catchment_avg,'(F9.5)'))

END FUNCTION objective_neutrons_kge_catchment_avg

END MODULE mo_objective_function
