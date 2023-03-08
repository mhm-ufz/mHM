!> \file mo_objective_function.f90
!> \brief   \copybrief mo_objective_function
!> \details \copydetails mo_objective_function

! ToDo: change comment for OF 15

!> \brief Objective Functions for Optimization of mHM.
!> \details This module provides a wrapper for several objective functions used to optimize mHM against various
!!       variables.
!!       If the objective is only regarding runoff move it to mRM/mo_mrm_objective_function_runoff.f90.
!!       If it contains besides runoff another variable like TWS implement it here.
!!
!!       All the objective functions are supposed to be minimized!
!!       - (10) SO: SM:       1.0 - KGE of catchment average soilmoisture
!!       - (11) SO: SM:       1.0 - Pattern dissimilarity (PD) of spatially distributed soil moisture
!!       - (12) SO: SM:       Sum of squared errors (SSE) of spatially distributed standard score (normalization) of soil moisture
!!       - (13) SO: SM:       1.0 - average temporal correlation of spatially distributed soil moisture
!!       - (15) SO: Q + TWS:  [1.0-KGE(Q)]*RMSE(domain_avg_TWS) - objective function using Q and domain average (standard score) TWS
!!       - (17) SO: N:        1.0 - KGE of spatio-temporal neutron data, catchment-average
!!       - (27) SO: ET:       1.0 - KGE of catchment average evapotranspiration
!> \changelog
!! - Oldrich Rakovec Oct 2015
!!   - added obj. func. 15 (objective_kge_q_rmse_tws) and extract_domain_avg_tws routine, former basin_avg
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!> \authors Juliane Mai
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_objective_function

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 moved all runoff only related objectives to mRM

  USE mo_kind, ONLY : i4, dp
  use mo_optimization_utils, only : eval_interface
  use mo_message, only : message, error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: objective
#ifdef MPI
  PUBLIC :: objective_master, objective_subprocess ! objective function wrapper for soil moisture only
#endif

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        objective

  !    PURPOSE
  !>       \brief Wrapper for objective functions.

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
  ! Stephan Thober Oct 2015 - moved all runoff related objective functions to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective(parameterset, eval, arg1, arg2, arg3)

    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : opti_function

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: objective

    real(dp), dimension(6) :: multiple_objective


    if (present(arg1) .or. present(arg2) .or. present(arg3)) then
      call error_message("Error mo_objective_function: Received unexpected argument, check optimization settings")
    end if

    ! set these to nan so compiler does not complain
    if (present(arg2)) then
      arg2 = nodata_dp
    end if
    if (present(arg3)) then
      arg3 = nodata_dp
    end if

    select case (opti_function)
    case (10)
      ! KGE of catchment average SM
      objective = objective_sm_kge_catchment_avg(parameterset, eval)
    case (11)
      ! pattern dissimilarity (PD) of SM fields
      objective = objective_sm_pd(parameterset, eval)
    case (12)
      ! sum of squared errors of standard_score SM
      objective = objective_sm_sse_standard_score(parameterset, eval)
    case (13)
      ! soil moisture correlation - temporal
      objective = objective_sm_corr(parameterset, eval)
    case (15)
      ! KGE for Q * RMSE for domain_avg TWS (standarized scored)
      objective = objective_kge_q_rmse_tws(parameterset, eval)
    case (17)
      ! KGE of catchment average SM
      objective = objective_neutrons_kge_catchment_avg(parameterset, eval)
    case (27)
      ! KGE of catchment average ET
      objective = objective_et_kge_catchment_avg(parameterset, eval)
    case (28)
      !  KGE for Q + SSE for SM (standarized scored)
      objective = objective_kge_q_sm_corr(parameterset, eval)
    case (29)
      !  KGE for Q + KGE of catchment average ET
      objective = objective_kge_q_et(parameterset, eval)
    case (30)
      ! KGE for Q * RMSE for domain_avg ET (standarized scored)
      objective = objective_kge_q_rmse_et(parameterset, eval)
    case (33)
      multiple_objective = objective_q_et_tws_kge_catchment_avg(parameterset, eval)
      objective = multiple_objective(1)
    case (34)
      ! KGE for Q * Absolute-Error for BFI
      objective = objective_kge_q_BFI(parameterset, eval)

    case default
      call error_message("Error objective: opti_function not implemented yet.")
    end select

  END FUNCTION objective

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_master

  !    PURPOSE
  !>       \brief Wrapper for objective functions.

  !>       \details The functions sends the parameterset to the subprocess, receives
  !>       the partial objective and calculates the final objective
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
  ! Stephan Thober Oct 2015 - moved all runoff related objective functions to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Maren Kaluza Jun 2019 - parallel version

#ifdef MPI
  FUNCTION objective_master(parameterset, eval, arg1, arg2, arg3)

    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : opti_function
    use mo_common_mpi_tools, only : distribute_parameterset
    use mo_common_variables, only : domainMeta
    use mo_string_utils, only : num2str
    use mpi_f08

    implicit none

    REAL(dp), DIMENSION(:), INTENT(IN) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: objective_master

    REAL(dp) :: partial_objective

    real(dp), dimension(6) :: multiple_partial_objective

    real(dp), dimension(6) :: multiple_master_objective

    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    integer(i4) :: iproc, nproc

    integer(i4) :: ierror

    type(MPI_Status) :: status


    if (present(arg1) .or. present(arg2) .or. present(arg3)) then
      call error_message("Error mo_objective_function: Received unexpected argument, check optimization settings")
    end if

    ! set these to nan so compiler does not complain
    if (present(arg2)) then
      arg2 = nodata_dp
    end if
    if (present(arg3)) then
      arg3 = nodata_dp
    end if

    call distribute_parameterset(parameterset)
    select case (opti_function)
    case (10 : 13, 17, 27 : 29)
      call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
      objective_master = 0.0_dp
      do iproc = 1, nproc - 1
        call MPI_Recv(partial_objective, 1, MPI_DOUBLE_PRECISION, iproc, 0, domainMeta%comMaster, status, ierror)
        objective_master = objective_master + partial_objective
      end do
      objective_master = objective_master**onesixth
    case (15)
      ! KGE for Q * RMSE for domain_avg TWS (standarized scored)
      call error_message("case 15, objective_kge_q_rmse_tws not implemented in parallel yet")
    case (30)
      ! KGE for Q * RMSE for domain_avg ET (standarized scored)
      ! objective_master = objective_kge_q_rmse_et(parameterset, eval)
      call message("case 30, objective_kge_q_rmse_et not implemented in parallel yet")
    case(33)
      call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
      objective_master = 0.0_dp
      multiple_master_objective(:) = 0.0_dp
      do iproc = 1, nproc - 1
        call MPI_Recv(multiple_partial_objective, 6, MPI_DOUBLE_PRECISION, iproc, 0, domainMeta%comMaster, status, ierror)
        multiple_master_objective = multiple_master_objective + multiple_partial_objective
      end do
      objective_master = objective_master + &
        (multiple_master_objective(1)+multiple_master_objective(2)+multiple_master_objective(3))
      objective_master = (objective_master/multiple_master_objective(4))**onesixth

    case default
      call error_message("Error objective_master: opti_function not implemented yet.")
    end select

    select case (opti_function)
    case(10)
      call message('    objective_sm_kge_catchment_avg = ', num2str(objective_master, '(F9.5)'))
    case(11)
      call message('    objective_sm_pd = ', num2str(objective_master, '(F9.5)'))
    case(12)
      call message('    objective_sm_sse_standard_score = ', num2str(objective_master, '(E12.5)'))
    case(13)
      call message('    objective_sm_corr = ', num2str(objective_master, '(F9.5)'))
    case(17)
      call message('    objective_neutrons_kge_catchment_avg = ', num2str(objective_master, '(F9.5)'))
    case(27)
      call message('    objective_et_kge_catchment_avg = ', num2str(objective_master, '(F9.5)'))
    case(28)
      call message('    objective_kge_q_sm_corr = ', num2str(objective_master, '(F9.5)'))
    case(29)
      call message('    objective_kge_q_et = ', num2str(objective_master, '(F9.5)'))
    case(33)
      call message('    objective_q_et_tws_kge_catchment_avg = ', num2str(objective_master, '(F9.5)'))
    case default
      call error_message("Error objective_master: opti_function not implemented yet, this part of the code should never execute.")
    end select

  END FUNCTION objective_master

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_subprocess

  !    PURPOSE
  !>       \brief Wrapper for objective functions.

  !>       \details The function receives the parameterset from the master
  !>       process, selects the objective function case defined in a namelist,
  !>       i.e. the global variable \e opti\_function.
  !>       It returns the partial objective function value for a specific parameter set.

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
  ! Stephan Thober Oct 2015 - moved all runoff related objective functions to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Maren Kaluza Jun 2019 - parallel version

  subroutine objective_subprocess(eval, arg1, arg2, arg3)

    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : opti_function
    use mo_common_mpi_tools, only : get_parameterset
    use mo_common_variables, only : domainMeta
    use mpi_f08

    implicit none

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp), optional, intent(in) :: arg1

    real(dp), optional, intent(out) :: arg2

    real(dp), optional, intent(out) :: arg3

    REAL(dp) :: partial_objective

    real(dp), dimension(6) :: multiple_partial_objective

    REAL(dp), DIMENSION(:), allocatable :: parameterset

    integer(i4) :: ierror

    type(MPI_Status) :: status

    logical :: do_obj_loop

    do ! a do loop without condition runs until exit
      call MPI_Recv(do_obj_loop, 1, MPI_LOGICAL, 0, 0, domainMeta%comMaster, status, ierror)

      if (.not. do_obj_loop) exit

      if (present(arg1) .or. present(arg2) .or. present(arg3)) then
        call error_message("Error mo_objective_function: Received unexpected argument, check optimization settings")
      end if

      ! set these to nan so compiler does not complain
      if (present(arg2)) then
        arg2 = nodata_dp
      end if
      if (present(arg3)) then
        arg3 = nodata_dp
      end if
      call get_parameterset(parameterset)
      select case (opti_function)
      case (10)
        ! KGE of catchment average SM
        partial_objective = objective_sm_kge_catchment_avg(parameterset, eval)
      case (11)
        ! pattern dissimilarity (PD) of SM fields
        partial_objective = objective_sm_pd(parameterset, eval)
      case (12)
        ! sum of squared errors of standard_score SM
        partial_objective = objective_sm_sse_standard_score(parameterset, eval)
      case (13)
        ! soil moisture correlation - temporal
        partial_objective = objective_sm_corr(parameterset, eval)
      case (15)
        ! KGE for Q * RMSE for domain_avg TWS (standarized scored)
        ! partial_objective = objective_kge_q_rmse_tws(parameterset, eval)
        call error_message("Error objective_subprocess: case 15 not supported with MPI.")
      case (17)
        ! KGE of catchment average SM
        partial_objective = objective_neutrons_kge_catchment_avg(parameterset, eval)
      case (27)
        ! KGE of catchment average ET
        partial_objective = objective_et_kge_catchment_avg(parameterset, eval)
      case (28)
        !  KGE for Q + SSE for SM (standarized scored)
        partial_objective = objective_kge_q_sm_corr(parameterset, eval)
      case (29)
        !  KGE for Q + KGE of catchment average ET
        partial_objective = objective_kge_q_et(parameterset, eval)
      case (30)
        ! KGE for Q * RMSE for domain_avg ET (standarized scored)
        ! partial_objective = objective_kge_q_rmse_et(parameterset, eval)
        call error_message("Error objective_subprocess: case 30 not supported with MPI.")
      case(33)
        multiple_partial_objective = objective_q_et_tws_kge_catchment_avg(parameterset, eval)
      case default
        call error_message("Error objective_subprocess: opti_function not implemented yet.")
      end select

      select case (opti_function)
      case (10 : 13, 17, 27 : 29)
        call MPI_Send(partial_objective,1, MPI_DOUBLE_PRECISION,0,0,domainMeta%comMaster,ierror)
      case(33)
        call MPI_Send(multiple_partial_objective, 6, MPI_DOUBLE_PRECISION,0,0,domainMeta%comMaster,ierror)
      case default
        call error_message("Error objective_subprocess: this part should not be executed -> error in the code.")
      end select

      deallocate(parameterset)
    end do

  END subroutine objective_subprocess

#endif
  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       soil mloisture (SM) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given domain \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date May 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_kge_catchment_avg(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, domainMeta
    use mo_errormeasures, only : KGE
    use mo_global_variables, only : L1_smObs
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_kge_catchment_avg

    ! domain loop counter
    integer(i4) :: iDomain

    ! time loop counter
    integer(i4) :: iTime

    ! number of time steps in simulated SM
    integer(i4) :: n_time_steps

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid timesteps
    real(dp) :: invalid_times
#ifndef MPI
    ! for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    ! spatial average of observed soil moisture
    real(dp), dimension(:), allocatable :: sm_catch_avg_domain

    ! spatial avergae of modeled  soil moisture
    real(dp), dimension(:), allocatable :: sm_opti_catch_avg_domain

    type(optidata_sim), dimension(:), allocatable :: smOptiSim

    ! mask for valid sm catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    allocate(smOptiSim(domainMeta%nDomains))
    call eval(parameterset, smOptiSim = smOptiSim)

    ! initialize some variables
    objective_sm_kge_catchment_avg = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! get domain information
      ncells1 = level1(iDomain)%ncells

      ! allocate
      allocate(mask_times              (size(smOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(sm_catch_avg_domain     (size(smOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(sm_opti_catch_avg_domain(size(smOptiSim(iDomain)%dataSim, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      sm_catch_avg_domain = nodata_dp
      sm_opti_catch_avg_domain = nodata_dp

      invalid_times = 0.0_dp
      ! calculate catchment average soil moisture
      n_time_steps = size(smOptiSim(iDomain)%dataSim, dim = 2)
      do iTime = 1, n_time_steps

        ! check for enough data points in timesteps for KGE calculation
        ! more then 10 percent avaiable in current field
        if (count(L1_smObs(iDomain)%maskObs(:, iTime)) .LE. (0.10_dp * real(nCells1, dp))) then
          invalid_times = invalid_times + 1.0_dp
          mask_times(iTime) = .FALSE.
          cycle
        end if
        sm_catch_avg_domain(iTime) = average(L1_smObs(iDomain)%dataObs(:, iTime), mask = L1_smObs(iDomain)%maskObs(:, iTime))
        sm_opti_catch_avg_domain(iTime) = average(smOptiSim(iDomain)%dataSim(:, iTime), mask = L1_smObs(iDomain)%maskObs(:, iTime))
      end do

      ! user information about invalid times
      if (invalid_times .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid timesteps (.LT. 10 valid data points).')
        call message('                          Fraction of invalid timesteps: ', &
                num2str(invalid_times / real(n_time_steps, dp), '(F4.2)'))
      end if


      ! calculate average soil moisture KGE over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg + &
              ((1.0_dp - KGE(sm_catch_avg_domain, sm_opti_catch_avg_domain, mask = mask_times)) / &
                        real(domainMeta%overallNumberOfDomains, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(sm_catch_avg_domain)
      deallocate(sm_opti_catch_avg_domain)
      call smOptiSim(iDomain)%destroy()
    end do
    deallocate(smOptiSim)

#ifndef MPI
    objective_sm_kge_catchment_avg = objective_sm_kge_catchment_avg**onesixth

    call message('    objective_sm_kge_catchment_avg = ', num2str(objective_sm_kge_catchment_avg, '(F9.5)'))
#endif


  END FUNCTION objective_sm_kge_catchment_avg

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_q_et_tws_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for et, tws and q.

  !>       \details The feature of this objective function is the
  !>                separation of the eval call into four
  !>                calls, each with another index list. The subroutine eval then only
  !>                uses the indices from that index list internally instead of having loops
  !>                over all domains. The integer array domainMeta%optidata decides which
  !>                indices to use. Therefore the array is split into disjunct subsets, and,
  !>                if chosen wisely in the namelist, also covers all domains.
  !>
  !>                With this the eval calls sum up in a way that for each domain eval was
  !>                called at most once, but for different opti_data.

  !    HISTORY
  !>       \authors Maren Kaluza

  !>       \date July 2019

  ! Modifications:

  FUNCTION objective_q_et_tws_kge_catchment_avg(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : domainMeta
    use mo_global_variables, only : L1_etObs, L1_twsaObs
    use mo_errormeasures, only : kge
    use mo_moment, only : average
    use mo_string_utils, only : num2str
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    !> the parameterset passed to the eval subroutine
    real(dp), dimension(:), intent(in) :: parameterset
    !> the eval subroutine called by this objective function
    procedure(eval_interface), INTENT(IN), POINTER :: eval
    !> the return value of the objective function. In this case it is
    !> an array to provide the possibility to weight the outcome accordingly
    real(dp), dimension(6) :: objective_q_et_tws_kge_catchment_avg

    !> domain loop counter
    integer(i4) :: iDomain

    !> counter for short loops
    integer(i4) :: i
#ifndef MPI
    !> for sixth root
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    !> modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff
    !> number of all gauges, aquired via runoff
    integer(i4) :: nGaugesTotal

    !> aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    !> measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    !> mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    !> kge_q(nGaugesTotal)
    real(dp) :: kge_q

    !> gauges counter
    integer(i4) :: gg, iCell

    !> number of q domains
    integer(i4) :: nQDomains

    !> number of et domains
    integer(i4) :: nEtDomains

    !> number of tws domains
    integer(i4) :: nTwsDomains

    !> number of TWS and ET domains (providing both)
    integer(i4) :: nEtTwsDomains

    !> index array of ET domains
    integer(i4), dimension(:), allocatable :: opti_domain_indices_ET

    !> index array of TWS domains
    integer(i4), dimension(:), allocatable :: opti_domain_indices_TWS

    !> index array of TWS and ET domains (providing both)
    integer(i4), dimension(:), allocatable :: opti_domain_indices_ET_TWS

    !> index array of ET domains
    integer(i4), dimension(:), allocatable :: opti_domain_indices_Q

    !> simulated et
    type(optidata_sim), dimension(:), allocatable :: etOptiSim

    !> simulated tws
    type(optidata_sim), dimension(:), allocatable :: twsOptiSim

    !> simulated twsa (anomaly)
    type(optidata_sim), dimension(:), allocatable :: twsaOptiSim

    real(dp) :: kge_tws

    real(dp) :: kge_et

    integer(i4) :: numberOfSummands


    ! initialize some variables
    objective_q_et_tws_kge_catchment_avg(:) = 0.0_dp
    kge_tws = 0.0_dp
    kge_et = 0.0_dp
    kge_q = 0.0_dp
    numberOfSummands = 0
    !--------------------------------------------
    ! ET & TWS
    !--------------------------------------------
    ! eval runs to get simulated output for et and tws
    ! before each eval call we generate an index list of the domains for which
    ! eval should be called. Read details for further information
    call init_indexarray_for_opti_data(domainMeta, 6, nEtTwsDomains, opti_domain_indices_ET_TWS)
    if (nEtTwsDomains > 0) then
      allocate( etOptiSim(domainMeta%nDomains))
      allocate(twsOptiSim(domainMeta%nDomains))
      allocate(twsaOptiSim(domainMeta%nDomains))
      call eval(parameterset, opti_domain_indices = opti_domain_indices_ET_TWS, &
                                                 twsOptiSim = twsOptiSim, etOptiSim = etOptiSim)
      ! for all domains that have ET and TWS
      do i = 1, size(opti_domain_indices_ET_TWS)
        iDomain = opti_domain_indices_ET_TWS(i)
        call convert_tws_to_twsa(twsOptiSim(iDomain), L1_twsaObs(iDomain), twsaOptiSim(iDomain))
        do iCell = 1, size(L1_etObs(iDomain)%maskObs(:, :), dim = 1)
          kge_et = kge_et + &
            (1.0_dp - KGE(L1_etObs(iDomain)%dataObs(iCell, :), etOptiSim(iDomain)%dataSim(iCell, :),&
                           mask = L1_etObs(iDomain)%maskObs(iCell, :)))**6
          numberOfSummands = numberOfSummands + 1
        end do
        do iCell = 1, size(L1_twsaObs(iDomain)%maskObs(:, :), dim = 1)
          kge_tws = kge_tws + &
            (1.0_dp - KGE(L1_twsaObs(iDomain)%dataObs(iCell, :), twsaOptiSim(iDomain)%dataSim(iCell, :),&
                           mask = L1_twsaObs(iDomain)%maskObs(iCell, :)))**6
          numberOfSummands = numberOfSummands + 1
        end do
        ! deallocate
        call etOptiSim(iDomain)%destroy()
        call twsOptiSim(iDomain)%destroy()
        call twsaOptiSim(iDomain)%destroy()
      end do
      deallocate(etOptiSim)
      deallocate(twsOptiSim)
      deallocate(twsaOptiSim)
     ! write(0,*) 'nEtTwsDomains, kge_tws', nEtTwsDomains, kge_tws
     ! write(0,*) 'nEtTwsDomains, kge_et', nEtTwsDomains, kge_et
    end if
    !--------------------------------------------
    ! TWS
    !--------------------------------------------
    ! eval runs to get simulated output for tws
    ! before each eval call we generate an index list of the domains for which
    ! eval should be called. Read details for further information
    call init_indexarray_for_opti_data(domainMeta, 3, nTwsDomains, opti_domain_indices_TWS)
    if (nTwsDomains > 0) then
      allocate(twsOptiSim(domainMeta%nDomains))
      allocate(twsaOptiSim(domainMeta%nDomains))
      call eval(parameterset, opti_domain_indices = opti_domain_indices_TWS, twsOptiSim = twsOptiSim)
      ! for all domains that have ET and TWS
      do i = 1, size(opti_domain_indices_TWS)
        iDomain = opti_domain_indices_TWS(i)
        call convert_tws_to_twsa(twsOptiSim(iDomain), L1_twsaObs(iDomain), twsaOptiSim(iDomain))
        do iCell = 1, size(L1_twsaObs(iDomain)%maskObs(:, :), dim = 1)
          kge_tws = kge_tws + &
            (1.0_dp - KGE(L1_twsaObs(iDomain)%dataObs(iCell, :), twsaOptiSim(iDomain)%dataSim(iCell, :),&
                           mask = L1_twsaObs(iDomain)%maskObs(iCell, :)))**6
          numberOfSummands = numberOfSummands + 1
        end do
        call twsOptiSim(iDomain)%destroy()
      end do
      deallocate(twsOptiSim)
    !  write(0,*) 'nTwsDomains, kge_tws', nTwsDomains, kge_tws
    end if
    objective_q_et_tws_kge_catchment_avg(2) = kge_tws

    !--------------------------------------------
    ! ET
    !--------------------------------------------
    ! eval runs to get simulated output for et
    ! before each eval call we generate an index list of the domains for which
    ! eval should be called. Read details for further information
    call init_indexarray_for_opti_data(domainMeta, 5, nEtDomains, opti_domain_indices_ET)
    if (nEtDomains > 0) then
      allocate(etOptiSim(domainMeta%nDomains))
      call eval(parameterset, opti_domain_indices = opti_domain_indices_ET, etOptiSim = etOptiSim)
      ! for all domains that have ET and TWS
      do i = 1, size(opti_domain_indices_ET)
        iDomain = opti_domain_indices_ET(i)
        do iCell = 1, size(L1_etObs(iDomain)%maskObs(:, :), dim = 1)
          kge_et = kge_et + &
            (1.0_dp - KGE(L1_etObs(iDomain)%dataObs(iCell, :), etOptiSim(iDomain)%dataSim(iCell, :),&
                           mask = L1_etObs(iDomain)%maskObs(iCell, :)))**6
          numberOfSummands = numberOfSummands + 1
        end do
        call etOptiSim(iDomain)%destroy()
      end do
      deallocate(etOptiSim)
     ! write(0,*) 'nEtDomains, kge_et', nEtDomains, kge_et
    end if
    objective_q_et_tws_kge_catchment_avg(3) = kge_et

    !--------------------------------------------
    ! RUNOFF
    !--------------------------------------------
    ! eval runs to get simulated output for runoff
    ! before the eval call we generate an index list of the domains for which
    ! eval should be called. Read details for further information
    ! ToDo:  The arrays for qTin, qTout, will be rewritten in the other calls when
    ! Q is not called last. Change that for more flexibility
    call init_indexarray_for_opti_data(domainMeta, 1, nQDomains, opti_domain_indices_Q)

    if (nQDomains > 0) then
      call eval(parameterset, opti_domain_indices = opti_domain_indices_Q, runoff = runoff)
      nGaugesTotal = size(runoff, dim = 2)
      do gg = 1, nGaugesTotal

        ! extract runoff
        call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

        kge_q = kge_q + &
              (1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask))**6
        numberOfSummands = numberOfSummands + 1
        deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
      end do
     ! write(0,*) 'nQDomains, kge_q', nQDomains, kge_q
    end if
    objective_q_et_tws_kge_catchment_avg(1) = kge_q

    objective_q_et_tws_kge_catchment_avg(4) = real(numberOfSummands, dp)


#ifndef MPI
    objective_q_et_tws_kge_catchment_avg(1) = ((kge_q+kge_et+kge_tws)/real(numberOfSummands, dp))**onesixth

    call message('    objective_q_et_tws_kge_catchment_avg = ', &
                      num2str(objective_q_et_tws_kge_catchment_avg(1), '(F9.5)'))
#endif

  END FUNCTION objective_q_et_tws_kge_catchment_avg

  ! ------------------------------------------------------------------

  !    NAME
  !        init_indexarray_for_opti_data

  !    PURPOSE
  !>       \brief creates an index array of the inidices of the domains eval
  !>              should MPI process.
  !
  !>       \details The data type domainMeta contains an array optidata of size
  !>                domainMeta%nDomains, telling us, which domains should be
  !>                optimized with which opti_data. This subroutine splits all
  !>                domains assigned to a process and returns an index list
  !>                corresponding to the value of domainMeta%optidata.
  !>
  !>                The index array opti_domain_indices can then be passed
  !>                as an optional argument to the eval subroutine. The
  !>                eval then instead of using loops over all domains only
  !>                uses the passed indices.
  !>
  !>                This subroutine also returns the size of that array since it
  !>                helps with the calculations of the optimization in the end.

  !    HISTORY
  !>       \authors Maren Kaluza

  !>       \date July 2019
  subroutine init_indexarray_for_opti_data(domainMeta, optidataOption, nOptiDomains, opti_domain_indices)
    use mo_common_types, only: domain_meta
    !> meta data for all domains assigned to that process
    type(domain_meta),                                intent(in)    :: domainMeta
    !> which opti data should be used in the eval called after calling this subroutine
    integer(i4),                                      intent(in)    :: optidataOption
    !> number of domains that will be optimized in the following eval call
    integer(i4),                                      intent(out)   :: nOptiDomains
    !> the indices of the domains that are to be processed in the following eval call
    integer(i4), dimension(:), allocatable,           intent(out)   :: opti_domain_indices

    !> domain loop counter
    integer(i4) :: iDomain, i

    if (allocated(opti_domain_indices)) deallocate(opti_domain_indices)
    ! count domains on MPI process that use optidata
    nOptiDomains = 0
    do iDomain = 1, domainMeta%nDomains
      if (domainMeta%optidata(iDomain) == optidataOption) nOptiDomains = nOptiDomains + 1
    end do
    ! write indices of these domains into an array
    if (nOptiDomains > 0) then
      allocate(opti_domain_indices(nOptiDomains))
      i = 0
      do iDomain = 1, domainMeta%nDomains
        if (domainMeta%optidata(iDomain) == optidataOption) then
          i = i + 1
          opti_domain_indices(i) = iDomain
        end if
      end do
    end if
  end subroutine init_indexarray_for_opti_data
  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_corr

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the Pearson correlation between observed and modeled soil
  !>       moisture on each grid cell \f$ j \f$ is compared
  !>       \f[ r_j = r^2(SM_{obs}^j, SM_{sim}^j) \f]
  !>       where
  !>       \f$ r^2\f$        = Pearson correlation coefficient,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       The observed data \f$ SM_{obs} \f$ are global in this module.

  !>       The the correlation is spatially averaged as
  !>       \f[ \phi_{i} = \frac{1}{K} \cdot \sum_{j=1}^K r_j \f]
  !>       where \f$ K \f$ denotes the number of valid cells in the study domain.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 }. \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_corr &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date March 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_corr(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_variables, only : level1, domainMeta
    use mo_global_variables, only : L1_smObs
    use mo_moment, only : correlation
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_corr

    ! domain loop counter
    integer(i4) :: iDomain

    ! cell loop counter
    integer(i4) :: iCell

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! domains wise objectives
    real(dp) :: objective_sm_corr_domain

#ifndef MPI
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    type(optidata_sim), dimension(:), allocatable :: smOptiSim


    allocate(smOptiSim(domainMeta%nDomains))
    call eval(parameterset, smOptiSim = smOptiSim)

    ! initialize some variables
    objective_sm_corr = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! init
      objective_sm_corr_domain = 0.0_dp
      ! get domain information
      ncells1 = level1(iDomain)%ncells

      invalid_cells = 0.0_dp
      ! temporal correlation is calculated on individual gridd cells

      do iCell = 1, size(L1_smObs(iDomain)%maskObs(:, :), dim = 1)

        ! check for enough data points in time for correlation
        if (count(L1_smObs(iDomain)%maskObs(iCell, :)) .LE. 0.10_dp * real(size(L1_smObs(iDomain)%dataObs(:, :), dim = 2), dp)) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if
        objective_sm_corr_domain = objective_sm_corr_domain + &
                correlation(L1_smObs(iDomain)%dataObs(iCell, :), smOptiSim(iDomain)%dataSim(iCell, :), &
                                                           mask = L1_smObs(iDomain)%maskObs(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invalid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if


      ! calculate average soil moisture correlation over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_sm_corr = objective_sm_corr + &
              ((1.0_dp - objective_sm_corr_domain / real(nCells1, dp)) / real(domainMeta%overallNumberOfDomains, dp))**6
    end do
#ifndef MPI
    objective_sm_corr = objective_sm_corr**onesixth

    call message('    objective_sm_corr = ', num2str(objective_sm_corr, '(F9.5)'))
#endif

  END FUNCTION objective_sm_corr

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_pd

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the Pattern Dissimilarity (PD) of observed and modeled soil
  !>       moisture fields is calculated - aim: matching spatial patters
  !>       \f[ E(t) = PD\left( SM_{obs}(t), SM_{sim}(t) \right) \f]
  !>       where
  !>       \f$ PD \f$        = pattern dissimilarity function,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       \f$ E(t)  \f$     = pattern dissimilarity at timestep \f$ t \f$.
  !>       The the pattern dissimilaity (E) is spatially averaged as
  !>       \f[ \phi_{i} = \frac{1}{T} \cdot \sum_{t=1}^T E_t \f]
  !>       where \f$ T \f$ denotes the number of time steps.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - \phi_{i})/N)^6 } . \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objecive_sm_pd &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date May 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_pd(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1, domainMeta
    use mo_global_variables, only : L1_smObs
    use mo_spatialsimilarity, only : PD
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    ! objective function value
    real(dp) :: objective_sm_pd

    ! domain loop counter
    integer(i4) :: iDomain

    ! time loop counter
    integer(i4) :: iTime

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! for sixth root
#ifndef MPI
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    ! matrices of SM from vectorized arrays
    real(dp), dimension(:, :), allocatable :: mat1, mat2

    ! pattern dissimilarity (pd) at every time step
    real(dp), dimension(:), allocatable :: pd_time_series

    ! simulated soil moisture
    type(optidata_sim), dimension(:), allocatable :: smOptiSim

    ! mask of valid cells at level1
    logical, dimension(:, :), allocatable :: mask1

    ! mask of valid sm cells
    logical, dimension(:, :), allocatable :: mask_sm

    ! mask for valid sm catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    allocate(smOptiSim(domainMeta%nDomains))
    call eval(parameterset, smOptiSim = smOptiSim)

    ! initialize some variables
    objective_sm_pd = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! get domain information
      mask1 = level1(iDomain)%mask
      ncols1 = level1(iDomain)%ncols
      nrows1 = level1(iDomain)%nrows

      ! allocate
      allocate(mask_times    (size(smOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(pd_time_series(size(smOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(mat1   (nrows1, ncols1))
      allocate(mat2   (nrows1, ncols1))
      allocate(mask_sm(nrows1, ncols1))

      ! initalize
      mask_times = .FALSE.
      pd_time_series = 0.0_dp

      ! calculate pattern similarity criterion
      do iTime = 1, size(smOptiSim(iDomain)%dataSim, dim = 2)
        mat1 = unpack(L1_smObs(iDomain)%dataObs(:, iTime), mask1, nodata_dp)
        mat2 = unpack(smOptiSim(iDomain)%dataSim(:, iTime), mask1, nodata_dp)
        mask_sm = unpack(L1_smObs(iDomain)%maskObs(:, iTime), mask1, .FALSE.)
        pd_time_series = PD(mat1, mat2, mask = mask_sm, valid = mask_times(itime))
      end do

      if (count(mask_times) > 0_i4) then
        ! calculate avergae PD over all domains with power law -domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
        objective_sm_pd = objective_sm_pd + &
                ((1.0_dp - sum(pd_time_series, mask = mask_times) / real(count(mask_times), dp)) / &
                                                    real(domainMeta%overallNumberOfDomains, dp))**6
      else
        call error_message('***ERROR: mo_objective_funtion: objective_sm_pd: No soil moisture observations available!')
      end if

      ! deallocate
      deallocate(mask_times)
      deallocate(pd_time_series)
      deallocate(mat1)
      deallocate(mat2)
      deallocate(mask_sm)
      call smOptiSim(iDomain)%destroy()
    end do
    deallocate(smOptiSim)

#ifndef MPI
    objective_sm_pd = objective_sm_pd**onesixth

    call message('    objective_sm_pd = ', num2str(objective_sm_pd, '(F9.5)'))
#endif

  END FUNCTION objective_sm_pd

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_sm_sse_standard_score

  !    PURPOSE
  !>       \brief Objective function for soil moisture.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore the sum of squared errors (SSE) of the standard score of observed and
  !>       modeled soil moisture is calculated. The standard score or normalization (anomaly)
  !>       make the objctive function bias insensitive and basically the dynamics of the soil moisture
  !>       is tried to capture by this obejective function.
  !>       \f[ phi_i = \sum_{j=1}^K \{ standard\_score( SM_{obs}(j) )- standard\_score(SM_{sim}(j)) \}^2 \f]
  !>       where
  !>       \f$  standard\_score \f$ = standard score function,
  !>       \f$ SM_{obs} \f$  = observed soil moisture,
  !>       \f$ SM_{sim}  \f$ = simulated soil moisture.
  !>       \f$ K  \f$ = valid elements in study domain.
  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum(\phi_{i}/N)^6 }.  \f]
  !>       The observed data L1_sm, L1_sm_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_sm_sse_standard_score &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date March 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_sm_sse_standard_score(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_variables, only : level1, domainMeta
    use mo_errormeasures, only : SSE
    use mo_global_variables, only : L1_smObs
    use mo_standard_score, only : standard_score
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_sm_sse_standard_score

    ! domain loop counter
    integer(i4) :: iDomain

    ! cell loop counter
    integer(i4) :: iCell

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! domains wise objectives
    real(dp) :: objective_sm_sse_standard_score_domain

#ifndef MPI
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    ! simulated soil moisture
    type(optidata_sim), dimension(:), allocatable :: smOptiSim


    call eval(parameterset, smOptiSim = smOptiSim)

    ! initialize some variables
    objective_sm_sse_standard_score = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! init
      objective_sm_sse_standard_score_domain = 0.0_dp
      ! get domain information
      nCells1 = level1(iDomain)%nCells

      invalid_cells = 0.0_dp
      ! standard_score signal is calculated on individual grid cells
      do iCell = 1, size(L1_smObs(iDomain)%maskObs(:, :), dim = 1)

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_smObs(iDomain)%maskObs(iCell, :)) .LE. (0.10_dp * real(size(L1_smObs(iDomain)%dataObs, dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if
        objective_sm_sse_standard_score_domain = objective_sm_sse_standard_score_domain + &
                SSE(standard_score(L1_smObs(iDomain)%dataObs(iCell, :), mask = L1_smObs(iDomain)%maskObs(iCell, :)), &
                        standard_score(smOptiSim(iDomain)%dataSim(iCell, :), mask = L1_smObs(iDomain)%maskObs(iCell, :)), &
                                                          mask = L1_smObs(iDomain)%maskObs(iCell, :))

      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invalid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture correlation over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_sm_sse_standard_score = objective_sm_sse_standard_score + &
              (objective_sm_sse_standard_score_domain / real(domainMeta%overallNumberOfDomains, dp))**6
    end do

#ifndef MPI
    objective_sm_sse_standard_score = objective_sm_sse_standard_score**onesixth

    call message('    objective_sm_sse_standard_score = ', num2str(objective_sm_sse_standard_score, '(E12.5)'))
#endif

  END FUNCTION objective_sm_sse_standard_score


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_rmse_tws

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and RMSE for domain_avg TWS (standarized scores)

  !>       \details Objective function of KGE for runoff and RMSE for domain_avg TWS (standarized scores)

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_rmse_tws &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Oldrich Rakovec, Rohini Kumar

  !>       \date Oct. 2015

  ! Modifications:
  ! Stephan Thober Oct 2015 - moved tws optimization from mo_mrm_objective_function_runoff here
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Maren Kaluza Oct 2019 - changed averaging function for tws, this will not produce the same output as before

  FUNCTION objective_kge_q_rmse_tws(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : domainMeta
    use mo_global_variables, only : L1_twsaObs
    use mo_errormeasures, only : rmse
    use mo_julian, only : caldat
    use mo_moment, only : mean
    use mo_standard_score, only : classified_standard_score
    use mo_string_utils, only : num2str
    use mo_temporal_aggregation, only : day2mon_average
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_rmse_tws

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    !> simulated tws
    type(optidata_sim), dimension(:), allocatable :: twsOptiSim

    ! domain counter, month counters
    integer(i4) :: domainID, iDomain, pp, mmm

    integer(i4) :: year, month, day

    real(dp), dimension(domainMeta%nDomains) :: initTime

    ! simulated tws
    real(dp), dimension(:), allocatable :: tws_catch_avg_domain

    ! measured tws
    real(dp), dimension(:), allocatable :: tws_opti_catch_avg_domain

    ! mask for measured tws
    logical, dimension(:), allocatable :: tws_obs_mask

    ! total number of months
    integer(i4) :: nMonths

    ! vector with months' classes
    integer(i4), dimension(:), allocatable :: month_classes

    ! monthly values anomaly time series
    real(dp), DIMENSION(:), allocatable :: tws_sim_m_anom, tws_obs_m_anom

    ! rmse_tws(domainMeta%nDomains)
    real(dp), dimension(:), allocatable :: rmse_tws

    ! obj. functions
    real(dp) :: rmse_tws_avg, kge_q_avg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! kge_q(nGaugesTotal)
    real(dp), dimension(:), allocatable :: kge_q

    ! gauges counter
    integer(i4) :: gg

    ! obtain hourly values of runoff and tws:
    allocate(twsOptiSim(domainMeta%nDomains))
    call eval(parameterset, runoff = runoff, twsOptiSim = twsOptiSim)

    !--------------------------------------------
    !! TWS
    !--------------------------------------------

    ! allocate
    allocate(rmse_tws(domainMeta%nDomains))
    rmse_tws(:) = nodata_dp

    do iDomain = 1, domainMeta%nDomains
      if (.not. (L1_twsaObs(iDomain)%timeStepInput == -2)) then
        call message('objective_kge_q_rmse_tws: current implementation of this subroutine only allows monthly timesteps')
      end if
      domainID = domainMeta%indices(iDomain)

      ! extract tws the same way as runoff using mrm
      ! Note that with the change from tws(iDomain, tt) to tws(tt, :) this
      ! will not work like before and also does maybe not make sense
      call create_domain_avg_tws(iDomain, twsOptiSim, tws_catch_avg_domain, tws_opti_catch_avg_domain, tws_obs_mask)

      ! check for potentially 2 years of data
      if (count(tws_obs_mask) .lt.  12 * 2) then
        call message('objective_kge_q_rmse_tws: Length of TWS data of domain ', trim(adjustl(num2str(domainID))), &
                ' less than 2 years: this is not recommended')
      end if

      ! get initial time of the evaluation period
      initTime(iDomain) = real(evalPer(iDomain)%julStart, dp)

      ! get calendar days, months, year
      call caldat(int(initTime(iDomain)), yy = year, mm = month, dd = day)

      nMonths = size(tws_obs_mask)

      allocate (month_classes(nMonths))
      allocate (tws_obs_m_anom(nMonths))
      allocate (tws_sim_m_anom(nMonths))

      month_classes(:) = 0
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

      ! calculate standard score
      tws_obs_m_anom = classified_standard_score(tws_opti_catch_avg_domain, month_classes, mask = tws_obs_mask)
      tws_sim_m_anom = classified_standard_score(tws_catch_avg_domain,      month_classes, mask = tws_obs_mask)
      rmse_tws(iDomain) = rmse(tws_sim_m_anom, tws_obs_m_anom, mask = tws_obs_mask)

      deallocate (month_classes)
      deallocate (tws_sim_m_anom)
      deallocate (tws_obs_m_anom)
      deallocate (tws_catch_avg_domain, tws_opti_catch_avg_domain, tws_obs_mask)

      call twsOptiSim(iDomain)%destroy()
    end do

    rmse_tws_avg = sum(rmse_tws(:), abs(rmse_tws - nodata_dp) .gt. eps_dp) / &
            real(count(abs(rmse_tws - nodata_dp) .gt. eps_dp), dp)
    deallocate(rmse_tws)
    deallocate(twsOptiSim)

    !--------------------------------------------
    !! RUNOFF
    !--------------------------------------------
    kge_q_avg = 0_dp
    nGaugesTotal = size(runoff, dim = 2)
    allocate(kge_q(nGaugesTotal))
    kge_q(:) = nodata_dp

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! check for potentially 2 years of data
      pp = count(runoff_agg .ge. 0.0_dp)
      if (pp .lt.  365 * 2) then
        call message('objective_kge_q_rmse_tws: The simulation at gauge ', trim(adjustl(num2str(gg))), &
        ' is not long enough. Please provide at least 730 days of data.')
      end if
      ! calculate KGE for each domain:
      kge_q(gg) = kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      deallocate (runoff_agg, runoff_obs, runoff_obs_mask)

    end do

    ! calculate average KGE value for runoff
    kge_q_avg = sum(kge_q(:), abs(kge_q - nodata_dp) .gt. eps_dp) / &
            real(count(abs(kge_q - nodata_dp) .gt. eps_dp), dp)
    deallocate(kge_q)

    !
    objective_kge_q_rmse_tws = rmse_tws_avg * (1._dp - kge_q_avg)

    call message('    objective_kge_q_rmse_tws = ', num2str(objective_kge_q_rmse_tws, '(F9.5)'))

  END FUNCTION objective_kge_q_rmse_tws

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_neutrons_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for neutrons.

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       neutrons (N) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment CORRELATION coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given domain \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_neutronsdata, L1_neutronsdata_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_neutrons_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine)

  !    HISTORY
  !>       \authors Martin Schroen

  !>       \date Jun 2015

  ! Modifications:
  ! Maren Kaluza Mar 2018 - changed format string to '(I10)'
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_neutrons_kge_catchment_avg(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : domainMeta
    use mo_errormeasures, only : KGE
    use mo_global_variables, only : L1_neutronsObs
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_neutrons_kge_catchment_avg

    ! domain loop counter
    integer(i4) :: iDomain

    ! time loop counter
    integer(i4) :: iTime

    ! for sixth root
#ifndef MPI
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    ! spatial average of observed neutrons
    real(dp), dimension(:), allocatable :: neutrons_catch_avg_domain

    ! spatial avergae of modeled  neutrons
    real(dp), dimension(:), allocatable :: neutrons_opti_catch_avg_domain

    ! simulated neutrons
    type(optidata_sim), dimension(:), allocatable :: neutronsOptiSim

    ! mask for valid neutrons catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    allocate(neutronsOptiSim(domainMeta%nDomains))
    call eval(parameterset, neutronsOptiSim = neutronsOptiSim)

    ! initialize some variables
    objective_neutrons_kge_catchment_avg = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! allocate
      allocate(mask_times                    (size(neutronsOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(neutrons_catch_avg_domain     (size(neutronsOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(neutrons_opti_catch_avg_domain(size(neutronsOptiSim(iDomain)%dataSim, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      neutrons_catch_avg_domain = nodata_dp
      neutrons_opti_catch_avg_domain = nodata_dp

      ! calculate catchment average soil moisture
      do iTime = 1, size(neutronsOptiSim(iDomain)%dataSim, dim = 2)

        ! check for enough data points in time for correlation
        if (all(.NOT. L1_neutronsObs(iDomain)%maskObs(:, iTime))) then
          call message('WARNING: neutrons data at time ', num2str(iTime, '(I10)'), ' is empty.')
          !call message('WARNING: objective_neutrons_kge_catchment_avg: ignored current time step since less than')
          !call message('         10 valid cells available in soil moisture observation')
          mask_times(iTime) = .FALSE.
          cycle
        end if
        neutrons_catch_avg_domain(iTime) = average(L1_neutronsObs(iDomain)%dataObs(:, iTime), &
                                            mask = L1_neutronsObs(iDomain)%maskObs(:, iTime))
        neutrons_opti_catch_avg_domain(iTime) = average(neutronsOptiSim(iDomain)%dataSim(:, iTime), &
                                         mask = L1_neutronsObs(iDomain)%maskObs(:, iTime))
      end do

      ! calculate average neutrons KGE over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg + &
        ((1.0_dp - KGE(neutrons_catch_avg_domain, neutrons_opti_catch_avg_domain, mask = mask_times)) / &
                                                                      real(domainMeta%overallNumberOfDomains, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(neutrons_catch_avg_domain)
      deallocate(neutrons_opti_catch_avg_domain)

      call neutronsOptiSim(iDomain)%destroy()
    end do
    deallocate(neutronsOptiSim)

#ifndef MPI
    objective_neutrons_kge_catchment_avg = objective_neutrons_kge_catchment_avg**onesixth

    call message('    objective_neutrons_kge_catchment_avg = ', num2str(objective_neutrons_kge_catchment_avg, '(F9.5)'))
#endif

  END FUNCTION objective_neutrons_kge_catchment_avg

  ! ------------------------------------------------------------------

  !    NAME
  !        objective_et_kge_catchment_avg

  !    PURPOSE
  !>       \brief Objective function for evpotranspirstion (et).

  !>       \details The objective function only depends on a parameter vector.
  !>       The model will be called with that parameter vector and
  !>       the model output is subsequently compared to observed data.

  !>       Therefore, the Kling-Gupta model efficiency \f$ KGE \f$ of the catchment average
  !>       evapotranspiration (et) is calculated
  !>       \f[ KGE = 1.0 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>       where
  !>       \f$ r \f$ = Pearson product-moment correlation coefficient
  !>       \f$ \alpha \f$ = ratio of simulated mean to observed mean SM
  !>       \f$ \beta  \f$ = ratio of similated standard deviation to observed standard deviation
  !>       is calculated and the objective function for a given domain \f$ i \f$ is
  !>       \f[ \phi_{i} = 1.0 - KGE_{i} \f]
  !>       \f$ \phi_{i} \f$ is the objective since we always apply minimization methods.
  !>       The minimal value of \f$ \phi_{i} \f$ is 0 for the optimal KGE of 1.0.

  !>       Finally, the overall objective function value \f$ OF \f$ is estimated based on the power-6
  !>       norm to combine the \f$ \phi_{i} \f$ from all domains \f$ N \f$.
  !>       \f[ OF = \sqrt[6]{\sum((1.0 - KGE_{i})/N)^6 }.  \f]
  !>       The observed data L1_et, L1_et_mask are global in this module.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_et_kge_catchment_avg &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date Feb 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_et_kge_catchment_avg(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : domainMeta
    use mo_errormeasures, only : KGE
    use mo_moment, only : average
    use mo_string_utils, only : num2str

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_et_kge_catchment_avg

    ! domain loop counter
    integer(i4) ::iDomain

    ! for sixth root
#ifndef MPI
    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp
#endif

    ! spatial average of observed et
    real(dp), dimension(:), allocatable :: et_catch_avg_domain

    ! spatial avergae of modeled  et
    real(dp), dimension(:), allocatable :: et_opti_catch_avg_domain

    !> simulated et
    type(optidata_sim), dimension(:), allocatable :: etOptiSim

    ! mask for valid et catchment avg time steps
    logical, dimension(:), allocatable :: mask_times


    call eval(parameterset, etOptiSim = etOptiSim)

    ! initialize some variables
    objective_et_kge_catchment_avg = 0.0_dp

    ! loop over domain - for applying power law later on
    allocate(etOptiSim(domainMeta%nDomains))
    do iDomain = 1, domainMeta%nDomains

      ! create et array input
      call create_domain_avg_et(iDomain, etOptiSim, et_catch_avg_domain, &
                                           et_opti_catch_avg_domain, mask_times)
      ! calculate average ET KGE over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6

      objective_et_kge_catchment_avg = objective_et_kge_catchment_avg + &
              ((1.0_dp - KGE(et_catch_avg_domain, et_opti_catch_avg_domain, mask = mask_times)) / &
                                        real(domainMeta%overallNumberOfDomains, dp))**6

      ! deallocate
      deallocate(mask_times)
      deallocate(et_catch_avg_domain)
      deallocate(et_opti_catch_avg_domain)
      call etOptiSim(iDomain)%destroy()
    end do
    deallocate(etOptiSim)

#ifndef MPI
    objective_et_kge_catchment_avg = objective_et_kge_catchment_avg**onesixth

    call message('    objective_et_kge_catchment_avg = ', num2str(objective_et_kge_catchment_avg, '(F9.5)'))
#endif

  END FUNCTION objective_et_kge_catchment_avg

  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_sm_corr

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and correlation for SM

  !>       \details Objective function of KGE for runoff and SSE for soil moisture (standarized scores).
  !>       Further details can be found in the documentation of objective functions
  !>       '14 - objective_multiple_gauges_kge_power6' and '13 - objective_sm_corr'.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_sse_sm &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Mar. 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_sm_corr(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_variables, only : level1, domainMeta
    use mo_global_variables, only : L1_smObs
    use mo_moment, only : correlation
    use mo_string_utils, only : num2str
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_sm_corr

    real(dp) :: objective_sm

    real(dp) :: objective_kge

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! domain loop counter
    integer(i4) :: iDomain

    ! cell loop counter
    integer(i4) :: iCell

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! domains wise objectives
    real(dp) :: objective_sm_domain

    ! simulated soil moisture
    type(optidata_sim), dimension(:), allocatable :: smOptiSim

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! run mHM
    allocate(smOptiSim(domainMeta%nDomains))
    call eval(parameterset, runoff = runoff, smOptiSim = smOptiSim)

    ! -----------------------------
    ! SOIL MOISTURE
    ! -----------------------------

    ! initialize some variables
    objective_sm = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! init
      objective_sm_domain = 0.0_dp
      ! get domain information
      nCells1 = level1(iDomain)%nCells

      ! correlation signal is calculated on individual grid cells
      invalid_cells = 0.0_dp
      do iCell = 1, size(L1_smObs(iDomain)%maskObs(:, :), dim = 1)

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_smObs(iDomain)%maskObs(iCell, :)) .LE. (0.10_dp * real(size(L1_smObs(iDomain)%dataObs, dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if

        ! calculate ojective function
        objective_sm_domain = objective_sm_domain + &
                correlation(L1_smObs(iDomain)%dataObs(iCell, :), smOptiSim(iDomain)%dataSim(iCell, :), &
                                                           mask = L1_smObs(iDomain)%maskObs(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_sm: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invalid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture objective over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_sm = objective_sm + &
              ((1.0_dp - objective_sm_domain / real(nCells1, dp)) / real(domainMeta%overallNumberOfDomains, dp))**6
      call smOptiSim(iDomain)%destroy()
    end do
    deallocate(smOptiSim)

    ! compromise solution - sixth root
    objective_sm = objective_sm**onesixth

    ! -----------------------------
    ! RUNOFF
    ! -----------------------------
    objective_kge = 0.0_dp
    nGaugesTotal = size(runoff, dim = 2)

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! KGE
      objective_kge = objective_kge + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6

    end do

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

    ! compromise solution - sixth root
    objective_kge = objective_kge**onesixth

    ! equal weighted compromise objective functions for discharge and soilmoisture
    ! ToDo: why do we use the sixth root of of objective_sm and objective_kge
    ! only to take the power to 6 here again, when we never need the
    ! intermediate results?
#ifdef MPI
    objective_kge_q_sm_corr = (objective_sm**6 + objective_kge**6)
#else
    objective_kge_q_sm_corr = (objective_sm**6 + objective_kge**6)**onesixth

    call message('    objective_kge_q_sm_corr = ', num2str(objective_kge_q_sm_corr, '(F9.5)'))
#endif
    !    print*, "1-SM 2-Q : ", 1.0_dp-objective_sm, 1.0_dp-objective_kge ! MZMZMZMZ

  END FUNCTION objective_kge_q_sm_corr


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_et

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and KGE for ET

  !>       \details Objective function of KGE for runoff and KGE for ET.
  !>       Further details can be found in the documentation of objective functions
  !>       '14 - objective_multiple_gauges_kge_power6'.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_et &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date July 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_et(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_variables, only : level1, domainMeta
    use mo_errormeasures, only : kge
    use mo_global_variables, only : L1_etObs
    use mo_string_utils, only : num2str
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_et

    real(dp) :: objective_et

    real(dp) :: objective_q

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! domain loop counter
    integer(i4) :: iDomain

    ! cell loop counter
    integer(i4) :: iCell

    ! ncells1 of level 1
    integer(i4) :: nCells1

    ! domains wise objectives
    real(dp) :: objective_et_domain

    !> simulated et
    type(optidata_sim), dimension(:), allocatable :: etOptiSim

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! run mHM
    allocate(etOptiSim(domainMeta%nDomains))
    call eval(parameterset, runoff = runoff, etOptiSim = etOptiSim)

    ! -----------------------------
    ! EVAPOTRANSPIRATION
    ! -----------------------------

    ! initialize some variables
    objective_et = 0.0_dp

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains

      ! init
      objective_et_domain = 0.0_dp
      ! get domain information
      nCells1 = level1(iDomain)%nCells

      ! correlation signal is calculated on individual grid cells
      invalid_cells = 0.0_dp
      do iCell = 1, size(L1_etObs(iDomain)%maskObs, dim = 1)

        ! check for enough data points in time for statistical calculations (e.g. mean, stddev)
        if (count(L1_etObs(iDomain)%maskObs(iCell, :)) .LE. &
                            (0.10_dp * real(size(L1_etObs(iDomain)%dataObs(:, :), dim = 2), dp))) then
          invalid_cells = invalid_cells + 1.0_dp
          cycle
        end if

        ! calculate ojective function
        objective_et_domain = objective_et_domain + &
                kge(L1_etObs(iDomain)%dataObs(iCell, :), etOptiSim(iDomain)%dataSim(iCell, :), &
                                                   mask = L1_etObs(iDomain)%maskObs(iCell, :))
      end do

      ! user information about invalid cells
      if (invalid_cells .GT. 0.5_dp) then
        call message('   WARNING: objective_et: Detected invalid cells in study area (.LT. 10 valid data points).')
        call message('                          Fraction of invalid cells: ', &
                num2str(invalid_cells / real(nCells1, dp), '(F4.2)'))
      end if

      ! calculate average soil moisture objective over all domains with power law
      ! domains are weighted equally ( 1 / real(domainMeta%overallNumberOfDomains,dp))**6
      objective_et = objective_et + &
              ((1.0_dp - objective_et_domain / real(nCells1, dp)) / real(domainMeta%overallNumberOfDomains, dp))**6
      call etOptiSim(iDomain)%destroy()
    end do
    deallocate(etOptiSim)

    ! compromise solution - sixth root
    objective_et = objective_et**onesixth

    ! -----------------------------
    ! RUNOFF
    ! -----------------------------
    objective_q = 0.0_dp
    nGaugesTotal = size(runoff, dim = 2)

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! KGE
      objective_q = objective_q + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6

    end do

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

    ! compromise solution - sixth root
    objective_q = objective_q**onesixth

    ! equal weighted compromise objective functions for discharge and soilmoisture
    ! ToDo: why do we use the sixth root of of objective_sm and objective_kge
    ! only to take the power to 6 here again, when we never need the
    ! intermediate results?
#ifdef MPI
    objective_kge_q_et = (objective_et**6 + objective_q**6)
#else
    objective_kge_q_et = (objective_et**6 + objective_q**6)**onesixth

    call message('    objective_kge_q_et = ', num2str(objective_kge_q_et, '(F9.5)'))
#endif
    !    print*, "1-SM 2-Q : ", 1.0_dp-objective_sm, 1.0_dp-objective_kge ! MZMZMZMZ

  END FUNCTION objective_kge_q_et


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_BFI

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and BFI absulute difference

  !>       \details Objective function of KGE for runoff and KGE for ET.
  !>       Further details can be found in the documentation of objective functions
  !>       '14 - objective_multiple_gauges_kge_power6'.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_BFI &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Sebastian Müller

  !>       \date Apr 2022

  FUNCTION objective_kge_q_BFI(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_variables, only : level1, domainMeta
    use mo_errormeasures, only : kge
    use mo_global_variables, only : BFI_obs
    use mo_string_utils, only : num2str
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_BFI
    real(dp) :: objective_BFI
    real(dp) :: objective_q

    real(dp), parameter :: onesixth = 1.0_dp / 6.0_dp

    ! number of invalid cells in catchment
    real(dp) :: invalid_cells

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    ! domain loop counter
    integer(i4) :: iDomain

    !> baseflow index for each domain
    real(dp), dimension(:), allocatable :: BFI

    ! counter
    integer(i4) :: gg, i

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    integer(i4), dimension(:), allocatable :: domain_ids, domain_ids_pack

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! run mHM
    allocate(BFI(domainMeta%nDomains))
    call eval(parameterset, runoff = runoff, BFI = BFI)

    ! -----------------------------
    ! BFI
    ! -----------------------------

    ! initialize some variables
    objective_BFI = 0.0_dp

    if ( any(BFI_obs < 0.0_dp) ) then
      allocate(domain_ids(domainMeta%nDomains))
      allocate(domain_ids_pack(count(BFI_obs < 0.0_dp)))
      domain_ids = [(i, i=1,size(domain_ids))]
      domain_ids_pack = pack(domain_ids, mask=(BFI_obs < 0.0_dp))
      call error_message( &
        "objective_kge_q_BFI: missing BFI values for domain ", &
        trim(adjustl(num2str(domain_ids_pack(1)))) &
      )
    end if

    ! loop over domain - for applying power law later on
    do iDomain = 1, domainMeta%nDomains
      objective_BFI = objective_BFI + abs(BFI(iDomain) - BFI_obs(iDomain)) / domainMeta%nDomains
    end do

    ! -----------------------------
    ! RUNOFF
    ! -----------------------------
    objective_q = 0.0_dp
    nGaugesTotal = size(runoff, dim = 2)

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! KGE
      objective_q = objective_q + &
              ((1.0_dp - kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)) / real(nGaugesTotal, dp))**6

    end do

    deallocate(runoff_agg, runoff_obs, runoff_obs_mask)

    ! compromise solution - sixth root
    objective_q = objective_q**onesixth

    objective_kge_q_BFI = (objective_BFI + 1._dp)*objective_q
    call message('    objective_kge_q_BFI = ', num2str(objective_kge_q_BFI, '(F9.5)'))

  END FUNCTION objective_kge_q_BFI


  ! -----------------------------------------------------------------

  !    NAME
  !        objective_kge_q_rmse_et

  !    PURPOSE
  !>       \brief Objective function of KGE for runoff and RMSE for domain_avg ET (standarized scores)

  !>       \details Objective function of KGE for runoff and RMSE for domain_avg ET (standarized scores)

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset"
  !>       \param[in] "procedure(eval_interface) :: eval"

  !    RETURN
  !>       \return real(dp) :: objective_kge_q_rmse_et &mdash; objective function value
  !>       (which will be e.g. minimized by an optimization routine like DDS)

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date July 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  FUNCTION objective_kge_q_rmse_et(parameterset, eval)

    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : eps_dp, nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : domainMeta
    use mo_errormeasures, only : rmse
    use mo_global_variables, only : L1_etObs
    use mo_julian, only : caldat
    use mo_moment, only : average, mean
    use mo_standard_score, only : classified_standard_score
    use mo_string_utils, only : num2str
    use mo_temporal_aggregation, only : day2mon_average
    use mo_errormeasures, only : kge
    use mo_mrm_objective_function_runoff, only : extract_runoff

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset

    procedure(eval_interface), INTENT(IN), POINTER :: eval

    real(dp) :: objective_kge_q_rmse_et

    ! modelled runoff for a given parameter set
    ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:, :) :: runoff

    !> simulated et
    type(optidata_sim), dimension(:), allocatable :: etOptiSim

    ! time loop counter
    integer(i4) :: iTime

    ! domain counter, month counters
    integer(i4) :: iDomain, pp, mmm

    integer(i4) :: year, month, day

    real(dp), dimension(domainMeta%nDomains) :: initTime

    ! total number of months
    integer(i4) :: nMonths

    ! vector with months' classes
    integer(i4), dimension(:), allocatable :: month_classes

    ! monthly values original time series
    real(dp), dimension(:), allocatable :: et_sim_m, et_obs_m

    ! monthly values anomaly time series
    real(dp), dimension(:), allocatable :: et_sim_m_anom, et_obs_m_anom

    logical, dimension(:), allocatable :: et_obs_m_mask

    ! kge_q(nGaugesTotal)
    real(dp), dimension(:), allocatable :: rmse_et

    ! obj. functions
    real(dp) :: rmse_et_avg, kge_q_avg

    ! spatial average of observed et
    real(dp), dimension(:), allocatable :: et_catch_avg_domain

    ! spatial avergae of modeled  et
    real(dp), dimension(:), allocatable :: et_opti_catch_avg_domain

    ! mask for valid et catchment avg time steps
    logical, dimension(:), allocatable :: mask_times

    ! rmse_et(domainMeta%nDomains)
    real(dp), dimension(:), allocatable :: kge_q

    ! gauges counter
    integer(i4) :: gg

    integer(i4) :: nGaugesTotal

    ! aggregated simulated runoff
    real(dp), dimension(:), allocatable :: runoff_agg

    ! measured runoff
    real(dp), dimension(:), allocatable :: runoff_obs

    ! mask for measured runoff
    logical, dimension(:), allocatable :: runoff_obs_mask

    ! obtain simulation values of runoff (hourly) and ET
    ! for ET only valid cells (domains concatenated)
    ! et_opti: aggregate ET to needed time step for optimization (see timeStep_et_input)
    allocate(etOptiSim(domainMeta%nDomains))
    call eval(parameterset, runoff = runoff, etOptiSim = etOptiSim)

    !--------------------------------------------
    !! EVAPOTRANSPIRATION
    !--------------------------------------------

    ! allocate
    allocate(rmse_et(domainMeta%nDomains))
    rmse_et(:) = nodata_dp

    do iDomain = 1, domainMeta%nDomains

      ! allocate
      allocate(mask_times              (size(etOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(et_catch_avg_domain     (size(etOptiSim(iDomain)%dataSim, dim = 2)))
      allocate(et_opti_catch_avg_domain(size(etOptiSim(iDomain)%dataSim, dim = 2)))

      ! initalize
      mask_times = .TRUE.
      et_catch_avg_domain = nodata_dp
      et_opti_catch_avg_domain = nodata_dp

      ! calculate catchment average evapotranspiration
      do iTime = 1, size(etOptiSim(iDomain)%dataSim, dim = 2)
        ! check for enough data points in time for correlation
        if (all(.NOT. L1_etObs(iDomain)%maskObs(:, iTime))) then
          !write (*,*) 'WARNING: et data at time ', iTime, ' is empty.'
          !call message('WARNING: objective_et_kge_catchment_avg: ignored current time step since less than')
          !call message('         10 valid cells available in evapotranspiration observation')
          mask_times(iTime) = .FALSE.
          cycle
        end if
        ! spatial average of observed ET
        et_catch_avg_domain(iTime) = average(L1_etObs(iDomain)%dataObs(:, iTime), &
                                      mask = L1_etObs(iDomain)%maskObs(:, iTime))
        ! spatial avergae of modeled ET
        et_opti_catch_avg_domain(iTime) = average(etOptiSim(iDomain)%dataSim(:, iTime), &
                                            mask = L1_etObs(iDomain)%maskObs(:, iTime))
      end do

      ! get initial time of the evaluation period
      initTime(iDomain) = real(evalPer(iDomain)%julStart, dp)

      ! get calendar days, months, year
      call caldat(int(initTime(iDomain)), yy = year, mm = month, dd = day)

      ! if evapotranspiration input daily
      select case(L1_etObs(iDomain)%timeStepInput)
        ! JBJBJB
        ! daily: aggregate to monthly mean
      case(-1)
        ! calculate monthly averages from daily values of the model
        call day2mon_average(et_opti_catch_avg_domain, year, month, day, et_sim_m, misval = nodata_dp)
        ! calculate monthly averages from daily values of the observations
        call day2mon_average(et_catch_avg_domain, year, month, day, et_obs_m, misval = nodata_dp)
        !
        ! monthly: proceed without action
      case(-2)
        ! simulation
        allocate(et_sim_m(size(etOptiSim(iDomain)%dataSim, dim = 2)))
        et_sim_m = et_opti_catch_avg_domain
        ! observation
        allocate(et_obs_m(size(etOptiSim(iDomain)%dataSim, dim = 2)))
        et_obs_m = et_catch_avg_domain

        ! yearly: ERROR stop program
      case(-3)
        call error_message('***ERROR: objective_kge_q_rmse_et: time step of evapotranspiration yearly.')
      end select
      ! remove mean from modelled time series
      et_sim_m(:) = et_sim_m(:) - mean(et_sim_m(:))
      ! remove mean from observed time series
      et_obs_m(:) = et_obs_m(:) - mean(et_obs_m(:))
      ! get number of months for given domain
      nMonths = size(et_obs_m)
      allocate (month_classes(nMonths))
      allocate (et_obs_m_mask(nMonths))
      allocate (et_obs_m_anom(nMonths))
      allocate (et_sim_m_anom(nMonths))

      month_classes(:) = 0
      et_obs_m_mask(:) = .TRUE.
      et_obs_m_anom(:) = nodata_dp
      et_sim_m_anom(:) = nodata_dp

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
      ! double check missing data
      ! define mask for missing data in observations (there are always data for simulations)
      where(abs(et_obs_m - nodata_dp) .lt. eps_dp) et_obs_m_mask = .FALSE.

      ! calculate standard score
      et_obs_m_anom = classified_standard_score(et_obs_m, month_classes, mask = et_obs_m_mask)
      et_sim_m_anom = classified_standard_score(et_sim_m, month_classes, mask = et_obs_m_mask)
      rmse_et(iDomain) = rmse(et_sim_m_anom, et_obs_m_anom, mask = et_obs_m_mask)

      deallocate (month_classes)
      deallocate (et_obs_m)
      deallocate (et_sim_m)
      deallocate (et_obs_m_mask)
      deallocate (et_sim_m_anom)
      deallocate (et_obs_m_anom)
      !end if

      call etOptiSim(iDomain)%destroy()
    end do

    rmse_et_avg = sum(rmse_et(:), abs(rmse_et - nodata_dp) .gt. eps_dp) / &
            real(count(abs(rmse_et - nodata_dp) .gt. eps_dp), dp)
    deallocate(rmse_et)
    deallocate(etOptiSim)

    !--------------------------------------------
    !! RUNOFF
    !--------------------------------------------
    kge_q_avg = 0_dp
    nGaugesTotal = size(runoff, dim = 2)
    allocate(kge_q(nGaugesTotal))
    kge_q(:) = nodata_dp

    do gg = 1, nGaugesTotal

      ! extract runoff
      call extract_runoff(gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask)

      ! check for whether to proceed with this domain or not
      ! potentially 3 years of data
      !pp = count(runoff_agg .ge. 0.0_dp )
      !if (pp .lt.  365*3 ) then
      !    deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
      !    cycle
      ! else
      ! calculate KGE for each domain:
      kge_q(gg) = kge(runoff_obs, runoff_agg, mask = runoff_obs_mask)
      deallocate (runoff_agg, runoff_obs, runoff_obs_mask)
      ! end if

    end do

    ! calculate average KGE value for runoff
    kge_q_avg = sum(kge_q(:), abs(kge_q - nodata_dp) .gt. eps_dp) / &
            real(count(abs(kge_q - nodata_dp) .gt. eps_dp), dp)
    deallocate(kge_q)

    !
    objective_kge_q_rmse_et = rmse_et_avg * (1._dp - kge_q_avg)

    call message('    objective_kge_q_rmse_et = ', num2str(objective_kge_q_rmse_et, '(F9.5)'))

  END FUNCTION objective_kge_q_rmse_et

  subroutine create_domain_avg_tws(iDomain, twsOptiSim, tws_catch_avg_domain, &
                                           tws_opti_catch_avg_domain, mask_times)
    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_global_variables, only : L1_twsaObs
    use mo_moment, only : average
    ! current domain Id
    integer(i4), intent(in) :: iDomain

    ! simulated tws
    type(optidata_sim), dimension(:), intent(in) :: twsOptiSim

    ! aggregated simulated
    real(dp), dimension(:), allocatable, intent(out) :: tws_catch_avg_domain

    ! extracted measured
    real(dp), dimension(:), allocatable, intent(out) :: tws_opti_catch_avg_domain

    ! mask of no data values
    logical, dimension(:), allocatable, intent(out) :: mask_times

    ! local
    ! time loop counter
    integer(i4) :: iTime

    ! allocate
    allocate(mask_times               (size(twsOptiSim(iDomain)%dataSim, dim = 2)))
    allocate(tws_catch_avg_domain     (size(twsOptiSim(iDomain)%dataSim, dim = 2)))
    allocate(tws_opti_catch_avg_domain(size(twsOptiSim(iDomain)%dataSim, dim = 2)))

    ! initalize
    mask_times = .TRUE.
    tws_catch_avg_domain = nodata_dp
    tws_opti_catch_avg_domain = nodata_dp

    ! calculate catchment average evapotranspiration
    do iTime = 1, size(twsOptiSim(iDomain)%dataSim, dim = 2)

      ! check for enough data points in time for correlation
      if (all(.NOT. L1_twsaObs(iDomain)%maskObs(:, iTime))) then
        !write (*,*) 'WARNING: et data at time ', iTime, ' is empty.'
        !call message('WARNING: objective_et_kge_catchment_avg: ignored current time step since less than')
        !call message('         10 valid cells available in evapotranspiration observation')
        mask_times(iTime) = .FALSE.
        cycle
      end if

      tws_catch_avg_domain(iTime) = average(L1_twsaObs(iDomain)%dataObs(:, iTime), &
                                     mask = L1_twsaObs(iDomain)%maskObs(:, iTime))
      tws_opti_catch_avg_domain(iTime) = average(twsOptiSim(iDomain)%dataSim(:, iTime), &
                                     mask = L1_twsaObs(iDomain)%maskObs(:, iTime))
    end do

  end subroutine create_domain_avg_tws

  subroutine create_domain_avg_et(iDomain, etOptiSim, et_catch_avg_domain, &
                                           et_opti_catch_avg_domain, mask_times)
    use mo_optimization_types, only : optidata_sim
    use mo_common_constants, only : nodata_dp
    use mo_global_variables, only : L1_etObs
    use mo_moment, only : average
    ! current domain Id
    integer(i4), intent(in) :: iDomain

    ! simulated et
    type(optidata_sim), dimension(:), intent(in) :: etOptiSim

    ! aggregated simulated
    real(dp), dimension(:), allocatable, intent(out) :: et_catch_avg_domain

    ! extracted measured
    real(dp), dimension(:), allocatable, intent(out) :: et_opti_catch_avg_domain

    ! mask of no data values
    logical, dimension(:), allocatable, intent(out) :: mask_times

    ! local
    ! time loop counter
    integer(i4) :: iTime

    ! allocate
    allocate(mask_times              (size(etOptiSim(iDomain)%dataSim, dim = 2)))
    allocate(et_catch_avg_domain     (size(etOptiSim(iDomain)%dataSim, dim = 2)))
    allocate(et_opti_catch_avg_domain(size(etOptiSim(iDomain)%dataSim, dim = 2)))

    ! initalize
    mask_times = .TRUE.
    et_catch_avg_domain = nodata_dp
    et_opti_catch_avg_domain = nodata_dp

    ! calculate catchment average evapotranspiration
    do iTime = 1, size(etOptiSim(iDomain)%dataSim, dim = 2)

      ! check for enough data points in time for correlation
      if (all(.NOT. L1_etObs(iDomain)%maskObs(:, iTime))) then
        !write (*,*) 'WARNING: et data at time ', iTime, ' is empty.'
        !call message('WARNING: objective_et_kge_catchment_avg: ignored current time step since less than')
        !call message('         10 valid cells available in evapotranspiration observation')
        mask_times(iTime) = .FALSE.
        cycle
      end if

      et_catch_avg_domain(iTime) = average(L1_etObs(iDomain)%dataObs(:, iTime), &
                                    mask = L1_etObs(iDomain)%maskObs(:, iTime))
      et_opti_catch_avg_domain(iTime) = average(etOptiSim(iDomain)%dataSim(:, iTime), &
                                          mask = L1_etObs(iDomain)%maskObs(:, iTime))
    end do

  end subroutine create_domain_avg_et

  subroutine convert_tws_to_twsa(twsOptiSim, L1_twsaObs, twsaOptiSim)
    use mo_optimization_types, only : optidata_sim, optidata
    use mo_moment, only : average
    ! simulated tws
    type(optidata_sim), intent(in)    :: twsOptiSim
    ! observed twsa
    type(optidata),     intent(in)    :: L1_twsaObs
    ! simulated twsa
    type(optidata_sim), intent(inout) :: twsaOptiSim

    ! local
    integer(i4) :: iCell
    real(dp)    :: twsa_av_cell

    allocate(twsaOptiSim%dataSim(size(twsOptiSim%dataSim(:, :), dim = 1), size(twsOptiSim%dataSim(:, :), dim = 2)))

    do iCell = 1, size(twsOptiSim%dataSim(:, :), dim = 1)
      twsa_av_cell = average(twsOptiSim%dataSim(iCell, :), mask = L1_twsaObs%maskObs(iCell, :))
      twsaOptiSim%dataSim(iCell, :) = twsOptiSim%dataSim(iCell, :) - twsa_av_cell
    end do

  end subroutine convert_tws_to_twsa

END MODULE mo_objective_function
