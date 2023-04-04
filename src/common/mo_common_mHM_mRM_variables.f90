!> \file mo_common_mHM_mRM_variables.f90
!> \brief \copybrief mo_common_mhm_mrm_variables
!> \details \copydetails mo_common_mhm_mrm_variables

!> \brief Provides structures needed by mHM, mRM and/or mpr.
!> \details Provides the global structure period that is used by both mHM and mRM.
!> \changelog
!! - Stephan Thober  Nov 2016
!!   - moved processdescription from mo_global_variables to here
!! - Robert Schweppe Dec 2017
!!   - merged more duplicated variables from mhm and mrm global variables
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!> \authors Stephan Thober
!> \date Sep 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_mHM_mRM_variables

  use mo_kind, only : i4, i8, dp
  use mo_common_types, only: period
  implicit none

  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of main part
  ! -------------------------------------------------------------------
  integer(i4) :: mrm_coupling_mode !-1 = no mrm (mHM only)
  !                                ! 0 = stand alone (mRM only)
  !                                ! 1 = general coupling to a model (not used)
  !                                ! 2 = specific coupling to mHM (mHM and mRM)

  integer(i4), public :: timeStep                   ! [h] simulation time step (= TS) in [h]
  real(dp), public :: c2TSTu            !       Unit transformation = timeStep/24
  real(dp), dimension(:), allocatable, public :: resolutionRouting          ! [m or degree] resolution of routing - Level 11
  logical, public :: read_restart               ! flag
  logical, public :: mrm_read_river_network     ! flag
  logical, public :: read_old_style_restart_bounds     ! flag
  logical, public :: restart_reset_fluxes_states     !< flag to reset fluxes and states read from restart to default values


  type(period), dimension(:), allocatable, public :: warmPer     ! time period for warming
  type(period), dimension(:), allocatable, public :: evalPer     ! time period for model evaluation
  type(period), dimension(:), allocatable, public :: simPer      ! warmPer + evalPer
  integer(i4), dimension(:), allocatable, public :: warmingDays ! number of days for warm up period
  integer(i4), dimension(:, :), allocatable, public :: LCyearId            ! Mapping of landcover scenes (1, 2,..) for each domain

  ! ------------------------------------------------------------------
  ! CONSTANT
  ! ------------------------------------------------------------------
  integer(i4), public :: nTstepDay          !       Number of time intervals per day

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nDomains
  character(256), dimension(:), allocatable, public :: mhmFileRestartIn! Directory where input of restart is read from
  character(256), dimension(:), allocatable, public :: mrmFileRestartIn! Directory where input of restart is read from

  ! -------------------------------------------------------------------
  ! OPTIMIZATION
  ! -------------------------------------------------------------------
  integer(i4), public :: opti_method         ! Optimization algorithm:
  !                                                                       ! 1 - DDS
  !                                                                       ! 2 - Simulated Annealing
  !                                                                       ! 3 - SCE
  integer(i4), public :: opti_function       ! Objective function:
  !                                                                       ! 1 - 1.0-NSE
  !                                                                       ! 2 - 1.0-lnNSE
  !                                                                       ! 3 - 1.0-0.5*(NSE+lnNSE)
  logical, public :: optimize            ! Optimization   (.true. ) or
  !                                                                       ! Evaluation run (.false.)
  logical, public :: optimize_restart    ! Optimization will be restarted from
  !                                                                       ! mo_<opti_method>.restart file (.true.)
  ! settings for optimization algorithms:
  integer(i8), public :: seed                ! seed used for optimization
  !                                                                       ! default: -9 --> system time
  integer(i4), public :: nIterations         ! number of iterations for optimization
  real(dp), public :: dds_r               ! DDS: perturbation rate
  !                                                                       !      default: 0.2
  real(dp), public :: sa_temp             ! SA:  initial temperature
  !                                                                       !      default: -9.0 --> estimated
  integer(i4), public :: sce_ngs             ! SCE: # of complexes
  !                                                                       !      default: 2
  integer(i4), public :: sce_npg             ! SCE: # of points per complex
  !                                                                       !      default: -9 --> 2n+1
  integer(i4), public :: sce_nps             ! SCE: # of points per subcomplex
  !                                                                       !      default: -9 --> n+1
  logical, public :: mcmc_opti           ! MCMC: Optimization (.true. ) or
  !                                                                       !       Only parameter uncertainty (.false.)
  integer(i4), public, parameter :: nerror_model = 2    ! # possible parameters in error model
  !                                                                       !       e.g. for opti_function=8: 2
  real(dp), public, dimension(nerror_model) :: mcmc_error_params   !       Parameters of error model if mcmc_opti=.false.
  !                                                                       !       e.g. for opti_function=8: 0.01, 0.3

end module mo_common_mHM_mRM_variables
