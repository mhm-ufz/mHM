!> \file mo_common_variables.f90

!> \brief Provides structures needed by mHM and mRM.

!> \details Provides the global structure period that is used
!>     by both mHM and mRM.

!> \author Stephan Thober
!> \date Sep 2015
!  Modified Stephan Thober, Nov 2016 - moved processdescription from mo_global_variables to here
module mo_common_variables
  use mo_kind, only: i4, i8, dp
  implicit none

  ! -------------------------------------------------------------------
  ! PROCESSES description
  ! -------------------------------------------------------------------
  integer(i4), parameter,                public :: nProcesses = 10 ! Number of possible processes to consider
  !                                                                !   process 1 :: interception
  !                                                                !   process 2 :: snow
  !                                                                !   process 3 :: soilmoisture
  !                                                                !   process 4 :: sealed area direct runoff
  !                                                                !   process 5 :: potential evapotranspiration
  !                                                                !   process 6 :: interflow
  !                                                                !   process 7 :: percolation
  !                                                                !   process 8 :: routing
  !                                                                !   process 9 :: baseflow
  !                                                                !   process 10:: neutrons
  integer(i4), dimension(nProcesses, 3), public :: processMatrix   ! Info about which process runs in which option and
  !                                                                ! number of parameters necessary for this option
  !                                                                !   col1: process_switch
  !                                                                !   col2: no. of parameters
  !                                                                !   col3: cum. no. of parameters

  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type period
      integer(i4) :: dStart      ! first day
      integer(i4) :: mStart      ! first month
      integer(i4) :: yStart      ! first year
      integer(i4) :: dEnd        ! last  day
      integer(i4) :: mEnd        ! last  month
      integer(i4) :: yEnd        ! last  year
      integer(i4) :: julStart    ! first julian day 
      integer(i4) :: julEnd      ! last  julian day 
      integer(i4) :: nObs        ! total number of observations
   end type period
   ! -------------------------------------------------------------------
   ! OPTIMIZATION
   ! -------------------------------------------------------------------
   integer(i4), public                              :: opti_method         ! Optimization algorithm:
   !                                                                       ! 1 - DDS
   !                                                                       ! 2 - Simulated Annealing
   !                                                                       ! 3 - SCE
   integer(i4), public                              :: opti_function       ! Objective function:
   !                                                                       ! 1 - 1.0-NSE
   !                                                                       ! 2 - 1.0-lnNSE
   !                                                                       ! 3 - 1.0-0.5*(NSE+lnNSE)
   logical,     public                              :: optimize            ! Optimization   (.true. ) or
   !                                                                       ! Evaluation run (.false.)
   logical,     public                              :: optimize_restart    ! Optimization will be restarted from
   !                                                                       ! mo_<opti_method>.restart file (.true.)
   ! settings for optimization algorithms: 
   integer(i8), public                              :: seed                ! seed used for optimization
   !                                                                       ! default: -9 --> system time 
   integer(i4), public                              :: nIterations         ! number of iterations for optimization
   real(dp),    public                              :: dds_r               ! DDS: perturbation rate
   !                                                                       !      default: 0.2
   real(dp),    public                              :: sa_temp             ! SA:  initial temperature
   !                                                                       !      default: -9.0 --> estimated
   integer(i4), public                              :: sce_ngs             ! SCE: # of complexes
   !                                                                       !      default: 2
   integer(i4), public                              :: sce_npg             ! SCE: # of points per complex
   !                                                                       !      default: -9 --> 2n+1
   integer(i4), public                              :: sce_nps             ! SCE: # of points per subcomplex
   !                                                                       !      default: -9 --> n+1
   logical,     public                              :: mcmc_opti           ! MCMC: Optimization (.true. ) or
   !                                                                       !       Only parameter uncertainty (.false.)
   integer(i4), public, parameter                   :: nerror_model = 2    !       # possible parameters in error model
   !                                                                       !       e.g. for opti_function=8: 2
   real(dp),    public, dimension(nerror_model)     :: mcmc_error_params   !       Parameters of error model if mcmc_opti=.false.
   !                                                                       !       e.g. for opti_function=8: 0.01, 0.3
   
   ! -------------------------------------------------------------------
   ! PARAMETERS
   ! -------------------------------------------------------------------
   real(dp),       dimension(:,:), allocatable, public :: global_parameters      ! Matrix of global parameters (former: gamma)
   !                                                                             !   col1: min,  col2: max, col3: initial, 
   !                                                                             !   col4: flag, col5: scaling
   character(256), dimension(:), allocatable,   public :: global_parameters_name ! Matrix of global parameters (former: gamma)
   !                                                                             !   col1: names
   ! -------------------------------------------------------------------
   ! ALMA convention
   ! -------------------------------------------------------------------
   logical :: ALMA_convention ! flag for ALMA convention
   !                          ! see http://www.lmd.jussieu.fr/~polcher/ALMA/convention_3.html
   !                          ! .True.: ALMA convention is used for Input/Output
   !                          ! .False.: default mHM units are used
   !                          ! CAUTION: only Qall is considered at the moment
   
end module mo_common_variables
