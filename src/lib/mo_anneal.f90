MODULE mo_anneal

  ! This module is minimizing a cost function via Simulated Annealing
  ! and is part of the UFZ CHS Fortran library.


  ! Written  Juliane Mai, Mar 2012 : module implementation
  ! Modified Juliane Mai, May 2012 : anneal: sp version
  !          Juliane Mai, May 2012 : anneal: documentation
  !          Juliane Mai, May 2012 : GetTemperature: sp and dp version
  !          Juliane Mai, Jun 2012 : weighted parameter selection
  !          Juliane Mai, Aug 2012 : - function anneal instead of subroutine
  !                                  - using new module get_timeseed as default for seeding
  !                                  - new optional for minimization or maximization
  !                                  - fixed parameter ranges possible instead of interface range
  !          Juliane Mai, Nov 2012 : history of achieved objective function values as optional out
  !                                  only in anneal but not anneal_valid
  !          Juliane Mai, Jan 2013 : - including DDS features in anneal, i.e. reflection at parameter boundaries,
  !                                    different parameter selection modes (one, all, neighborhood), and
  !                                    different parameter pertubation modes (flexible r=dR (anneal version) or
  !                                    constant r=0.2 (dds version))
  !                                  - remove sp versions
  !                                  - fixed and flexible parameter ranges are now in one function
  !                                    using optional arguments
  !                                  - undef_funcval instead of anneal_valid function
  !          Juliane Mai, Feb 2013 : - xor4096 optionals combined in save_state

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Juliane Mai.

  ! Copyright 2012-13 Juliane Mai

  USE mo_kind,    ONLY: i4, i8, dp
  USE mo_utils,   ONLY: le, ge
  USE mo_xor4096, ONLY: get_timeseed, xor4096, xor4096g, n_save_state

  IMPLICIT NONE

  PUBLIC :: anneal                  ! Minimization of a cost function via Simaulated Annealing
  PUBLIC :: GetTemperature          ! Function which returns the optimal initial temperature for
  !                                 ! a given acceptance ratio and initial parameter set

  ! ------------------------------------------------------------------

  !     NAME
  !>        \brief anneal

  !     PURPOSE
  !>        \details Optimizes a user provided cost function using the Simulated Annealing strategy.


  !     CALLING SEQUENCE
  !         para = (/ 1.0_dp , 2.0_dp /)
  !         prange(1,:) = (/ 0.0_dp, 10.0_dp /)
  !         prange(2,:) = (/ 0.0_dp, 50.0_dp /)
  !
  !         user defined function cost_dp which calculates the cost function value for a
  !         parameter set (the interface given below has to be used for this function!)

  !         parabest = anneal(cost_dp, para, prange)

  !         see also test folder for a detailed example

  !     INTENT(IN)
  !>        \param[in]  "INTERFACE                       :: cost_dp"         interface calculating the
  !>                                                                         cost function at a given point

  !>        \param[in]  "REAL(DP),    DIMENSION(:)       :: para"            initial parameter set

  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in]  "REAL(DP), DIMENSION(size(para),2), optional :: prange"
  !>                                                                         lower and upper bound per parameter
  !             OR
  !>        \param[in]  "INTERFACE, optional             :: prange_func"     interface calculating the
  !>                                                                         feasible range for a parameter
  !>                                                                         at a certain point, if ranges are variable

  !>        \param[in]  "REAL(DP), optional              :: temp"            initial temperature \n
  !>                                                                         DEFAULT: Get_Temperature

  !>        \param[in]  "REAL(DP), optional              :: DT"              geometrical decreement of temperature \n
  !>                                                                         0.7<DT<0.999 \n
  !>                                                                         DEFAULT: 0.9_dp

  !>        \param[in]  "INTEGER(I4), optional           :: nITERmax"        maximal number of iterations \n
  !>                                                                         will be increased by 10% if stopping criteria of
  !>                                                                         acceptance ratio or epsilon decreement of cost
  !>                                                                         function is not fullfilled \n
  !>                                                                         DEFAULT: 1000_i4

  !>        \param[in]  "INTEGER(I4), optional           :: LEN"             Length of Markov Chain \n
  !>                                                                         DEFAULT: MAX(250_i4, size(para,1))

  !>        \param[in]  "INTEGER(I4), optional           :: nST"             Number of consecutive LEN steps \n
  !>                                                                         DEFAULT: 5_i4

  !>        \param[in]  "REAL(DP), optional              :: eps"             Stopping criteria of epsilon decreement of
  !>                                                                         cost function \n
  !>                                                                         DEFAULT: 0.0001_dp

  !>        \param[in]  "REAL(DP), optional              :: acc"             Stopping criteria for Acceptance Ratio \n
  !>                                                                         acc    <= 0.1_dp \n
  !>                                                                         DEFAULT: 0.1_dp

  !>        \param[in]  "INTEGER(I4/I8), DIMENSION(3), optional      :: seeds" 
  !>                                                                         Seeds of random numbers used for random parameter
  !>                                                                         set generation \n
  !>                                                                         DEFAULT: dependent on current time

  !>        \param[in]  "LOGICAL, optional                           :: printflag" 
  !>                                                                         If .true. detailed command line output is written \n
  !>                                                                         DEFAULT: .false.

  !>        \param[in]  "LOGICAL, DIMENSION(size(para)), optional  :: maskpara"
  !>                                                                         maskpara(i) = .true.  --> parameter is optimized \n
  !>                                                                         maskpara(i) = .false. --> parameter is discarded
  !>                                                                                                   from optimiztaion \n
  !>                                                                         DEFAULT: .true.

  !>        \param[in]  "REAL(DP), DIMENSION(size(para)), optional :: weight" 
  !>                                                                         vector of weights per parameter: \n
  !>                                                                         gives the frequency of parameter to be chosen for
  !>                                                                         optimization (will be scaled to a CDF internally) \n
  !>                                                                         eg. [1,2,1] --> parameter 2 is chosen twice as
  !>                                                                         often as parameter 1 and 2 \n
  !>                                                                         DEFAULT: weight = 1.0_dp

  !>        \param[in]  "INTEGER(I4), optional           :: changeParaMode"  which and how many param.are changed in one step \n
  !>                                                                         1 = one parameter \n
  !>                                                                         2 = all parameter \n
  !>                                                                         3 = neighborhood parameter \n
  !>                                                                         DEFAULT: 1_i4

  !>        \param[in]  "LOGICAL, optional               :: reflectionFlag"  if new parameter values are Gaussian 
  !>                                                                         distributed and reflected (.true.) or
  !>                                                                         uniform in range (.false.) \n
  !>                                                                         DEFAULT: .false.

  !>        \param[in]  "LOGICAL, optional               :: pertubFlexFlag"  if pertubation of Gaussian distributed 
  !>                                                                         parameter values is constant 
  !>                                                                         at 0.2 (.false.) or 
  !>                                                                         depends on dR (.true.) \n
  !>                                                                         DEFAULT: .true.

  !>        \param[in]  "LOGICAL, optional               :: maxit"           maximizing (.true.) or minimizing (.false.)
  !>                                                                         a function \n
  !>                                                                         DEFAULT: .false. (minimization)

  !>        \param[in]  "REAL(DP), optional              :: undef_funcval"   objective function value defining invalid
  !>                                                                         model output, e.g. -999.0_dp \n

  !>        \param[in]  "character(len=*) , optional     :: tmp_file"        file with temporal output

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out]  "REAL(DP), optional             :: funcbest"        minimized value of cost function

  !>        \param[out]  "REAL(DP), DIMENSION(:,:), ALLOCATABLE, optional
  !>                                                     :: history"         returns a vector of achieved objective
  !                                                                          after ith model evaluation

  !     RETURN                                                            
  !>        \return real(dp) :: parabest(size(para))  &mdash;  Parameter set minimizing the cost function.

  !     RESTRICTIONS
  !>        \note Either fixed parameter range (prange) OR flexible parameter range (function interface prange_func)
  !>              has to be given in calling sequence. \n
  !>              Only double precision version available. \n
  !>              If single precision is needed not only DP has to be replaced by SP
  !>              but also I8 of save_state (random number variables) has to be replaced by I4. \n
  !>              ParaChangeMode > 1 is not applied in GetTemperature. 
  !>              For Temperature estimation always only one single parameter is changed (ParaChangeMode=1)
  !>              which should give theoretically always the best estimate. \n
  !>              Cost and prange_func are user defined functions. See interface definition.

  !     EXAMPLE
  !         see test/test_mo_anneal/

  !     LITERATURE
  !         (1) S. Kirkpatrick, C. D. Gelatt, and M. P. Vecchi.
  !             Optimization by simulated annealing.
  !             Science, 220:671-680, 1983.

  !         (2) N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, and E. Teller.
  !             Equation of state calculations by fast computing machines.
  !             J. Chem. Phys., 21:1087-1092, June 1953.
  !
  !         (3) B. A. Tolson, and C. A. Shoemaker.
  !             Dynamically dimensioned search algorithm for computationally efficient watershed model calibration.
  !             WRR, 43(1), W01413, 2007.
  !
  !     HISTORY
  !        Written  Samaniego,   Jan   2000 : Created
  !                 Samaniego,   Mar   2003 : Re-heating
  !                 Juliane Mai, March 2012 : modular version
  !        Modified Juliane Mai, May   2012 : sp version
  !                 Juliane Mai, May   2012 : documentation

  ! ------------------------------------------------------------------

  INTERFACE anneal
     MODULE PROCEDURE anneal_dp
  END INTERFACE anneal

  ! ------------------------------------------------------------------

  !     NAME
  !>        \brief GetTemperature

  !     PURPOSE
  !>        \details Determines an initial temperature for Simulated Annealing achieving
  !         certain acceptance ratio.


  !     CALLING SEQUENCE
  !         para = (/ 1.0_dp , 2.0_dp /)
  !         acc_goal   = 0.95_dp
  !         prange(1,:) = (/ 0.0_dp, 10.0_dp /)
  !         prange(2,:) = (/ 0.0_dp, 50.0_dp /)
  !
  !         user defined function cost_dp which calculates the cost function value for a
  !         parameter set (the interface given below has to be used for this function!)

  !         temp = GetTemperature(para, cost_dp, acc_goal, prange)

  !         see also test folder for a detailed example

  !     INTENT(IN)
  !>        \param[in] "REAL(DP),    DIMENSION(:)   :: paraset"   initial (valid) parameter set
  !>        \param[in] "INTERFACE                   :: cost_dp"   interface calculating the
  !>                                                              cost function at a given point
  !>        \param[in] "REAL(DP)                    :: acc_goal"  Acceptance Ratio which has to be achieved

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "REAL(DP), DIMENSION(size(para),2), optional :: prange"
  !>                                                                         lower and upper bound per parameter
  !             OR
  !>        \param[in] "INTERFACE, optional              :: prange_func"     interface calculating the
  !>                                                                         feasible range for a parameter
  !>                                                                         at a certain point, if ranges are variable

  !>        \param[in] "INTEGER(I4), optional            :: samplesize"      number of iterations the estimation of temperature
  !>                                                                         is based on \n
  !>                                                                         DEFAULT: Max(20_i4*n,250_i4)

  !>        \param[in] "LOGICAL, DIMENSION(size(para)), optional    :: maskpara"
  !>                                                                         maskpara(i) = .true.  --> parameter is optimized \n
  !>                                                                         maskpara(i) = .false. --> parameter is discarded
  !>                                                                                                   from optimiztaion \n
  !>                                                                         DEFAULT: .true.

  !>        \param[in] "INTEGER(I4/I8), DIMENSION(2), optional      :: seeds"
  !>                                                                         Seeds of random numbers used for random parameter
  !>                                                                         set generation \n
  !>                                                                         DEFAULT: dependent on current time

  !>        \param[in] "LOGICAL, optional                           :: printflag"
  !>                                                                         If .true. detailed command line output is written
  !>                                                                         DEFAULT: .false.

  !>        \param[in] "REAL(DP), DIMENSION(size(para,1)), optional :: weight"
  !>                                                                         vector of weights per parameter \n
  !>                                                                         gives the frequency of parameter to be chosen for
  !>                                                                         optimization (will be scaled to a CDF internally) \n
  !>                                                                         eg. [1,2,1] --> parameter 2 is chosen twice as
  !>                                                                                         often as parameter 1 and 2 \n
  !                                                                          DEFAULT: weight = 1.0_dp

  !>        \param[in] "LOGICAL, optional                :: maxit"            minimizing (.false.) or maximizing (.true.) 
  !>                                                                          a function \n
  !>                                                                          DEFAULT: .false. (minimization)

  !>        \param[in] "REAL(DP), optional               :: undef_funcval"    objective function value defining invalid
  !>                                                                          model output, e.g. -999.0_dp 

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: temperature  &mdash;  Temperature achieving a certain
  !>                                                  acceptance ratio in Simulated Annealing

  !     RESTRICTIONS
  !>        \note Either fixed parameter range (prange) OR flexible parameter range (function interface prange_func)
  !>              has to be given in calling sequence. \n

  !>              Only double precision version available.
  !>              If single precision is needed not only DP has to be replaced by SP
  !>              but also I8 of save_state (random number variables) has to be replaced by I4. \n

  !>              ParaChangeMode > 1 is not applied in GetTemperature.
  !>              For Temperature estimation always only one single parameter is changed (ParaChangeMode=1)
  !>              which should give theoretically always the best estimate. \n

  !>              Cost and prange_func are user defined functions. See interface definition.

  !     EXAMPLE
  !         see test/test_mo_anneal/

  !     LITERATURE
  !         (1) Walid Ben-Ameur
  !             Compututing the Initial Temperature of Simulated Annealing
  !             Comput. Opt. and App. (2004).

  !     HISTORY
  !        Written  Juliane Mai,   May   2012

  ! ------------------------------------------------------------------

  INTERFACE GetTemperature
     MODULE PROCEDURE GetTemperature_dp
  END INTERFACE GetTemperature

  PRIVATE

  INTERFACE Generate_Neighborhood_weight
     MODULE PROCEDURE Generate_Neighborhood_weight_dp
  END INTERFACE Generate_Neighborhood_weight

  ! ------------------------------------------------------------------

CONTAINS

  FUNCTION anneal_dp(  cost, para,                &   ! obligatory
       prange, prange_func,                       &   ! optional IN: exactly one of both
       temp, Dt, nITERmax, Len, nST, eps, acc,    &   ! optional IN
       seeds, printflag, maskpara, weight,        &   ! optional IN
       changeParaMode, reflectionFlag,            &   ! optional IN
       pertubFlexFlag, maxit, undef_funcval,      &   ! optional IN
       tmp_file,                                  &   ! optional IN
       funcbest, history                          &   ! optional OUT
       ) &
       result(parabest)

    IMPLICIT NONE

    INTERFACE
       FUNCTION cost(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind
         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
         REAL(DP)                            :: cost
       END FUNCTION cost
    END INTERFACE

    INTERFACE
       SUBROUTINE prange_func(paraset,iPar,rangePar)
         ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
         use mo_kind
         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
         INTEGER(I4),            INTENT(IN)  :: iPar
         REAL(DP), DIMENSION(2), INTENT(OUT) :: rangePar
       END SUBROUTINE prange_func
    END INTERFACE
    optional :: prange_func

    REAL(DP),    DIMENSION(:),                          INTENT(IN)    :: para
    !                                                                  ! initial parameter
    REAL(DP),    DIMENSION(size(para,1))                              :: parabest
    !                                                                  ! parameter set minimizing objective

    ! optionals
    REAL(DP),    OPTIONAL, DIMENSION(size(para,1),2),   INTENT(IN)  :: prange
    !                                                                  ! lower and upper limit per parameter
    REAL(DP),    OPTIONAL,                              INTENT(IN)  :: temp
    !                                                                  ! starting temperature
    !                                                                  ! (DEFAULT: Get_Temperature)
    REAL(DP),    OPTIONAL,                              INTENT(IN)  :: Dt
    !                                                                  ! geometrical decreement, 0.7<DT<0.999
    !                                                                  ! (DEFAULT: 0.9)
    INTEGER(I4), OPTIONAL,                              INTENT(IN)  :: nITERmax
    !                                                                  ! maximal number of iterations
    !                                                                  ! (DEFAULT: 1000)
    INTEGER(I4), OPTIONAL,                              INTENT(IN)  :: Len
    !                                                                  ! Length of Markov Chain,
    !                                                                  ! DEFAULT: max(250, size(para,1))
    INTEGER(I4), OPTIONAL,                              INTENT(IN)  :: nST
    !                                                                  ! Number of consecutive LEN steps
    !                                                                  ! (DEFAULT: 5)
    REAL(DP),    OPTIONAL,                              INTENT(IN)  :: eps
    !                                                                  ! epsilon decreement of cost function
    !                                                                  ! (DEFAULT: 0.01)
    REAL(DP),    OPTIONAL,                              INTENT(IN)  :: acc
    !                                                                  ! Acceptance Ratio, <0.1 stopping criteria
    !                                                                  ! (DEFAULT: 0.1)
    INTEGER(I8), OPTIONAL,                              INTENT(IN)  :: seeds(3)
    !                                                                  ! Seeds of random numbers
    !                                                                  ! (DEFAULT: Get_timeseed)
    LOGICAL,     OPTIONAL,                              INTENT(IN)  :: printflag
    !                                                                  ! If command line output is written (.true.)
    !                                                                  ! (DEFAULT: .false.)
    LOGICAL,     OPTIONAL, DIMENSION(size(para,1)),     INTENT(IN)  :: maskpara
    !                                                                  ! true if parameter will be optimized
    !                                                                  ! false if parameter is discarded in optimization
    !                                                                  ! (DEFAULT: .true.)
    REAL(DP),    OPTIONAL, DIMENSION(size(para,1)),     INTENT(IN)  :: weight
    !                                                                  ! vector of weights per parameter
    !                                                                  ! gives the frequency of parameter to be
    !                                                                  ! chosen for optimization
    !                                                                  ! (DEFAULT: uniform)
    INTEGER(I4), OPTIONAL,                              INTENT(IN)  :: changeParaMode
    !                                                                  ! which and how many param. are changed in a step
    !                                                                  ! 1 = one parameter
    !                                                                  ! 2 = all parameter
    !                                                                  ! 3 = neighborhood parameter
    !                                                                  ! (DEFAULT: 1_i4)
    LOGICAL,     OPTIONAL,                              INTENT(IN)  :: reflectionFlag
    !                                                                  ! if new parameter values are selected normal
    !                                                                  ! distributed and reflected (.true.) or
    !                                                                  ! uniform in range (.false.)
    !                                                                  ! (DEFAULT: .false.)
    LOGICAL,     OPTIONAL,                              INTENT(IN)  :: pertubFlexFlag
    !                                                                  ! if pertubation of normal distributed parameter
    !                                                                  ! values is constant 0.2 (.false.) or
    !                                                                  ! depends on dR (.true.)
    !                                                                  ! (DEFAULT: .true.)
    LOGICAL,     OPTIONAL,                              INTENT(IN)  :: maxit
    !                                                                  ! Maximization or minimization of function
    !                                                                  ! maximization = .true., minimization = .false.
    !                                                                  ! (DEFAULT: .false.)
    REAL(DP),    OPTIONAL,                              INTENT(IN)  :: undef_funcval
    !                                                                  ! objective function value occuring if
    !                                                                  ! parameter set leads to  invalid model results,
    !                                                                  ! e.g. -9999.0_dp
    !                                                                  ! (DEFAULT: not present)
    CHARACTER(LEN=*), OPTIONAL,                         INTENT(IN)  :: tmp_file    
    !                                                                  ! file for temporal output
    REAL(DP),    OPTIONAL,                              INTENT(OUT) :: funcbest
    !                                                                  ! minimized value of cost function
    !                                                                  ! (DEFAULT: not present)
    REAL(DP),    OPTIONAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: history
    !                                                                  ! returns a vector of achieved objective
    !                                                                  ! after ith model evaluation
    !                                                                  ! (DEFAULT: not present)

    ! local variables
    INTEGER(I4)                              :: n              ! Number of parameters
    integer(I4)                              :: iPar, par      ! counter for parameter
    real(DP), DIMENSION(2)                   :: iParRange      ! parameter's range
    REAL(DP)                                 :: T_in           ! Temperature
    REAL(DP)                                 :: DT_IN          ! Temperature decreement
    INTEGER(I4)                              :: nITERmax_in    ! maximal number of iterations
    INTEGER(I4)                              :: LEN_IN         ! Length of Markov Chain
    INTEGER(I4)                              :: nST_in         ! Number of consecutive LEN steps
    REAL(DP)                                 :: eps_in         ! epsilon decreement of cost function
    REAL(DP)                                 :: acc_in         ! Acceptance Ratio, <0.1 stopping criteria
    INTEGER(I8)                              :: seeds_in(3)    ! Seeds of random numbers
    LOGICAL                                  :: printflag_in   ! If command line output is written
    LOGICAL                                  :: coststatus     ! Checks status of cost function value,
    !                                                          ! i.e. is parameter set is feasible
    LOGICAL,     DIMENSION(size(para,1))     :: maskpara_in    ! true if parameter will be optimized
    INTEGER(I4), DIMENSION(:), ALLOCATABLE   :: truepara       ! indexes of parameters to be optimized
    REAL(DP),    DIMENSION(size(para,1))     :: weight_in      ! CDF of parameter chosen for optimization
    REAL(DP),    DIMENSION(size(para,1))     :: weightUni      ! uniform CDF of parameter chosen for optimization
    REAL(DP),    DIMENSION(size(para,1))     :: weightGrad     ! linear combination of weight and weightUni
    INTEGER(I4)                              :: changeParaMode_inin  ! 1=one, 2=all, 3=neighborhood
    LOGICAL                                  :: reflectionFlag_inin  ! true=gaussian distributed and reflected,
    !                                                                ! false=uniform distributed in range
    LOGICAL                                  :: pertubFlexFlag_inin  ! variance of normal distributed parameter values
    !                                                                ! .true. = depends on dR, .false. = constant 0.2
    REAL(DP)                                 :: maxit_in             ! maximization = -1._dp, minimization = 1._dp
    REAL(DP)                                 :: costbest             ! minimized value of cost function
    REAL(DP),    DIMENSION(:,:), ALLOCATABLE :: history_out          ! best cost function value after k iterations
    REAL(DP),    DIMENSION(:,:), ALLOCATABLE :: history_out_tmp      ! helper vector

    type paramLim
       real(DP)                                 :: min         ! minimum value
       real(DP)                                 :: max         ! maximum value
       real(DP)                                 :: new         ! new state value
       real(DP)                                 :: old         ! old state value
       real(DP)                                 :: best        ! best value found
       real(DP)                                 :: dMult       ! sensitivity multiplier for parameter search
    end type paramLim
    type (paramLim), dimension (size(para,1)), target :: gamma ! Parameter

    ! for random numbers
    real(DP)                                :: RN1, RN2, RN3     ! Random numbers
    integer(I8), dimension(n_save_state)    :: save_state_1      ! optional arguments for restarting RN stream 1
    integer(I8), dimension(n_save_state)    :: save_state_2      ! optional arguments for restarting RN stream 2
    integer(I8), dimension(n_save_state)    :: save_state_3      ! optional arguments for restarting RN stream 3
    ! for dds parameter selection
    logical, dimension(size(para,1))        :: neighborhood      ! selected parameter in neighborhood
    real(dp)                                :: pertubationR      ! neighborhood pertubation size parameter

    ! for SA
    integer(I4)                              :: idummy,i
    real(DP)                                 :: NormPhi
    real(DP)                                 :: ac_ratio, pa
    integer(i4), dimension(:,:), allocatable :: iPos_iNeg_history      ! stores iPos and iNeg of nST last Markov Chains
    real(DP)                                 :: fo, fn, df, fBest
    real(DP)                                 :: rho, fInc, fbb, dr
    real(DP)                                 :: T0, DT0
    real(DP),  parameter                     :: small = -700._DP
    integer(I4)                              :: j, iter, kk
    integer(I4)                              :: Ipos, Ineg
    integer(I4)                              :: iConL, iConR, iConF
    integer(I4)                              :: iTotalCounter          ! includes reheating for final conditions
    integer(I4)                              :: iTotalCounterR         ! counter of interations in one reheating
    logical                                  :: iStop
    logical                                  :: ldummy

    ! CHECKING OPTIONALS

    n = size(para,1)

    ! either range or rangfunc has to be given
    if ( present(prange) .eqv. present(prange_func) ) then
       stop 'anneal: Either range or prange_func has to be given'
    end if

    if (present(Dt)) then
       if ( (Dt .lt. 0.7_dp) .or. (Dt .gt. 0.999_dp) ) then
          stop 'Input argument DT must lie between 0.7 and 0.999'
       else
          DT_IN = Dt
       end if
    else
       DT_IN = 0.9_dp
    endif

    if (present(nITERmax)) then
       if (nITERmax .lt. 1_I4) then
          stop 'Input argument nITERmax must be greater 0'
       else
          nITERmax_in = nITERmax
       end if
    else
       nITERmax_in = 1000_i4
    endif

    if (present(Len)) then
       if (Len .lt. Max(20_i4*n,250_i4)) then
          print*, 'WARNING: Input argument LEN should be greater than Max(250,20*N), N=number of parameters'
          LEN_IN = Len
       else
          LEN_IN = Len
       end if
    else
       LEN_IN = Max(20_i4*n,250_i4)
    endif

    idummy = nITERmax_in/LEN_IN+1_i4
    allocate(history_out(idummy,2))

    if (present(nST)) then
       if (nST .lt. 0_i4) then
          stop 'Input argument nST must be greater than 0'
       else
          nST_in = nST
       end if
    else
       nST_in = 5_i4
    endif

    allocate(iPos_iNeg_history(nST_in,2))
    iPos_iNeg_history = 0_i4

    if (present(eps)) then
       if ( le(eps, 0.0_dp) ) then
          stop 'Input argument eps must be greater than 0'
       else
          eps_in = eps
       end if
    else
       eps_in = 0.0001_dp
    endif

    if (present(acc)) then
       if ( le(acc, 0.0_dp)  .or. ge(acc, 1.0_dp) ) then
          stop 'Input argument acc must lie between 0.0 and 1.0'
       else
          acc_in = acc
       end if
    else
       acc_in = 0.1_dp
    endif

    if (present(seeds)) then
       seeds_in = seeds
    else
       ! Seeds depend on actual time
       call get_timeseed(seeds_in)
    endif

    if (present(printflag)) then
       printflag_in = printflag
    else
       printflag_in = .false.
    endif

    if (present(maskpara)) then
       if (count(maskpara) .eq. 0_i4) then
          stop 'Input argument maskpara: At least one element has to be true'
       else
          maskpara_in = maskpara
       end if
    else
       maskpara_in = .true.
    endif

    if (present(maxit)) then
       ldummy = maxit
       if (maxit) then
          maxit_in = -1._dp
       else
          maxit_in = 1._dp
       end if
    else
       ldummy = .false.
       maxit_in = 1._dp
    endif

    if (present(temp)) then
       if ( (temp .lt. 0.0_dp) ) then
          stop 'Input argument temp must be greater then zero'
       else
          T_in = temp
       end if
    else
       if (present(prange_func)) then
          if (present(undef_funcval)) then
             T_in = GetTemperature( para, cost, 0.95_dp,  prange_func=prange_func, &
                  maskpara=maskpara_in, samplesize=2_i4*LEN_IN, &
                  seeds=seeds_in(1:2), printflag=printflag_in, &
                  maxit=ldummy,undef_funcval=undef_funcval)
          else
             T_in = GetTemperature( para, cost, 0.95_dp,  prange_func=prange_func, &
                  maskpara=maskpara_in, samplesize=2_i4*LEN_IN, &
                  seeds=seeds_in(1:2), printflag=printflag_in, &
                  maxit=ldummy)
          end if
       else
          if (present(undef_funcval)) then
             T_in = GetTemperature( para, cost, 0.95_dp, prange=prange, &
                  maskpara=maskpara_in, samplesize=2_i4*LEN_IN, &
                  seeds=seeds_in(1:2), printflag=printflag_in, &
                  maxit=ldummy,undef_funcval=undef_funcval)
          else
             T_in = GetTemperature( para, cost, 0.95_dp, prange=prange, &
                  maskpara=maskpara_in, samplesize=2_i4*LEN_IN, &
                  seeds=seeds_in(1:2), printflag=printflag_in, &
                  maxit=ldummy)
          end if
       end if
    end if

    if (present(changeParaMode)) then
       changeParaMode_inin = changeParaMode
    else
       changeParaMode_inin = 1_i4
    end if

    if (present(reflectionFlag)) then
       reflectionFlag_inin = reflectionFlag
    else
       reflectionFlag_inin = .false.
    end if

    if (present(pertubFlexFlag)) then
       pertubFlexFlag_inin = pertubFlexFlag
    else
       pertubFlexFlag_inin = .true.
    end if

    ! Temporal file writing
    if(present(tmp_file)) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', status = 'unknown')
       write(999,*) '# settings :: general'
       write(999,*) '# nIterations    iseed'
       write(999,*) nITERmax_in, seeds_in
       write(999,*) '# settings :: anneal specific'
       write(999,*) '# sa_tmp'
       write(999,*) T_in
       write(999,*) '# iter   bestf   (bestx(j),j=1,nopt)'
       close(999)
    end if

    ! INITIALIZATION

    allocate ( truepara(count(maskpara_in)) )
    idummy = 0_i4
    do i=1,N
       if ( maskpara_in(i) ) then
          idummy = idummy+1_i4
          truepara(idummy) = i
       end if
    end do

    if (printflag_in) then
       print*, 'Following parameters will be optimized: ',truepara
    end if

    weight_in = 0.0_dp
    weightUni = 0.0_dp
    if (present(weight)) then
       where ( maskpara_in(:) )
          weight_in(:)    = weight(:)
          weightUni(:) = 1.0_dp
       end where
    else
       where ( maskpara_in(:) )
          weight_in(:)    = 1.0_dp
          weightUni(:) = 1.0_dp
       end where
    endif
    ! scaling the weights
    weight_in    = weight_in/sum(weight_in)
    weightUni = weightUni/sum(weightUni)
    ! cummulating the weights
    do i=2,n
       weight_in(i)    = weight_in(i)    + weight_in(i-1)
       weightUni(i) = weightUni(i) + weightUni(i-1)
    end do

    call xor4096 (seeds_in(1), RN1,    save_state=save_state_1)
    call xor4096 (seeds_in(2), RN2,    save_state=save_state_2)
    call xor4096g(seeds_in(3), RN3,    save_state=save_state_3)
    seeds_in = 0_i8

    ! Start Simulated Annealing routine
    gamma(:)%dmult = 1.0_dp
    gamma(:)%new = para(:)
    gamma(:)%old = para(:)
    gamma(:)%best = para(:)
    NormPhi = -9999.9_dp
    T0=      T_in
    DT0=     DT_IN

    ! Generate and evaluate the initial solution state
    fo = cost(gamma(:)%old) * maxit_in
    if ( abs(fo) .lt. tiny(0.0_dp) ) fo = 0.0000001_dp * maxit_in

    file_write: if (present(tmp_file)) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append', recl=(n+2)*30)
       if (.not. ldummy) then
          write(999,*) '0', fo, gamma(:)%old
       else
          write(999,*) '0', -fo, gamma(:)%old
       end if
       close(999)
    end if file_write

    ! initialize counters /var for new SA
    iConL=           0_i4
    iConR=           1_i4
    iConF=           0_i4
    iTotalCounter=   0_i4
    iTotalCounterR = 0_i4
    fn =             0.0_DP
    ac_ratio =       1.0_DP
    iStop=           .TRUE.
    iter =           0_i4
    ! restoring initial T param.
    DT_IN=     DT0
    ! Storing the best solution so far
    NormPhi = fo * maxit_in
    fo      = fo/NormPhi
    fBest   = 1.0_dp * maxit_in

    if (printflag_in) then
       print '(A15,E15.7,A4,E15.7,A4)',     ' start NSe   = ', fBest*maxit_in , '  ( ',fBest*normPhi*maxit_in,' )  '
       print '(I8, 2I5, 4E15.7)',   1_i4, 0_i4, 0_i4, 1._dp, T_in, fo*normPhi*maxit_in, fBest*normPhi*maxit_in
    end if

    ! ****************** Stop Criterium  *******************
    ! Repeat until the % reduction of nST consecutive Markov
    ! chains (LEN) of the objective function (f) <= epsilon
    ! ******************************************************
    loopTest: do while (iStop)
       iter = iter + 1_i4
       Ipos= 0_i4
       Ineg= 0_i4
       fbb=  fBest
       !LEN_IN = int( real(iter,dp)*sqrt(real(iter,dp)) / nITER + 1.5*real(N,dp),i4)
       ! Repeat LEN times with feasible solution
       j=1
       loopLEN: do while (j .le. LEN_IN)

          iTotalCounterR =  iTotalCounterR + 1_i4
          iTotalCounter = iTotalCounter + 1_i4

          ! Generate a random subsequent state and evaluate its objective function
          ! (1) Generate new parameter set
          dR=(1.0_DP - real(iTotalCounterR,dp) / real(nITERmax_in,dp))**2.0_DP
          if (  .not. ( ge(dR , 0.05_DP) ) .and. &
               (iTotalCounterR <= int(real(nIterMax_in,dp)/3._dp*4_dp,i4))) then
             dR = 0.05_DP
          end if

          if (pertubFlexFlag_inin) then
             pertubationR = max(dR/5.0_dp,0.05_dp)
          else
             pertubationR = 0.2_dp
          end if

          ! gradual weights:
          !        acc_ratio close to 1 --> weight
          !        acc_ratio close to 0 --> weight uniform
          ! gradual weighting is linear combination of these two weights
          weightGrad(:) = weightUni(:)+(ac_ratio)*(weight_in(:)-weightUni(:))

          select case(changeParaMode_inin)
          case(1_i4)  ! only one parameter is changed
             ! (1a) Change one parameter traditionally
             call xor4096(seeds_in(2),RN2, save_state=save_state_2)
             iPar=1_i4
             !
             do while (weightGrad(iPar) .lt. RN2)
                iPar = iPar + 1_i4
             end do
             if (present(prange_func)) then
                call prange_func(gamma(:)%old, iPar, iParRange)
             else
                iParRange(1) = prange(iPar,1)
                iParRange(2) = prange(iPar,2)
             endif
             gamma(iPar)%min = iParRange(1)
             gamma(iPar)%max = iParRange(2)
             if (reflectionFlag_inin) then
                call xor4096g(seeds_in(3),RN3, save_state=save_state_3)
                gamma(iPar)%new = parGen_dds_dp( gamma(iPar)%old, pertubationR, &
                     gamma(iPar)%min, gamma(iPar)%max,RN3)
             else
                call xor4096(seeds_in(2),RN2, save_state=save_state_2)
                gamma(iPar)%new = parGen_anneal_dp( gamma(iPar)%old, dR, &
                     gamma(iPar)%min, gamma(iPar)%max,RN2)
             end if
          case(2_i4)  ! all parameter are changed
             do par=1, size(truepara)
                iPar = truepara(par)
                if (present(prange_func)) then
                   call prange_func(gamma(:)%old, iPar, iParRange)
                else
                   iParRange(1) = prange(iPar,1)
                   iParRange(2) = prange(iPar,2)
                endif
                gamma(iPar)%min = iParRange(1)
                gamma(iPar)%max = iParRange(2)
                if (reflectionFlag_inin) then
                   call xor4096g(seeds_in(3),RN3, save_state=save_state_3)
                   gamma(iPar)%new = parGen_dds_dp( gamma(iPar)%old, pertubationR, &
                        gamma(iPar)%min, gamma(iPar)%max,RN3)
                else
                   call xor4096(seeds_in(2),RN2, save_state=save_state_2)
                   gamma(iPar)%new = parGen_anneal_dp( gamma(iPar)%old, dR, &
                        gamma(iPar)%min, gamma(iPar)%max,RN2)
                end if
             end do
          case(3_i4)  ! parameter in neighborhood are changed
             ! Generate new neighborhood
             call generate_neighborhood_weight_dp( truepara, weightGrad, save_state_1, &
                  iTotalCounter, nIterMax_in, neighborhood)
             !
             ! change parameter in neighborhood
             do iPar=1, n
                if (neighborhood(iPar)) then
                   ! find range of parameter
                   if (present(prange_func)) then
                      call prange_func(gamma(:)%old, iPar, iParRange)
                   else
                      iParRange(1) = prange(iPar,1)
                      iParRange(2) = prange(iPar,2)
                   endif
                   gamma(iPar)%min = iParRange(1)
                   gamma(iPar)%max = iParRange(2)
                   !
                   if (reflectionFlag_inin) then
                      ! generate gaussian distributed new parameter value which is reflected if out of bound
                      call xor4096g(seeds_in(3),RN3, save_state=save_state_3)
                      gamma(iPar)%new = parGen_dds_dp( gamma(iPar)%old, pertubationR, &
                           gamma(iPar)%min, gamma(iPar)%max,RN3)
                   else
                      ! generate new parameter value uniform distributed in range (no reflection)
                      call xor4096(seeds_in(2),RN2, save_state=save_state_2)
                      gamma(iPar)%new = parGen_anneal_dp( gamma(iPar)%old, dR, &
                           gamma(iPar)%min, gamma(iPar)%max,RN2)
                   end if
                end if
             end do
          end select

          ! (2) Calculate new objective function value
          fn = cost(gamma(:)%new) * maxit_in
          coststatus = .true.
          if (present(undef_funcval)) then
             if ( abs(fn*maxit_in-undef_funcval) .lt. tiny(1.0_dp) ) then
                coststatus = .false.
             end if
          end if

          feasible: if (coststatus) then   ! feasible parameter set
             fn = fn/normPhi

             ! Change in cost function value
             df = fn-fo

             ! analyze change in the objective function: df
             if (df < 0.0_DP) then

                ! accept the new state
                Ipos=Ipos+1_i4
                fo = fn
                gamma(:)%old   = gamma(:)%new

                ! keep best solution
                if (fo < fBest) then
                   fBest =  fo
                   gamma(:)%best   = gamma(:)%new
                endif
             else
                if ( df >  eps_in ) then
                   rho=-df/T_in
                   if (rho < small) then
                      pa=0.0_DP
                   else
                      pa=EXP(rho)
                   end if
  	           !
                   call xor4096(seeds_in(1), RN1, save_state=save_state_1)
                   !
                   if (pa > RN1) then
                      ! accept new state with certain probability
                      Ineg=Ineg+1_i4
                      fo = fn
                      ! save old state
                      gamma(:)%old   = gamma(:)%new
                   end if
                end if
             end if
             j=j+1
          else
             ! function value was not valid
             iTotalCounterR =  iTotalCounterR - 1_i4
             iTotalCounter = iTotalCounter - 1_i4
          end if feasible !valid parameter set
       end do loopLEN

       ! estimate acceptance ratio
       ac_ratio=(real(Ipos,dp) + real(Ineg,dp))/real(LEN_IN,dp)
       if (modulo(iTotalCounter,LEN_IN*nST_in)/LEN_IN .gt. 0_i4) then
          idummy = modulo(iTotalCounter,LEN_IN*nST_in)/LEN_IN
       else
          idummy = nST_in
       end if
       iPos_iNeg_history(idummy,:) = (/ iPos, iNeg /)

       ! store current best in history vector
       history_out(iTotalCounter/LEN_IN,1) = real(iTotalCounter,dp)
       history_out(iTotalCounter/LEN_IN,2) = maxit_in * fBest*normPhi

       file_write2: if (present(tmp_file)) then
          open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append',recl=(n+2)*30)
          write(999,*) iTotalCounter, maxit_in * fBest*normPhi, gamma(:)%best
          close(999)
       end if file_write2

       if (printflag_in) then
          print '(I8, 2I5, E15.7, 3E15.7)', ITotalCounter, Ipos, Ineg, ac_ratio, T_in, &
               fo*NormPhi*maxit_in, fBest*NormPhi*maxit_in
       end if

       ! Cooling schedule
       if (fbb .gt. tiny(1.0_dp)) then
          fInc= (fbb-fBest)/fbb
       else
          fInc = 1.0_dp
       end if
       if (fInc < 0.00000001_dp) then
          iConF= iConF+1_i4
       else
          iConF=0_i4
       end if

       if ((ac_ratio < 0.15_dp) .and. (iConF > 5_i4) .and. (iConR <= -3_i4)) then     ! - iConR  no reheating
          ! Re-heating
          if (printflag_in) then
             print *, 'Re-heating: ', iConR
          end if
          iConR          = iConR+1_i4
          T_in           = T0/2._DP
          iter           = 0_i4       ! for LEN
          iTotalCounterR = 0_i4       ! for dR
          ! start from current best
          gamma(:)%old   = gamma(:)%best
       else
          ! Update Temperature (geometrical decrement)
          if (ac_ratio < 0.4_dp)  then
             DT_IN=0.995_dp
          else
             DT_IN=DT0
          end if
          T_in = T_in*DT_IN
       end if

       ! Stop Criteria: consecutive MC with marginal decrements and acceptance ratio
       if (fInc < eps_in) then
          iConL= iConL+1_i4
       else
          iConL=0_i4
       end if
       if ( (iConL > nST_in) .and. (ac_ratio < acc_in) .and. (iConR > 2_i4) .and. &
            (sum(iPos_iNeg_history(:,:)) .lt. 1_i4 ) ) then
          iStop=.FALSE.
       end if
       ! a way out maximum
       if ( (iTotalCounter > nITERmax_in) .and. & !(ac_ratio > acc_in) .and. &
            (sum(iPos_iNeg_history(:,:)) .ge. 1_i4 )) then

          ! store achieved history so far
          allocate(history_out_tmp(size(history_out,1),size(history_out,2)))
          history_out_tmp = history_out
          deallocate(history_out)

          nITERmax_in = max(nITERmax_in+LEN_IN,int(real(nITERmax_in,dp)* 1.10_dp,i4))
          if (printflag_in) then
             print *,                       'nITERmax changed to =', nITERmax_in
          end if

          ! increase lenght of history vecor
          idummy = nITERmax_in/LEN_IN+1_i4
          allocate(history_out(idummy,2))

          do j=1,size(history_out_tmp,1)
             history_out(j,:) = history_out_tmp(j,:)
          end do

          deallocate(history_out_tmp)

       end if
       ! STOP condition
       if ( iTotalCounter > nITERmax_in ) then
          iStop=.FALSE.
       end if

    end do loopTest

    ! calculate cost function again (only for check and return values)
    parabest = gamma(:)%best
    costbest = cost(parabest) * maxit_in
    if (present (funcbest)) then
       funcbest = costbest * maxit_in
    end if
    fo = costbest/NormPhi

    if (printflag_in) then
       print *, '   '
       print '(A15,E15.7)',                   ' end cost    = ', maxit_in * fBest*normPhi
       print *,           'end parameter: '
       do kk = 1,N
          print '(A10,I3,A3, E15.7)' ,    '    para #',kk,' = ', gamma(kk)%best
       end do

       print *, 'Final check:    ', (fo - fBest)
    end if

    if (present(history)) then
       allocate(history(size(history_out,1)-1,size(history_out,2)))
       history(:,:) = history_out(1:size(history_out,1)-1,:)
    end if

    deallocate(truepara)
    deallocate(history_out)

    return

  END FUNCTION anneal_dp

  real(DP) function GetTemperature_dp( paraset, cost, acc_goal, &    ! obligatory
       prange, prange_func,                                     &    ! optional IN: exactly one of both
       samplesize, maskpara, seeds, printflag,                  &    ! optional IN
       weight, maxit, undef_funcval                             &    ! optional IN
       )
    use mo_kind, only: dp, i4, i8
    implicit none

    real(dp),               dimension(:),                 intent(in)    :: paraset
    !                                                                      ! a valid parameter set of the model
    real(dp),                                             intent(in)    :: acc_goal
    !                                                                      ! acceptance ratio to achieve
    real(dp),    optional,  dimension(size(paraset,1),2), intent(in)    :: prange
    !                                                                      ! lower and upper limit per parameter
    integer(i4), optional,                                intent(in)    :: samplesize
    !                                                                      ! size of random set the
    !                                                                      ! acc_estimate is based on
    !                                                                      ! DEFAULT: Max(250, 20*Number paras)
    logical,     optional,  dimension(size(paraset,1)),   intent(in)    :: maskpara
    !                                                                      ! true if parameter will be optimized
    !                                                                      ! false if parameter is discarded in
    !                                                                      ! optimization
    !                                                                      ! DEFAULT: .true.
    integer(i8), optional,  dimension(2),                 intent(in)    :: seeds
    !                                                                      ! Seeds of random numbers
    !                                                                      ! DEFAULT: time dependent
    logical,     optional,                                intent(in)    :: printflag
    !                                                                      ! .true. if detailed temperature
    !                                                                      ! estimation is printed
    !                                                                      ! DEFAULT: .false.
    REAL(DP),    OPTIONAL,  DIMENSION(size(paraset,1)),   INTENT(IN)    :: weight
    !                                                                      ! vector of weights per parameter
    !                                                                      ! gives the frequency of parameter to
    !                                                                      ! be chosen for optimization
    !                                                                      ! DEFAULT: equal weighting
    LOGICAL,     OPTIONAL,                                INTENT(IN)    :: maxit
    !                                                                      ! Maxim. or minim. of function
    !                                                                      ! maximization = .true.,
    !                                                                      ! minimization = .false.
    !                                                                      ! DEFAULT: .false.
    REAL(DP),    OPTIONAL,                                INTENT(IN)    :: undef_funcval
    !                                                                      ! objective function value occuring
    !                                                                      ! if parameter set leads to
    !                                                                      ! invalid model results, e.g. -999.0_dp
    !                                                                      ! DEFAULT: not present

    INTERFACE
       FUNCTION cost(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind
         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
         REAL(DP)                            :: cost
       END FUNCTION cost
    END INTERFACE

    INTERFACE
       SUBROUTINE prange_func(paraset,iPar,rangePar)
         ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
         use mo_kind
         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
         INTEGER(I4),            INTENT(IN)  :: iPar
         REAL(DP), DIMENSION(2), INTENT(OUT) :: rangePar
       END SUBROUTINE prange_func
    END INTERFACE
    OPTIONAL :: prange_func
    !
    ! Local variables

    integer(i4)                           :: n
    integer(i4)                           :: samplesize_in
    integer(i4)                           :: idummy, i, j
    integer(i4)                           :: iPar
    real(dp),      dimension(2)           :: iParRange
    real(dp)                              :: NormPhi
    real(dp)                              :: fo, fn, dr
    real(dp)                              :: T

    LOGICAL                                   :: coststatus     ! true if model output is valid
    LOGICAL,     DIMENSION(size(paraset,1))   :: maskpara_in    ! true if parameter will be optimized
    INTEGER(I4), DIMENSION(:), ALLOCATABLE    :: truepara       ! indexes of parameters to be optimized
    LOGICAL                                   :: printflag_in   ! if detailed estimation of temperature is printed
    REAL(DP),    DIMENSION(size(paraset,1))   :: weight_in      ! CDF of parameter to chose for optimization
    REAL(DP)                                  :: maxit_in       ! Maximization or minimization of function:
    !                                                           ! -1 = maxim, 1 = minim

    type paramLim
       real(DP)                                 :: min                 ! minimum value
       real(DP)                                 :: max                 ! maximum value
       real(DP)                                 :: new                 ! new state value
       real(DP)                                 :: old                 ! old state value
       real(DP)                                 :: best                ! best value found
       real(DP)                                 :: dMult               ! sensitivity multiplier
       !                                                               ! for parameter search
    end type paramLim
    type (paramLim), dimension (size(paraset,1)), target   :: gamma    ! Parameter

    ! for random numbers
    INTEGER(I8)                             :: seeds_in(2)
    real(DP)                                :: RN1, RN2     ! Random numbers
    integer(I8), dimension(n_save_state)    :: save_state_1      ! optional arguments for restarting RN stream 1
    integer(I8), dimension(n_save_state)    :: save_state_2      ! optional arguments for restarting RN stream 2

    ! for initial temperature estimate
    real(DP)                              :: acc_estim  ! estimate of acceptance probability
    ! depends on temperature
    ! goal: find a temperature such that acc_estim ~ 0.9
    real(DP), dimension(:,:), allocatable :: Energy     ! dim(LEN,2) cost function values before (:,1) and
    !                                                   ! after (:,2) transition

    n = size(paraset,1)

    ! either range or rangfunc has to be given
    if ( present(prange) .eqv. present(prange_func) ) then
       stop 'anneal: Either range or prange_func has to be given'
    end if

    if (present(samplesize)) then
       if (samplesize .lt. Max(20_i4*n,250_i4)) then
          !stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
          print*, 'WARNING (GetTemperature): '
          print*, 'Input argument samplesize should be greater than Max(250,20*N), N=number of parameters'
          samplesize_in = samplesize
       else
          samplesize_in = samplesize
       end if
    else
       samplesize_in = Max(20_i4*n,250_i4)
    endif

    if (present(maskpara)) then
       if (count(maskpara) .eq. 0_i4) then
          stop 'Input argument maskpara: At least one element has to be true'
       else
          maskpara_in = maskpara
       end if
    else
       maskpara_in = .true.
    endif

    if (present(seeds)) then
       seeds_in = seeds
    else
       ! Seeds depend on actual time
       call get_timeseed(seeds_in)
       print*,'temp: seeds(1)=', seeds_in(1)
    endif

    if (present(printflag)) then
       printflag_in = printflag
    else
       printflag_in = .false.
    endif

    if (present(maxit)) then
       if (maxit) then
          maxit_in = -1._dp
       else
          maxit_in = 1._dp
       end if
    else
       maxit_in = 1._dp
    endif

    allocate(Energy(samplesize_in,2))

    allocate ( truepara(count(maskpara_in)) )
    idummy = 0_i4
    do i=1,N
       if ( maskpara_in(i) ) then
          idummy = idummy+1_i4
          truepara(idummy) = i
       end if
    end do

    weight_in    = 0.0_dp
    if (present(weight)) then
       where ( maskpara_in(:) )
          weight_in(:) = weight(:)
       end where
    else
       where ( maskpara_in(:) )
          weight_in(:) = 1.0_dp
       end where
    endif
    ! scaling the weights
    weight_in = weight_in/sum(weight_in)
    ! cummulating the weights
    do i=2,n
       weight_in(i) = weight_in(i) + weight_in(i-1)
    end do

    ! Setting up the RNG
    ! (1) Seeds depend on actual time or on input seeds
    ! (2) Initialize the streams
    call xor4096(seeds_in(1), RN1, save_state=save_state_1)
    call xor4096(seeds_in(2), RN2, save_state=save_state_2)
    seeds_in = 0_i8
    ! (3) Now ready for calling

    gamma(:)%dmult = 1.0_dp
    gamma(:)%new   = paraset(:)
    gamma(:)%old   = paraset(:)
    gamma(:)%best  = paraset(:)
    NormPhi        = -9999.9_dp

    fo =  cost(paraset) * maxit_in
    if ( abs(fo) .lt. tiny(0.0_dp) ) fo = 0.0000001_dp * maxit_in
    NormPhi = fo
    fo      = fo/NormPhi * maxit_in

    j=1
    loopSamplesize: do while (j .le. samplesize_in)
       ! Generate a random subsequent state and evaluate its objective function
       ! (1)  Generate new parameter set
       dR=1.0_DP

       ! (1a) Select parameter to be changed
       call xor4096(seeds_in(1),RN1, save_state=save_state_1)
       iPar=1_i4
       do while (weight_in(iPar) .lt. RN1)
          iPar = iPar + 1_i4
       end do

       ! (1b) Generate new value of selected parameter
       if ( present(prange_func) ) then
          call prange_func(gamma(:)%old, iPar, iParRange )
       else
          iParRange(1) = prange(iPar,1)
          iParRange(2) = prange(iPar,2)
       end if
       gamma(iPar)%min = iParRange(1)
       gamma(iPar)%max = iParRange(2)
       call xor4096(seeds_in(2),RN2, save_state=save_state_2)
       gamma(iPar)%new = parGen_anneal_dp( gamma(iPar)%old, gamma(iPar)%dMult*dR, &
            gamma(iPar)%min, gamma(iPar)%max, RN2)
       !
       ! (2)  Calculate new objective function value and normalize it
       fn = cost(gamma(:)%new) * maxit_in
       coststatus = .true.
       if (present(undef_funcval)) then
          if ( abs(fn*maxit_in-undef_funcval) .lt. tiny(1.0_dp) ) then
             coststatus = .false.
          end if
       end if

       feasible: if (coststatus) then   ! feasible parameter set
          fn = fn/normPhi
          ! Save Energy states of feasible transitions
          ! for adaption of optimal initial temperature
          ! Walid Ben-Ameur: "Comput. the Initial Temperature of Sim. Annealing"
          ! Comput. Opt. and App. 2004
          Energy(j,2) = fn     ! E_max_t
          Energy(j,1) = fo     ! E_min_t
          j=j+1
       end if feasible !valid parameter set
    end do loopSamplesize

    ! estimation of the acceptance probability based on the random set ||<Samplesize>||
    ! only if actual temperature (T) equals initial temperature (temp)
    T = maxval(Energy)

    acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
    if (printflag_in) then
       print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
    end if
    Do While ( (acc_estim .lt. 1.0_dp) .and. (abs(acc_estim - acc_goal) .gt. 0.0001_dp))
       T = T * (Log(acc_estim)/Log(acc_goal))**(0.5_dp) ! **(1.0/p)  with p=1.0
       if ( all(T .gt. Energy(:,1)/709._dp) .and. all(T .gt. Energy(:,2)/709._dp) ) then
          acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
          if (printflag_in) then
             print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
          end if
       else
          T = T/(Log(acc_estim)/Log(acc_goal))**(0.5_dp)
          exit
       end if
    end do
    GetTemperature_dp = T

  end function GetTemperature_dp

  !***************************************************
  !*               PRIVATE FUNCTIONS                 *
  !***************************************************

  real(DP) function parGen_anneal_dp(old, dMax, oMin, oMax,RN)
    use mo_kind, only: dp,i4
    implicit none

    real(DP),    intent(IN)   :: old, dMax, oMin,oMax
    real(DP),    intent(IN)   :: RN

    ! local variables
    real(DP)                  :: oi, ox, delta, dMaxScal
    integer(I4)               :: iDigit,iszero             ! were intent IN before
    !
    iszero = 1_i4
    iDigit = 8_i4

    oi=oMin
    ox=oMax
    ! scaling (new)
    dMaxScal = dMax * ABS(oMax - oMin)

    if( ge(oi , ox) ) then
       parGen_anneal_dp = (oi+ox)/2.0_DP
    else
       if ( (old - dMaxScal) > oi)  oi = old - dMaxScal
       if ( (old + dMaxScal) < ox)  ox = old + dMaxScal
       delta = (oi + RN * (ox - oi) ) - old
       if ( delta > dMaxScal) delta =  dMaxScal
       if ( delta <-dMaxScal) delta = -dMaxScal
       parGen_anneal_dp = old + dChange_dp(delta,iDigit,isZero)
    endif

    ! Parameter is bounded between Max and Min.
    ! Correction from Kumar and Luis

    if (parGen_anneal_dp < oMin) then
       parGen_anneal_dp = oMin
    elseif(parGen_anneal_dp > oMax)then
       parGen_anneal_dp = oMax
    end if
  end function parGen_anneal_dp

  real(DP) function parGen_dds_dp(old, perturb, oMin, oMax, RN)
    ! PURPOSE:
    !    Perturb variable using a standard normal random number RN
    !    and reflecting at variable bounds if necessary
    !    Tolson et al. (2007): STEP 4
    !
    implicit none

    real(DP),    intent(IN)   :: old, perturb, oMin, oMax
    real(DP),    intent(IN)   :: RN

    ! generating new value
    parGen_dds_dp = old + perturb * (oMax - oMin) * RN

    ! reflect one time and set to boundary value
    if (parGen_dds_dp .lt. oMin) then
       parGen_dds_dp = oMin + ( oMin - parGen_dds_dp )
       if (parGen_dds_dp .gt. oMax) parGen_dds_dp = oMin
    end if

    if (parGen_dds_dp .gt. oMax) then
       parGen_dds_dp = oMax - ( parGen_dds_dp - oMax )
       if (parGen_dds_dp .lt. oMin) parGen_dds_dp = oMax
    end if

  end function parGen_dds_dp

  real(DP) function  dChange_dp(delta,iDigit,isZero)
    use mo_kind, only: i4,i8, dp
    implicit none

    integer(I4), intent(IN)   :: iDigit,isZero
    real(DP),    intent(IN)   :: delta

    ! local variables
    integer(I4)               :: ioszt
    integer(I8)               :: iDelta

    ioszt = 10**iDigit
    iDelta = int( delta * real(ioszt,dp),i8 )
    if(isZero == 1_i4) then
       if(iDelta == 0_i8) then
          if(delta < 0_i4) then
             iDelta = -1_i8
          else
             iDelta = 1_i8
          endif
       endif
    endif
    dChange_dp=real(iDelta,dp)/real(ioszt,dp)
  end function  dChange_dp

  subroutine generate_neighborhood_weight_dp(truepara, cum_weight, save_state_xor, iTotalCounter, &
       nITERmax, neighborhood)
    ! PURPOSE:
    !    generates a new neighborhood
    !    Tolson et al. (2007): STEP 3
    !
    integer(i4), dimension(:),                intent(in)    :: truepara
    real(dp),    dimension(:),                intent(in)    :: cum_weight
    integer(i8), dimension(n_save_state),     intent(inout) :: save_state_xor
    integer(i4),                              intent(in)    :: iTotalCounter
    integer(i4),                              intent(in)    :: nITERmax
    logical,     dimension(size(cum_weight)), intent(out)   :: neighborhood

    ! local variables
    integer(i4)                           :: ipar, npar
    integer(i4)                           :: iSize, size_neighbor
    real(dp)                              :: prob
    real(dp)                              :: rn
    real(dp)                              :: weight, norm
    real(dp), dimension(size(cum_weight)) :: local_cum_weight

    prob = 1.0_dp - Log(real(iTotalCounter,dp)) / Log(real(nITERmax,dp))
    neighborhood = .false.

    npar = size(cum_weight)
    local_cum_weight = cum_weight

    ! How many parameters will be selected for neighborhood?
    size_neighbor = 0_i4
    do ipar=1, size(truepara)
       call xor4096(0_i8, rn, save_state=save_state_xor)
       if (rn < prob) then
          size_neighbor = size_neighbor + 1_i4
       end if
    end do
    ! at least one...
    size_neighbor = max(1_i4, size_neighbor)

    ! Which parameter will be used for neighborhood?
    do iSize = 1, size_neighbor
       ! (1) generate RN
       call xor4096(0_i8, rn, save_state=save_state_xor)
       !
       ! (2) find location <iPar> in cummulative distribution function
       iPar = 1_i4
       do while (local_cum_weight(iPar) .lt. rn)
          iPar = iPar + 1_i4
       end do
       !
       ! (3) add parameter to neighborhood
       neighborhood(iPar) = .true.
       !
       ! (4) recalculate cummulative distribution function
       ! (4.1) Which weight had iPar?
       if (iPar .gt. 1_i4) then
          weight = local_cum_weight(iPar) - local_cum_weight(iPar-1)
       else
          weight = local_cum_weight(1)
       end if
       ! (4.2) Substract this weight from cummulative array starting at iPar
       local_cum_weight(iPar:nPar) = local_cum_weight(iPar:nPar) - weight
       ! (4.3) Renormalize cummulative weight to one
       if ( count(neighborhood) .lt. size(truepara) ) then
          norm = 1.0_dp/local_cum_weight(nPar)
       else
          norm = 1.0_dp
       end if
       local_cum_weight(:) = local_cum_weight(:) * norm
       !
    end do

  end subroutine generate_neighborhood_weight_dp

END MODULE mo_anneal
