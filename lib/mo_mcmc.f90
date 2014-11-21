!> \file mo_mcmc.f90

!> \brief This module is Monte Carlo Markov Chain sampling of a posterior parameter distribution.

!> \details This module is Monte Carlo Markov Chain sampling of a posterior parameter distribution.

!> \authors Maren Goehler, Juliane Mai
!> \date Aug 2012

MODULE mo_mcmc

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Juliane Mai.

  ! Copyright 2012 Juliane Mai, Maren Goehler

  USE mo_kind,    only: i4, i8, dp
  USE mo_xor4096, only: xor4096, xor4096g, get_timeseed, n_save_state
  USE mo_append,  only: append
  USE mo_moment,  only: stddev
  !$ USE omp_lib,    only: OMP_GET_NUM_THREADS
  use mo_ncwrite, only: dump_netcdf

  IMPLICIT NONE

  PUBLIC :: mcmc          ! Sample posterior parameter distribution

  !-----------------------------------------------------------------------------------------------
  !

  ! NAME
  !>        \brief This module is Monte Carlo Markov Chain sampling of a posterior parameter distribution
  !

  ! PURPOSE
  !>        \details Sample posterior parameter distribution with Metropolis Algorithm.\n
  !>                 This sampling is performed in two steps, i.e. the burn-in phase for adjusting 
  !>                 model dependent parameters for the second step which is the proper 
  !>                 sampling using the Metropolis Hastings Algorithm.\n\n
  !
  !>     <b>1. BURN IN PHASE: FIND THE OPTIMAL STEP SIZE </b>\n\n
  !
  !>          <b><i>Purpose:</i></b>\n
  !>            Find optimal stepsize for each parameter such that the
  !>            acceptance ratio converges to a value around 0.3.\n
  !>
  !>          <b><i>Important variables:</i></b>\n
  !>
  !>            <b> Variable           |  Description                                                     </b> 
  !>            ---------------------- | -----------------------------------------------------------------------
  !>            burnin_iter            | length of markov chain performed to calculate acceptance ratio\n
  !>            acceptance ratio       | ratio between accepted jumps and all trials (LEN)\n
  !>            acceptance multiplier  | stepsize of a parameter is multiplied with this value when jump is 
  !>                                     accepted (initial : 1.01)
  !>            rejection multiplier   | stepsize of a parameter is multiplied with this value when jump is rejected 
  !>                                     (initial : 0.99 and will never be changed)
  !>            stepsize               | a new parameter value is chosen based on a uniform distribution 
  !>                                     pnew_i = pold_i + Unif(-stepsize_i, stepsize_i) (initial : stepsize_i = 1.0 for all i)
  !
  !>          \n
  !>          <b><i>Algorithm:</i></b>\n
  !>          <ol>
  !>          <li>start a new markov chain of length burnin_iter with initial parameter set is the OPTIMAL one</i>\n
  !>                - select a set of parameters to change:\n
  !>                  * accurate --> choose one parameter,\n
  !>                  * comput. efficient --> choose all parameters,\n
  !>                  * moderate accurate & efficient --> choose half of the parameters\n
  !>                - change parameter(s) based on their stepsize\n
  !>                - decide whether changed parameter set is accepted or rejected:\n
  !>                  * \f$ \mathrm{odds Ratio} = \frac{\mathrm{likelihood}(p_{new})}{\mathrm{likelihood}(p_{old})} \f$\n
  !>                  * random number r = Uniform[0,1]\n
  !>                  * odds Ratio > 0 --> positive accept\n
  !>                    odds Ratio > r --> negative accept\n
  !>                    odds Ratio < r --> reject\n
  !>                - adapt stepsize of parameters changed:\n
  !>                  * accepted step: stepsize_i = stepsize_i * accptance multiplier\n
  !>                  * rejected step: stepsize_i = stepsize_i * rejection multiplier\n
  !>                - if step is accepted: for all changed parameter(s) change stepsize\n
  !
  !>          <li>calculate acceptance ratio of the Markov Chain</i>\n
  !
  !>          <li>adjust acceptance multiplier acc_mult and
  !>                store good ratios in history list\n
  !>                - acceptance ratio < 0.23 --> acc_mult = acc_mult * 0.99\n
  !>                                              delete history list\n
  !>                - acceptance ratio > 0.44 --> acc_mult = acc_mult * 1.01\n
  !>                                              delete history list\n
  !>                - 0.23 < acceptance ratio < 0.44\n
  !>                                              add acceptance ratio to history list\n
  !
  !>          <li>check if already 10 values are stored in history list and
  !>                if they have converged to a value above 0.3\n
  !>                ( mean above 0.3 and variance less \f$\sqrt{1/12*0.05^2}\f$ = Variance
  !>                  of uniform [acc_ratio +/- 2.5%] )\n
  !>                - if check is positive abort and save stepsizes\n
  !>                  else goto (1)\n\n
  !>          </ol>
  !
  !>
  !>     <b>2. MONTE CARLO MARKOV CHAIN: SAMPLE POSTERIOR DISTRIBUTION OF PARAMETER</b>\n\n
  !>          <b><i>Purpose:</i></b>\n
  !>            use the previous adapted stepsizes and perform ONE monte carlo markov chain\n
  !>            the accepted parameter sets show the posterior distribution of parameters\n
  !>
  !>          <b><i>Important variables:</i></b>\n
  !>            <b> Variable           |  Description                                                     </b> 
  !>            ---------------------- | -----------------------------------------------------------------------
  !>            iter_mcmc              | length of the markov chain (>> iter_burnin)\n
  !>            stepsize               | a new parameter value is chosen based on a uniform distribution 
  !>                                     pnew_i = pold_i + Unif(-stepsize_i, stepsize_i) use stepsizes of the burn-in (1)\n
  !
  !>          \n
  !>          <b><i>Algorithm:</i></b>\n
  !>          <ol>
  !>          <li> select a set of parameters to change\n
  !>               * accurate --> choose one parameter,\n
  !>               * comput. efficient --> choose all parameters,\n
  !>               * moderate accurate & efficient --> choose half of the parameters\n
  !>          <li> change parameter(s) based on their stepsize\n
  !>          <li> decide whether changed parameter set is accepted or rejected:\n
  !>               * \f$ \mathrm{odds Ratio} = \frac{\mathrm{likelihood}(p_{new})}{\mathrm{likelihood}(p_{old})} \f$\n
  !>               * random number r = Uniform[0,1]\n
  !>               * odds Ratio > 0 --> positive accept\n
  !>                 odds Ratio > r --> negative accept\n
  !>                 odds Ratio < r --> reject\n
  !>          <li> if step is accepted: save parameter set\n
  !>          <li> goto (1)\n
  !>          </ol>

  !     CALLING SEQUENCE
  !         call mcmc(  likelihood, stddev_function, para, rangePar, mcmc_paras, burnin_paras, &
  !                     seed_in=seeds, printflag_in=printflag, maskpara_in=maskpara,           &
  !                     tmp_file=tmp_file, loglike_in=loglike,                                 &
  !                     ParaSelectMode_in=ParaSelectMode,                                      &
  !                     iter_burnin_in=iter_burnin, iter_mcmc_in=iter_mcmc,                    &
  !                     chains_in=chains, stepsize_in=stepsize)

  !     INTENT(IN)
  !>        \param[in]  "real(dp) :: likelihood(x,sigma,stddev_new,likeli_new)"              
  !>                                                                   Interface Function which calculates likelihood 
  !>                                                                   of given parameter set x and given standard deviation sigma 
  !>                                                                   and returns optionally the standard deviation stddev_new
  !>                                                                   of the errors using x and 
  !>                                                                   likelihood likeli_new using stddev_new

  !>        \param[in]  "real(dp) :: para(:)"                          Inital parameter set (should be GOOD approximation 
  !>                                                                   of best parameter set)

  !>        \param[in]  "real(dp) :: rangePar(size(para),2)"           Min/max range of parameters

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out]  "real(dp), allocatable :: mcmc_paras(:,:)"    Parameter sets sampled in proper MCMC part of algorithm
  !>        \param[out]  "real(dp), allocatable :: burnin_paras(:,:)"  Parameter sets sampled during burn-in part of algorithm

  !     INTENT(IN), OPTIONAL
  !>        \param[in]  "integer(i8),      optional :: seed_in"        User seed to initialise the random number generator \n
  !>                                                                   (default: none --> initialized with timeseed)

  !>        \param[in]  "logical,          optional :: printflag_in"   Print of output on command line \n
  !>                                                                   (default: .False.)

  !>        \param[in]  "logical,          optional :: maskpara_in(size(para))"   
  !>                                                                   Parameter will be sampled (.True.) or not (.False.) \n
  !>                                                                   (default: .True.)

  !>        \param[in]  "character(len=*), optional :: tmp_file"       filename for temporal data saving: \n
  !>                                                                   every iter_mcmc_in iterations parameter sets are 
  !>                                                                   appended to this file \n
  !>                                                                   the number of the chain will be prepended to filename \n
  !>                                                                   output format: netcdf \n
  !>                                                                   (default: no file writing)

  !>        \param[in]  "logical,          optional :: loglike_in"     true if loglikelihood function is given instead of 
  !>                                                                   likelihood function \n
  !>                                                                   (default: .false.)

  !>        \param[in]  "integer(i4),      optional :: ParaSelectMode_in"  
  !>                                                                   How many parameters will be changed at once?
  !>                                                                   - half of the parameter  --> 1_i4
  !>                                                                   - only one parameter     --> 2_i4
  !>                                                                   - all parameter          --> 3_i4
  !>                                                                   (default: 2_i4)

  !>        \param[in]  "integer(i4),      optional :: iter_burnin_in" Length of Markov chains of initial burn-in part \n
  !>                                                                   (default: Max(250, 200*count(maskpara)) )

  !>        \param[in]  "integer(i4),      optional :: iter_mcmc_in"   Length of Markov chains of proper MCMC part \n
  !>                                                                   (default: 1000 * count(maskpara) )

  !>        \param[in]  "integer(i4),      optional :: chains_in"      number of parallel mcmc chains \n
  !>                                                                   (default: 5_i4)

  !>        \param[in]  "real(dp), DIMENSION(size(para,1)), optional :: stepsize_in"
  !>                                                                   stepsize for each parameter \n
  !>                                                                   if given burn-in is discarded \n
  !>                                                                   (default: none --> adjusted in burn-in)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         \note Likelihood has to be defined as a function interface
  !         -> see also example in test directory

  !     EXAMPLE
  !         -> see also example in test directory

  !     LITERATURE
  !         Gelman et. al (1995)
  !             Baysian Data Analyis, Chapman & Hall.
  !         Gelman et. al (1997)
  !             Efficient Metropolis Jumping Rules, Baysian Statistics 5, pp. 599-607
  !         Tautenhahn et. al (2012)
  !             Beyond distance-invariant survival in inverse recruitment modeling:
  !             A case study in Siberian Pinus sylvestris forests, Ecological Modelling, 233,
  !             90-103. doi:10.1016/j.ecolmodel.2012.03.009.

  !     HISTORY
  !         Written  Goehler,    Sep. 2012 : Created using copy of Simulated Annealing
  !                                          - constant temperature T
  !                                          - burn-in for stepsize adaption
  !                                          - acceptance/ rejection multiplier
  !         Modified Mai,        Sep. 2012 : Cleaning code and introduce likelihood
  !                                          - likelihood instead of objective function
  !                                          - odds ratio
  !                                          - convergence criteria for burn-in
  !                                          - different modes of parameter selection
  !                                          OpenMP for chains of MCMC
  !                                          optional file for temporal output
  !                  Mai,        Nov. 2012 : Temporary file writing as NetCDF
  !                  Mai,        Aug. 2013 : New likelihood interface to reduce number of function evaluations
  !                                          Only one seed has to be given
  INTERFACE mcmc
     MODULE PROCEDURE mcmc_dp
  END INTERFACE mcmc

  !-----------------------------------------------------------------------------------------------

  PRIVATE

  !-----------------------------------------------------------------------------------------------

CONTAINS

  !-----------------------------------------------------------------------------------------------

  SUBROUTINE mcmc_dp(likelihood, para, rangePar, &   ! obligatory IN
       mcmc_paras, burnin_paras,                 &   ! obligatory OUT
       seed_in, printflag_in, maskpara_in,       &   ! optional IN
       tmp_file,                                 &   ! optional IN : filename for temporal output of 
       !                                             !               MCMC parasets
       loglike_in,                               &   ! optional IN : true if loglikelihood is given
       ParaSelectMode_in,                        &   ! optional IN : (=1) half, (=2) one, (=3) all
       iter_burnin_in,                           &   ! optional IN : markov length of (1) burn-in
       iter_mcmc_in,                             &   ! optional IN : markov length of (2) mcmc
       chains_in,                                &   ! optional IN : number of parallel chains of MCMC
       stepsize_in)                                  ! optional_IN : stepsize for each param. (no burn-in)

    IMPLICIT NONE

    INTERFACE
       FUNCTION likelihood(paraset,sigma,stddev_new,likeli_new)
         ! calculates the likelihood function at a certain parameter set paraset
         use mo_kind
         REAL(DP), DIMENSION(:), INTENT(IN)            :: paraset
         REAL(DP),               INTENT(IN)            :: sigma
         REAL(DP),               INTENT(OUT), OPTIONAL :: stddev_new       ! standard deviation of errors using paraset
         REAL(DP),               INTENT(OUT), OPTIONAL :: likeli_new       ! likelihood using stddev_new, 
         !                                                                 ! i.e. using new parameter set
         REAL(DP)                                      :: likelihood
       END FUNCTION likelihood
    END INTERFACE

    REAL(DP),    DIMENSION(:,:), INTENT(IN)    :: rangePar           ! range for each parameter
    REAL(DP),    DIMENSION(:),   INTENT(IN)    :: para               ! initial parameter i
    INTEGER(I8), OPTIONAL,       INTENT(IN)    :: seed_in           ! Seeds of random numbers: dim1=chains, dim2=3
    !                                                                ! (DEFAULT: Get_timeseed)
    LOGICAL,     OPTIONAL,       INTENT(IN)    :: printflag_in       ! If command line output is written (.true.)
    !                                                                ! (DEFAULT: .false.)
    LOGICAL,     OPTIONAL, DIMENSION(size(para,1)), &
         INTENT(IN)    :: maskpara_in                                ! true if parameter will be optimized
    !                                                                ! false if parameter is discarded in optimization
    !                                                                ! DEFAULT = .true.
    LOGICAL,     OPTIONAL,       INTENT(IN)    :: loglike_in         ! true if loglikelihood is given instead of likelihood
    !                                                                ! DEFAULT: .false.
    INTEGER(I4), OPTIONAL,       INTENT(IN)    :: ParaSelectMode_in  ! how many parameters changed per iteration:
    !                                                                !   1: half of the parameters
    !                                                                !   2: only one (DEFAULT)
    !                                                                !   3: all
    INTEGER(I4), OPTIONAL,       INTENT(IN)    :: iter_burnin_in     ! # iterations before checking ac_ratio
    !                                                                ! DEFAULT= MIN(250, 20_i4*n)
    INTEGER(I4), OPTIONAL,       INTENT(IN)    :: iter_mcmc_in       ! # iterations per chain for posterior sampling
    !                                                                ! DEFAULT= 1000_i4*n
    INTEGER(I4), OPTIONAL,       INTENT(IN)    :: chains_in          ! number of parallel mcmc chains
    !                                                                ! DEFAULT= 5_i4
    REAL(DP),    OPTIONAL, DIMENSION(size(para,1)), &
         INTENT(IN)    :: stepsize_in                                ! stepsize for each parameter
    !                                                                ! if given burn-in is discarded
    CHARACTER(len=*), OPTIONAL,  INTENT(IN)    :: tmp_file           ! filename for temporal data saving: every iter_mcmc_in 
    !                                                                ! iterations parameter sets are appended to this file
    !                                                                ! the number of the chain will be prepended to filename
    !                                                                ! DEFAULT= no file writing
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)    :: mcmc_paras
    !                                                                ! array to save para values of MCMC runs,
    !                                                                ! dim1=sets of all chains, dim2=paras
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE, INTENT(OUT)  :: burnin_paras
    !                                                                ! array to save para values of Burn in run


    ! local variables
    INTEGER(I4)                            :: n           ! Number of parameters
    LOGICAL                                :: printflag   ! If command line output is written
    LOGICAL,     DIMENSION(size(para,1))   :: maskpara    ! true if parameter will be optimized
    INTEGER(I4), DIMENSION(:), ALLOCATABLE :: truepara    ! indexes of parameters to be optimized
    LOGICAL                                :: loglike     ! if loglikelihood is given

    ! for random numbers
    INTEGER(I8), dimension(:,:), allocatable  :: seeds             ! Seeds of random numbers: dim1=chains, dim2=3
    REAL(DP)                                  :: RN1, RN2, RN3     ! Random numbers
    integer(I8), dimension(:,:), allocatable  :: save_state_1      ! optional argument for restarting RN stream 1: 
    integer(I8), dimension(:,:), allocatable  :: save_state_2      ! optional argument for restarting RN stream 2: 
    integer(I8), dimension(:,:), allocatable  :: save_state_3      ! optional argument for restarting RN stream 3: 
    !                                                              !     dim1=chains, dim2=n_save_state

    ! Dummies
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: tmp
    integer(I4)                             :: idummy
    logical                                 :: oddsSwitch1, oddsSwitch2
    character(100)                          :: str
    character(200)                          :: outputfile
    integer(i4)                             :: slash_pos
    integer(i4)                             :: len_filename
    character(200)                          :: filename
    character(200)                          :: path 

    ! FOR BURN-IN AND MCMC
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE        :: mcmc_paras_3d     ! array to save para values of MCMC runs,
    !                                                                   ! dim1=sets, dim2=paras, dim3=chain
    integer(I4)                                    :: i
    integer(I4), dimension(:), allocatable         :: Ipos, Ineg
    logical                                        :: iStop
    INTEGER(I4)                                    :: ParaSelectMode    ! how many parameters are changed per jump
    INTEGER(I4)                                    :: iter_burnin       ! number fo iterations before checking ac_ratio
    INTEGER(I4)                                    :: iter_mcmc         ! number fo iterations for posterior sampling
    INTEGER(I4)                                    :: chains            ! number of parallel mcmc chains
    INTEGER(I4)                                    :: chain             ! counter for chains, 1...chains
    REAL(DP)                                       :: accMult           ! acceptance multiplier for stepsize
    REAL(DP)                                       :: rejMult           ! rejection multiplier for stepsize
    LOGICAL,DIMENSION(size(para,1))                :: ChangePara        ! logical array to switch if parameter is changed
    INTEGER(I4)                                    :: trial             ! number of trials for burn-in
    REAL(DP),DIMENSION(size(para,1))               :: stepsize          ! stepsize adjusted by burn-in and used by mcmc
    REAL(DP),DIMENSION(size(para,1))               :: paraold           ! old parameter set
    REAL(DP),DIMENSION(size(para,1))               :: paranew           ! new parameter set
    REAL(DP),DIMENSION(size(para,1))               :: parabest          ! best parameter set overall
    REAL(DP)                                       :: stddev_data       ! approximated stddev of data for best paraset found
    LOGICAL                                        :: parabestChanged   ! if better parameter set was found during burn-in
    REAL(DP)                                       :: likeliold         ! likelihood of old parameter set
    REAL(DP)                                       :: likelinew         ! likelihood of new parameter set
    REAL(DP)                                       :: likelibest        ! likelihood of best parameter set overall
    REAL(DP)                                       :: stddev_new        ! standard deviation of errors with current parameter set
    REAL(DP)                                       :: likeli_new        ! likelihood using stddev_new instead of stddev_data
    INTEGER(I4)                                    :: markov            ! counter for markov chain
    REAL(DP), DIMENSION(:,:), ALLOCATABLE          :: burnin_paras_part ! accepted parameter sets of one MC in burn-in
    REAL(DP)                                       :: oddsRatio         ! ratio of likelihoods = likelinew/likeliold
    REAL(DP), dimension(:), allocatable            :: accRatio          ! acceptance ratio = (pos/neg) accepted/iter_burnin
    REAL(DP)                                       :: accratio_stddev   ! stddev of accRatios in history
    REAL(DP), DIMENSION(:), ALLOCATABLE            :: history_accratio  ! history of 'good' acceptance ratios
    INTEGER(I4)                                    :: accratio_n        ! number of 'good' acceptance ratios
    !REAL(DP), DIMENSION(:,:), ALLOCATABLE          :: history_stepsize  ! history of 'good' stepsizes

    ! for checking convergence of MCMC
    real(dp), allocatable, dimension(:)   :: sqrtR
    real(dp), allocatable, dimension(:)   :: vDotJ
    real(dp), allocatable, dimension(:)   :: s2
    real(dp)                              :: vDotDot
    real(dp)                              :: B
    real(dp)                              :: W
    integer(i4)                           :: n_start, n_end, iPar
    LOGICAL                               :: converged         ! if MCMC already converged

    ! for OMP
    !$  integer(i4)                           :: n_threads

    ! CHECKING OPTIONALS

    if (present(loglike_in)) then
       loglike = loglike_in
    else
       loglike = .false.
    end if

    if (present(maskpara_in)) then
       if (count(maskpara_in) .eq. 0_i4) then
          stop 'Input argument maskpara: At least one element has to be true'
       else
          maskpara = maskpara_in
       end if
    else
       maskpara = .true.
    endif

    allocate ( truepara(count(maskpara)) )
    idummy = 0_i4
    do i=1,size(para,1)
       if ( maskpara(i) ) then
          idummy = idummy+1_i4
          truepara(idummy) = i
       end if
    end do

    n = size(truepara,1)

    if (present(ParaSelectMode_in)) then
       ParaSelectMode = ParaSelectMode_in
    else
       ! change only one parameter per jump
       ParaSelectMode = 2_i4
    end if

    ! after how many iterations do we compute ac_ratio???
    if (present(iter_burnin_in)) then
       if (iter_burnin_in .le. 0_i4)  then
          stop 'Input argument iter_burn_in must be greater than  0!'
       else
          iter_burnin = iter_burnin_in
       end if
    else
       iter_burnin = Max(250_i4, 1000_i4*n)
    endif

    ! how many iterations ('jumps') are performed in MCMC
    ! iter_mcmc_in is handled later properly (after acceptance ratio of burn_in is known)
    if (present(iter_mcmc_in)) then
       if (iter_mcmc_in .le. 0_i4)  then
          stop 'Input argument iter_mcmc must be greater than  0!'
       else
          iter_mcmc = iter_mcmc_in
       end if
    else
       iter_mcmc = 1000_i4 * n
    endif

    if (present(chains_in)) then
       if (chains_in .lt. 2_i4)  then
          stop 'Input argument chains must be at least 2!'
       end if
       chains = chains_in
    else
       chains = 5_i4
    endif

    if (present(stepsize_in)) then
       stepsize   = stepsize_in
    end if

    if (present(printflag_in)) then
       printflag = printflag_in
    else
       printflag = .false.
    endif

    !$  write(*,*) '--------------------------------------------------'
    !$  write(*,*) ' This program is parallel.'
    !$OMP parallel
    !$    n_threads = OMP_GET_NUM_THREADS()
    !$OMP end parallel
    !$  write(*,*) ' ',chains,' MCMC chains will run in ',n_threads,' threads'
    !$  write(*,*) '--------------------------------------------------'

    if (printflag) then
       print*, 'Following parameters will be sampled with MCMC: ',truepara
    end if

    allocate(seeds(chains,3))
    allocate(save_state_1(chains,n_save_state))
    allocate(save_state_2(chains,n_save_state))
    allocate(save_state_3(chains,n_save_state))
    allocate(Ipos(chains), Ineg(chains), accRatio(chains))
    allocate(vDotJ(chains), s2(chains))
    allocate(sqrtR(size(para)))

    if (present(seed_in)) then
       seeds(1,:) = (/ seed_in, seed_in + 1000_i8, seed_in + 2000_i8 /)
    else
       ! Seeds depend on actual time
       call get_timeseed(seeds(1,:))
    endif
    do chain=2,chains
       seeds(chain,:) = seeds(chain-1_i4,:) + 3000_i8
    end do

    do chain=1,chains
       call xor4096(seeds(chain,1),  RN1, save_state=save_state_1(chain,:))
       call xor4096(seeds(chain,2),  RN2, save_state=save_state_2(chain,:))
       call xor4096g(seeds(chain,3), RN3, save_state=save_state_3(chain,:))
    end do
    seeds = 0_i8

    parabest   = para
    ! print*, parabest

    ! initialize likelihood and sigma
    likelibest  = likelihood(parabest,1.0_dp,stddev_new=stddev_new,likeli_new=likeli_new)
    likelibest  = likeli_new
    stddev_data = stddev_new

    !----------------------------------------------------------------------
    ! (1) BURN IN
    !----------------------------------------------------------------------

    if (.not. present(stepsize_in)) then
       if (printflag) then
          print*, ''
          print*, '--------------------------------------------------'
          print*, 'Starting Burn-In   (iter_burnin = ',iter_burnin,')'
          print*, '--------------------------------------------------'
          print*, ''
       end if

       ! INITIALIZATION

       ! probably too large, but large enough to store values of one markovchain
       allocate(burnin_paras_part(iter_burnin,size(para)))

       parabestChanged = .true.

       ! -------------------------------------------------------------------------------------
       ! restarts if a better parameter set was found in between
       ! -------------------------------------------------------------------------------------
       betterParaFound: do while (parabestChanged)

          if (allocated(burnin_paras))     deallocate(burnin_paras)
          if (allocated(history_accRatio)) deallocate(history_accRatio)

          parabestChanged = .false.
          stepsize   = 1.0_dp
          trial      = 1_i4
          iStop      = .false.
          accMult    = 1.01_dp
          rejMult    = 0.99_dp
          !stddev_data = stddev_function(parabest)
          !likelibest  = likelihood(parabest,stddev_data,stddev_new=stddev_new,likeli_new=likeli_new)
          !likelibest  = likeli_new
          !stddev_data = stddev_new

          if (printflag) then
             print*, ' '
             print*, 'Restart Burn-In with new approximation of std.dev. of data: '
             print*, '   parabest   = ', parabest
             print*, '   stddev     = ', stddev_data
             print*, '   likelihood = ', likelibest
             print*, '   ssq        = ', -2.0_dp * stddev_data**2 * likelibest
             print*, ' '
          end if


          ! ----------------------------------------------------------------------------------
          ! repeats until convergence of acceptance ratio or better parameter set was found
          ! ----------------------------------------------------------------------------------
          convergedBURNIN: do while ( (.not. iStop) .and. (.not. parabestChanged) )

             Ipos = 0_i4   ! positive accepted
             Ineg = 0_i4   ! negative accepted
             paraold   = parabest
             likeliold = likelibest !likelihood(paraold,stddev_data)

             ! -------------------------------------------------------------------------------
             ! do a short-cut MCMC
             ! -------------------------------------------------------------------------------
             markovchain: do markov=1, iter_burnin

                ! (A) Generate new parameter set
                ChangePara = .false.
                paranew    = paraold
                ! using RN from chain #1
                call GenerateNewParameterset_dp(ParaSelectMode, paraold, truepara, rangePar, stepsize, &
                     save_state_2(1,:),&
                     save_state_3(1,:),&
                     paranew,ChangePara)

                ! (B) new likelihood
                likelinew = likelihood(paranew,stddev_data,stddev_new=stddev_new,likeli_new=likeli_new)

                oddsSwitch1 = .false.
                if (loglike) then
                   oddsRatio = likelinew-likeliold
                   if (oddsRatio .gt. 0.0_dp) oddsSwitch1 = .true.
                else
                   oddsRatio = likelinew/likeliold
                   if (oddsRatio .gt. 1.0_dp) oddsSwitch1 = .true.
                end if
! print*, 'oddsRatio = ',oddsRatio

                ! (C) Accept or Reject?
                If (oddsSwitch1) then

                   ! positive accept
                   Ipos(1)   = Ipos(1) + 1_i4
                   paraold   = paranew
                   likeliold = likelinew
                   where (changePara)
                      stepsize = stepsize * accMult
                   end where
                   if (likelinew .gt. likelibest) then
                      parabest        = paranew
                      likeliold       = likeli_new !JM
                      ! Here the sigma is reset!
                      likelibest      = likeli_new
                      stddev_data     = stddev_new
                      !
                      parabestChanged = .true.
                      print*, ''
                      print*, 'best para changed: ',paranew
                      print*, 'likelihood new:    ',likelinew
                      print*, ''
                   end if
                   burnin_paras_part(Ipos(1)+Ineg(1),:) = paranew(:)
! print*, 'positive accept'

                else

                   call xor4096(seeds(1,1), RN1, save_state=save_state_1(1,:))
                   oddsSwitch2 = .false.
                   if (loglike) then
                      if (oddsRatio .lt. -700.0_dp) oddsRatio = -700.0_dp     ! to avoid underflow
                      if (rn1 .lt. exp(oddsRatio)) oddsSwitch2 = .true.
                   else
                      if (rn1 .lt. oddsRatio)      oddsSwitch2 = .true.
                   end if

                   If ( oddsSwitch2 ) then

                      ! negative accept
                      Ineg(1)   = Ineg(1) + 1_i4
                      paraold   = paranew
                      likeliold = likelinew
                      ! stddev_data     = stddev_new !JM
                      where (changePara)
                         stepsize = stepsize * accMult
                      end where
                      burnin_paras_part(Ipos(1)+Ineg(1),:) = paranew(:) 
! print*, 'negative accept'

                   else

                      ! reject
                      where (changePara)
                         stepsize = stepsize * rejMult
                      end where
! print*, 'reject'

                   end if

                end if

! print*, ''

             end do markovchain

             accRatio(1) = real(Ipos(1) + Ineg(1),dp) / real(iter_burnin,dp)
             if (printflag) then
                write(str,'(A,I03,A)') '(A7,I4,A15,F5.3,A17,',size(para,1),'(E9.2,1X),A1)'
                write(*,str) 'trial #',trial,'   acc_ratio = ',accRatio(1),'     (stepsize = ',stepsize,')'
             end if

             if (Ipos(1)+Ineg(1) .gt. 0_i4) then
                call append(burnin_paras, burnin_paras_part(1:Ipos(1)+Ineg(1),:))
             end if

             ! adjust acceptance multiplier
             if (accRatio(1) .lt. 0.234_dp) accMult = accMult * 0.99_dp
             if (accRatio(1) .gt. 0.441_dp) accMult = accMult * 1.01_dp

             ! store good accRatios in history and delete complete history if bad one appears
             if (accRatio(1) .lt. 0.234_dp .or. accRatio(1) .gt. 0.441_dp) then
                if( allocated(history_accRatio) ) deallocate(history_accRatio)
             else
                call append(history_accRatio, accRatio(1))
             end if

             ! if in history more than 10 values, check for mean and variance
             if ( allocated(history_accRatio) ) then
                accRatio_n = size(history_accRatio,1)
                if ( accRatio_n .ge. 10_i4 ) then
                   idummy = accRatio_n-9_i4
                   accRatio_stddev = stddev( history_accRatio(idummy:accRatio_n) )

                   ! Check of Convergence
                   if ( (accRatio_stddev .lt. Sqrt( 1._dp/12._dp * 0.05_dp**2 )) ) then
                      iStop = .true.
                      if (printflag) then
                         print*, ''
                         print*, 'STOP BURN-IN with accaptence ratio of ', history_accRatio(accRatio_n)
                         print*, 'final stepsize: ',stepsize
                         if (parabestChanged) then
                            print*, 'better parameter set was found: ',parabest
                            print*, 'with likelihood: ',likelibest
                         end if
                      end if
                   end if

                end if
             end if

             trial = trial + 1_i4

          end do convergedBURNIN

       end do betterParaFound

    end if

    !----------------------------------------------------------------------
    ! (2) MCMC
    !----------------------------------------------------------------------

    if (printflag) then
       print*, ''
       print*, '--------------------------------------------------'
       print*, 'Starting MCMC  (chains = ',chains,',   iter_mcmc = ',iter_mcmc,')'
       print*, '--------------------------------------------------'
       print*, ''
    end if
    Ipos(:) = 0_i4   ! positive accepted
    Ineg(:) = 0_i4   ! negative accepted

    ! if all parameters converged: Sqrt(R_i) < 1.1 (see Gelman et. al: Baysian Data Analysis, p. 331ff
    converged = .False.

    convergedMCMC: do while (.not. converged)

       if (.not. allocated(mcmc_paras_3d)) then
          ! first mcmc iterations for all chains
          ! probably too large, but will be resized at the end
          allocate(mcmc_paras_3d(iter_mcmc,size(para),chains))
       else
          ! later mcmc iterations for all chains: iter_mcmc was increased
          idummy = Minval(Ipos+Ineg)

          ! resize mcmc_paras_3d
          allocate(tmp(idummy,size(para),chains))
          tmp(:,:,:) = mcmc_paras_3d(1:idummy,:,:)
          deallocate(mcmc_paras_3d)
          allocate(mcmc_paras_3d(iter_mcmc,size(para),chains))
          mcmc_paras_3d(1:idummy,:,:) = tmp(:,:,:)
          deallocate(tmp)
       end if

       !$OMP parallel default(shared) &
       !$OMP private(chain, paraold, paranew, likeliold, likelinew, oddsSwitch1, oddsSwitch2, RN1, oddsRatio, ChangePara)
       !$OMP do
       parallelchain: do chain=1, chains

          if (Ipos(chain)+Ineg(chain) .eq. 0_i4) then
             paraold = parabest
          else
             paraold = mcmc_paras_3d(Ipos(chain)+Ineg(chain),:,chain)
          end if
          likeliold = likelihood(paraold,stddev_data)

          markovchainMCMC: do

             ! (A) Generate new parameter set
             ChangePara = .false.
             paranew    = paraold

             call GenerateNewParameterset_dp(ParaSelectMode, paraold, truepara, rangePar, stepsize, &
                  save_state_2(chain,:),&
                  save_state_3(chain,:),&
                  paranew,ChangePara)

             ! (B) new likelihood
             likelinew = likelihood(paranew,stddev_data)
             oddsSwitch1 = .false.
             if (loglike) then
                oddsRatio = likelinew-likeliold
                if (oddsRatio .gt. 0.0_dp) oddsSwitch1 = .true.
             else
                oddsRatio = likelinew/likeliold
                if (oddsRatio .gt. 1.0_dp) oddsSwitch1 = .true.
             end if

             ! (C) Accept or Reject?
             If (oddsSwitch1) then

                ! positive accept
                Ipos(chain) = Ipos(chain) + 1_i4
                paraold     = paranew
                likeliold   = likelinew
                mcmc_paras_3d(Ipos(chain)+Ineg(chain),:,chain) = paranew(:)

                if ( printflag ) then
                   if ( (Ipos(chain)+Ineg(chain) .gt. 0_i4 ) .and. &
                        (mod(Ipos(chain)+Ineg(chain),10000_i4)) .eq. 0 ) then
                      print*, 'Chain ',chain,': Done ',Ipos(chain)+Ineg(chain),' samples ...'
                   end if
                end if

                if (likelinew .gt. likelibest) then
                   parabest   = paranew
                   likelibest = likelinew
                   if (printflag) then
                      print*, ''
                      print*, 'best para changed: ',paranew
                      print*, 'likelihood new:    ',likelinew
                      print*, ''
                   end if
                end if

             else

                call xor4096(seeds(chain,1), RN1, save_state=save_state_1(chain,:))
                oddsSwitch2 = .false.
                if (loglike) then
                   if (rn1 .lt. exp(oddsRatio)) oddsSwitch2 = .true.
                else
                   if (rn1 .lt. oddsRatio)      oddsSwitch2 = .true.
                end if

                If ( oddsSwitch2 ) then

                   ! negative accept
                   Ineg(chain) = Ineg(chain) + 1_i4
                   paraold     = paranew
                   likeliold   = likelinew
                   mcmc_paras_3d(Ipos(chain)+Ineg(chain),:,chain) = paranew(:)

                   if ( printflag ) then
                      if ( (Ipos(chain)+Ineg(chain) .gt. 0_i4 ) .and. &
                           (mod(Ipos(chain)+Ineg(chain),10000_i4)) .eq. 0 ) then
                         print*, 'Chain ',chain,': Done ',Ipos(chain)+Ineg(chain),' samples ...'
                      end if
                   end if

                end if

             end if

             ! enough samples were found
             if ( (Ipos(chain)+Ineg(chain) .gt. 0_i4 ) .and. &
                  (mod(Ipos(chain)+Ineg(chain),iter_mcmc) .eq. 0_i4 ) ) exit

          end do markovchainMCMC

       end do parallelchain
       !$OMP end do
       !$OMP end parallel

       ! write parameter sets to temporal file
       if (present(tmp_file)) then
          ! splitting into path and filename
          slash_pos    = index(tmp_file, '/', .true.)
          len_filename = len_trim(tmp_file)
          path         = tmp_file(1:slash_pos)
          filename     = tmp_file(slash_pos+1:len_filename)
          !
          do chain=1,chains
             write(str,*) chain
             write(outputfile,*) trim(adjustl(path)), trim(adjustl(str)), '_' , trim(adjustl(filename))
             if (present(iter_mcmc_in)) then
                allocate(tmp(iter_mcmc_in,size(para,1),1))
                tmp(:,:,1) = mcmc_paras_3d(iter_mcmc-iter_mcmc_in+1_i4:iter_mcmc,:,chain)
                if (iter_mcmc .ne. iter_mcmc_in) then
                   ! append
                   call dump_netcdf(trim(adjustl(outputfile)), tmp, append=.true.)
                else
                   ! first time of writing
                   call dump_netcdf(trim(adjustl(outputfile)), tmp)
                end if
                deallocate(tmp)
             else
                allocate(tmp(1000_i4*n,size(para,1),1))
                tmp(:,:,1) = mcmc_paras_3d(iter_mcmc-(1000_i4*n)+1_i4:iter_mcmc,:,chain)
                if (iter_mcmc .ne. 1000_i4*n) then
                   ! append
                   call dump_netcdf(trim(adjustl(outputfile)), tmp, append=.true.)
                else
                   ! first time of writing
                   call dump_netcdf(trim(adjustl(outputfile)), tmp)
                end if
                deallocate(tmp)
             end if
          end do
       end if

       ! test for convergence: Gelman et. al: Baysian Data Analysis, p. 331ff
       !    sqrt(R) = sqrt( ( (n-1)/n * W + 1/n * B ) / W ) < 1.1 for each parameter
       !    n ... last half of the sequence
       !    W ... within sequence variance
       !    B ... between chain variance

       if ( printflag ) then
          print*, ' '
          print*, 'Checking for convergence ....'
       end if
       sqrtR   = 0.0_dp
       n_end   = minval(Ipos+Ineg)
       n_start = n_end / 2_i4 + 1_i4

       do iPar=1, size(truepara)

          ! Between chain variances
          do chain=1, chains
             vDotJ(chain) = 1.0_dp / real(n_end-n_start+1_i4,dp) * &
                  sum(mcmc_paras_3d(n_start:n_end,truepara(iPar),chain))
          end do
          vDotDot = 1.0_dp / real(chains,dp) * sum(vDotJ)
          B = real(n_end -n_start + 1_i4,dp) / real(chains-1_i4,dp) * &
               sum((vDotJ - vDotDot)*(vDotJ - vDotDot))

          ! Within chain variances
          do chain=1, chains
             s2(chain) = 1.0_dp / real(n_end-n_start,dp) * &
                  sum((mcmc_paras_3d(n_start:n_end,truepara(iPar),chain)-vDotJ(chain))**2)
          end do
          W = 1.0_dp / real(chains,dp) * Sum(s2)

          ! ratio sqrtR
          if ( (w .lt. tiny(1.0_dp)) .and. (b .lt. tiny(1.0_dp))) then
             ! Mathematica says that this is the limit, if w and b both go to zero
             sqrtR(truepara(iPar)) = sqrt(real(n_end-n_start,dp) / real(n_end-n_start+1_i4,dp))
          else
             sqrtR(truepara(iPar)) = real(n_end-n_start,dp) / real(n_end-n_start+1_i4,dp) * W + &
                  1.0_dp / real(n_end-n_start+1_i4,dp) * B
             sqrtR(truepara(iPar)) = sqrt( sqrtR(truepara(iPar)) / W )
          end if
       end do
       if (printflag) then
          do i=1, size(para)
             if ( sqrtR(i) .gt. 1.1_dp ) then
                print*, '   sqrtR para #',i,' : ', sqrtR(i),'  <-- FAILED'
             else
                print*, '   sqrtR para #',i,' : ', sqrtR(i)
             end if
          end do
       end if
       if ( all(sqrtR .lt. 1.1_dp) ) then
          converged = .true.
          if ( printflag ) then
             print*, '   --> converged (all less than 1.1)'
             print*, ' '
          end if
       else
          ! increase number of iterations
          if (present(iter_mcmc_in)) then
             iter_mcmc = iter_mcmc + iter_mcmc_in
          else
             iter_mcmc = iter_mcmc + 1000_i4 * n
          endif

          if ( printflag ) then
             print*, '   --> not converged (not all less than 1.1)'
             print*, '       increasing iterations to ',iter_mcmc
             print*, ' '
          end if
       end if

    end do convergedMCMC

    ! reshape of mcmc_paras_3d: return only 2d matrix mcmc_paras
    allocate(mcmc_paras(size(mcmc_paras_3d,1)*size(mcmc_paras_3d,3),size(mcmc_paras_3d,2)))
    do chain=1, chains
       mcmc_paras( (chain-1_i4)*size(mcmc_paras_3d,1)+1_i4 : chain*size(mcmc_paras_3d,1),:) = mcmc_paras_3d(:,:,chain)
    end do

    RETURN
  END SUBROUTINE mcmc_dp

  !***************************************************
  !*               PRIVATE FUNCTIONS                 *
  !***************************************************

  !
  real(DP) function parGen_dp(old, dMax, oMin, oMax,RN,inbound)
    use mo_kind, only: dp
    implicit none
    !
    REAL(DP),    intent(IN)      :: old, dMax, oMin,oMax
    REAL(DP),    intent(IN)      :: RN
    LOGICAL,OPTIONAL,INTENT(OUT) :: inbound

    ! local
    LOGICAL   :: inbound2
    !
    inbound2 = .true.
    parGen_dp = old + (rn*2.0_dp*dMax - dMax)

    if (parGen_dp < oMin) then
       parGen_dp = oMin
       inbound2 = .false.
    elseif(parGen_dp > oMax)then
       parGen_dp = oMax
       inbound2 = .false.
    end if

    if (present(inbound)) inbound = inbound2

  end function parGen_dp

  ! ************************************************************************
  ! normalized Parameter Change
  ! ************************************************************************

  real(DP) function parGenNorm_dp(old, dMax, oMin, oMax,RN,inbound)
    use mo_kind, only: dp
    implicit none
    !
    REAL(DP),    intent(IN)      :: old, dMax, oMin,oMax
    REAL(DP),    intent(IN)      :: RN
    LOGICAL,OPTIONAL,INTENT(OUT) :: inbound

    ! LOCAL VARIABLES
    REAL(DP)                     :: RN_scale   ! scaled RN
    REAL(DP)                     :: old_scaled ! for normalization
    LOGICAL                      :: inbound2

    inbound2 =.true.

    !(1) normalize
    old_scaled  = ( old - oMin) / ( oMax - oMin)

    !(2) scale RN = [-dMax, dMax]
    !RN_scale = (RN * 2.0_dp * dMax) - dMax
    !(2) scale RN distrib. N(0, dMax)
    RN_scale = (RN * dMax)

    !(3) generate new parameter value + check inbound
    parGenNorm_dp = old_scaled + RN_scale
    !    Parameter is bounded between Max and Min.
    ! if (parGenNorm_dp < 0.0_dp) then
    !    parGenNorm_dp = 0.0_dp
    !    inbound2 = .false.
    ! elseif(parGenNorm_dp > 1.0_dp)then
    !    parGenNorm_dp = 1.0_dp
    !    inbound2 = .false.
    ! end if

    if (present(inbound)) inbound = inbound2

    !(4) convert back to original range
    parGenNorm_dp = parGenNorm_dp * ( oMax - oMin) + oMin

  end function parGenNorm_dp

  recursive subroutine GenerateNewParameterset_dp(ParaSelectMode, paraold, truepara, rangePar, stepsize, &
       save_state_2,save_state_3,paranew,ChangePara)

    integer(i4),                           intent(in)    :: ParaSelectMode
    real(dp),    dimension(:),             intent(in)    :: paraold
    integer(i4), dimension(:),             intent(in)    :: truepara
    real(dp),    dimension(:,:),           intent(in)    :: rangePar
    real(dp),    dimension(:),             intent(in)    :: stepsize
    integer(I8), dimension(n_save_state),  intent(inout) :: save_state_2  ! optional argument for restarting RN stream 2
    integer(I8), dimension(n_save_state),  intent(inout) :: save_state_3  ! optional argument for restarting RN stream 3
    real(dp),    dimension(size(paraold)), intent(out)   :: paranew
    logical,     dimension(size(paraold)), intent(out)   :: ChangePara

    ! local variables
    REAL(DP)       :: RN2, RN3
    logical        :: inbound
    integer(i4)    :: i, iPar

    paranew = paraold
    ChangePara = .False.

    select case(ParaSelectMode)
    case(1_i4) ! change half of the parameters
       do i=1,size(truepara)
          call xor4096(0_i8,RN2, save_state=save_state_2)
          if ( RN2 > 0.5_dp ) then
             iPar = truepara(i)
             inbound = .false.
             call xor4096(0_i8,RN3, save_state=save_state_3)
             paranew(iPar) = parGenNorm_dp( paraold(iPar), stepsize(iPar), rangePar(iPar,1), rangePar(iPar,2),RN3,inbound)
             ChangePara(iPar) = .True.
          end if
       end do
       ! if no parameter was selected, recall and change one
       if (count(ChangePara) .eq. 0_i4) then
          call xor4096(0_i8,RN2, save_state=save_state_2)
          iPar = truepara(int(( RN2 * real(size(truepara),dp))  + 1.0_dp, i4 ))
          inbound = .false.
          call xor4096(0_i8,RN3, save_state=save_state_3)
          paranew(iPar) = parGenNorm_dp( paraold(iPar), stepsize(iPar), rangePar(iPar,1), rangePar(iPar,2),RN3,inbound)
          ChangePara(iPar) = .True.
       end if

    case(2_i4)     ! change only one
       call xor4096(0_i8,RN2, save_state=save_state_2)
       iPar = truepara(int(( RN2 * real(size(truepara),dp))  + 1.0_dp, i4 ))
       inbound = .false.
       call xor4096g(0_i8,RN3, save_state=save_state_3)
       paranew(iPar) = parGenNorm_dp( paraold(iPar), stepsize(iPar), rangePar(iPar,1), rangePar(iPar,2),RN3,inbound)
 ! print*, 'p_old(',iPar,') = ',paraold(iPar)
 ! print*, 'p_new(',iPar,') = ',paranew(iPar)
       ChangePara(iPar) = .True.

    case(3_i4)    ! change all
       do i=1,size(truepara)
          iPar = truepara(i)
          inbound = .false.
          do while( .not. inbound)
             call xor4096g(0_i8,RN2, save_state=save_state_3)
             paranew(iPar) = parGenNorm_dp( paraold(iPar), stepsize(iPar), rangePar(iPar,1), rangePar(iPar,2),RN3,inbound)
          end do
          ChangePara(iPar) = .True.
       end do

    end select

  end subroutine GenerateNewParameterset_dp

END MODULE mo_mcmc
