!> \file mo_sce.f90

!> \brief Shuffled Complex Evolution optimization algorithm.

!> \details Optimization algorithm using Shuffled Complex Evolution strategy. 
!>          Original version 2.1 of Qingyun Duan (1992) rewritten in Fortran 90.

!> \authors Juliane Mai
!> \date Feb 2013

MODULE mo_sce

  ! This module is the Shuffled Complex Evolution optimization algorithm.

  ! Written Juliane Mai, Feb 2013

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Juliane Mai or Matthias Cuntz.

  ! Copyright 2011-2014 Juliane Mai, Matthias Cuntz

  IMPLICIT NONE

  PUBLIC :: sce      ! sce optimization

  ! ------------------------------------------------------------------

  !     NAME
  !         sce

  !     PURPOSE
  !>        \brief Shuffled Complex Evolution (SCE) algorithm for global optimization.
  !
  !>        \details  Shuffled Complex Evolution method for global optimization\n
  !>                  -- version 2.1\n
  !>                  \n
  !>                  by Qingyun Duan\n
  !>                  Department of Hydrology & Water Resources\n
  !>                  University of Arizona, Tucson, AZ 85721\n
  !>                  (602) 621-9360, email: duan@hwr.arizona.edu\n
  !>                  \n
  !>                  Written by Qingyun Duan,   Oct 1990.\n
  !>                  Revised by Qingyun Duan,   Aug 1991.\n
  !>                  Revised by Qingyun Duan,   Apr 1992.\n
  !>                  \n
  !>                  Re-written by Juliane Mai, Feb 2013.\n
  !>                  \n
  !>                  Statement by Qingyun Duan:\n
  !>                  ----------------------------------\n
  !>                  \n
  !>                     This general purpose global optimization program is developed at
  !>                     the Department of Hydrology & Water Resources of the University
  !>                     of Arizona.  Further information regarding the SCE-UA method can
  !>                     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
  !>                     at the address and phone number listed above.  We request all
  !>                     users of this program make proper reference to the paper entitled
  !>                     'Effective and Efficient Global Optimization for Conceptual
  !>                     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
  !>                     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.\n
  !>                  \n
  !>                 The function to be minimized is the first argument of DDS and must be defined as \n
  !>                 \code
  !>                     function functn(p)
  !>                          use mo_kind, only: dp
  !>                          implicit none
  !>                          real(dp), dimension(:), intent(in) :: p
  !>                          real(dp) :: functn
  !>                     end function functn
  !>                 \endcode
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp) :: functn(p)"               Function on which to search the optimum
  !>        \param[in] "real(dp) :: pini(:)"                 inital value of decision variables
  !>        \param[in] "real(dp) :: prange(size(pini),2)"    Min/max range of decision variables
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i8), optional :: mymaxn"    max no. of trials allowed before optimization is terminated\n
  !>                                                             DEFAULT: 1000_i8
  !>        \param[in] "logical,     optional :: mymaxit"   maximization (.true.) or minimization (.false.) of function\n
  !>                                                             DEFAULT: false
  !>        \param[in] "integer(i4), optional :: mykstop"   number of shuffling loops in which the criterion value must
  !>                                                             change by given percentage before optimiz. is terminated\n
  !>                                                             DEFAULT: 10_i4
  !>        \param[in] "real(dp),    optional :: mypcento"  percentage by which the criterion value must change in
  !>                                                             given number of shuffling loops\n
  !>                                                             DEFAULT: 0.0001_dp
  !>        \param[in] "real(dp),    optional :: mypeps"    optimization is terminated if volume of complex has 
  !>                                                             converged to given percentage of feasible space\n
  !>                                                             DEFAULT: 0.001_dp
  !>        \param[in] "integer(i8), optional :: myseed"    initial random seed\n
  !>                                                             DEFAULT: get_timeseed
  !>        \param[in] "integer(i4), optional :: myngs"     number of complexes in the initial population\n
  !>                                                             DEFAULT: 2_i4
  !>        \param[in] "integer(i4), optional :: mynpg"     number of points in each complex\n
  !>                                                             DEFAULT: 2*n+1
  !>        \param[in] "integer(i4), optional :: mynps"     number of points in a sub-complex\n
  !>                                                             DEFAULT: n+1
  !>        \param[in] "integer(i4), optional :: mynspl"    number of evolution steps allowed for each complex before
  !>                                                             complex shuffling\n
  !>                                                             DEFAULT: 2*n+1
  !>        \param[in] "integer(i4), optional :: mymings"   minimum number of complexes required, if the number of 
  !>                                                             complexes is allowed to reduce as the 
  !>                                                             optimization proceeds\n
  !>                                                             DEFAULT: ngs = number of complexes in initial population
  !>        \param[in] "integer(i4), optional :: myiniflg"  flag on whether to include the initial point in population\n
  !>                                                             0, not included\n
  !>                                                             1, included (DEFAULT)
  !>        \param[in] "integer(i4), optional :: myprint"   flag for controlling print-out after each shuffling loop\n
  !>                                                             0, print information on the best point of the population\n
  !>                                                             1, print information on every point of the population\n
  !>                                                             2, no printing (DEFAULT)
  !>                                                             3, same as 0 but print progress '.' on every function call
  !>                                                             4, same as 1 but print progress '.' on every function call
  !>        \param[in] "real(dp),    optional :: myalpha"   parameter for reflection  of points in complex\n
  !>                                                             DEFAULT: 0.8_dp
  !>        \param[in] "real(dp),    optional :: mybeta"    parameter for contraction of points in complex\n
  !>                                                             DEFAULT: 0.45_dp
  !>        \param[in]  "character(len=*), optional  :: tmp_file"    file for temporal output: write results after evolution loop\n
  !>                                                                    # of headlines: 7\n
  !>                                                                    format: '# nloop   icall   ngs1   bestf   worstf ... \n
  !>                                                                             ... gnrng   (bestx(j),j=1,nn)'
  !>        \param[in]  "character(len=*), optional  :: popul_file"  file for temporal output: writes whole population \n
  !>                                                                    # of headlines: 1 \n
  !>                                                                    format: #_evolution_loop, xf(i), (x(i,j),j=1,nn)\n
  !>                                                                    total number of lines written <= neval <= mymaxn\n
  !>        \param[in]  "logical, optional  :: popul_file_append"    if true, append to existing population file (default: false)\n
  !>        \param[in]  "logical, optional  :: parallel"    sce runs in parallel (true) or not (false)
  !>                                                             parallel sce should only be used if model/ objective 
  !>                                                             is not parallel
  !>                                                             DEAFULT: .false.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp),    optional              :: bestf"       the best value of the function.
  !>        \param[out] "integer(i8), optional              :: neval"       number of function evaluations needed.
  !>        \param[out] "real(dp),    optional, allocatable :: history(:)"  the history of best function values,
  !>                                                                        history(neval)=bestf
  !
  !     RETURN
  !>        \return real(dp) :: bestx(size(pini))  &mdash;  The parameters of the point which is estimated 
  !>                                                        to minimize/maximize the function.
  !
  !     RESTRICTIONS
  !>       \note No mask of parameters implemented yet.
  !>             SCE is OpenMP enabled on the loop over the complexes.
  !>             OMP_NUM_THREADS > 1 does not give reproducible results even when seeded!
  !
  !     EXAMPLE
  !         use mo_opt_functions, only: griewank
  !
  !         prange(:,1) = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  !         prange(:,2) = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  !         pini        = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
  !                            -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  !
  !         opt = sce(griewank, pini, prange)
  !
  !         -> see also example in test directory

  !     LITERATURE  
  !         Duan, Q., S. Sorooshian, and V.K. Gupta - 
  !             Effective and Efficient Global Optimization for Conceptual Rainfall-runoff Models
  !             Water Resources Research, Vol 28(4), pp.1015-1031, 1992.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified Juliane Mai, Matthias Cuntz, Jul 2013 - OpenMP
  !                                                        - NaN and Inf in objective function
  !                  Juliane Mai,                 Oct 2013 - added peps as optional argument
  !                                                        - allow for masked parameters
  !                                                        - write population to file --> popul_file
  !                                                        - write intermediate results to file --> tmp_file
  !                                                        - flag parallel introduced
  !                  Matthias Cuntz,              Nov 2013 - progress dots
  !                                                        - use iso_fortran_env
  !                                                        - treat functn=NaN as worse function value in cce
  !                  Matthias Cuntz,              May 2014 - sort -> orderpack
  !                  Matthias Cuntz,              May 2014 - popul_file_append
  !                  Matthias Cuntz,              May 2014 - sort with NaNs

  ! ------------------------------------------------------------------

  PRIVATE

  ! chkcst   ! check if a trial point satisfies all constraints.
  ! comp     ! reduces the number of complexes and points in a population
  ! parstt   ! checks for parameter convergence
  ! getpnt   ! generated a new point within its range
  ! cce      ! generates new point(s) from sub-complex

  ! ------------------------------------------------------------------

CONTAINS

  function sce(functn,pini,prange,                        & ! IN
       mymaxn,mymaxit,mykstop,mypcento,mypeps,myseed,     & ! Optional IN
       myngs,mynpg,mynps,mynspl,mymings,myiniflg,myprint, & ! Optional IN
       mymask,myalpha, mybeta,                            & ! Optional IN
       tmp_file, popul_file,                              & ! Optional IN
       popul_file_append,                                 & ! Optional IN
       parallel,                                          & ! OPTIONAL IN
       bestf,neval,history                                & ! Optional OUT
       ) result(bestx)

    use mo_kind,         only: i4, i8, dp
    use mo_orderpack,    only: sort
    use mo_string_utils, only: num2str, compress
    use mo_xor4096,      only: get_timeseed, n_save_state, xor4096, xor4096g
    !$ use omp_lib,      only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    use iso_fortran_env, only: output_unit, error_unit
#ifndef GFORTRAN
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, IEEE_QUIET_NAN
#endif

    implicit none
    ! 
    interface
       function functn(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind, only: dp
         real(dp), dimension(:),    intent(in)       :: paraset
         real(dp)                                    :: functn
       end function functn
    end interface
    real(dp), dimension(:),   intent(in)             :: pini        ! initial parameter set
    real(dp), dimension(:,:), intent(in)             :: prange      ! lower and upper bound on parameters
    !
    ! algorithmic control & convergence parameter 
    integer(i8), optional,               intent(in)  :: mymaxn      ! max no. of trials allowed before optimization is terminated
    !                                                               !     DEFAULT: 1000_i8
    logical,     optional,               intent(in)  :: mymaxit     ! maximization (.true.) or minimization (.false.) of function
    !                                                               !     DEFAULT: false
    integer(i4), optional,               intent(in)  :: mykstop     ! number of shuffling loops in which the criterion value must
    !                                                               !     change by given percentage before optimiz. is terminated
    !                                                               !     DEFAULT: 10_i4
    real(dp),    optional,               intent(in)  :: mypcento    ! percentage by which the criterion value must change in
    !                                                               !     given number of shuffling loops
    !                                                               !     DEFAULT: 0.0001_dp
    real(dp),    optional,               intent(in)  :: mypeps      ! optimization is terminated if volume of complex has 
    !                                                               ! converged to given percentage of feasible space
    !                                                               !     DEFAULT: 0.001_dp
    integer(i8), optional,               intent(in)  :: myseed      ! initial random seed
    !                                                               !     DEFAULT: get_timeseed
    integer(i4), optional,               intent(in)  :: myngs       ! number of complexes in the initial population
    !                                                               !     DEFAULT: 2_i4
    integer(i4), optional,               intent(in)  :: mynpg       ! number of points in each complex
    !                                                               !     DEFAULT: 2*n+1
    integer(i4), optional,               intent(in)  :: mynps       ! number of points in a sub-complex
    !                                                               !     DEFAULT: n+1
    integer(i4), optional,               intent(in)  :: mynspl      ! number of evolution steps allowed for each complex before
    !                                                               !     complex shuffling
    !                                                               !     DEFAULT: 2*n+1
    integer(i4), optional,               intent(in)  :: mymings     ! minimum number of complexes required, if the number of 
    !                                                               !     complexes is allowed to reduce as the 
    !                                                               !     optimization proceeds
    !                                                               !     DEFAULT: ngs = number of complexes in initial population
    integer(i4), optional,               intent(in)  :: myiniflg    ! flag on whether to include the initial point in population
    !                                                               !     0, not included
    !                                                               !     1, included (DEFAULT)
    integer(i4), optional,               intent(in)  :: myprint     ! flag for controlling print-out after each shuffling loop
    !                                                               !     0, print information on the best point of the population
    !                                                               !     1, print information on every point of the population
    !                                                               !     2, no printing (DEFAULT)
    logical,     optional,               intent(in), &
         dimension(size(pini,1))                     :: mymask      ! parameter included in optimization (true) or discarded (false)
    !                                                               !     DEFAULT: .true.
    real(dp),    optional,               intent(in)  :: myalpha     ! parameter for reflection  of points in complex
    !                                                               !     DEFAULT: 0.8_dp
    real(dp),    optional,               intent(in)  :: mybeta      ! parameter for contraction of points in complex
    !                                                               !     DEFAULT: 0.45_dp
    character(len=*), optional,          intent(in)  :: tmp_file    ! file for temporal output: write results after evolution loop
    !                                                               !     # of headlines: 7
    !                                                               !     format: '# nloop   icall   ngs1   bestf   worstf ...
    !                                                               !              ... gnrng   (bestx(j),j=1,nn)'
    character(len=*), optional,          intent(in)  :: popul_file  ! file for temporal output: writes whole population
    !                                                               !     # of headlines: 1
    !                                                               !     format: #_evolution_loop, xf(i), (x(i,j),j=1,nn)
    !                                                               !     total number of lines written <= neval <= mymaxn
    logical,          optional,          intent(in)  :: popul_file_append ! if true, append to popul_file
    logical,     optional,               intent(in)  :: parallel    ! sce runs in parallel (true) or not (false)
    !                                                               !     parallel sce should only be used if model/ objective 
    !                                                               !     is not parallel
    !                                                               !     DEAFULT: .false.
    real(dp),    optional,               intent(out) :: bestf       ! function value of bestx(.)
    integer(i8), optional,               intent(out) :: neval       ! number of function evaluations
    real(dp),    optional, &
         dimension(:), allocatable,      intent(out) :: history     ! history of best function values after each iteration
    real(dp), dimension(size(pini,1))                :: bestx       ! best point at current shuffling loop (is RETURN)
    !
    ! optionals transfer
    real(dp), dimension(size(pini,1))                :: bl          ! lower bound on parameters
    real(dp), dimension(size(pini,1))                :: bu          ! upper bound on parameters
    integer(i8)                                      :: maxn        ! max no. of trials allowed before optimization is terminated
    logical                                          :: maxit       ! minimization (false) or maximization (true)
    integer(i4)                                      :: kstop       ! number of shuffling loops in which the criterion value
    !                                                               !     must change
    real(dp)                                         :: pcento      ! percentage by which the criterion value must change
    real(dp)                                         :: peps        ! optimization is terminated if volume of complex has 
    !                                                               ! converged to given percentage of feasible space
    integer(i4)                                      :: ngs         ! number of complexes in the initial population
    integer(i4)                                      :: npg         ! number of points in each complex
    integer(i4)                                      :: nps         ! number of points in a sub-complex
    integer(i4)                                      :: nspl        ! number of evolution steps allowed for each complex before 
    !                                                               !     complex shuffling
    integer(i4)                                      :: mings       ! minimum number of complexes required
    integer(i4)                                      :: iniflg      ! flag on whether to include the initial point in population
    integer(i4)                                      :: iprint      ! flag for controlling print-out after each shuffling loop
    logical                                          :: idot        ! controls progress report .
    logical, dimension(size(pini,1))                 :: maskpara    ! mask(i) = .true.  --> parameter i will be optimized 
    !                                                               ! mask(i) = .false. --> parameter i will not be optimized 
    real(dp)                                         :: alpha       ! parameter for reflection  of points in complex
    real(dp)                                         :: beta        ! parameter for contraction of points in complex
    real(dp), dimension(:), allocatable              :: history_tmp ! history of best function values after each iteration
    real(dp), dimension(:), allocatable              :: ihistory_tmp ! local history for OpenMP
    real(dp)                                         :: bestf_tmp   ! function value of bestx(.)
    logical                                          :: parall      ! if sce is used parallel or not
    !
    ! local variables
    integer(i4)                                      :: nopt        ! number of parameters to be optimized
    integer(i4)                                      :: nn          ! total number of parameters 
    integer(i4)                                      :: npt         ! total number of points in initial population (npt=ngs*npg)
    real(dp)                                         :: fpini       ! function value at initial point
    real(dp),    dimension(:,:), allocatable         :: x           ! coordinates of points in the population 
    real(dp),    dimension(:),   allocatable         :: xf          ! function values of x                    
    real(dp),    dimension(size(pini,1))             :: xx          ! coordinates of a single point in x
    real(dp),    dimension(:,:), allocatable         :: cx          ! coordinates of points in a complex      
    real(dp),    dimension(:),   allocatable         :: cf          ! function values of cx
    real(dp),    dimension(:,:), allocatable         :: s           ! coordinates of points in the current sub-complex !simplex
    real(dp),    dimension(:),   allocatable         :: sf          ! function values of s
    real(dp),    dimension(:),   allocatable         :: worstx      ! worst point at current shuffling loop
    real(dp)                                         :: worstf      ! function value of worstx(.)
    real(dp),    dimension(:),   allocatable         :: xnstd       ! standard deviation of parameters in the population
    real(dp)                                         :: gnrng       ! normalized geometric mean of parameter ranges
    integer(i4), dimension(:),   allocatable         :: lcs         ! indices locating position of s(.,.) in x(.,.)  
    real(dp),    dimension(:),   allocatable         :: bound       ! length of range of ith variable being optimized
    real(dp),    dimension(:),   allocatable         :: unit        ! ?????
    integer(i4)                                      :: ngs1        ! number of complexes in current population
    integer(i4)                                      :: ngs2        ! number of complexes in last population
    real(dp),    dimension(:),   allocatable         :: criter      ! vector containing the best criterion values of the last
    !                                                               !     (kstop+1) shuffling loops
    integer(i4)                                      :: ipcnvg      ! flag indicating whether parameter convergence is reached
    !                                                               !      (i.e., check if gnrng is less than 0.001)
    !                                                               !      0, parameter convergence not satisfied
    !                                                               !      1, parameter convergence satisfied
    integer(i4)                                      :: nloop
    integer(i4)                                      :: loop
    integer(i4)                                      :: ii, jj, kk, ll
    integer(i4)                                      :: lpos
    logical                                          :: lpos_ok     ! for selction of points based on triangular 
    !                                                               ! probability distribution
    integer(i4)                                      :: npt1
    integer(i8)                                      :: icall       ! counter for function evaluations
    integer(i8)                                      :: icall_merk, iicall ! local icall for OpenMP
    integer(i4)                                      :: igs
    integer(i4)                                      :: k1, k2
    real(dp)                                         :: denomi      ! for checking improvement of last steps
    real(dp)                                         :: timeou      ! for checking improvement of last steps
    real(dp)                                         :: rtiny       ! for checking improvement of last steps
    character(4), dimension(:),  allocatable         :: xname       ! parameter names: "p1", "p2", "p3", ...
    character(512)                                   :: format_str1 ! format string
    character(512)                                   :: format_str2 ! format string
    ! for random numbers
    real(dp)                                         :: rand        ! random number
    integer(i8), dimension(:),   allocatable         :: iseed       ! initial random seed  
    integer(i8), dimension(:,:), allocatable         :: save_state_unif    ! save state of uniform stream
    integer(i8), dimension(:,:), allocatable         :: save_state_gauss   ! save state of gaussian stream
    real(dp),    dimension(:),   allocatable         :: rand_tmp    ! random number
    integer(i4)                                      :: ithread, n_threads   ! OMP or not
    integer(i8), dimension(n_save_state)             :: itmp
    real(dp),    dimension(:,:), allocatable         :: xtmp          ! tmp array for complex reduction
    real(dp),    dimension(:),   allocatable         :: ftmp          !            %
    real(dp)                                         :: large         ! for treating NaNs
    logical                                          :: ipopul_file_append
    integer(i4)                                      :: nonan         ! # of non-NaN in history_tmp
    real(dp), dimension(:), allocatable              :: htmp          ! tmp storage for history_tmp
#ifdef GFORTRAN
    real(dp)                                         :: NaN           ! NaN value
#endif

    if (present(parallel)) then
       parall = parallel
    else
       parall = .false.
    end if

    if (parall) then
       ! OpenMP or not
       n_threads = 1
       !$  write(output_unit,*) '--------------------------------------------------'
       !$OMP parallel
       !$    n_threads = OMP_GET_NUM_THREADS()
       !$OMP end parallel
       !$  write(output_unit,*) '   SCE is parellel with ',n_threads,' threads'
       !$  write(output_unit,*) '--------------------------------------------------'
    else
       n_threads = 1
    end if

    ! One random number chain per OpenMP thread
    allocate(rand_tmp(n_threads))
    allocate(iseed(n_threads))
    allocate(save_state_unif(n_threads,n_save_state))
    allocate(save_state_gauss(n_threads,n_save_state))

    if (present(mymask)) then
       if (.not. any(mymask)) stop 'mo_sce: all parameters are masked --> none will be optimized' 
       maskpara = mymask
    else
       maskpara = .true.
    end if

    ! number of parameters to optimize
    nopt = count(maskpara,dim=1)
    ! total number of parameters
    nn   = size(pini,1)

    ! input checking
    if (size(prange,dim=1) .ne. size(pini,1)) then
       stop 'mo_sce: prange has not matching rows'
    end if
    if (size(prange,dim=2) .ne. 2) then
       stop 'mo_sce: two colums expected for prange'
    end if
    bl(:)   = prange(:,1)
    bu(:)   = prange(:,2)
    do ii=1, nn
       if( ((bu(ii)-bl(ii)) .lt. tiny(1.0_dp) ) .and. maskpara(ii) ) then
          write(error_unit,*) 'para #',ii,'  :: range = ( ',bl(ii),' , ',bu(ii),' )'
          stop 'mo_sce: inconsistent or too small parameter range'
       end if
    end do
    
    !
    ! optionals checking
    if (present(mymaxn)) then
       if (mymaxn .lt. 2_i4) stop 'mo_sce: maxn has to be at least 2'
       maxn = mymaxn
    else
       maxn = 1000_i8
    end if
    if (present(mymaxit)) then
       maxit = mymaxit
    else
       maxit = .false.
    end if
    if (present(mykstop)) then
       if (mykstop .lt. 1_i4) stop 'mo_sce: kstop has to be at least 1'
       kstop = mykstop
    else
       kstop = 10_i4
    end if
    if (present(mypcento)) then
       if (mypcento .lt. 0_dp) stop 'mo_sce: pcento should be positive'
       pcento = mypcento
    else
       pcento = 0.0001_dp
    end if
    if (present(mypeps)) then
       if (mypeps .lt. 0_dp) stop 'mo_sce: peps should be positive'
       peps = mypeps
    else
       peps = 0.001_dp
    end if
    if (present(myseed)) then
       if (myseed .lt. 1_i8) stop 'mo_sce: seed should be non-negative'
       forall(ii=1:n_threads) iseed(ii) = myseed + (ii-1)*1000_i8
    else
       call get_timeseed(iseed)
    end if
    if (present(myngs)) then
       if (myngs .lt. 1_i4) stop 'mo_sce: ngs has to be at least 1'
       ngs = myngs
    else
       ngs = 2_i4
    end if
    if (present(mynpg)) then
       if (mynpg .lt. 3_i4) stop 'mo_sce: npg has to be at least 3'
       npg = mynpg
    else
       npg = 2*nopt + 1
    end if
    if (present(mynps)) then
       if (mynps .lt. 2_i4) stop 'mo_sce: nps has to be at least 2'
       nps = mynps
    else
       nps = nopt + 1_i4
    end if
    if (present(mynspl)) then
       if (mynspl .lt. 3_i4) stop 'mo_sce: nspl has to be at least 3'
       nspl = mynspl
    else
       nspl = 2*nopt + 1
    end if
    if (present(mymings)) then
       if (mymings .lt. 1_i4) stop 'mo_sce: mings has to be at least 1'
       mings = mymings
    else
       mings = ngs  ! no reduction of complexes
    end if
    if (present(myiniflg)) then
       if ( (myiniflg .ne. 1_i4) .and. (myiniflg .ne. 0_i4)) stop 'mo_sce: iniflg has to be 0 or 1'
       iniflg = myiniflg
    else
       iniflg = 1_i4
    end if
    idot = .false.
    if (present(myprint)) then
       if ( (myprint .lt. 0_i4) .or. (myprint .gt. 4_i4)) stop 'mo_sce: iprint has to be between 0 and 4'
       if (myprint > 2_i4) then
          iprint = myprint-3
          idot = .true.
       else
          iprint = myprint
       endif
    else
       iprint = 2_i4  ! no printing
       iprint = 0_i4
    end if
    if (present(myalpha)) then
       alpha = myalpha
    else
       alpha = 0.8_dp
    end if
    if (present(mybeta)) then
       beta = mybeta
    else
       beta = 0.45_dp
    end if

    if (present(popul_file_append)) then
       ipopul_file_append = popul_file_append
    else
       ipopul_file_append = .false.
    endif

    if(present(tmp_file)) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', status = 'unknown')
       write(999,*) '# settings :: general'
       write(999,*) '# nIterations    seed'
       write(999,*) maxn, iseed
       write(999,*) '# settings :: sce specific'
       write(999,*) '# sce_ngs    sce_npg    sce_nps'
       write(999,*) ngs, npg, nps
       write(999,*) '# nloop   icall   ngs1   bestf   worstf   gnrng   (bestx(j),j=1,nn)'
       close(999)
    end if

    if (present(popul_file) .and. (.not. ipopul_file_append)) then
       open(unit=999,file=trim(adjustl(popul_file)), action='write', status = 'unknown')
       write(999,*) '#   xf(i)   (x(i,j),j=1,nn)'
       close(999)
    end if

    ! allocation of arrays
    allocate(x(ngs*npg,nn))
    allocate(xf(ngs*npg))
    allocate(worstx(nn))
    allocate(xnstd(nn))
    allocate(bound(nn))
    allocate(unit(nn))
    allocate(criter(kstop+1))  
    allocate(xname(nn))
    allocate(history_tmp(maxn+3*ngs*nspl))
    allocate(xtmp(npg,nn)) 
    allocate(ftmp(npg))
    if (maxit) then
       large = -huge(1.0_dp)
    else
       large = huge(1.0_dp)
    endif
    criter(:) = large
    !
    !  initialize variables
    do ii=1,nn
       xname(ii) = compress('p'//num2str(ii))
    end do

    nloop = 0_i4
    loop  = 0_i4
    igs   = 0_i4
    !
    !  compute the total number of points in initial population
    npt = ngs * npg
    ngs1 = ngs
    npt1 = npt
    !
    if (iprint .lt. 2) then
       write(output_unit,*) '=================================================='
       write(output_unit,*) 'ENTER THE SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH'
       write(output_unit,*) '=================================================='
    end if
    !
    !  Print seed
    if (iprint .lt. 2) then
       write (*,*) ' Seeds used : ',iseed
    end if
    !  initialize random number generator: Uniform
    call xor4096(iseed, rand_tmp, save_state=save_state_unif)
    !  initialize random number generator: Gaussian
    iseed = iseed + 1000_i8
    call xor4096g(iseed, rand_tmp, save_state=save_state_gauss)
    iseed = 0_i8
    !
    !  compute the bound for parameters being optimized
    do ii = 1, nn
       bound(ii) = bu(ii) - bl(ii)
       unit(ii) = 1.0_dp
    end do
    !--------------------------------------------------
    !  compute the function value of the initial point
    !--------------------------------------------------
    ! function evaluation will be counted later...
    if (idot) write(output_unit,'(A1)') '.'
    if (.not. maxit) then
       fpini = functn(pini)
       history_tmp(1) = fpini
    else
       fpini = -functn(pini)
       history_tmp(1) = -fpini
    end if

    !  print the initial point and its criterion value
    bestx = pini
    bestf_tmp = fpini
    if (iprint .lt. 2) then
       write(output_unit,*) ''
       write(output_unit,*) ''
       write(output_unit,*) '*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***'
       call write_best_final()
    end if
    file_write_initial: if (present(tmp_file)) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append',recl=(nn+6)*30)
       if (.not. maxit) then
          write(999,*) 0,1,ngs1,fpini,fpini,1.0_dp, pini
       else
          write(999,*) 0,1,ngs1,-fpini,-fpini,1.0_dp, pini
       end if
       close(999)
    end if file_write_initial
    !
    !  generate an initial set of npt1 points in the parameter space
    !  if iniflg is equal to 1, set x(1,.) to initial point pini(.)
    if (iniflg .eq. 1) then
       do ii = 1, nn
          x(1,ii) = pini(ii)
       end do
       xf(1) = fpini
       !
       !  else, generate a point randomly and set it equal to x(1,.)
    else
       itmp = save_state_unif(1,:)
       call getpnt(1,bl(1:nn),bu(1:nn),unit(1:nn),pini(1:nn),maskpara,itmp,xx)
       save_state_unif(1,:) = itmp
       do ii=1, nn
          x(1,ii) = xx(ii)
       end do
       if (idot) write(output_unit,'(A1)') '.'
       if (.not. maxit) then
          xf(1) = functn(xx)
       else
          xf(1) = -functn(xx)
       end if
    end if
    !
    ! count function evaluation of the first point
    icall = 1_i8
    ! if (icall .ge. maxn) return 
    !
    if (parall) then
       
       ! -----------------------------------------------------------------------
       ! Parallel version of complex-loop
       ! -----------------------------------------------------------------------
       !  generate npt1-1 random points distributed uniformly in the parameter
       !  space, and compute the corresponding function values
       ithread = 1
       !$OMP parallel default(shared) private(ii,jj,ithread,xx)
       !$OMP do
       do ii = 2, npt1
          !$ ithread = OMP_GET_THREAD_NUM() + 1
          itmp = save_state_unif(ithread,:)
          call getpnt(1,bl(1:nn),bu(1:nn),unit(1:nn),pini(1:nn),maskpara,itmp,xx)
          save_state_unif(ithread,:) = itmp
          do jj = 1, nn
             x(ii,jj) = xx(jj)
          end do
          if (idot) write(output_unit,'(A1)') '.'
          if (.not. maxit) then
             xf(ii) = functn(xx)
             history_tmp(ii) = xf(ii) ! min(history_tmp(ii-1),xf(ii)) --> will be sorted later
          else
             xf(ii) = -functn(xx)
             history_tmp(ii) = -xf(ii) ! max(history_tmp(ii-1),-xf(ii)) --> will be sorted later
          end if
       end do
       !$OMP end do
       !$OMP end parallel

    else

       ! -----------------------------------------------------------------------
       ! Non-Parallel version of complex-loop
       ! -----------------------------------------------------------------------
       !  generate npt1-1 random points distributed uniformly in the parameter
       !  space, and compute the corresponding function values
       ithread = 1

       do ii = 2, npt1
          call getpnt(1,bl(1:nn),bu(1:nn),unit(1:nn),pini(1:nn),maskpara,save_state_unif(ithread,:),xx)
          do jj = 1, nn
             x(ii,jj) = xx(jj)
          end do
          if (idot) write(output_unit,'(A1)') '.'
          if (.not. maxit) then
             xf(ii) = functn(xx)
             icall = icall + 1_i8
             history_tmp(icall) = xf(ii) ! min(history_tmp(icall-1),xf(ii)) --> will be sorted later
          else
             xf(ii) = -functn(xx)
             icall = icall + 1_i8
             history_tmp(icall) = -xf(ii) ! max(history_tmp(icall-1),-xf(ii)) --> will be sorted later
          end if
          
          if (icall .ge. maxn) then
             npt1 = ii
             exit
          end if
       end do

    end if

    icall = int(npt1,i8)
    !
    !  arrange the points in order of increasing function value
    if (maxit) then
       large = minval(xf(1:npt1))
       large = merge(0.9_dp*large, 1.1_dp*large, large>0._dp)
    else
       large = maxval(xf(1:npt1))
       large = merge(1.1_dp*large, 0.9_dp*large, large>0._dp)
    endif
#ifndef GFORTRAN
    xf(1:npt1) = merge(xf(1:npt1), large, ieee_is_finite(xf(1:npt1))) ! NaN and Infinite
    ! sort does not work with NaNs
    ! -> get history_tmp w/o NaN, sort it, and set the rest to NaN
    nonan = size(pack(history_tmp(1:npt1), mask=ieee_is_finite(history_tmp(1:npt1))))
    if (nonan /= npt1) then
       allocate(htmp(nonan))
       htmp(1:nonan) = pack(history_tmp(1:npt1), mask=ieee_is_finite(history_tmp(1:npt1)))
       call sort(htmp(1:nonan))
       history_tmp(1:nonan) = htmp(1:nonan)
       history_tmp(nonan+1:npt1) = ieee_value(1.0_dp, IEEE_QUIET_NAN)
       deallocate(htmp)
    else
       call sort(history_tmp(1:npt1))
    endif
#else
#ifdef GFORTRAN41
    xf(1:npt1) = merge(xf(1:npt1), large, xf(1:npt1) == xf(1:npt1))   ! only NaN
    ! sort does not work with NaNs
    ! -> get the NaN value
    ! -> get history_tmp w/o NaN, sort it, and set the rest to NaN
    nonan = size(pack(history_tmp(1:npt1), mask=(xf(1:npt1)==xf(1:npt1))))
    if (nonan /= npt1) then
       NaN = large
       do ii=1, npt1
          if (xf(ii) /= xf(ii)) then
             NaN = xf(ii)
             exit
          endif
       end do
       allocate(htmp(nonan))
       htmp(1:nonan) = pack(history_tmp(1:npt1), mask=(xf(1:npt1)==xf(1:npt1)))
       call sort(htmp(1:nonan))
       history_tmp(1:nonan) = htmp(1:nonan)
       history_tmp(nonan+1:npt1) = NaN
       deallocate(htmp)
    else
       call sort(history_tmp(1:npt1))
    endif
#else
    xf(1:npt1) = merge(xf(1:npt1), large, .not. isnan(xf(1:npt1)))   ! only NaN
    ! sort does not work with NaNs
    ! -> get the NaN value
    ! -> get history_tmp w/o NaN, sort it, and set the rest to NaN
    nonan = size(pack(history_tmp(1:npt1), mask=(.not. isnan(xf(1:npt1)))))
    if (nonan /= npt1) then
       NaN = large
       do ii=1, npt1
          if (isnan(xf(ii))) then
             NaN = xf(ii)
             exit
          endif
       end do
       allocate(htmp(nonan))
       htmp(1:nonan) = pack(history_tmp(1:npt1), mask=(.not. isnan(xf(1:npt1))))
       call sort(htmp(1:nonan))
       history_tmp(1:nonan) = htmp(1:nonan)
       history_tmp(nonan+1:npt1) = NaN
       deallocate(htmp)
    else
       call sort(history_tmp(1:npt1))
    endif
#endif
#endif
    call sort_matrix(x(1:npt1,1:nn),xf(1:npt1))
    !
    !  record the best and worst points
    do ii = 1, nn
       bestx(ii) = x(1,ii)
       worstx(ii) = x(npt1,ii)
    end do
    bestf_tmp = xf(1)
    worstf = xf(npt1)
    !
    !  compute the parameter range for the initial population
    call parstt(x(1:npt1,1:nn),bound,peps,maskpara,xnstd,gnrng,ipcnvg)
    !
    ! write currently best x and xf to temporal file
    if (present(tmp_file)) then
       call write_best_intermediate(.true.)
    end if

    ! write population to file
    if (present(popul_file)) then
       call write_population(.true.)
    end if
    !
    !  print the results for the initial population
    print: if (iprint .lt. 2) then
       ! write currently best x and xf to screen
       call write_best_intermediate(.false.)

       ! write population on screen
       if (iprint .eq. 1) then
          call write_population(.false.)
       end if
    end if print
    !
    ! Maximum number of function evaluations reached?
    if (icall .ge. maxn) then 
       if (iprint .lt. 2) then
          ! maximum trials reached
          call write_termination_case(1)
          call write_best_final()
       end if
       call set_optional()
       ! -----------------------
       ! Abort subroutine
       ! -----------------------
       return
    end if
    !
    ! Feasible parameter space converged?
    if (ipcnvg .eq. 1) then 
       if (iprint .lt. 2) then
          ! converged because feasible parameter space small
          call write_termination_case(3)
          call write_best_final()
       end if
       call set_optional()
       ! -----------------------
       ! Abort subroutine
       ! -----------------------
       return
    end if
    !
    !  begin the main loop 
    mainloop:do while (icall .lt. maxn)
       nloop = nloop + 1
       !
       if (iprint .lt. 2) then
          write(output_unit, *) ''
          write(output_unit,'(A28,I4)') ' ***  Evolution Loop Number ',nloop
       end if
       !
       !  begin loop on complexes
       !     <beta> loop from duan(1993)

       if (parall) then

          ! -----------------------------------------------------------------------
          ! Parallel version of complex-loop
          ! -----------------------------------------------------------------------

          ithread = 1
          !$OMP parallel default(shared) &
          !$OMP private(igs, loop, ithread, kk, k1, k2, jj, lpos_ok, lpos, rand, cx, cf, lcs, s, sf) &
          !$OMP private(icall_merk, iicall, ihistory_tmp, large)
          allocate(cx(npg,nn)) 
          allocate(cf(npg))
          allocate(lcs(nps))
          allocate(s(nps,nn))
          allocate(sf(nps))
          allocate(ihistory_tmp(maxn+3*ngs*nspl))
          !$OMP do
          comploop_omp: do igs = 1, ngs1
             !$ ithread = OMP_GET_THREAD_NUM() + 1
             !
             !  assign points into complexes
             do k1 = 1, npg
                k2 = (k1-1) * ngs1 + igs
                do jj = 1, nn
                   cx(k1,jj) = x(k2,jj)
                end do
                cf(k1) = xf(k2)
             end do
             !
             !  begin inner loop - random selection of sub-complexes 
             !     <alpha> loop from duan(1993)
             subcomploop_omp: do loop = 1, nspl
                !
                !  choose a sub-complex (nps points) according to a linear
                !  probability distribution
                !
                ! if number of points in subcomplex (nps) = number of points in complex (npg)
                ! --> select all
                if (nps .eq. npg) then
                   do kk = 1, nps
                      lcs(kk) = kk
                   end do
                else
                   call xor4096(iseed(ithread), rand, save_state=save_state_unif(ithread,:))
                   lcs(1) = 1 + int(real(npg,dp) + 0.5_dp &
                        - sqrt( (real(npg,dp)+.5_dp)**2 - real(npg,dp) * real(npg+1,dp) * rand ))
                   do kk = 2, nps
                      lpos_ok = .false.
                      do while (.not. lpos_ok)
                         lpos_ok = .true.
                         call xor4096(iseed(ithread), rand, save_state=save_state_unif(ithread,:))
                         lpos = 1 + int(real(npg,dp) + 0.5_dp &
                              - sqrt( (real(npg,dp)+.5_dp)**2 - real(npg,dp) * real(npg+1,dp) * rand ))
                         ! test if point already chosen: returns lpos_ok=false if any point already exists
                         do k1 = 1, kk-1
                            if (lpos .eq. lcs(k1)) lpos_ok = (lpos_ok .and. .false.)
                         end do
                      end do
                      lcs(kk) = lpos
                   end do
                   !
                   !  arrange the sub-complex in order of increasing function value
                   call sort(lcs(1:nps))
                end if
                !
                !  create the sub-complex arrays
                do kk = 1, nps
                   do jj = 1, nn
                      s(kk,jj) = cx(lcs(kk),jj)
                   end do
                   sf(kk) = cf(lcs(kk))
                end do
                !
                ! remember largest for treating of NaNs
                if (maxit) then
                   large = minval(cf(1:npg))
                   large = merge(0.9_dp*large, 1.1_dp*large, large>0._dp)
                else
                   large = maxval(cf(1:npg))
                   large = merge(1.1_dp*large, 0.9_dp*large, large>0._dp)
                endif
                !
                !  use the sub-complex to generate new point(s)
                icall_merk = icall
                iicall     = icall
                ihistory_tmp = history_tmp
                call cce(s(1:nps,1:nn),sf(1:nps),bl(1:nn),bu(1:nn),maskpara,xnstd(1:nn),  &
                     iicall,maxn,maxit,save_state_gauss(ithread,:),functn,alpha,beta,ihistory_tmp,idot)
                history_tmp(icall+1:icall+(iicall-icall_merk)) = ihistory_tmp(icall_merk+1:iicall)
                icall = icall + (iicall-icall_merk)
                !
                !  if the sub-complex is accepted, replace the new sub-complex
                !  into the complex
                do kk = 1, nps
                   do jj = 1, nn
                      cx(lcs(kk),jj) = s(kk,jj)
                   end do
                   cf(lcs(kk)) = sf(kk)
                end do
                !
                !  sort the points
#ifndef GFORTRAN
                cf(1:npg) = merge(cf(1:npg), large, ieee_is_finite(cf(1:npg))) ! NaN and Infinite
#else
#ifdef GFORTRAN41
                cf(1:npg) = merge(cf(1:npg), large, cf(1:npg) == cf(1:npg))   ! only NaN
#else
                cf(1:npg) = merge(cf(1:npg), large, .not. isnan(cf(1:npg)))   ! only NaN
#endif
#endif
                call sort_matrix(cx(1:npg,1:nn),cf(1:npg))
                !
                ! !  if maximum number of runs exceeded, break out of the loop
                ! if (icall .ge. maxn) exit
                !
             end do subcomploop_omp ! <alpha loop>
             !
             !  replace the new complex into original array x(.,.)
             do k1 = 1, npg
                k2 = (k1-1) * ngs1 + igs
                do jj = 1, nn
                   x(k2,jj) = cx(k1,jj)
                end do
                xf(k2) = cf(k1)
             end do
             ! if (icall .ge. maxn) exit
             !
             !  end loop on complexes
          end do comploop_omp  ! <beta loop> 
          !$OMP end do
          deallocate(cx)
          deallocate(cf)
          deallocate(lcs)
          deallocate(s)
          deallocate(sf)
          deallocate(ihistory_tmp)
          !$OMP end parallel

       else

          ! -----------------------------------------------------------------------
          ! Non-parallel version of complex-loop
          ! -----------------------------------------------------------------------

          allocate(cx(npg,nn)) 
          allocate(cf(npg))
          allocate(lcs(nps))
          allocate(s(nps,nn))
          allocate(sf(nps))

          ithread = 1

          comploop: do igs = 1, ngs1
             !
             !  assign points into complexes
             do k1 = 1, npg
                k2 = (k1-1) * ngs1 + igs
                do jj = 1, nn
                   cx(k1,jj) = x(k2,jj)
                end do
                cf(k1) = xf(k2)
             end do
             !
             !  begin inner loop - random selection of sub-complexes 
             !     <alpha> loop from duan(1993)
             subcomploop: do loop = 1, nspl
                !
                !  choose a sub-complex (nps points) according to a linear
                !  probability distribution
                !
                ! if number of points in subcomplex (nps) = number of points in complex (npg)
                ! --> select all
                if (nps .eq. npg) then
                   do kk = 1, nps
                      lcs(kk) = kk
                   end do
                else
                   call xor4096(iseed(ithread), rand, save_state=save_state_unif(ithread,:))
                   lcs(1) = 1 + int(real(npg,dp) + 0.5_dp &
                        - sqrt( (real(npg,dp)+.5_dp)**2 - real(npg,dp) * real(npg+1,dp) * rand ))
                   do kk = 2, nps
                      lpos_ok = .false.
                      do while (.not. lpos_ok)
                         lpos_ok = .true.
                         call xor4096(iseed(ithread), rand, save_state=save_state_unif(ithread,:))
                         lpos = 1 + int(real(npg,dp) + 0.5_dp &
                              - sqrt( (real(npg,dp)+.5_dp)**2 - real(npg,dp) * real(npg+1,dp) * rand ))
                         ! test if point already chosen: returns lpos_ok=false if any point already exists
                         do k1 = 1, kk-1
                            if (lpos .eq. lcs(k1)) lpos_ok = (lpos_ok .and. .false.)
                         end do
                      end do
                      lcs(kk) = lpos
                   end do
                   !
                   !  arrange the sub-complex in order of increasing function value
                   call sort(lcs(1:nps))
                end if
                !
                !  create the sub-complex arrays
                do kk = 1, nps
                   do jj = 1, nn
                      s(kk,jj) = cx(lcs(kk),jj)
                   end do
                   sf(kk) = cf(lcs(kk))
                end do
                !
                ! remember largest for treating of NaNs
                if (maxit) then
                   large = minval(cf(1:npg))
                   large = merge(0.9_dp*large, 1.1_dp*large, large>0._dp)
                else
                   large = maxval(cf(1:npg))
                   large = merge(1.1_dp*large, 0.9_dp*large, large>0._dp)
                endif
                !
                !  use the sub-complex to generate new point(s)
                call cce(s(1:nps,1:nn),sf(1:nps),bl(1:nn),bu(1:nn),maskpara,xnstd(1:nn),  &
                     icall,maxn,maxit,save_state_gauss(ithread,:),functn, alpha,beta,history_tmp,idot)
                !
                !  if the sub-complex is accepted, replace the new sub-complex
                !  into the complex
                do kk = 1, nps
                   do jj = 1, nn
                      cx(lcs(kk),jj) = s(kk,jj)
                   end do
                   cf(lcs(kk)) = sf(kk)
                end do
                !
                !  sort the points
#ifndef GFORTRAN
                cf(1:npg) = merge(cf(1:npg), large, ieee_is_finite(cf(1:npg))) ! NaN and Infinite
#else
#ifdef GFORTRAN41
                cf(1:npg) = merge(cf(1:npg), large, cf(1:npg) == cf(1:npg))   ! only NaN
#else
                cf(1:npg) = merge(cf(1:npg), large, .not. isnan(cf(1:npg)))   ! only NaN
#endif
#endif
                call sort_matrix(cx(1:npg,1:nn),cf(1:npg))
                !
                !  if maximum number of runs exceeded, break out of the loop
                if (icall .ge. maxn) exit
                !
             end do subcomploop ! <alpha loop>
             !
             !  replace the new complex into original array x(.,.)
             do k1 = 1, npg
                k2 = (k1-1) * ngs1 + igs
                do jj = 1, nn
                   x(k2,jj) = cx(k1,jj)
                end do
                xf(k2) = cf(k1)
             end do
             if (icall .ge. maxn) exit
             !
             !  end loop on complexes
          end do comploop  ! <beta loop> 

          deallocate(cx)
          deallocate(cf)
          deallocate(lcs)
          deallocate(s)
          deallocate(sf)

       end if ! end parallel/non-parallel region
       !
       !  re-sort the points
       call sort_matrix(x(1:npt1,1:nn),xf(1:npt1))
       !
       !  record the best and worst points
       do jj = 1, nn
          bestx(jj) = x(1,jj)
          worstx(jj) = x(npt1,jj)
       end do
       bestf_tmp = xf(1)
       worstf    = xf(npt1)
       !
       !  test the population for parameter convergence
       call parstt(x(1:npt1,1:nn),bound,peps,maskpara,xnstd,gnrng,ipcnvg)
       !
       ! write currently best x and xf to temporal file
       if (present(tmp_file)) then
          call write_best_intermediate(.true.)
       end if
       !
       ! write population to file
       if (present(popul_file)) then
          call write_population(.true.)
       end if
       !
       !  print the results for current population
       if (iprint .lt. 2) then
          call write_best_intermediate(.false.)

          if (iprint .eq. 1) then
             call write_population(.false.)
          end if

       end if
       !
       !  test if maximum number of function evaluations exceeded
       ! Maximum number of function evaluations reached?
       if (icall .ge. maxn) then
          if (iprint .lt. 2) then
             ! maximum trials reached
             call write_termination_case(1)
             call write_best_final()
          end if
          call set_optional()
          ! -----------------------
          ! Abort subroutine
          ! -----------------------
          return
       end if
       !
       !  compute the count on successive loops w/o function improvement
       criter(kstop+1) = bestf_tmp
       if (nloop .gt. kstop) then 
          denomi = dabs(criter((kstop+1)-kstop) + criter(kstop+1)) / 2.
          rtiny = 1.0/(1.0**15)
          denomi = max(denomi, rtiny)
          timeou = dabs(criter((kstop+1)-kstop) - criter(kstop+1)) / denomi
          if (timeou .lt. pcento) then 
             if (iprint .lt. 2) then
                ! criterion value has not changed during last loops
                call write_termination_case(2)
                call write_best_final()
             end if
             call set_optional()
             ! -----------------------
             ! Abort subroutine
             ! -----------------------
             return
          end if
       end if
       do ll = 1, kstop
          criter(ll) = criter(ll+1)
       end do
       !
       !  if population is converged into a sufficiently small space
       if (ipcnvg .eq. 1) then 
          if (iprint .lt. 2) then
             ! converged because feasible parameter space small
             call write_termination_case(3)
             call write_best_final()
          end if
          call set_optional()
          ! -----------------------
          ! Abort subroutine
          ! -----------------------
          return
       end if
       !
       !  none of the stopping criteria is satisfied, continue search
       !
       !  check for complex number reduction
       if (ngs1 .gt. mings) then
          ngs2 = ngs1
          ngs1 = ngs1 - 1
          npt1 = ngs1 * npg
          call comp(ngs2,npg,x(1:ngs2*npg,1:nn),xf(1:ngs2*npg),xtmp(1:ngs2*npg,1:nn),ftmp(1:ngs2*npg))
       end if
    end do mainloop

    deallocate(xtmp)
    deallocate(ftmp)

    call set_optional()

  contains

    subroutine write_best_intermediate(to_file)

      implicit none

      logical, intent(in) :: to_file

      if (to_file) then
         open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append',recl=(nn+6)*30)
         if (.not. maxit) then
            write(999,*) nloop,icall,ngs1,bestf_tmp,worstf,gnrng, (bestx(jj),jj=1,nn)
         else
            write(999,*) nloop,icall,ngs1,-bestf_tmp,-worstf,gnrng, (bestx(jj),jj=1,nn)
         end if
         close(999)
      else
         write(format_str1,'(A13,I3,A8)') '( A49, ',nn,'(6X,A4))'
         write(format_str2,'(A26,I3,A8)') '(I5,1X,I5,3X,I5,3(E22.14), ',nn,'(G10.3))'
         if (nloop == 0) then
            write(output_unit,*) ''
            write(output_unit,'(A44)') ' *** PRINT THE RESULTS OF THE SCE SEARCH ***'
            write(output_unit,*) ''
            write(output_unit,format_str1) ' LOOP TRIALS COMPLXS  BEST F   WORST F   PAR RNG ',(xname(jj),jj=1,nn)
         end if
         if (.not. maxit) then
            write(output_unit,format_str2) nloop,icall,ngs1,bestf_tmp,worstf,gnrng, (bestx(jj),jj=1,nn)
         else
            write(output_unit,format_str2) nloop,icall,ngs1,-bestf_tmp,-worstf,gnrng, (bestx(jj),jj=1,nn)
         end if
      end if

    end subroutine write_best_intermediate

    subroutine write_best_final()

      implicit none

      write(format_str1,'(A13,I3,A8)') '( A10, ',nn,'(6X,A4))'
      write(format_str2,'(A14,I3,A8)') '(E22.14, ',nn,'(G10.3))'

      write(output_unit,format_str1) 'CRITERION ',(xname(jj),jj=1,nn)
      if (.not. maxit) then
         write(output_unit,format_str2) bestf_tmp,(bestx(jj),jj=1,nn)
      else
         write(output_unit,format_str2) -bestf_tmp,(bestx(jj),jj=1,nn)
      end if

    end subroutine write_best_final

    subroutine write_population(to_file)

      logical, intent(in) :: to_file

      if (to_file) then
         write(format_str2,'(A13,I3,A9)') '(I4, E22.14, ',nn,'(E22.14))'
         open(unit=999,file=trim(adjustl(popul_file)), action='write', position='append',recl=(nn+2)*30)
         if (.not. maxit) then
            do ii = 1, npt1
               write(999,*) nloop, xf(ii), (x(ii,jj),jj=1,nn)
            end do
         else
            do ii = 1, npt1
               write(999,*) nloop, -xf(ii), (x(ii,jj),jj=1,nn)
            end do
         end if
         close(999)
      else
         write(format_str2,'(A12,I3,A8)') '(I4, E22.14, ',nn,'(G10.3))'
         write(output_unit,*) ''
         write(output_unit,'(A22,I3)') '   POPULATION AT LOOP ',nloop
         write(output_unit,'(A27)')    '---------------------------'
         if (.not. maxit) then
            do ii = 1, npt1
               write(output_unit,format_str2) nloop, xf(ii), (x(ii,jj),jj=1,nn)
            end do
         else
            do ii = 1, npt1
               write(output_unit,format_str2) nloop, -xf(ii), (x(ii,jj),jj=1,nn)
            end do
         end if
      end if

    end subroutine write_population

    subroutine write_termination_case(case)

      integer(i4), intent(in) :: case

      select case (case)
      case (1) ! maximal number of iterations reached
         write(output_unit,*) ''
         write(output_unit,'(A46,A39,I7,A46,I4,A12,I4,A19,I4,A4)') &
              '*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE', &
              ' LIMIT ON THE MAXIMUM NUMBER OF TRIALS ', maxn,   &
              ' EXCEEDED.  SEARCH WAS STOPPED AT SUB-COMPLEX ', &
              loop, ' OF COMPLEX ', igs,' IN SHUFFLING LOOP ', nloop, ' ***'
         !
      case(2) ! objective not changed during last evolution loops
         write(output_unit,*) ''
         write(output_unit,'(A72,F8.4,A12,I3,A20)') '*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION VALUE HAS NOT CHANGED ', &
              pcento*100._dp,' PERCENT IN ', kstop, ' SHUFFLING LOOPS ***'
         !
      case(3) ! complexes converged
         write(output_unit,*) ''
         write(output_unit,'(A50,A20,F8.4,A34)') '*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION', &
              ' HAS CONVERGED INTO ', gnrng*100._dp, ' PERCENT OF THE FEASIBLE SPACE ***'
         !
      case default
         write(error_unit,*) 'This termination case is not implemented!'
         stop
      end select

      write(output_unit,*) ''
      write(output_unit,'(A66)') '*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS CRITERION VALUE ***'

    end subroutine write_termination_case

    subroutine set_optional()

      implicit none

      if (present(neval))                   neval = icall
      if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
      if (present(bestf) .and. maxit)       bestf = -bestf_tmp
      if (present(history)) then
         allocate(history(icall))
         history(:) = history_tmp(1:icall)
      end if

    end subroutine set_optional

  end function sce


  subroutine parstt(x,bound,peps,mask,xnstd,gnrng,ipcnvg)   
    !
    !  subroutine checking for parameter convergence
    !
    use mo_kind, only: i4, dp

    implicit none  

    real(dp),    dimension(:,:),           intent(in)  :: x      ! points in population, cols=nn, rows=npt1 (!)
    real(dp),    dimension(:),             intent(in)  :: bound  ! difference of upper and lower limit per parameter
    real(dp),                              intent(in)  :: peps   ! optimization is terminated if volume of complex has 
    !                                                            ! converged to given percentage of feasible space
    logical,     dimension(:),            intent(in)   :: mask   ! mask of parameters
    real(dp),    dimension(size(bound,1)), intent(out) :: xnstd  ! std. deviation of points in population per parameter
    real(dp),                              intent(out) :: gnrng  ! fraction of feasible space covered by complexes
    integer(i4),                           intent(out) :: ipcnvg ! 1 : population converged into sufficiently small space
    !                                                            ! 0 : not converged
    !
    ! local variables
    integer(i4)                         :: nn      ! number of parameters
    integer(i4)                         :: npt1    ! number of points in current population
    integer(i4)                         :: ii, kk
    real(dp)                            :: xsum1   ! sum           of all values per parameter
    real(dp)                            :: xsum2   ! sum of square of all values per parameter
    real(dp)                            :: gsum    ! sum of all (range scaled) currently covered parameter ranges: 
    !                                              !    (max-min)/bound
    real(dp), dimension(size(bound,1))  :: xmin    ! minimal value per parameter
    real(dp), dimension(size(bound,1))  :: xmax    ! maximal value per parameter
    real(dp), dimension(size(bound,1))  :: xmean   ! mean    value per parameter
    real(dp),                 parameter :: delta = tiny(1.0_dp)
    !
    nn   = size(x,2)
    npt1 = size(x,1)

    !  compute maximum, minimum and standard deviation of parameter values
    gsum = 0._dp
    do kk = 1, nn
       if (mask(kk)) then
          xmax(kk) = -huge(1.0_dp)
          xmin(kk) =  huge(1.0_dp)
          xsum1 = 0._dp
          xsum2 = 0._dp
          do ii = 1, npt1
             xmax(kk) = dmax1(x(ii,kk), xmax(kk))
             xmin(kk) = dmin1(x(ii,kk), xmin(kk))
             xsum1 = xsum1 + x(ii,kk)
             xsum2 = xsum2 + x(ii,kk)*x(ii,kk)
          end do
          xmean(kk) = xsum1 / real(npt1,dp)
          xnstd(kk) = (xsum2 / real(npt1,dp) - xmean(kk)*xmean(kk))
          if (xnstd(kk) .le. delta) then
             xnstd(kk) = delta
          end if
          xnstd(kk) = sqrt(xnstd(kk))
          xnstd(kk) = xnstd(kk) / bound(kk)
          gsum = gsum + log( delta + (xmax(kk)-xmin(kk))/bound(kk) )
       end if
    end do
    gnrng = exp(gsum/real(nn,dp))
    !
    !  check if normalized standard deviation of parameter is <= eps
    ipcnvg = 0_i4
    if (gnrng .le. peps) then
       ipcnvg = 1_i4
    end if

  end subroutine parstt

  subroutine comp(ngs2,npg,a,af,b,bf)   
    !
    !  This subroutine reduce 
    !  input matrix a(n,ngs2*npg) to matrix b(n,ngs1*npg) and 
    !  input vector af(ngs2*npg)  to vector bf(ngs1*npg)
    !
    !  Needed for complex reduction: 
    !     Number of complexes decreases from ngs2 to ngs1
    !     Therefore, number of points in population decreases (npt+npg) to npt
    !
    use mo_kind, only: i4, dp

    implicit none  

    integer(i4),                              intent(in)    :: ngs2  ! OLD number of complexes = ngs1+1
    integer(i4),                              intent(in)    :: npg   ! number of points in each complex
    real(dp), dimension(:,:),                 intent(inout) :: a     ! OLD points,          ncols=n, nrows=ngs2*npg
    real(dp), dimension(:),                   intent(inout) :: af    ! OLD function values,          nrows=ngs2*npg
    real(dp), dimension(size(a,1),size(a,2)), intent(out)   :: b     ! NEW points,          ncols=n, nrows=ngs1*npg=npt 
    real(dp), dimension(size(af)),            intent(out)   :: bf    ! NEW function values,          nrows=ngs1*npg=npt
    !
    ! local variables
    integer(i4) :: n      ! cols: number of parameters: nopt
    integer(i4) :: npt    ! rows: number of points in NEW population: npt1
    integer(i4) :: ngs1   ! NEW number of complexes
    integer(i4) :: i, j
    integer(i4) :: igs    ! current new complex
    integer(i4) :: ipg    ! current new point in complex igs
    integer(i4) :: k1, k2

    n = size(a,dim=2)
    ! NEW number of complexes
    ngs1 = ngs2 - 1
    ! number of points in NEW population
    npt  = ngs1 * npg

    do igs=1, ngs1
       do ipg=1, npg
          k1=(ipg-1)*ngs2 + igs
          k2=(ipg-1)*ngs1 + igs
          do i=1, n
             b(k2,i) = a(k1,i)
          end do
          bf(k2) = af(k1)
       end do
    end do
    !
    do j=1, npt
       do i=1, n
          a(j,i) = b(j,i)
       end do
       af(j) = bf(j)
    end do
    !
  end subroutine comp

  subroutine sort_matrix(rb,ra)
    !
    use mo_kind,      only: i4, dp
    use mo_orderpack, only: sort_index

    implicit none    

    real(dp), dimension(:),   intent(inout) :: ra    ! rows: number of points in population --> npt1
    real(dp), dimension(:,:), intent(inout) :: rb    ! rows: number of points in population --> npt1
    !                                                ! cols: number of parameters --> nn
    ! local variables
    integer(i4)                        :: n          ! number of points in population --> npt1
    integer(i4)                        :: m          ! number of parameters --> nn
    integer(i4)                        :: i
    integer(i4), dimension(size(rb,1)) :: iwk
    !
    n = size(rb,1)
    m = size(rb,2)

    ! indexes of sorted reference vector
    iwk(:) = sort_index(ra(1:n))

    ! sort reference vector
    ra(1:n) = ra(iwk)

    ! sort each column of second array
    do i=1, m
       rb(1:n,i) = rb(iwk,i)
    end do

  end subroutine sort_matrix

  subroutine chkcst(x,bl,bu,mask,ibound)
    !
    !     This subroutine check if the trial point satisfies all
    !     constraints.
    !
    !     ibound - violation indicator
    !            = 0  no violation
    !            = 1  violation
    !     nn = number of optimizing variables
    !     ii = the ii'th variable of the arrays x, bl, and bu
    !
    !     Note: removed checking for implicit constraints (was anyway empty)
    !
    use mo_kind, only: i4, dp
    implicit none

    real(dp), dimension(:),    intent(in)  :: x      ! trial point dim=(nn)
    real(dp), dimension(:),    intent(in)  :: bl     ! lower bound dim=(nn)
    real(dp), dimension(:),    intent(in)  :: bu     ! upper bound dim=(nn)
    logical,  dimension(:),    intent(in)  :: mask   ! parameter mask
    integer(i4),               intent(out) :: ibound ! violation indicator

    ! local variables
    integer(i4) :: ii
    integer(i4) :: nn
    !
    ibound = 0_i4
    nn   = size(x,1)

    do ii=1, nn
       if (mask(ii)) then
          if ( (x(ii) .lt. bl(ii)) .or. (x(ii) .gt. bu(ii)) ) then
             ! This constraint is violated   
             ibound = 1_i4
             return
          end if
       end if
    end do

  end subroutine chkcst

  subroutine getpnt(idist,bl,bu,std,xi,mask,save_state,x)
    !
    !     This subroutine generates a new point within feasible region
    !
    !     Note: checking of implicit constraints removed
    !
    use mo_kind, only: i4, i8, dp
    use mo_xor4096, only: n_save_state, xor4096, xor4096g

    implicit none 

    integer(i4),                          intent(in)    :: idist            ! idist = probability flag
    !                                                                       !       = 1 - uniform distribution
    !                                                                       !       = 2 - Gaussian distribution
    real(dp),    dimension(:),            intent(in)    :: bl               ! lower bound
    real(dp),    dimension(:),            intent(in)    :: bu               ! upper bound
    real(dp),    dimension(:),            intent(in)    :: std              ! standard deviation of probability distribution
    real(dp),    dimension(:),            intent(in)    :: xi               ! focal point
    logical,     dimension(:),            intent(in)    :: mask             ! mask of parameters
    integer(i8), dimension(n_save_state), intent(inout) :: save_state       ! save state of random number stream
    !                                                                       ! --> stream has to be according to idist
    real(dp),    dimension(size(xi,1)),   intent(out)   :: x                ! new point
    ! 
    ! local variables
    integer(i4) :: nn    ! number of parameters
    integer(i4) :: jj
    integer(i4) :: ibound   ! 0=point in bound, 1=point out of bound
    real(dp)    :: rand
    !
    nn = size(xi,1)
    do jj=1, nn
       if (mask(jj)) then
          ibound = 1
          do while (ibound .eq. 1)
             if (idist .eq. 1) then
                call xor4096(0_i8, rand, save_state=save_state)
                x(jj) = bl(jj) + std(jj) * rand * (bu(jj) - bl(jj))
             end if
             if (idist .eq. 2) then
                call xor4096g(0_i8, rand, save_state=save_state)  
                x(jj) = xi(jj) + std(jj) * rand * (bu(jj) - bl(jj))
             end if
             !
             !     Check explicit constraints
             !      
             call chkcst((/x(jj)/),(/bl(jj)/),(/bu(jj)/),(/mask(jj)/),ibound)
          end do
       else
          x(jj) = xi(jj)
       end if
    end do

  end subroutine getpnt

  subroutine cce(s,sf,bl,bu,maskpara,xnstd,icall,maxn,maxit,save_state_gauss, functn, &
       alpha,beta,history,idot)
    !
    !  algorithm generate a new point(s) from a sub-complex
    !
    !  Note: new intent IN    variables for flexible reflection & contraction: alpha, beta
    !        new intent INOUT variable for history of objective function values
    !
    use mo_kind,         only: i4, i8, dp
    use mo_xor4096,      only: n_save_state
    use iso_fortran_env, only: output_unit
#ifndef GFORTRAN
    use ieee_arithmetic, only: ieee_is_finite
#endif

    implicit none    

    real(dp),    dimension(:,:),          intent(inout) :: s                ! points in sub-complex, cols=nn, rows=nps
    real(dp),    dimension(:),            intent(inout) :: sf               ! objective function value of points in 
    !                                                                       ! sub-complex, rows=nps
    real(dp),    dimension(:),            intent(in)    :: bl               ! lower bound per parameter
    real(dp),    dimension(:),            intent(in)    :: bu               ! upper bound per parameter
    logical,     dimension(:),            intent(in)    :: maskpara         ! mask of parameters
    real(dp),    dimension(:),            intent(in)    :: xnstd            ! standard deviation of points in 
    !                                                                       ! sub-complex per parameter
    integer(i8),                          intent(inout) :: icall            ! number of function evaluations
    integer(i8),                          intent(in)    :: maxn             ! maximal number of function evaluations allowed
    logical,                              intent(in)    :: maxit            ! minimization (true) or maximization (false)
    integer(i8), dimension(n_save_state), intent(inout) :: save_state_gauss ! seed for random number generation
    interface
       function functn(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind
         real(dp), dimension(:), intent(in)    :: paraset
         real(dp)                              :: functn
       end function functn
    end interface
    real(dp),                    intent(in)    :: alpha   ! parameter for reflection steps
    real(dp),                    intent(in)    :: beta    ! parameter for contraction steps
    real(dp), dimension(maxn),   intent(inout) :: history ! history of best function value
    logical,                     intent(in)    :: idot    ! if true: progress report with '.'
    !
    ! local variables
    integer(i4)                                :: nn      ! number of parameters
    integer(i4)                                :: nps     ! number of points in a sub-complex
    integer(i4)                                :: i, j
    integer(i4)                                :: n       ! nps  = number of points in sub-complex 
    integer(i4)                                :: m       ! nn = number of parameters
    integer(i4)                                :: ibound  ! ibound = flag indicating if constraints are violated
    !                                                     !        = 1   yes
    !                                                     !        = 0   no
    real(dp), dimension(size(s,2))             :: sw      ! the worst point of the simplex
    real(dp), dimension(size(s,2))             :: sb      ! the best point of the simplex
    real(dp), dimension(size(s,2))             :: ce      ! the centroid of the simplex excluding wo
    real(dp), dimension(size(s,2))             :: snew    ! new point generated from the simplex
    real(dp)                                   :: fw      ! function value of the worst point
    real(dp)                                   :: fnew    ! function value of the new point

    nps  = size(s,1)
    nn   = size(s,2)
    !
    ! equivalence of variables for readabilty of code
    n = nps
    m = nn
    !
    ! identify the worst point wo of the sub-complex s
    ! compute the centroid ce of the remaining points
    ! compute step, the vector between wo and ce
    ! identify the worst function value fw
    do j = 1, m
       sb(j) = s(1,j)
       sw(j) = s(n,j)
       ce(j) = 0.0_dp
       do i = 1, n-1
          ce(j) = ce(j) + s(i,j)
       end do
       ce(j) = ce(j)/real(n-1,dp)
    end do
    fw = sf(n)
    !
    ! compute the new point snew
    !
    ! first try a reflection step
    do j = 1, m
       if (maskpara(j)) then
          snew(j) = ce(j) + alpha * (ce(j) - sw(j))
       else
          snew(j) = s(1,j)
       end if
    end do
    !
    ! check if snew satisfies all constraints
    call chkcst(snew(1:nn),bl(1:nn),bu(1:nn),maskpara(1:nn),ibound)
    !
    ! snew is outside the bound,
    ! choose a point at random within feasible region according to
    ! a normal distribution with best point of the sub-complex
    ! as mean and standard deviation of the population as std
    if (ibound .eq. 1) then
       call getpnt(2,bl(1:nn),bu(1:nn),xnstd(1:nn),sb(1:nn),maskpara(1:nn),save_state_gauss,snew)
    end if
    !
    ! compute the function value at snew
    if (idot) write(output_unit,'(A1)') '.'
    if (.not. maxit) then
       fnew = functn(snew)
       icall = icall + 1_i8
       history(icall) = min(history(icall-1),fnew)
    else
       fnew = -functn(snew)
       icall = icall + 1_i8
       history(icall) = max(history(icall-1),-fnew)
    end if
    !
    ! maximum numbers of function evaluations reached
    if (icall .ge. maxn) return
    !
    ! compare fnew with the worst function value fw
    ! fnew is greater than fw, so try a contraction step
#ifndef GFORTRAN
    if ((fnew .gt. fw) .or. (.not. ieee_is_finite(fnew))) then
#else
#ifdef GFORTRAN41
    if ((fnew .gt. fw) .or. (fnew .ne. fnew)) then
#else
    if ((fnew .gt. fw) .or. isnan(fnew)) then
#endif
#endif
       do j = 1, m
          if (maskpara(j)) then
             snew(j) = ce(j) - beta * (ce(j) - sw(j))
          else
             snew(j) = s(1,j)
          end if
       end do
       !
       ! compute the function value of the contracted point
       if (idot) write(output_unit,'(A1)') '.'
       if (.not. maxit) then
          fnew = functn(snew)
          icall = icall + 1_i8
          history(icall) = min(history(icall-1),fnew)
       else
          fnew = -functn(snew)
          icall = icall + 1_i8
          history(icall) = max(history(icall-1),-fnew)
       end if
       !
       ! maximum numbers of function evaluations reached
       if (icall .ge. maxn) return
       !
       ! compare fnew to the worst value fw
       ! if fnew is less than or equal to fw, then accept the point and return
#ifndef GFORTRAN
       if ((fnew .gt. fw) .or. (.not. ieee_is_finite(fnew))) then
#else
#ifdef GFORTRAN41
       if ((fnew .gt. fw) .or. (fnew .ne. fnew)) then
#else
       if ((fnew .gt. fw) .or. isnan(fnew)) then
#endif
#endif
          !
          ! if both reflection and contraction fail, choose another point
          ! according to a normal distribution with best point of the sub-complex
          ! as mean and standard deviation of the population as std
          call getpnt(2,bl(1:nn),bu(1:nn),xnstd(1:nn),sb(1:nn),maskpara(1:nn),save_state_gauss,snew)
          !
          ! compute the function value at the random point
          if (idot) write(output_unit,'(A1)') '.'
          if (.not. maxit) then
             fnew = functn(snew)
             icall = icall + 1_i8
             history(icall) = min(history(icall-1),fnew)
          else
             fnew = -functn(snew)
             icall = icall + 1_i8
             history(icall) = max(history(icall-1),-fnew)
          end if
          !
          ! maximum numbers of function evaluations reached
          if (icall .ge. maxn) return
          !
          ! successful mutation
          ! replace the worst point by the new point
          do j = 1, m
             s(n,j) = snew(j)
          end do
          sf(n) = fnew
       end if
    end if
    !
    ! replace the worst point by the new point
    do j = 1, m
       s(n,j) = snew(j)
    end do
    sf(n) = fnew

  end subroutine cce

END MODULE mo_sce
