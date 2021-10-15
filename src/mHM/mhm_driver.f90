!>       \file mhm_driver.f90

!>       \brief Distributed precipitation-runoff model mHM

!>       \details This is the main driver of mHM, which calls
!>       one instance of mHM for a multiple domains and a given period.
!>       \image html  mhm5-logo.png "Typical mHM cell"
!>       \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>       \authors Luis Samaniego & Rohini Kumar (UFZ)

!>       \date Jun 2018

!>       \version \htmlinclude version.txt \latexinclude version.txt

!>       \copyright (c) \f$2005 - \the\year{}\f$, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!!       All rights reserved.
!!
!!       This code is a property of:
!!
!!       ----------------------------------------------------------
!!
!!       Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ
!!       Registered Office: Leipzig
!!       Registration Office: Amtsgericht Leipzig
!!       Trade Register: Nr. B 4703
!!       Chairman of the Supervisory Board: MinDirig Wilfried Kraus
!!       Scientific Director: Prof. Dr. Georg Teutsch
!!       Administrative Director: Dr. Heike Grassmann
!!
!!       ----------------------------------------------------------
!!
!!       NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!!       EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!!       OF THIS SOFTWARE. If software is modified to produce
!!       derivative works, such modified software should be
!!       clearly marked, so as not to confuse it with the version
!!       available from UFZ.  This code can be used for research
!!       purposes ONLY provided that the following sources are
!!       acknowledged:
!!
!!       Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!!       parameter regionalization of a grid-based hydrologic
!!       model at the mesoscale.  Water Resour. Res., 46,
!!       W05523, doi:10.1029/2008WR007327.
!!
!!       Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!!       of distributed hydrologic model parameterization on water
!!       fluxes at multiple scales and locations, Water Resour. Res.,
!!       49, doi:10.1029/2012WR012195.
!!
!!       For commercial applications you have to consult the
!!       authorities of the UFZ.

! Modifications:
! Stephan Thober                Nov 2013 - added read in of latitude longitude fields
! Matthias Zink                 Mar 2013 - edited screen output for gauges added inflow gauges
! Matthias Cuntz & Juliane Mai  Mar 2014 - Likelihood Kavetski uses 2 more parameters for the
!                                          error model global_parameters -> local_parameters
! Rohini Kumar                  Apr 2014 - implementation of the mHM run on a single cell configuration
!                                          also for the routing mode.
!                                        - run mHM at the input data level i.e. L0 grid
! Rohini Kumar                  May 2014 - model run on a regular lat-lon grid or on a regular X-Y coordinate system
! Stephan Thober                May 2014 - moved read meteo forcings to mo_mhm_eval
! Matthias Cuntz & Juliane Mai  Nov 2014 - LAI input from daily, monthly or yearly files
! Matthias Zink                 Mar 2015 - added optional soil mositure read in for calibration
! Luis Samaniego                Jul 2015 - added temporal directories for optimization
! Stephan Thober                Aug 2015 - removed routing related variables
! Stephan Thober                Oct 2015 - reorganized optimization (now compatible with mRM)
! Oldrich Rakovec, Rohini Kumar Oct 2015 - added reading of domain averaged TWS and objective function 15
!                                          for simultaneous calibration based on runoff and TWS
! Rohini Kumar                  Mar 2016 - options to handle different soil databases modified MPR to included
!                                          soil horizon specific properties/parameters
! Stephan Thober                Nov 2016 - implemented adaptive timestep for routing
! Rohini Kumar                  Dec 2016 - options to read (monthly mean) LAI fields
! Robert Schweppe               Jun 2018 - refactoring and reformatting
! Maren Kaluza                  Oct 2019 - TWS to data structure
! M.C. Demirel, Simon Stisen    Jun 2020 - New Soil Moisture Process: Feddes and FC dependency on root fraction coefficient processCase(3) = 4
PROGRAM mhm_driver

  use mo_file, only: &
          file_namelist_mhm, unamelist_mhm, &      ! filename of namelist: main setup
          file_namelist_mhm_param, unamelist_mhm_param    ! filename of namelist: mhm model parameter
  use mo_finish, only: finish                         ! Finish with style
  use mo_global_variables, only: &
          timestep_model_inputs, & !frequency of input read
          L1_twsaObs, &
          L1_etObs, &
          L1_neutronsObs, &
          L1_smObs
  use mo_optimization_types, only: &
          optidata ! type for opti data
  use mo_common_variables, only: &
          optimize, opti_function, &                                   ! optimization on/off and optimization method
          write_restart, &      ! restart writing flags
          mhmFileRestartOut, &
          dirConfigOut, &
          domainMeta, &
#ifdef MPI
          comm, &
#endif
          processMatrix, &      ! domain information,  processMatrix
          global_parameters, global_parameters_name      ! mhm parameters (gamma) and their clear names
  use mo_kind, only: i4, dp                         ! number precision
  use mo_message, only : error_message, message          ! For print out
  use mo_meteo_forcings, only: prepare_meteo_forcings_data
  use mo_mhm_eval, only: mhm_eval
  use mo_read_optional_data, only: readOptidataObs ! read optional observed data
  use mo_common_read_config, only: common_read_config, &       ! Read main configuration files
                                    check_optimization_settings ! Read main configuration files
  use mo_mhm_read_config, only: mhm_read_config                    ! Read main configuration files
  use mo_restart, only: write_restart_files
  use mo_startup, only: mhm_initialize
  use mo_string_utils, only: num2str             ! String magic
  use mo_timer, only: &
          timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  use mo_write_ascii, only: &
          write_configfile, &      ! Writing Configuration file
          write_optifile, &      ! Writing optimized parameter set and objective
          write_optinamelist     ! Writing optimized parameter set to a namelist
  use mo_objective_function, only: &
#ifdef MPI
          objective_subprocess, &
          objective_master, &
#endif
          objective                 ! objective functions and likelihoods
  use mo_optimization, only: optimization
  use mo_mrm_objective_function_runoff, only: &
#ifdef MPI
          single_objective_runoff_master, &
          single_objective_runoff_subprocess, &
#endif
          single_objective_runoff
  use mo_mrm_init, only: mrm_init, mrm_configuration
  use mo_mrm_write, only : mrm_write

  !$ use omp_lib, only: OMP_GET_NUM_THREADS           ! OpenMP routines
#ifdef MPI
  use mpi_f08
#endif

  use mo_check, only: check_dir
  use mo_mhm_cli, only: parse_command_line
  use mo_mhm_messages, only: startup_message, domain_dir_check_message, finish_message

  IMPLICIT NONE

  ! local
  !$ integer(i4)                        :: n_threads        ! OpenMP number of parallel threads
  integer(i4) :: domainID, iDomain               ! Counters
  integer(i4) :: itimer           ! Current timer number
  real(dp) :: funcbest         ! best objective function achivied during optimization
  logical, dimension(:), allocatable :: maskpara ! true  = parameter will be optimized, = parameter(i,4) = 1
  !                                              ! false = parameter will not be optimized = parameter(i,4) = 0
  procedure(mhm_eval), pointer :: eval
  procedure(objective), pointer :: obj_func

  logical :: ReadLatLon

#ifdef MPI
  integer             :: ierror
  integer(i4)         :: nproc
  integer(i4)         :: rank, oldrank
  logical :: compiled_with_mpi = .true.

  ! Initialize MPI
  call MPI_Init(ierror)
  call MPI_Comm_dup(MPI_COMM_WORLD, comm, ierror)
  ! find number of processes nproc
  call MPI_Comm_size(comm, nproc, ierror)
  ! find the number the process is referred to, called rank
  call MPI_Comm_rank(comm, rank, ierror)
  oldrank = rank
  write(*,*) 'MPI!, comm', rank, nproc
#else
  logical :: compiled_with_mpi = .false.
#endif

  ! parse command line arguments
  call parse_command_line()
  ! startup message
  call startup_message(compiled_with_mpi)

  ! read configs
  call common_read_config(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)
#ifdef MPI
  call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
  ! find the number the process is referred to, called rank
  call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
#endif
  call mhm_read_config(file_namelist_mhm, unamelist_mhm)
  call mrm_configuration(file_namelist_mhm, unamelist_mhm, &
          file_namelist_mhm_param, unamelist_mhm_param, ReadLatLon)
  call check_optimization_settings()

  ! Message about input directories
  call domain_dir_check_message()

  ! Start timings
  call timers_init

  ! --------------------------------------------------------------------------
  ! READ AND INITIALIZE
  ! --------------------------------------------------------------------------
  itimer = 1
#ifdef MPI
  ! ComLocal is a communicator, i.e. a group of processes assigned to the same
  ! domain, with a master and subprocesses. Only the master processes of these
  ! groups need to read the data. The master process with rank 0 only
  ! coordinates the other processes and does not need to read the data.
  if (rank > 0 .and. domainMeta%isMasterInComLocal) then
#endif
  call message()

  ! read data for every domain
  itimer = itimer + 1
  call message('  Initialize domains ...')
  call timer_start(itimer)
  call mhm_initialize(global_parameters(:, 3), global_parameters_name)
  call timer_stop(itimer)
  call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  if (processMatrix(8, 1) > 0) call mrm_init(file_namelist_mhm, unamelist_mhm, &
          file_namelist_mhm_param, unamelist_mhm_param, ReadLatLon)

  itimer = itimer + 1
  call message('  Read forcing and optional data ...')
  call timer_start(itimer)

  do iDomain = 1, domainMeta%nDomains
    domainID = domainMeta%indices(iDomain)
    ! read meteorology now, if optimization is switched on
    ! meteorological forcings (reading, upscaling or downscaling)
    if (timestep_model_inputs(iDomain) .eq. 0_i4) then
      call prepare_meteo_forcings_data(iDomain, 1)
    end if

    ! read optional optional data if necessary
    if (optimize) then
      select case (opti_function)
        case(10 : 13, 28)
          ! read optional spatio-temporal soil mositure data
          call readOptidataObs(iDomain, domainID, L1_smObs(iDomain))
        case(17)
          ! read optional spatio-temporal neutrons data
          call readOptidataObs(iDomain, domainID, L1_neutronsObs(iDomain))
        case(27, 29, 30)
          ! read optional spatio-temporal evapotranspiration data
          call readOptidataObs(iDomain, domainID, L1_etObs(iDomain))
        case(15)
          ! read optional spatio-temporal tws data
          call readOptidataObs(iDomain, domainID, L1_twsaObs(iDomain))
        case(33)
          ! read optional spatio-temporal evapotranspiration data
          if (domainMeta%optidata(iDomain) == 0 .or. domainMeta%optidata(iDomain) == 5 .or. &
            domainMeta%optidata(iDomain) == 6 ) then
            call readOptidataObs(iDomain, domainID, L1_etObs(iDomain))
          end if
          ! read optional spatio-temporal tws data
          if (domainMeta%optidata(iDomain) == 0 .or. domainMeta%optidata(iDomain) == 3 .or. &
            domainMeta%optidata(iDomain) == 6 ) then
            call readOptidataObs(iDomain, domainID, L1_twsaObs(iDomain))
          end if
      end select
    end if

  end do
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

  !this call may be moved to another position as it writes the master config out file for all domains
  call write_configfile()

#ifdef MPI
  end if
#endif
  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  itimer = itimer + 1
  call message()
  if (optimize) then
    eval => mhm_eval

    select case(opti_function)
     case(1 : 9, 14, 31 : 32)
      ! call optimization against only runoff (no other variables)
      obj_func => single_objective_runoff
#ifdef MPI
      if (rank == 0 .and. domainMeta%isMasterInComLocal) then
        obj_func => single_objective_runoff_master
        call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
      else if (domainMeta%isMasterInComLocal) then
        ! In case of a master process from ComLocal, i.e. a master of a group of
        ! processes that are assigned to a single domain, this process calls the
        ! objective subroutine directly. The master over all processes collects
        ! the data and runs the dds/sce/other opti method.
        call single_objective_runoff_subprocess(eval)
      end if
#else
      call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
#endif
     case(10 : 13, 15, 17, 27, 28, 29, 30, 33)
      ! call optimization for other variables
      obj_func => objective
#ifdef MPI
      if (rank == 0 .and. domainMeta%isMasterInComLocal) then
        obj_func => objective_master
        call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
      else if (domainMeta%isMasterInComLocal) then
        ! In case of a master process from ComLocal, i.e. a master of a group of
        ! processes that are assigned to a single domain, this process calls the
        ! objective subroutine directly. The master over all processes collects
        ! the data and runs the dds/sce/other opti method.
        call objective_subprocess(eval)
      end if
#else
      call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
#endif
    case default
      call error_message('***ERROR: mhm_driver: The given objective function number ', &
              trim(adjustl(num2str(opti_function))), ' in mhm.nml is not valid!')
    end select
#ifdef MPI
  if (rank == 0 .and. domainMeta%isMasterInComLocal) then
#endif
    ! write a file with final objective function and the best parameter set
    call write_optifile(funcbest, global_parameters(:, 3), global_parameters_name(:))
    ! write a file with final best parameter set in a namlist format
    call write_optinamelist(processMatrix, global_parameters, maskpara, global_parameters_name(:))
    deallocate(maskpara)
#ifdef MPI
  end if
#endif
  else

#ifdef MPI
    if (rank > 0 .and. domainMeta%isMasterInComLocal) then
#endif
      ! --------------------------------------------------------------------------
      ! call mHM
      ! get runoff timeseries if possible (i.e. when domainMeta%doRouting,
      ! processMatrix(8,1) > 0)
      ! get other model outputs  (i.e. gridded fields of model output)
      ! --------------------------------------------------------------------------
      call message('  Run mHM')
      call timer_start(itimer)
      call mhm_eval(global_parameters(:, 3))
      call timer_stop(itimer)
      call message('    in ', trim(num2str(timer_get(itimer), '(F12.3)')), ' seconds.')
#ifdef MPI
    endif
#endif

  end if

#ifdef MPI
  if (rank > 0 .and. domainMeta%isMasterInComLocal) then
#endif
  ! --------------------------------------------------------------------------
  ! WRITE RESTART files
  ! --------------------------------------------------------------------------
  if (write_restart  .AND. (.NOT. optimize)) then
    itimer = itimer + 1
    call message()
    call message('  Write restart file')
    call timer_start(itimer)
    call write_restart_files(mhmFileRestartOut)
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  end if

  ! --------------------------------------------------------------------------
  ! WRITE RUNOFF (INCLUDING RESTART FILES, has to be called after mHM restart
  ! files are written)
  ! --------------------------------------------------------------------------
  if (processMatrix(8, 1) > 0) call mrm_write()

#ifdef MPI
  end if
#endif

  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call finish_message()

#ifdef MPI
  ! find number of processes nproc
  call MPI_Comm_size(comm, nproc, ierror)
  call MPI_Comm_rank(comm, rank, ierror)
  write(*,*) 'MPI finished', rank, nproc
  call MPI_Finalize(ierror)
#endif

END PROGRAM mhm_driver
