!> \file    mo_mhm_interface.f90
!> \brief   \copybrief mo_mhm_interface
!> \details \copydetails mo_mhm_interface

!> \brief   Module providing interfaces for mHM.
!> \details Interfaces to control the mHM workflow from outside (init, run, get infos, etc.).
!> \authors Sebastian Mueller
!> \version 0.1
!> \date    Oct 2021
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
module mo_mhm_interface

  use mo_kind, only: i4, dp
  use mo_message, only: message, error_message
  use mo_string_utils, only: num2str

#ifdef MPI
  use mpi_f08
#endif

  implicit none

  private

  public :: mhm_interface_init
  public :: mhm_interface_get_parameter
  public :: mhm_interface_get_parameter_number
  public :: mhm_interface_run
  public :: mhm_interface_run_optimization
  public :: mhm_interface_finalize

contains

  !> \brief initialize mHM from given namelist paths.
  subroutine mhm_interface_init(namelist_mhm, namelist_mhm_param, namelist_mhm_output, namelist_mrm_output, cwd)
    use mo_file, only: &
      file_namelist_mhm, &
      unamelist_mhm, &
      file_namelist_mhm_param, &
      unamelist_mhm_param, &
      file_defOutput
    use mo_mrm_file, only: mrm_file_defOutput => file_defOutput
    use mo_common_read_config, only: &
      common_read_config
    use mo_common_mHM_mRM_read_config, only : &
      common_mHM_mRM_read_config, &
      check_optimization_settings
    use mo_mpr_read_config, only: mpr_read_config
    use mo_mhm_read_config, only: mhm_read_config
    use mo_read_wrapper, only : read_data
    use mo_mrm_init, only: &
      mrm_init, &
      mrm_configuration
    use mo_common_variables, only: &
      level0, &
      level1, &
      itimer, &
      domainMeta, &
      processMatrix
    use mo_common_mHM_mRM_variables, only : &
      timeStep, &
      simPer, &
      optimize, &
      opti_function, &
      mrm_coupling_mode, &
      read_restart
    use mo_mhm_messages, only: &
      startup_message, &
      domain_dir_check_message
    use mo_timer, only: &
      timers_init, &
      timer_start, &
      timer_stop, &
      timer_get
    use mo_startup, only: mhm_initialize
    use mo_file, only: &
      unamelist_mhm, &
      unamelist_mhm_param
    use mo_global_variables, only: &
      couple_cfg, &
      meteo_handler, &
      L1_twsaObs, &
      L1_etObs, &
      L1_neutronsObs, &
      L1_smObs, &
      BFI_calc
    use mo_read_optional_data, only: readOptidataObs
    use mo_write_ascii, only: write_configfile
    use mo_mhm_bfi, only: calculate_BFI
    use mo_os, only: change_dir

    implicit none

    character(*), optional, intent(in) :: namelist_mhm !< path to mHM configuration namelist
    character(*), optional, intent(in) :: namelist_mhm_param !< path to mHM parameter namelist
    character(*), optional, intent(in) :: namelist_mhm_output !< path to mHM output namelist
    character(*), optional, intent(in) :: namelist_mrm_output !< path to mRM output namelist
    character(*), optional, intent(in) :: cwd !< desired working directory

    integer(i4) :: domainID, iDomain

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank
#endif

    ! reset nml paths if wanted
    if (present(namelist_mhm)) file_namelist_mhm = namelist_mhm
    if (present(namelist_mhm_param)) file_namelist_mhm_param = namelist_mhm_param
    if (present(namelist_mhm_output)) file_defOutput = namelist_mhm_output
    if (present(namelist_mrm_output)) mrm_file_defOutput = namelist_mrm_output
    ! change working directory
    if (present(cwd)) call change_dir(cwd)

    ! startup message
    call startup_message()

    ! coupling configuration
    call couple_cfg%read_config(file_namelist_mhm, unamelist_mhm)
    ! read configs
    call common_read_config(file_namelist_mhm, unamelist_mhm)
    call mpr_read_config(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)
    call common_mHM_mRM_read_config(file_namelist_mhm, unamelist_mhm)
    call mhm_read_config(file_namelist_mhm, unamelist_mhm)
    call couple_cfg%check(domainMeta, optimize)
    call meteo_handler%config(file_namelist_mhm, unamelist_mhm, optimize, domainMeta, processMatrix, timestep, couple_cfg)
    mrm_coupling_mode = 2_i4 ! TODO: this shouldn't be needed
    call mrm_configuration(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)
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
    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
    ! ComLocal is a communicator, i.e. a group of processes assigned to the same
    ! domain, with a master and subprocesses. Only the master processes of these
    ! groups need to read the data. The master process with rank 0 only
    ! coordinates the other processes and does not need to read the data.
    if (rank > 0 .and. domainMeta%isMasterInComLocal) then
#endif
    call message()

    if (.not. read_restart) then
      call message('  Read data ...')
      call timer_start(itimer)
      ! for DEM, slope, ... define nGvar local
      ! read_data has a domain loop inside
      call read_data(simPer)
      call timer_stop(itimer)
      call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    end if

    ! read data for every domain
    itimer = itimer + 1
    call message('  Initialize domains ...')
    call timer_start(itimer)
    call mhm_initialize()
    call meteo_handler%init_level2(level0, level1)
    call timer_stop(itimer)
    call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    if (processMatrix(8, 1) > 0) &
        call mrm_init(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)

    itimer = itimer + 1
    call message('  Read forcing and optional data ...')
    call timer_start(itimer)

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      ! read meteorology now, if it should be loaded in one go
      if (meteo_handler%single_read(iDomain)) call meteo_handler%prepare_data(1, iDomain, level1, simPer)

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

    ! calculate observed BFI if wanted
    if ( optimize .and. opti_function==34 .and. BFI_calc ) call calculate_BFI()

    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

    !this call may be moved to another position as it writes the master config out file for all domains
    call write_configfile(meteo_handler%dirPrecipitation, meteo_handler%dirReferenceET, meteo_handler%dirTemperature)

#ifdef MPI
    end if
#endif

  end subroutine mhm_interface_init

  !> \brief Get current global parameter value of mHM.
  subroutine mhm_interface_get_parameter(para)
    use mo_common_variables, only: global_parameters

    implicit none

    real(dp), dimension(:), allocatable, intent(out) :: para !< global parameter values of mHM

    allocate(para(size(global_parameters, dim=1)))

    para = global_parameters(:, 3)

  end subroutine mhm_interface_get_parameter

  !> \brief Get number of current global parameter value of mHM.
  subroutine mhm_interface_get_parameter_number(n)
    use mo_common_variables, only: global_parameters

    implicit none

    integer(i4), intent(out) :: n !< number of global parameter values of mHM

    n = size(global_parameters, dim=1)

  end subroutine mhm_interface_get_parameter_number

  !> \brief Run mHM with current settings.
  subroutine mhm_interface_run()
    use mo_common_variables, only: &
#ifdef MPI
      domainMeta, &
#endif
      itimer, &
      global_parameters
    use mo_timer, only: &
      timer_start, &
      timer_stop, &
      timer_get
    use mo_mhm_eval, only: mhm_eval

    implicit none

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank

    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
#endif

    itimer = itimer + 1

    ! --------------------------------------------------------------------------
    ! call mHM
    ! get runoff timeseries if possible (i.e. when domainMeta%doRouting,
    ! processMatrix(8,1) > 0)
    ! get other model outputs  (i.e. gridded fields of model output)
    ! --------------------------------------------------------------------------

#ifdef MPI
    if (rank > 0 .and. domainMeta%isMasterInComLocal) then
#endif

    call message('  Run mHM')
    call timer_start(itimer)
    call mhm_eval(global_parameters(:, 3))
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F12.3)')), ' seconds.')

#ifdef MPI
    endif
#endif

  end subroutine mhm_interface_run

  !> \brief Run mHM optimization with current settings.
  subroutine mhm_interface_run_optimization()
    use mo_common_variables, only: &
#ifdef MPI
      domainMeta, &
#endif
      itimer, &
      dirConfigOut, &
      global_parameters,&
      global_parameters_name, &
      processMatrix
    use mo_common_mHM_mRM_variables, only : &
      opti_function
    use mo_timer, only: &
      timer_start, &
      timer_stop, &
      timer_get
    use mo_mhm_eval, only: mhm_eval
    use mo_objective_function, only: &
#ifdef MPI
      objective_subprocess, &
      objective_master, &
#endif
      objective
    use mo_optimization, only: optimization
    use mo_mrm_objective_function_runoff, only: &
#ifdef MPI
      single_objective_runoff_master, &
      single_objective_runoff_subprocess, &
#endif
      single_objective_runoff
    use mo_write_ascii, only: &
      write_optifile, &      ! Writing optimized parameter set and objective
      write_optinamelist     ! Writing optimized parameter set to a namelist

    implicit none

    procedure(mhm_eval), pointer :: eval
    procedure(objective), pointer :: obj_func

    real(dp) :: funcbest         ! best objective function achivied during optimization
    logical, dimension(:), allocatable :: maskpara ! true  = parameter will be optimized, = parameter(i,4) = 1
    !                                              ! false = parameter will not be optimized = parameter(i,4) = 0

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank

    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
#endif

    itimer = itimer + 1
    call message('  Run mHM optimization')
    call timer_start(itimer)

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

      case(10 : 13, 15, 17, 27, 28, 29, 30, 33, 34)
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

    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F12.3)')), ' seconds.')

  end subroutine mhm_interface_run_optimization

  !> \brief Write mHM restart.
  subroutine mhm_interface_finalize()
    use mo_common_variables, only: &
#ifdef MPI
      domainMeta, &
#endif
      itimer, &
      mhmFileRestartOut, &
      write_restart, &
      processMatrix
    use mo_common_mHM_mRM_variables, only : &
      optimize
    use mo_timer, only: &
      timer_start, &
      timer_stop, &
      timer_get
    use mo_restart, only: write_restart_files
    use mo_mrm_write, only : mrm_write
    use mo_mhm_messages, only: finish_message
    use mo_clean_up, only: deallocate_global_variables

    implicit none

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank

    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
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

    call finish_message()

    ! clean up all allocated variables
    call deallocate_global_variables()

  end subroutine mhm_interface_finalize

end module mo_mhm_interface
