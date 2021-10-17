!> \file    mo_mhm_interface.f90
!> \brief   Module providing interfaces for mHM.
!> \details \copydetails mo_mhm_interface

!> \brief   Module providing interfaces for mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Oct 2021
!> \details Interfaces to control the mHM workflow from outside (init, run, get infos, etc.).
module mo_mhm_interface

  use mo_kind, only: i4
  use mo_message, only: message
  use mo_string_utils, only: num2str, separator

  implicit none

  private

  public :: mhm_interface_init

contains

  !> \brief initialize mHM from given namelist paths.
  subroutine mhm_interface_init(namelist_mhm, namelist_mhm_param, namelist_mhm_output, namelist_mrm_output)
    use mo_file, only: &
      file_namelist_mhm, &
      unamelist_mhm, &
      file_namelist_mhm_param, &
      unamelist_mhm_param, &
      file_defOutput
    use mo_mrm_file, only: mrm_file_defOutput => file_defOutput
    use mo_common_read_config, only: &
      common_read_config, &
      check_optimization_settings
    use mo_mhm_read_config, only: mhm_read_config
    use mo_mrm_init, only: &
      mrm_init, &
      mrm_configuration
    use mo_common_variables, only: &
      itimer, &
      domainMeta, &
      global_parameters, &
      global_parameters_name, &
      processMatrix, &
      optimize, &
      opti_function
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
      timestep_model_inputs, &
      L1_twsaObs, &
      L1_etObs, &
      L1_neutronsObs, &
      L1_smObs
    use mo_meteo_forcings, only: prepare_meteo_forcings_data
    use mo_read_optional_data, only: readOptidataObs
    use mo_write_ascii, only: write_configfile

#ifdef MPI
    use mpi_f08
#endif

    implicit none

    character(*), optional, intent(in) :: namelist_mhm !< path to mHM configuration namelist
    character(*), optional, intent(in) :: namelist_mhm_param !< path to mHM parameter namelist
    character(*), optional, intent(in) :: namelist_mhm_output !< path to mHM output namelist
    character(*), optional, intent(in) :: namelist_mrm_output !< path to mRM output namelist

    integer(i4) :: domainID, iDomain

    ! MPI variables
    integer             :: ierror
    integer(i4)         :: nproc, rank, oldrank

    ! reset nml paths if wanted
    if (present(namelist_mhm)) file_namelist_mhm = namelist_mhm
    if (present(namelist_mhm_param)) file_namelist_mhm_param = namelist_mhm_param
    if (present(namelist_mhm_output)) file_defOutput = namelist_mhm_output
    if (present(namelist_mrm_output)) mrm_file_defOutput = namelist_mrm_output

    ! startup message
    call startup_message()

    ! read configs
    call common_read_config(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)
#ifdef MPI
    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(domainMeta%comMaster, rank, ierror)
#endif
    call mhm_read_config(file_namelist_mhm, unamelist_mhm)
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
    if (processMatrix(8, 1) > 0) &
        call mrm_init(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)

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
  end subroutine mhm_interface_init

end module mo_mhm_interface
