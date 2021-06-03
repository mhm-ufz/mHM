!>       \file mo_common_read_config.f90

!>       \brief Reading of main model configurations.

!>       \details This routine reads the configurations of namelists commonly used by mHM, mRM and MPR

!>       \authors Matthias Zink

!>       \date Dec 2012

! Modifications:

MODULE mo_common_read_config

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: common_read_config, set_land_cover_scenes_id, common_check_resolution, check_optimization_settings
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        common_read_config

  !    PURPOSE
  !>       \brief Read main configurations commonly used by mHM and mRM

  !>       \details Read the main configurations commonly used by mHM and mRM, namely:
  !>       project_description, directories_general, mainconfig, processSelection, LCover

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist" name of file
  !>       \param[in] "integer :: unamelist"          id of file

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Dec  2018 - refactoring and restructuring

  subroutine common_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_common_constants, only : maxNLcovers, maxNoDomains
    use mo_common_variables, only : Conventions, LC_year_end, LC_year_start, LCfilename, contact, &
                                    dirCommonFiles, dirConfigOut, dirLCover, dirMorpho, dirOut, &
                                    mhmFileRestartOut, mrmFileRestartOut, &
                                    fileLatLon, history, mHM_details, domainMeta, nLandCoverPeriods, &
                                    nProcesses, nuniqueL0Domains, processMatrix, project_details, resolutionHydrology, &
                                    setup_description, simulation_type, write_restart, &
                                    dds_r, mhmFileRestartIn, mrmFileRestartIn, evalPer,&
                                    mcmc_error_params, mcmc_opti, nIterations, &
                                    opti_function, opti_method, optimize, optimize_restart, &
                                    read_restart, mrm_read_river_network, resolutionRouting, sa_temp, &
                                    sce_ngs, sce_npg, sce_nps, seed, &
                                    warmPer, warmingDays
    use mo_common_datetime_type, only: LCyearId, simPer, timestep, nTStepDay, period
    use mo_julian, only : caldat, julday
    use mo_message, only : message
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str
    use mo_grid, only : iFlag_coordinate_sys


    implicit none

    ! name of file
    character(*), intent(in) :: file_namelist
    ! id of file
    integer, intent(in) :: unamelist
    ! name of file
    character(*), intent(in) :: file_namelist_param
    ! id of file
    integer, intent(in) :: unamelist_param

    ! Choosen process description number
    integer(i4), dimension(nProcesses) :: processCase

    character(256), dimension(maxNoDomains) :: dir_Morpho

    character(256), dimension(maxNoDomains) :: mhm_file_RestartOut

    character(256), dimension(maxNoDomains) :: mrm_file_RestartOut

    character(256), dimension(maxNoDomains) :: dir_LCover

    character(256), dimension(maxNoDomains) :: dir_Out

    character(256), dimension(maxNoDomains) :: file_LatLon

    real(dp), dimension(maxNoDomains) :: resolution_Hydrology

    integer(i4), dimension(maxNoDomains) :: L0Domain

    integer(i4), dimension(maxNoDomains) :: read_opt_domain_data

    ! starting year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart

    ! ending year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd

    ! filename of Lcover file
    character(256), dimension(maxNLCovers) :: LCoverfName

    integer(i4) :: i, newDomainID, domainID, iDomain, nDomains

    ! flag to advance nuniqueL0Domain counter
    logical :: addCounter

    integer(i4) :: jday

    integer(i4), dimension(maxNoDomains) :: warming_Days

    type(period), dimension(maxNoDomains) :: eval_Per

    real(dp), dimension(maxNoDomains) :: resolution_Routing

    character(256), dimension(maxNoDomains) :: mhm_file_RestartIn
    character(256), dimension(maxNoDomains) :: mrm_file_RestartIn

    ! define namelists
    ! namelist directories
    namelist /project_description/ project_details, setup_description, simulation_type, &
            Conventions, contact, mHM_details, history
    namelist /directories_general/ dirConfigOut, dirCommonFiles, &
            dir_Morpho, dir_LCover, &
            dir_Out, mhm_file_RestartOut, mrm_file_RestartOut, &
            file_LatLon
    ! namelist spatial & temporal resolution, optimization information
    namelist /mainconfig/ iFlag_coordinate_sys, resolution_Hydrology, nDomains, L0Domain, write_restart, &
            read_opt_domain_data
    ! namelist process selection
    namelist /processSelection/ processCase

    ! namelist for land cover scenes
    namelist/LCover/nLandCoverPeriods, LCoverYearStart, LCoverYearEnd, LCoverfName
    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig_mhm_mrm/ timestep, resolution_Routing, optimize, &
            optimize_restart, opti_method, opti_function, &
            read_restart, mrm_read_river_network, mhm_file_RestartIn, mrm_file_RestartIn
    ! namelist for optimization settings
    namelist /Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, &
            sce_npg, sce_nps, mcmc_opti, mcmc_error_params
    ! namelist for time settings
    namelist /time_periods/ warming_Days, eval_Per

    ! set default values for optional arguments
    mrm_read_river_network = .false.

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist specifying the project description
    !===============================================================
    call position_nml('project_description', unamelist)
    read(unamelist, nml = project_description)

    !===============================================================
    !  Read namelist specifying the model configuration
    !===============================================================
    call position_nml('mainconfig', unamelist)
    read(unamelist, nml = mainconfig)

    call init_domain_variable(nDomains, read_opt_domain_data(1:nDomains), domainMeta)

    if (nDomains .GT. maxNoDomains) then
      call message()
      call message('***ERROR: Number of domains is resticted to ', trim(num2str(maxNoDomains)), '!')
      stop 1
    end if

    ! allocate patharray sizes
    allocate(resolutionHydrology(domainMeta%nDomains))
    allocate(mhmFileRestartOut(domainMeta%nDomains))
    allocate(mrmFileRestartOut(domainMeta%nDomains))
    allocate(dirOut(domainMeta%nDomains))

    ! TODO: MPR this block will go
    allocate(dirMorpho(domainMeta%nDomains))
    allocate(dirLCover(domainMeta%nDomains))
    allocate(fileLatLon(domainMeta%nDomains))

    nuniqueL0Domains = 0_i4
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      resolutionHydrology(iDomain) = resolution_Hydrology(domainID)
      ! if a domain uses the same L0 data as a previous one, write
      ! the index into domainMeta%L0DataFrom
      newDomainID = L0Domain(domainID)
      domainMeta%L0DataFrom(iDomain) = iDomain
      !
      addCounter = .True.
      do i = 1, iDomain - 1
        if (newDomainID == domainMeta%indices(i)) then
          domainMeta%L0DataFrom(iDomain) = i
          addCounter = .False.
        end if
      end do
      if (addCounter) nuniqueL0Domains = nuniqueL0Domains + 1_i4
    end do

    ! check for possible options
    if(.NOT. (iFlag_coordinate_sys == 0 .OR. iFlag_coordinate_sys == 1)) then
      call message()
      call message('***ERROR: coordinate system for the model run should be 0 or 1')
      stop 1
    end if

    ! TODO: MPR this will go
    !===============================================================
    ! Read land cover
    !===============================================================
    call position_nml('LCover', unamelist)
    read(unamelist, nml = LCover)
    ! put land cover scenes to corresponding file name and LuT
    ! this is done already here for MPR, which does not check for the time periods
    allocate(LCfilename(nLandCoverPeriods))
    allocate(LC_year_start(nLandCoverPeriods))
    allocate(LC_year_end(nLandCoverPeriods))
    LCfilename(:) = LCoverfName(1 : nLandCoverPeriods)
    LC_year_start(:) = LCoverYearStart(1 : nLandCoverPeriods)
    LC_year_end(:) = LCoverYearEnd(1 : nLandCoverPeriods)

    !===============================================================
    ! Read namelist for mainpaths
    !===============================================================
    call position_nml('directories_general', unamelist)
    read(unamelist, nml = directories_general)

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      domainMeta%optidata(iDomain) = read_opt_domain_data(domainID)
      mhmFileRestartOut(iDomain)   = mhm_file_RestartOut(domainID)
      mrmFileRestartOut(iDomain)   = mrm_file_RestartOut(domainID)
      dirOut(iDomain)              = dir_Out(domainID)
      ! TODO: MPR this will go
      dirMorpho(iDomain)           = dir_Morpho(domainID)
      dirLCover(iDomain)           = dir_LCover(domainID)
      fileLatLon(iDomain)          = file_LatLon(domainID)
    end do

    ! ToDo: add test if opti_function matches at least one domainMeta%optidata
    ! as soon as common and common_mRM_mHM are merged, if that is the plan


    !===============================================================
    ! Read process selection list
    !===============================================================
    ! init the processCase matrix to 0 to be backward compatible
    ! if cases were added later (then there would be no values if not init here)
    processCase = 0_i4
    call position_nml('processselection', unamelist)
    read(unamelist, nml = processSelection)

    processMatrix = 0_i4
    processMatrix(:, 1) = processCase
    if (processMatrix(8, 1) == 0) then
      domainMeta%doRouting(:) = .FALSE.
    else
      domainMeta%doRouting(:) = .TRUE.
    end if

    call position_nml('mainconfig_mhm_mrm', unamelist)
    read(unamelist, nml = mainconfig_mhm_mrm)
    ! consistency between read_restart and mrm_read_river_network
    if (read_restart) then
       if (.not. mrm_read_river_network) then
          call message('***WARNING: mrm_read_river_network is set to .true. because read_restart is .true.')
       end if
       mrm_read_river_network = .true.
    end if

    allocate(resolutionRouting(domainMeta%nDomains))
    allocate(mhmFileRestartIn(domainMeta%nDomains))
    allocate(mrmFileRestartIn(domainMeta%nDomains))
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      mhmFileRestartIn(iDomain) = mhm_file_RestartIn(domainID)
      mrmFileRestartIn(iDomain) = mrm_file_RestartIn(domainID)
      resolutionRouting(iDomain) = resolution_Routing(domainID)
    end do

    ! check for optimize and read restart
    if ((read_restart) .and. (optimize)) then
      call message()
      call message('***ERROR: cannot read states from restart file when optimizing')
      stop 1
    end if

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      if (processMatrix(8, 1) > 0 .and. domainMeta%optidata(iDomain) > 1 .and. optimize) then
        domainMeta%doRouting(iDomain) = .FALSE.
        call message('Warning: although defined in namelist, routing is switched off for domain', trim(num2str(domainID)))
        call message('         since the calibration of Q is not possible with the chosen opti input')
      end if
    end do

    !===============================================================
    !  INIT !!! (merged from mo_startup and mo_mrm_read_config)
    !===============================================================
    ! transformation of time units & constants
    if (mod(24, timeStep) > 0) then
      call message('mo_startup: timeStep must be a divisor of 24: ', num2str(timeStep))
      stop 1
    end if
    nTStepDay = 24_i4 / timeStep            ! # of time steps per day

    ! allocate time periods
    allocate(simPer(domainMeta%nDomains))
    allocate(evalPer(domainMeta%nDomains))
    allocate(warmingDays(domainMeta%nDomains))
    allocate(warmPer(domainMeta%nDomains))

    !===============================================================
    !  read simulation time periods incl. warming days
    !===============================================================
    call position_nml('time_periods', unamelist)
    read(unamelist, nml = time_periods)
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      warmingDays(iDomain) = warming_Days(domainID)
      call evalPer(iDomain)%init(eval_Per(domainID)%dStart, eval_Per(domainID)%mStart, eval_Per(domainID)%yStart, &
                                 eval_Per(domainID)%dEnd, eval_Per(domainID)%mEnd, eval_Per(domainID)%yEnd)
    end do
    ! evalPer = eval_Per(1 : domainMeta%nDomains)

    !===============================================================
    !  determine simulation time period incl. warming days for each
    !  domain
    !===============================================================
    do iDomain = 1, domainMeta%nDomains
      ! julian days for evaluation period
      jday = julday(dd = evalPer(iDomain)%dStart, mm = evalPer(iDomain)%mStart, yy = evalPer(iDomain)%yStart)
      evalPer(iDomain)%julStart = jday

      jday = julday(dd = evalPer(iDomain)%dEnd, mm = evalPer(iDomain)%mEnd, yy = evalPer(iDomain)%yEnd)
      evalPer(iDomain)%julEnd = jday

      ! determine warming period
      warmPer(iDomain)%julStart = evalPer(iDomain)%julStart - warmingDays(iDomain)
      warmPer(iDomain)%julEnd = evalPer(iDomain)%julStart - 1

      call caldat(warmPer(iDomain)%julStart, dd = warmPer(iDomain)%dStart, mm = warmPer(iDomain)%mStart, &
                  yy = warmPer(iDomain)%yStart)
      call caldat(warmPer(iDomain)%julEnd, dd = warmPer(iDomain)%dEnd, mm = warmPer(iDomain)%mEnd, &
                  yy = warmPer(iDomain)%yEnd)

      ! simulation Period = warming Period + evaluation Period
      simPer(iDomain)%dStart = warmPer(iDomain)%dStart
      simPer(iDomain)%mStart = warmPer(iDomain)%mStart
      simPer(iDomain)%yStart = warmPer(iDomain)%yStart
      simPer(iDomain)%julStart = warmPer(iDomain)%julStart
      simPer(iDomain)%dEnd = evalPer(iDomain)%dEnd
      simPer(iDomain)%mEnd = evalPer(iDomain)%mEnd
      simPer(iDomain)%yEnd = evalPer(iDomain)%yEnd
      simPer(iDomain)%julEnd = evalPer(iDomain)%julEnd
    end do

    ! TODO: MPR this will be exchanged
    ! !===============================================================
    ! ! Read land cover
    ! !===============================================================
    ! call position_nml('LCover', unamelist)
    ! read(unamelist, nml = LCover)
    ! allocate(LCyearId(minval(simPer(1:domainMeta%nDomains)%yStart):maxval(simPer(1:domainMeta%nDomains)%yEnd), domainMeta%nDomains))
    call set_land_cover_scenes_id(simPer, LCyearId)

    !===============================================================
    ! Settings for Optimization
    !===============================================================
    ! namelist for Optimization settings
    call position_nml('Optimization', unamelist)
    read(unamelist, nml = Optimization)
    ! checking of settings and default value initialization moved to new subroutine
    ! because global_parameters need to be set, which is not the case right now

    call close_nml(unamelist)

    call read_mhm_parameters(file_namelist, unamelist, file_namelist_param, unamelist_param)

  end subroutine common_read_config

    ! ------------------------------------------------------------------

  !    NAME
  !        set_land_cover_scenes_id

  !    PURPOSE
  !>       \brief Read main configurations commonly used by mHM, mRM and MPR

  !>       \details Read the main configurations commonly used by mHM, mRM and MPR, namely:
  !>       project_description, directories_general, mainconfig, processSelection, LCover

  !    INTENT(IN)
  !>       \param[in] "type(period), dimension(:) :: sim_Per"

  !    INTENT(INOUT)
  !>       \param[inout] "integer(i4), dimension(:, :) :: LCyear_Id"
  !>       \param[inout] "character(256), dimension(:) :: LCfilename"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Dec  2018 - refactoring and restructuring

  subroutine set_land_cover_scenes_id(sim_Per, LCyear_Id)

    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : LC_year_end, LC_year_start, domainMeta, nLandCoverPeriods
    use mo_common_datetime_type, only: period
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    type(period), dimension(:), intent(in) :: sim_Per

    integer(i4), dimension(:, :), allocatable, intent(inout) :: LCyear_Id

    integer(i4) :: ii, iDomain


    ! countercheck if land cover covers simulation period
    if (LC_year_start(1) .GT. minval(sim_Per(1 : domainMeta%nDomains)%yStart)) then
      call message()
      call message('***ERROR: Land cover for warming period is missing!')
      call message('   SimStart   : ', trim(num2str(minval(sim_Per(1 : domainMeta%nDomains)%yStart))))
      call message('   LCoverStart: ', trim(num2str(LC_year_start(1))))
      stop 1
    end if
    if (LC_year_end(nLandCoverPeriods) .LT. maxval(sim_Per(1 : domainMeta%nDomains)%yEnd)) then
      call message()
      call message('***ERROR: Land cover period shorter than modelling period!')
      call message('   SimEnd   : ', trim(num2str(maxval(sim_Per(1 : domainMeta%nDomains)%yEnd))))
      call message('   LCoverEnd: ', trim(num2str(LC_year_end(nLandCoverPeriods))))
      stop 1
    end if
    !
    allocate(LCyear_Id(minval(sim_Per(1 : domainMeta%nDomains)%yStart) : maxval(sim_Per(1 : domainMeta%nDomains)%yEnd), &
                   domainMeta%nDomains))
    LCyear_Id = nodata_i4
    do iDomain = 1, domainMeta%nDomains
      do ii = 1, nLandCoverPeriods
        ! land cover before model period or land cover after model period
        if ((LC_year_end(ii) .LT. sim_Per(iDomain)%yStart) .OR. &
                (LC_year_start(ii) .GT. sim_Per(iDomain)%yEnd)) then
          cycle
          ! land cover period fully covers model period
        else if ((LC_year_start(ii) .LE. sim_Per(iDomain)%yStart) .AND. &
                (LC_year_end(ii) .GE. sim_Per(iDomain)%yEnd)) then
          LCyear_Id(sim_Per(iDomain)%yStart : sim_Per(iDomain)%yEnd, iDomain) = ii
          exit
          ! land cover period covers beginning of model period
        else if ((LC_year_start(ii) .LE. sim_Per(iDomain)%yStart) .AND. &
                (LC_year_end(ii) .LT. sim_Per(iDomain)%yEnd)) then
          LCyear_Id(sim_Per(iDomain)%yStart : LC_year_end(ii), iDomain) = ii
          ! land cover period covers end of model period
        else if ((LC_year_start(ii) .GT. sim_Per(iDomain)%yStart) .AND. &
                (LC_year_end(ii) .GE. sim_Per(iDomain)%yEnd)) then
          LCyear_Id(LC_year_start(ii) : sim_Per(iDomain)%yEnd, iDomain) = ii
          ! land cover period covers part of model_period
        else
          LCyear_Id(LC_year_start(ii) : LC_year_end(ii), iDomain) = ii
        end if
      end do
    end do


  end subroutine set_land_cover_scenes_id

!< author: Maren Kaluza
!< date: September 2019
!< summary: Initialization of the domain variable for all domain loops and if activated for parallelization

!< In case of MPI parallelization domainMeta%overAllNumberOfDomains is a
!< variable where the number of domains from the namelist is stored. By this
!< every process knows the total number of domains. Then, in a loop the
!< domains are distributed onto the processes. There is a master process
!< and several subprocesses. The master process only reads the confings in the
!< mHM driver.
!<
!< The subprocesses get a number of domains. domainMeta%nDomain refers
!< to the number of domains assigned to a specific process. It is a local
!< variable and therefore has a different value for each process.
!<
!< In case more domains are there than processes, currently the domains
!< are distributed round robin, i.e. like cards in a card game.
!<
!< In case less domains than processes exist, all remaining processes
!< are assigned to the routing domains round robin. In that case the
!< local communicator is of interest: It is a group of processes assigned
!< to a routing domain again with a master process
!< (domainMeta%isMasterInComLocal) and subprocesses. This communicator can
!< in future be passed to the routing parallelization.
  subroutine init_domain_variable(nDomains, optiData, domainMeta)
    use mo_common_variables, only: domain_meta
#ifdef MPI
    use mo_common_variables, only: comm
    use mpi_f08
#endif
    integer(i4),       intent(in)    :: nDomains
    integer(i4), dimension(:), intent(in) :: optiData
    type(domain_meta), intent(inout) :: domainMeta

    integer             :: ierror
    integer(i4)         :: nproc
    integer(i4)         :: rank
    integer(i4)         :: iDomain
    integer(i4)         :: colDomain, colMasters

    domainMeta%overallNumberOfDomains = nDomains
#ifdef MPI
    ! find number of processes nproc
    call MPI_Comm_size(comm, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(comm, rank, ierror)
    if (nproc < 2) then
      stop 'at least 2 processes are required'
    end if
    ! if there are more processes than domains
    if (nproc > domainMeta%overallNumberOfDomains + 1) then
      domainMeta%nDomains = 0
      ! master reads only metadata of all domains
      if (rank == 0) then
        call init_domain_variable_for_master(domainMeta, colMasters, colDomain)
      ! all other nodes read metadata but also data of assigned domains
      else
        ! currently each domain gets one process except it is a routing domain.
        ! in that case the remaining processes are distributed round robin to
        ! the routing domains.
        call distribute_processes_to_domains_according_to_role(optiData, rank, &
                                               domainMeta, colMasters, colDomain)
      end if
      ! two communicators are created, i.e. groups of processes that talk about
      ! a certain topic:
      ! comMaster is the communicator of all processes that need to read all
      ! data. These are the processes that are masters in the comLocal plus the
      ! master over all processes. comLocal is a communicator for a group of
      ! processes assigned to the same domain.
      call MPI_Comm_split(comm, colMasters, rank, domainMeta%comMaster, ierror)
      call MPI_Comm_split(comm, colDomain, rank, domainMeta%comLocal, ierror)
      call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    else
      ! in case of more domains than processes, distribute domains round robin
      ! onto the processes
      call MPI_Comm_dup(comm, domainMeta%comMaster, ierror)
      domainMeta%isMasterInComLocal = .true.
      domainMeta%nDomains = 0
      ! master reads only metadata of all domains
      if (rank == 0) then
        domainMeta%nDomains = domainMeta%overallNumberOfDomains
        call domainMeta%allocate_domains()
        do iDomain = 1, domainMeta%nDomains
          domainMeta%indices(iDomain) = iDomain
        end do
      ! all other nodes read metadata but also data of assigned domains
      else
        call distributeDomainsRoundRobin(nproc, rank, domainMeta)
      end if
    end if ! round robin
#else
    domainMeta%nDomains = nDomains
    call domainMeta%allocate_domains()
    do iDomain = 1, domainMeta%nDomains
      domainMeta%indices(iDomain) = iDomain
    end do
#endif

  end subroutine init_domain_variable

#ifdef MPI
  subroutine init_domain_variable_for_master(domainMeta, colMasters, colDomain)
    use mo_common_variables, only: domain_meta
    type(domain_meta), intent(inout) :: domainMeta
    integer(i4),       intent(out)   :: colMasters
    integer(i4),       intent(out)   :: colDomain
    !local
    integer(i4) :: iDomain

    domainMeta%nDomains = domainMeta%overallNumberOfDomains
    call domainMeta%allocate_domains()
    do iDomain = 1, domainMeta%nDomains
      domainMeta%indices(iDomain) = iDomain
    end do
    colMasters = 1
    colDomain = 0
    domainMeta%isMasterInComLocal = .true.

  end subroutine init_domain_variable_for_master


  subroutine distributeDomainsRoundRobin(nproc, rank, domainMeta)
    use mo_common_variables, only: domain_meta
    integer(i4),       intent(in)    :: nproc
    integer(i4),       intent(in)    :: rank
    type(domain_meta), intent(inout) :: domainMeta

    integer(i4) :: iDomain, iProcDomain

    do iDomain = 1 , domainMeta%overallNumberOfDomains
      if (rank == (modulo(iDomain + nproc - 2, (nproc - 1)) + 1)) then
        domainMeta%nDomains = domainMeta%nDomains + 1
      end if
    end do
    call domainMeta%allocate_domains()
    iProcDomain = 0
    do iDomain = 1 , domainMeta%overallNumberOfDomains
      if (rank == (modulo(iDomain + nproc - 2, (nproc - 1)) + 1)) then
        iProcDomain = iProcDomain + 1
        domainMeta%indices(iProcDomain) = iDomain
      end if
    end do
  end subroutine distributeDomainsRoundRobin

  subroutine distribute_processes_to_domains_according_to_role(optiData, rank, &
                                               domainMeta, colMasters, colDomain)
    use mo_common_variables, only: domain_meta
    integer(i4), dimension(:), intent(in)    :: optiData
    integer(i4),               intent(in)    :: rank
    type(domain_meta),         intent(inout) :: domainMeta
    integer(i4),               intent(out)   :: colMasters
    integer(i4),               intent(out)   :: colDomain

    ! local
    integer(i4) :: nDomainsAll, nTreeDomains, i, iDomain
    integer(i4), dimension(:), allocatable :: treeDomainList

    nDomainsAll = domainMeta%overallNumberOfDomains
    nTreeDomains = 0
    do iDomain = 1, nDomainsAll
      !ToDo: should also routing but not opti domains be counted?
      if (optiData(iDomain) == 1) then
        nTreeDomains = nTreeDomains + 1
      end if
    end do
    allocate(treeDomainList(nTreeDomains))
    i = 0
    do iDomain = 1, nDomainsAll
      if (optiData(iDomain) == 1) then
        i = i + 1
        treeDomainList(i) = iDomain
      end if
    end do
    if (rank < nDomainsAll + 1) then
      colMasters = 1
      colDomain = rank
      domainMeta%isMasterInComLocal = .true.
      domainMeta%nDomains = 1
      call domainMeta%allocate_domains()
      domainMeta%indices(1) = rank
    else
      colMasters = 0
      if (nTreeDomains > 0) then
        colDomain = treeDomainList(mod(rank, nTreeDomains) + 1)
      else
        colDomain = 1
      end if
      domainMeta%isMasterInComLocal = .false.
      domainMeta%nDomains = 1
      call domainMeta%allocate_domains()
      ! ToDo : temporary solution, this should either not read data at all
      ! or data corresponding to the master process
      domainMeta%indices(1) = 1
    end if
    deallocate(treeDomainList)
  end subroutine
#endif

    subroutine check_optimization_settings

    use mo_common_variables, only : dds_r, nIterations, sce_ngs, sce_npg, sce_nps, global_parameters
    use mo_message, only : message

    implicit none

    integer(i4) :: n_true_pars


    ! check and set default values
    if (nIterations .le. 0_i4) then
      call message('Number of iterations for Optimization (nIterations) must be greater than zero')
      stop 1
    end if
    if (dds_r .lt. 0.0_dp .or. dds_r .gt. 1.0_dp) then
      call message('dds_r must be between 0.0 and 1.0')
      stop 1
    end if
    if (sce_ngs .lt. 1_i4) then
      call message ('number of complexes in SCE (sce_ngs) must be at least 1')
      stop 1
    end if
    ! number of points in each complex: default = 2n+1
    if (sce_npg .lt. 0_i4) then
      n_true_pars = count(nint(global_parameters(:, 4)) .eq. 1)
      sce_npg = 2 * n_true_pars + 1_i4
    end if
    ! number of points in each sub-complex: default = n+1
    if (sce_nps .lt. 0_i4) then
      n_true_pars = count(nint(global_parameters(:, 4)) .eq. 1)
      sce_nps = n_true_pars + 1_i4
    end if
    if (sce_npg .lt. sce_nps) then
      call message ('number of points per complex (sce_npg) must be greater or')
      call message ('equal number of points per sub-complex (sce_nps)')
      stop 1
    end if

  end subroutine check_optimization_settings

  subroutine common_check_resolution(do_message, allow_subgrid_routing)

    use mo_common_variables, only : resolutionRouting
    use mo_common_variables, only : domainMeta, resolutionHydrology
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    logical, intent(in) :: do_message

    logical, intent(in) :: allow_subgrid_routing

    integer(i4) :: iDomain, domainID

    ! conversion factor L11 to L1
    real(dp) :: cellFactorRbyH


    !===============================================================
    ! check matching of resolutions: hydrology, forcing and routing
    !===============================================================
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      cellFactorRbyH = resolutionRouting(iDomain) / resolutionHydrology(iDomain)
      if (do_message) then
        call message()
        call message('domain ', trim(adjustl(num2str(domainID))), ': ')
        call message('resolution Hydrology (domain ', trim(adjustl(num2str(domainID))), ')     = ', &
                trim(adjustl(num2str(resolutionHydrology(iDomain)))))
        call message('resolution Routing (domain ', trim(adjustl(num2str(domainID))), ')       = ', &
                trim(adjustl(num2str(resolutionRouting(iDomain)))))
      end if
      !
      if(nint(cellFactorRbyH * 100.0_dp) .eq. 100) then
        if (do_message) then
          call message()
          call message('Resolution of routing and hydrological modeling are equal!')
        end if

      else if ((nint(cellFactorRbyH * 100.0_dp) .gt. 100) .and. .not.allow_subgrid_routing) then
        if(nint(mod(cellFactorRbyH, 2.0_dp) * 100.0_dp) .ne. 0) then
          call message()
          call message('***ERROR: Resolution of routing is not a multiple of hydrological model resolution!')
          call message('   FILE: mhm.nml, namelist: mainconfig, variable: resolutionRouting')
          STOP
        end if
        !
        if (do_message) then
          call message()
          call message('Resolution of routing is bigger than hydrological model resolution by ', &
                  trim(adjustl(num2str(nint(cellFactorRbyH)))), ' times !')
        end if
      end if
      !
    end do

  end subroutine common_check_resolution

  subroutine read_mhm_parameters(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_append, only : append
    use mo_common_constants, only : eps_dp, maxNoDomains, nColPars, nodata_dp
    use mo_common_functions, only : in_bound
    use mo_common_variables, only : global_parameters, global_parameters_name, domainMeta, processMatrix, dummy_global_parameters, &
          dummy_global_parameters_name
    use mo_message, only : message
    use mo_mpr_constants, only : maxGeoUnit, &
                                 maxNoSoilHorizons
    use mo_mpr_global_variables, only : HorizonDepth_mHM, dirgridded_LAI, fracSealed_cityArea, iFlag_soilDB, &
                                        inputFormat_gridded_LAI, nGeoUnits, nSoilHorizons_mHM, tillageDepth
    use mo_common_datetime_type, only : timeStep_LAI_input
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str
    use mo_utils, only : EQ
    use mo_global_variables, only: soilHorizonBoundaries, nSoilHorizons

    implicit none

    character(*), intent(in) :: file_namelist

    integer, intent(in) :: unamelist

    character(*), intent(in) :: file_namelist_param

    integer, intent(in) :: unamelist_param

    integer(i4) :: ii

    ! depth of the single horizons
    real(dp), dimension(maxNoSoilHorizons) :: soil_Depth

    ! directory of gridded LAI data
    ! used when timeStep_LAI_input<0
    character(256), dimension(maxNoDomains) :: dir_gridded_LAI

    character(256) :: dummy

    ! space holder for routing parameters
    real(dp), dimension(5, nColPars) :: dummy_2d_dp

    ! space holder for routing parameters
    real(dp), dimension(1, nColPars) :: dummy_2d_dp_2

    real(dp), dimension(nColPars) :: canopyInterceptionFactor

    real(dp), dimension(nColPars) :: snowThresholdTemperature

    real(dp), dimension(nColPars) :: degreeDayFactor_forest

    real(dp), dimension(nColPars) :: degreeDayFactor_impervious

    real(dp), dimension(nColPars) :: degreeDayFactor_pervious

    real(dp), dimension(nColPars) :: increaseDegreeDayFactorByPrecip

    real(dp), dimension(nColPars) :: maxDegreeDayFactor_forest

    real(dp), dimension(nColPars) :: maxDegreeDayFactor_impervious

    real(dp), dimension(nColPars) :: maxDegreeDayFactor_pervious

    real(dp), dimension(nColPars) :: orgMatterContent_forest

    real(dp), dimension(nColPars) :: orgMatterContent_impervious

    real(dp), dimension(nColPars) :: orgMatterContent_pervious

    real(dp), dimension(nColPars) :: PTF_lower66_5_constant

    real(dp), dimension(nColPars) :: PTF_lower66_5_clay

    real(dp), dimension(nColPars) :: PTF_lower66_5_Db

    real(dp), dimension(nColPars) :: PTF_higher66_5_constant

    real(dp), dimension(nColPars) :: PTF_higher66_5_clay

    real(dp), dimension(nColPars) :: PTF_higher66_5_Db

    real(dp), dimension(nColPars) :: infiltrationShapeFactor

    real(dp), dimension(nColPars) :: PTF_Ks_constant

    real(dp), dimension(nColPars) :: PTF_Ks_sand

    real(dp), dimension(nColPars) :: PTF_Ks_clay

    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope

    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest

    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious

    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious

    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1

    real(dp), dimension(nColPars) :: FCmin_glob

    real(dp), dimension(nColPars) :: FCdelta_glob

    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand

    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay

    real(dp), dimension(nColPars) :: imperviousStorageCapacity

    real(dp), dimension(nColPars) :: PET_a_forest

    real(dp), dimension(nColPars) :: PET_a_impervious

    real(dp), dimension(nColPars) :: PET_a_pervious

    real(dp), dimension(nColPars) :: PET_b

    real(dp), dimension(nColPars) :: PET_c

    real(dp), dimension(nColPars) :: minCorrectionFactorPET

    real(dp), dimension(nColPars) :: maxCorrectionFactorPET

    real(dp), dimension(nColPars) :: aspectThresholdPET

    real(dp), dimension(nColPars) :: minCorrectionFactorPET_HS

    real(dp), dimension(nColPars) :: maxCorrectionFactorPET_HS

    real(dp), dimension(nColPars) :: aspectThresholdPET_HS

    real(dp), dimension(nColPars) :: HargreavesSamaniCoeff

    real(dp), dimension(nColPars) :: PriestleyTaylorCoeff

    real(dp), dimension(nColPars) :: PriestleyTaylorLAIcorr

    real(dp), dimension(nColPars) :: canopyheight_forest

    real(dp), dimension(nColPars) :: canopyheight_impervious

    real(dp), dimension(nColPars) :: canopyheight_pervious

    real(dp), dimension(nColPars) :: displacementheight_coeff

    real(dp), dimension(nColPars) :: roughnesslength_momentum_coeff

    real(dp), dimension(nColPars) :: roughnesslength_heat_coeff

    real(dp), dimension(nColPars) :: stomatal_resistance

    real(dp), dimension(nColPars) :: interflowStorageCapacityFactor

    real(dp), dimension(nColPars) :: interflowRecession_slope

    real(dp), dimension(nColPars) :: fastInterflowRecession_forest

    real(dp), dimension(nColPars) :: slowInterflowRecession_Ks

    real(dp), dimension(nColPars) :: exponentSlowInterflow

    real(dp), dimension(nColPars) :: rechargeCoefficient

    real(dp), dimension(nColPars) :: rechargeFactor_karstic

    real(dp), dimension(nColPars) :: gain_loss_GWreservoir_karstic

    real(dp), dimension(maxGeoUnit, nColPars) :: GeoParam

    real(dp), dimension(nColPars) :: Desilets_N0

    real(dp), dimension(nColPars) :: COSMIC_N0

    real(dp), dimension(nColPars) :: COSMIC_N1

    real(dp), dimension(nColPars) :: COSMIC_N2

    real(dp), dimension(nColPars) :: COSMIC_alpha0

    real(dp), dimension(nColPars) :: COSMIC_alpha1

    real(dp), dimension(nColPars) :: COSMIC_L30

    real(dp), dimension(nColPars) :: COSMIC_L31

    integer(i4) :: iDomain, domainID


    ! namelist directories
    namelist /directories_MPR/ dir_gridded_LAI
    ! namelist soil database
    namelist /soildata/ iFlag_soilDB, tillageDepth, nSoilHorizons_mHM, soil_Depth
    ! namelist for LAI related data
    namelist /LAI_data_information/ inputFormat_gridded_LAI, timeStep_LAI_input
    ! namelist for land cover scenes
    namelist /LCover_MPR/ fracSealed_cityArea

    ! namelist parameters
    namelist /mhm_parameters/ canopyInterceptionFactor, snowThresholdTemperature, degreeDayFactor_forest, &
            degreeDayFactor_impervious, &
            degreeDayFactor_pervious, increaseDegreeDayFactorByPrecip, maxDegreeDayFactor_forest, &
            maxDegreeDayFactor_impervious, maxDegreeDayFactor_pervious, orgMatterContent_forest, &
            orgMatterContent_impervious, orgMatterContent_pervious, &
            PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant, &
            PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant, &
            PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope, &
            rootFractionCoefficient_forest, rootFractionCoefficient_impervious, &
            rootFractionCoefficient_pervious, infiltrationShapeFactor, jarvis_sm_threshold_c1, &
            FCmin_glob, FCdelta_glob, &
            rootFractionCoefficient_sand, rootFractionCoefficient_clay, imperviousStorageCapacity, &
            PET_a_forest, PET_a_impervious, PET_a_pervious, PET_b, PET_c, minCorrectionFactorPET, &
            maxCorrectionFactorPET, aspectThresholdPET, &
            minCorrectionFactorPET_HS, maxCorrectionFactorPET_HS, aspectThresholdPET_HS, HargreavesSamaniCoeff, &
            PriestleyTaylorCoeff, PriestleyTaylorLAIcorr, canopyheight_forest, canopyheight_impervious, &
            canopyheight_pervious, displacementheight_coeff, &
            roughnesslength_momentum_coeff, roughnesslength_heat_coeff, stomatal_resistance, &
            interflowStorageCapacityFactor, interflowRecession_slope, fastInterflowRecession_forest, &
            slowInterflowRecession_Ks, exponentSlowInterflow, &
            rechargeCoefficient, rechargeFactor_karstic, gain_loss_GWreservoir_karstic, &
            Desilets_N0, COSMIC_N0, COSMIC_N1, COSMIC_N2, COSMIC_alpha0, COSMIC_alpha1, COSMIC_L30, COSMIC_L31, &
            GeoParam

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    soil_Depth = 0.0_dp
    dummy_2d_dp = nodata_dp
    dummy_2d_dp_2 = nodata_dp
    GeoParam = nodata_dp

    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist for LCover
    !===============================================================
    call position_nml('LCover_MPR', unamelist)
    read(unamelist, nml = LCover_MPR)

    !===============================================================
    ! Read soil layering information
    !===============================================================
    call position_nml('soildata', unamelist)
    read(unamelist, nml = soildata)

    nSoilHorizons = nSoilHorizons_mHM
    allocate(soilHorizonBoundaries(nSoilHorizons))
    allocate(HorizonDepth_mHM(nSoilHorizons))
    soilHorizonBoundaries(:) = 0.0_dp
    ! last layer is reset to 0 in MPR in case of iFlag_soilDB is 0
    soilHorizonBoundaries(1 : nSoilHorizons) = soil_Depth(1 : nSoilHorizons)

    ! counter checks -- soil horizons
    if (nSoilHorizons .GT. maxNoSoilHorizons) then
      call message()
      call message('***ERROR: Number of soil horizons is resticted to ', trim(num2str(maxNoSoilHorizons)), '!')
      stop
    end if

    ! the default is the HorizonDepths are all set up to last
    ! as is the default for option-1 where horizon specific information are taken into consideration
    if(iFlag_soilDB .eq. 0) then
      ! classical mhm soil database
      soilHorizonBoundaries(nSoilHorizons) = 0.0_dp
    else if(iFlag_soilDB .ne. 1) then
      call message()
      call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
      stop
    end if
    ! TODO: MPR remove this duplications
    HorizonDepth_mHM = soilHorizonBoundaries
    ! some consistency checks for the specification of the tillage depth
    if(iFlag_soilDB .eq. 1) then
      if(count(abs(soilHorizonBoundaries(:) - tillageDepth) .lt. eps_dp)  .eq. 0) then
        call message()
        call message('***ERROR: Soil tillage depth must conform with one of the specified horizon (lower) depth.')
        stop
      end if
    end if

    !===============================================================
    ! Read LAI related information
    !===============================================================
    call position_nml('LAI_data_information', unamelist)
    read(unamelist, nml = LAI_data_information)

    if (timeStep_LAI_input .ne. 0) then
      !===============================================================
      !  Read namelist for main directories
      !===============================================================
      call position_nml('directories_MPR', unamelist)
      read(unamelist, nml = directories_MPR)

      allocate(dirgridded_LAI(domainMeta%nDomains))
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)
        dirgridded_LAI(iDomain) = dir_gridded_LAI(domainID)
      end do

      if (timeStep_LAI_input .GT. 1) then
        call message()
        call message('***ERROR: option for selected timeStep_LAI_input not coded yet')
        stop
      end if
    end if

    call close_nml(unamelist)

    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    call open_nml(file_namelist_param, unamelist_param, quiet = .true.)
    ! decide which parameters to read depending on specified processes
    call position_nml('mhm_parameters', unamelist_param)
    read(unamelist_param, nml = mhm_parameters)

    ! Process 1 - interception
    select case (processMatrix(1, 1))
      ! 1 - maximum Interception
    case(1)
      processMatrix(1, 2) = 1_i4
      processMatrix(1, 3) = 1_i4
      call append(global_parameters, reshape(canopyInterceptionFactor, [1, nColPars]))

      call append(global_parameters_name, [  &
              'canopyInterceptionFactor'])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "interception" does not exist!')
      stop 1
    end select

    ! Process 2 - snow
    select case (processMatrix(2, 1))
      ! 1 - degree-day approach
    case(1)
      processMatrix(2, 2) = 8_i4
      processMatrix(2, 3) = sum(processMatrix(1 : 2, 2))
      call append(global_parameters, reshape(snowThresholdTemperature, [1, nColPars]))
      call append(global_parameters, reshape(degreeDayFactor_forest, [1, nColPars]))
      call append(global_parameters, reshape(degreeDayFactor_impervious, [1, nColPars]))
      call append(global_parameters, reshape(degreeDayFactor_pervious, [1, nColPars]))
      call append(global_parameters, reshape(increaseDegreeDayFactorByPrecip, [1, nColPars]))
      call append(global_parameters, reshape(maxDegreeDayFactor_forest, [1, nColPars]))
      call append(global_parameters, reshape(maxDegreeDayFactor_impervious, [1, nColPars]))
      call append(global_parameters, reshape(maxDegreeDayFactor_pervious, [1, nColPars]))

      call append(global_parameters_name, [  &
                      'snowThresholdTemperature       ', &
                      'degreeDayFactor_forest         ', &
                      'degreeDayFactor_impervious     ', &
                      'degreeDayFactor_pervious       ', &
                      'increaseDegreeDayFactorByPrecip', &
                      'maxDegreeDayFactor_forest      ', &
                      'maxDegreeDayFactor_impervious  ', &
                      'maxDegreeDayFactor_pervious    '])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "snow" does not exist!')
      stop 1
    end select

    ! Process 3 - soilmoisture
    call append(global_parameters, reshape(orgMatterContent_forest, [1, nColPars]))
    call append(global_parameters, reshape(orgMatterContent_impervious, [1, nColPars]))
    call append(global_parameters, reshape(orgMatterContent_pervious, [1, nColPars]))
    call append(global_parameters, reshape(PTF_lower66_5_constant, [1, nColPars]))
    call append(global_parameters, reshape(PTF_lower66_5_clay, [1, nColPars]))
    call append(global_parameters, reshape(PTF_lower66_5_Db, [1, nColPars]))
    call append(global_parameters, reshape(PTF_higher66_5_constant, [1, nColPars]))
    call append(global_parameters, reshape(PTF_higher66_5_clay, [1, nColPars]))
    call append(global_parameters, reshape(PTF_higher66_5_Db, [1, nColPars]))
    call append(global_parameters, reshape(PTF_Ks_constant, [1, nColPars]))
    call append(global_parameters, reshape(PTF_Ks_sand, [1, nColPars]))
    call append(global_parameters, reshape(PTF_Ks_clay, [1, nColPars]))
    call append(global_parameters, reshape(PTF_Ks_curveSlope, [1, nColPars]))
    call append(global_parameters, reshape(rootFractionCoefficient_forest, [1, nColPars]))
    call append(global_parameters, reshape(rootFractionCoefficient_impervious, [1, nColPars]))
    call append(global_parameters, reshape(rootFractionCoefficient_pervious, [1, nColPars]))
    call append(global_parameters, reshape(infiltrationShapeFactor, [1, nColPars]))

    call append(global_parameters_name, [     &
              'orgMatterContent_forest           ', &
                      'orgMatterContent_impervious       ', &
                      'orgMatterContent_pervious         ', &
                      'PTF_lower66_5_constant            ', &
                      'PTF_lower66_5_clay                ', &
                      'PTF_lower66_5_Db                  ', &
                      'PTF_higher66_5_constant           ', &
                      'PTF_higher66_5_clay               ', &
                      'PTF_higher66_5_Db                 ', &
                      'PTF_Ks_constant                   ', &
                      'PTF_Ks_sand                       ', &
                      'PTF_Ks_clay                       ', &
                      'PTF_Ks_curveSlope                 ', &
                      'rootFractionCoefficient_forest    ', &
                      'rootFractionCoefficient_impervious', &
                      'rootFractionCoefficient_pervious  ', &
                    'infiltrationShapeFactor           '])

    select case (processMatrix(3, 1))

      ! 1 - Feddes equation for PET reduction, bucket approach, Brooks-Corey like
    case(1)
      processMatrix(3, 2) = 17_i4
      processMatrix(3, 3) = sum(processMatrix(1:3, 2))

      call append(dummy_global_parameters, [ &
              rootFractionCoefficient_sand(3), &
              rootFractionCoefficient_clay(3), &
              jarvis_sm_threshold_c1(3), &
              FCmin_glob(3), &
              FCdelta_glob(3) &
      ])
      call append(dummy_global_parameters_name, [&
              'rootFractionCoefficient_sand      ', &
                      'rootFractionCoefficient_clay      ', &
                      'jarvis_sm_threshold_c1            ', &
                      'FCmin_glob                        ', &
                      'FCdelta_glob                      ' &
      ])
      ! 2- Jarvis equation for PET reduction, bucket approach, Brooks-Corey like
    case(2)
      processMatrix(3, 2) = 18_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(jarvis_sm_threshold_c1, [1, nColPars]))
      call append(global_parameters_name, ['jarvis_sm_threshold_c1            '])
      call append(dummy_global_parameters, [ &
              rootFractionCoefficient_sand(3), &
              rootFractionCoefficient_clay(3), &
              FCmin_glob(3), &
              FCdelta_glob(3) &
      ])
      call append(dummy_global_parameters_name, [&
              'rootFractionCoefficient_sand      ', &
                      'rootFractionCoefficient_clay      ', &
                      'FCmin_glob                        ', &
                      'FCdelta_glob                      ' &
      ])

      ! 3- Jarvis equation for ET reduction and FC dependency on root fraction coefficient
    case(3)
      processMatrix(3, 2) = 22_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(rootFractionCoefficient_sand, [1, nColPars]))
      call append(global_parameters, reshape(rootFractionCoefficient_clay, [1, nColPars]))
      call append(global_parameters, reshape(FCmin_glob, [1, nColPars]))
      call append(global_parameters, reshape(FCdelta_glob, [1, nColPars]))
      call append(global_parameters, reshape(jarvis_sm_threshold_c1, [1, nColPars]))

      call append(global_parameters_name, [ &
                      'rootFractionCoefficient_sand      ', &
                      'rootFractionCoefficient_clay      ', &
                      'FCmin_glob                        ', &
                      'FCdelta_glob                      ', &
                      'jarvis_sm_threshold_c1            ' &
])
      ! 4- Feddes equation for ET reduction and FC dependency on root fraction coefficient
    case(4)
      processMatrix(3, 2) = 21_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(rootFractionCoefficient_sand, [1, nColPars]))
      call append(global_parameters, reshape(rootFractionCoefficient_clay, [1, nColPars]))
      call append(global_parameters, reshape(FCmin_glob, [1, nColPars]))
      call append(global_parameters, reshape(FCdelta_glob, [1, nColPars]))

      call append(global_parameters_name, [ &
                      'rootFractionCoefficient_sand      ', &
                      'rootFractionCoefficient_clay      ', &
                      'FCmin_glob                        ', &
                      'FCdelta_glob                      ' &
      ])
      call append(dummy_global_parameters, [ &
              jarvis_sm_threshold_c1(3) &
      ])
      call append(dummy_global_parameters_name, [&
              'jarvis_sm_threshold_c1            ' &
      ])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "soilmoisture" does not exist!')
      stop
    end select

    ! Process 4 - sealed area directRunoff
    select case (processMatrix(4, 1))
      ! 1 - bucket exceedance approach
    case(1)
      processMatrix(4, 2) = 1_i4
      processMatrix(4, 3) = sum(processMatrix(1 : 4, 2))
      call append(global_parameters, reshape(imperviousStorageCapacity, [1, nColPars]))

      call append(global_parameters_name, ['imperviousStorageCapacity'])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "directRunoff" does not exist!')
      stop
    end select

    ! Process 5 - potential evapotranspiration (PET)
    select case (processMatrix(5, 1))
    case(-1) ! 0 - PET is input, correct PET by LAI
      processMatrix(5, 2) = 5_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(PET_a_forest, [1, nColPars]))
      call append(global_parameters, reshape(PET_a_impervious, [1, nColPars]))
      call append(global_parameters, reshape(PET_a_pervious, [1, nColPars]))
      call append(global_parameters, reshape(PET_b, [1, nColPars]))
      call append(global_parameters, reshape(PET_c, [1, nColPars]))

      call append(global_parameters_name, [ &
                      'PET_a_forest     ', &
                      'PET_a_impervious ', &
                      'PET_a_pervious   ', &
                      'PET_b            ', &
                      'PET_c            '])
      call append(dummy_global_parameters, [&
              HargreavesSamaniCoeff(3), &
                      minCorrectionFactorPET(3), &
                      maxCorrectionFactorPET(3), &
                      aspectThresholdPET(3), &
                      PriestleyTaylorCoeff(3), &
                      PriestleyTaylorLAIcorr(3), &
                      canopyheight_forest(3), &
                      canopyheight_impervious(3), &
                      canopyheight_pervious(3), &
                      displacementheight_coeff(3), &
                      roughnesslength_momentum_coeff(3), &
                      roughnesslength_heat_coeff(3), &
                      stomatal_resistance(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'HargreavesSamaniCoeff         ', &
                      'minCorrectionFactorPET        ', &
                      'maxCorrectionFactorPET        ', &
                      'aspectThresholdPET            ', &
                      'PriestleyTaylorCoeff          ', &
                      'PriestleyTaylorLAIcorr        ', &
                      'canopyheight_forest           ', &
                      'canopyheight_impervious       ', &
                      'canopyheight_pervious         ', &
                      'displacementheight_coeff      ', &
                      'roughnesslength_momentum_coeff', &
                      'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           ' &
              ])

    case(0) ! 0 - PET is input, correct PET by aspect
      processMatrix(5, 2) = 3_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(minCorrectionFactorPET, [1, nColPars]))
      call append(global_parameters, reshape(maxCorrectionFactorPET, [1, nColPars]))
      call append(global_parameters, reshape(aspectThresholdPET, [1, nColPars]))

      call append(global_parameters_name, [ &
                  'minCorrectionFactorPET ', &
                  'maxCorrectionFactorPET ', &
                      'aspectThresholdPET     '])
      call append(dummy_global_parameters, [&
              PET_a_forest(3), &
                      PET_a_impervious(3), &
                      PET_a_pervious(3), &
                      PET_b(3), &
                      PET_c(3), &
                      HargreavesSamaniCoeff(3), &
                      PriestleyTaylorCoeff(3), &
                      PriestleyTaylorLAIcorr(3), &
                      canopyheight_forest(3), &
                      canopyheight_impervious(3), &
                      canopyheight_pervious(3), &
                      displacementheight_coeff(3), &
                      roughnesslength_momentum_coeff(3), &
                      roughnesslength_heat_coeff(3), &
                      stomatal_resistance(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'PET_a_forest                  ', &
                      'PET_a_impervious              ', &
                      'PET_a_pervious                ', &
                      'PET_b                         ', &
                      'PET_c                         ', &
                      'HargreavesSamaniCoeff         ', &
                      'PriestleyTaylorCoeff          ', &
                      'PriestleyTaylorLAIcorr        ', &
                      'canopyheight_forest           ', &
                      'canopyheight_impervious       ', &
                      'canopyheight_pervious         ', &
                      'displacementheight_coeff      ', &
                      'roughnesslength_momentum_coeff', &
                      'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           ' &
              ])

    case(1) ! 1 - Hargreaves-Samani method (HarSam) - additional input needed: Tmin, Tmax
      processMatrix(5, 2) = 4_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(minCorrectionFactorPET_HS, [1, nColPars]))
      call append(global_parameters, reshape(maxCorrectionFactorPET_HS, [1, nColPars]))
      call append(global_parameters, reshape(aspectThresholdPET_HS, [1, nColPars]))
      call append(global_parameters, reshape(HargreavesSamaniCoeff, [1, nColPars]))
      call append(global_parameters_name, [ &
                   'minCorrectionFactorPET', &
                   'maxCorrectionFactorPET', &
                      'aspectThresholdPET    ', &
                      'HargreavesSamaniCoeff '])
      call append(dummy_global_parameters, [&
              PET_a_forest(3), &
                      PET_a_impervious(3), &
                      PET_a_pervious(3), &
                      PET_b(3), &
                      PET_c(3), &
                      PriestleyTaylorCoeff(3), &
                      PriestleyTaylorLAIcorr(3), &
                      canopyheight_forest(3), &
                      canopyheight_impervious(3), &
                      canopyheight_pervious(3), &
                      displacementheight_coeff(3), &
                      roughnesslength_momentum_coeff(3), &
                      roughnesslength_heat_coeff(3), &
                      stomatal_resistance(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'PET_a_forest                  ', &
                      'PET_a_impervious              ', &
                      'PET_a_pervious                ', &
                      'PET_b                         ', &
                      'PET_c                         ', &
                      'PriestleyTaylorCoeff          ', &
                      'PriestleyTaylorLAIcorr        ', &
                      'canopyheight_forest           ', &
                      'canopyheight_impervious       ', &
                      'canopyheight_pervious         ', &
                      'displacementheight_coeff      ', &
                      'roughnesslength_momentum_coeff', &
                      'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           ' &
              ])

    case(2) ! 2 - Priestley-Taylor method (PrieTay) - additional input needed: net_rad
      processMatrix(5, 2) = 2_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(PriestleyTaylorCoeff, [1, nColPars]))
      call append(global_parameters, reshape(PriestleyTaylorLAIcorr, [1, nColPars]))
      call append(global_parameters_name, [ &
                   'PriestleyTaylorCoeff  ', &
                      'PriestleyTaylorLAIcorr'])
      call append(dummy_global_parameters, [&
              PET_a_forest(3), &
                      PET_a_impervious(3), &
                      PET_a_pervious(3), &
                      PET_b(3), &
                      PET_c(3), &
                      HargreavesSamaniCoeff(3), &
                      minCorrectionFactorPET(3), &
                      maxCorrectionFactorPET(3), &
                      aspectThresholdPET(3), &
                      canopyheight_forest(3), &
                      canopyheight_impervious(3), &
                      canopyheight_pervious(3), &
                      displacementheight_coeff(3), &
                      roughnesslength_momentum_coeff(3), &
                      roughnesslength_heat_coeff(3), &
                      stomatal_resistance(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'PET_a_forest                  ', &
                      'PET_a_impervious              ', &
                      'PET_a_pervious                ', &
                      'PET_b                         ', &
                      'PET_c                         ', &
                      'HargreavesSamaniCoeff         ', &
                      'minCorrectionFactorPET        ', &
                      'maxCorrectionFactorPET        ', &
                      'aspectThresholdPET            ', &
                      'canopyheight_forest           ', &
                      'canopyheight_impervious       ', &
                      'canopyheight_pervious         ', &
                      'displacementheight_coeff      ', &
                      'roughnesslength_momentum_coeff', &
                      'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           ' &
              ])

    case(3) ! 3 - Penman-Monteith method - additional input needed: net_rad, abs. vapour pressue, windspeed
      processMatrix(5, 2) = 7_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))

      call append(global_parameters, reshape(canopyheight_forest, [1, nColPars]))
      call append(global_parameters, reshape(canopyheight_impervious, [1, nColPars]))
      call append(global_parameters, reshape(canopyheight_pervious, [1, nColPars]))
      call append(global_parameters, reshape(displacementheight_coeff, [1, nColPars]))
      call append(global_parameters, reshape(roughnesslength_momentum_coeff, [1, nColPars]))
      call append(global_parameters, reshape(roughnesslength_heat_coeff, [1, nColPars]))
      call append(global_parameters, reshape(stomatal_resistance, [1, nColPars]))

      call append(global_parameters_name, [ &
              'canopyheight_forest           ', &
                      'canopyheight_impervious       ', &
                      'canopyheight_pervious         ', &
           'displacementheight_coeff      ', &
           'roughnesslength_momentum_coeff', &
           'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           '])
      call append(dummy_global_parameters, [&
              PET_a_forest(3), &
                      PET_a_impervious(3), &
                      PET_a_pervious(3), &
                      PET_b(3), &
                      PET_c(3), &
                      HargreavesSamaniCoeff(3), &
                      minCorrectionFactorPET(3), &
                      maxCorrectionFactorPET(3), &
                      aspectThresholdPET(3), &
                      PriestleyTaylorCoeff(3), &
                      PriestleyTaylorLAIcorr(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'PET_a_forest          ', &
                      'PET_a_impervious      ', &
                      'PET_a_pervious        ', &
                      'PET_b                 ', &
                      'PET_c                 ', &
                      'HargreavesSamaniCoeff ', &
                      'minCorrectionFactorPET', &
                      'maxCorrectionFactorPET', &
                      'aspectThresholdPET    ', &
                      'PriestleyTaylorCoeff  ', &
                      'PriestleyTaylorLAIcorr' &
              ])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "actualET" does not exist!')
      stop
    end select


    ! Process 6 - interflow
    select case (processMatrix(6, 1))
      ! 1 - parallel soil reservoir approach
    case(1)
      processMatrix(6, 2) = 5_i4
      processMatrix(6, 3) = sum(processMatrix(1 : 6, 2))
      call append(global_parameters, reshape(interflowStorageCapacityFactor, [1, nColPars]))
      call append(global_parameters, reshape(interflowRecession_slope, [1, nColPars]))
      call append(global_parameters, reshape(fastInterflowRecession_forest, [1, nColPars]))
      call append(global_parameters, reshape(slowInterflowRecession_Ks, [1, nColPars]))
      call append(global_parameters, reshape(exponentSlowInterflow, [1, nColPars]))

      call append(global_parameters_name, [ &
           'interflowStorageCapacityFactor', &
           'interflowRecession_slope      ', &
           'fastInterflowRecession_forest ', &
           'slowInterflowRecession_Ks     ', &
                      'exponentSlowInterflow         '])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "interflow" does not exist!')
      stop
    end select

    ! Process 7 - percolation
    select case (processMatrix(7, 1))
      ! 1 - GW layer is assumed as bucket
    case(1)
      processMatrix(7, 2) = 3_i4
      processMatrix(7, 3) = sum(processMatrix(1 : 7, 2))
      call append(global_parameters, reshape(rechargeCoefficient, [1, nColPars]))
      call append(global_parameters, reshape(rechargeFactor_karstic, [1, nColPars]))
      call append(global_parameters, reshape(gain_loss_GWreservoir_karstic, [1, nColPars]))

      call append(global_parameters_name, [ &
              'rechargeCoefficient          ', &
              'rechargeFactor_karstic       ', &
                      'gain_loss_GWreservoir_karstic'])

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "percolation" does not exist!')
      stop
    end select

    ! Process 8 - routing
    select case (processMatrix(8, 1))
    case(0)
      ! 0 - deactivated
      call message()
      call message('***CAUTION: Routing is deativated! ')

      processMatrix(8, 2) = 0_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
    case(1)
      ! parameter values and names are set in mRM
      ! 1 - Muskingum approach
      processMatrix(8, 2) = 5_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      ! this is overwritten in read_mrm_routing_params
      call append(global_parameters, dummy_2d_dp)
      call append(global_parameters_name, ['dummy', 'dummy', 'dummy', 'dummy', 'dummy'])
    case(2)
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      ! this is overwritten in read_mrm_routing_params
      call append(global_parameters, dummy_2d_dp_2)
      call append(global_parameters_name, ['dummy'])
    case(3)
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      call append(global_parameters, dummy_2d_dp_2)
      call append(global_parameters_name, ['dummy'])
    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "routing" does not exist!')
      stop
    end select

    !===============================================================
    ! Geological formations
    !===============================================================
    dummy = dummy // ''   ! only to avoid warning

    ! Process 9 - geoparameter
    select case (processMatrix(9, 1))
    case(1)
      ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
      ! search number of geological parameters
      do ii = 1, size(GeoParam, 1) ! no while loop to avoid risk of endless loop
        if (EQ(GeoParam(ii, 1), nodata_dp)) then
          nGeoUnits = ii - 1
          exit
        end if
      end do

      ! for geology parameters
      processMatrix(9, 2) = nGeoUnits
      processMatrix(9, 3) = sum(processMatrix(1 : 9, 2))

      call append(global_parameters, GeoParam(1 : nGeoUnits, :))

      ! create names
      do ii = 1, nGeoUnits
        dummy = 'base_flow_geo_unit_'
        if (ii < 10) then
          dummy = trim(dummy)//'0'
        end if
        dummy = trim(dummy)//trim(adjustl(num2str(ii)))
        call append(global_parameters_name, [ trim(dummy) ])
      end do

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "geoparameter" does not exist!')
      stop
    end select

    ! Process 10 - neutrons
    !   0 - deactivated
    !   1 - inverse N0 based on Desilets et al. 2010
    !   2 - COSMIC forward operator by Shuttlworth et al. 2013
    if (processMatrix(10, 1) > 0) then

      processMatrix(10, 2) = 8_i4
      processMatrix(10, 3) = sum(processMatrix(1 : 10, 2))
      call append(global_parameters, reshape(Desilets_N0, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_N0, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_N1, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_N2, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_alpha0, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_alpha1, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_L30, [1, nColPars]))
      call append(global_parameters, reshape(COSMIC_L31, [1, nColPars]))

      call append(global_parameters_name, [  &
              'Desilets_N0   ', &
                      'COSMIC_N0     ', &
                      'COSMIC_N1     ', &
                      'COSMIC_N2     ', &
                      'COSMIC_alpha0 ', &
                      'COSMIC_alpha1 ', &
                      'COSMIC_L30    ', &
                      'COSMIC_L31    '])

    else
      call message(' INFO: Process (10, neutrons) is deactivated, so output will be suppressed.')
      ! this is done below, where nml_output is read
      processMatrix(10, 2) = 0_i4
      processMatrix(10, 3) = sum(processMatrix(1 : 10, 2))
      call append(dummy_global_parameters, [&
              Desilets_N0(3), &
                      COSMIC_N0(3), &
                      COSMIC_N1(3), &
                      COSMIC_N2(3), &
                      COSMIC_alpha0(3), &
                      COSMIC_alpha1(3), &
                      COSMIC_L30(3), &
                      COSMIC_L31(3) &
              ])
      call append(dummy_global_parameters_name, [&
              'Desilets_N0   ', &
                      'COSMIC_N0     ', &
                      'COSMIC_N1     ', &
                      'COSMIC_N2     ', &
                      'COSMIC_alpha0 ', &
                      'COSMIC_alpha1 ', &
                      'COSMIC_L30    ', &
                      'COSMIC_L31    ' &
              ])

    end if

    call close_nml(unamelist_param)
    ! check if parameter are in range
    if (.not. in_bound(global_parameters)) then
      call message('***ERROR: parameters in namelist "mhm_parameters" out of bound in ', &
              trim(adjustl(file_namelist_param)))
      stop 1
    end if

  end subroutine read_mhm_parameters

END MODULE mo_common_read_config
