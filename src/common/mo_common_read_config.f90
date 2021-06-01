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
  !TODO: MPR read_mhm_parameters needs to go here?!
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

  subroutine common_read_config(file_namelist, unamelist)

    use mo_common_constants, only : maxNLcovers, maxNoDomains
    use mo_common_variables, only : Conventions, LC_year_end, LC_year_start, LCfilename, contact, &
                                    dirCommonFiles, dirConfigOut, dirLCover, dirMorpho, dirOut, &
                                    mhmFileRestartOut, mrmFileRestartOut, &
                                    fileLatLon, history, mHM_details, domainMeta, nLcoverScene, &
                                    nProcesses, nuniqueL0Domains, processMatrix, project_details, resolutionHydrology, &
                                    setup_description, simulation_type, write_restart, &
                                    dds_r, mhmFileRestartIn, mrmFileRestartIn, evalPer,&
                                    mcmc_error_params, mcmc_opti, nIterations, &
                                    opti_function, opti_method, optimize, optimize_restart, &
                                    read_restart, mrm_read_river_network, resolutionRouting, sa_temp, &
                                    sce_ngs, sce_npg, sce_nps, seed, &
                                    warmPer, warmingDays
    use mo_common_variables, only : LCfilename, domainMeta, processMatrix
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
    namelist/LCover/nLcoverScene, LCoverYearStart, LCoverYearEnd, LCoverfName
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
    allocate(domainMeta%L0DataFrom(domainMeta%nDomains))
    allocate(domainMeta%optidata(domainMeta%nDomains))
    allocate(domainMeta%doRouting(domainMeta%nDomains))

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
    allocate(LCfilename(nLCoverScene))
    allocate(LC_year_start(nLCoverScene))
    allocate(LC_year_end(nLCoverScene))
    LCfilename(:) = LCoverfName(1 : nLCoverScene)
    LC_year_start(:) = LCoverYearStart(1 : nLCoverScene)
    LC_year_end(:) = LCoverYearEnd(1 : nLCoverScene)

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

    call read_mhm_parameters(file_namelist_param, unamelist_param)

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
    use mo_common_variables, only : LC_year_end, LC_year_start, domainMeta, nLcoverScene
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
    if (LC_year_end(nLCoverScene) .LT. maxval(sim_Per(1 : domainMeta%nDomains)%yEnd)) then
      call message()
      call message('***ERROR: Land cover period shorter than modelling period!')
      call message('   SimEnd   : ', trim(num2str(maxval(sim_Per(1 : domainMeta%nDomains)%yEnd))))
      call message('   LCoverEnd: ', trim(num2str(LC_year_end(nLCoverScene))))
      stop 1
    end if
    !
    allocate(LCyear_Id(minval(sim_Per(1 : domainMeta%nDomains)%yStart) : maxval(sim_Per(1 : domainMeta%nDomains)%yEnd), &
                   domainMeta%nDomains))
    LCyear_Id = nodata_i4
    do iDomain = 1, domainMeta%nDomains
      do ii = 1, nLCoverScene
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
        allocate(domainMeta%indices(domainMeta%nDomains))
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
    allocate(domainMeta%indices(domainMeta%nDomains))
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
    allocate(domainMeta%indices(domainMeta%nDomains))
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
    allocate(domainMeta%indices(domainMeta%nDomains))
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
      allocate(domainMeta%indices(domainMeta%nDomains))
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
      allocate(domainMeta%indices(domainMeta%nDomains))
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


END MODULE mo_common_read_config
