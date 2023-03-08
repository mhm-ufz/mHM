!> \file mo_common_read_config.f90
!> \brief   \copybrief mo_common_read_config
!> \details \copydetails mo_common_read_config

!> \brief Reading of main model configurations.
!> \details This routine reads the configurations of namelists commonly used by mHM, mRM and MPR
!> \authors Matthias Zink
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_common_read_config

  USE mo_kind, ONLY : i4, dp
  use mo_message, only: error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: common_read_config, set_land_cover_scenes_id

  ! ------------------------------------------------------------------

CONTAINS


  !> \brief Read main configurations commonly used by mHM, mRM and MPR
  !> \details Read the main configurations commonly used by mHM, mRM and MPR, namely:
  !! project_description, directories_general, mainconfig, processSelection, LCover
  !> \changelog
  !! - Robert Schweppe Dec  2018
  !!   - refactoring and restructuring
  !! - Sebastian Müller Mar 2023
  !!   - added check_L0Domain
  !> \authors Matthias Zink
  !> \date Dec 2012
  subroutine common_read_config(file_namelist, unamelist)

    use mo_common_constants, only : maxNLcovers, maxNoDomains
    use mo_common_variables, only : Conventions, LC_year_end, LC_year_start, LCfilename, contact, &
                                    dirCommonFiles, dirConfigOut, dirLCover, dirMorpho, dirOut, &
                                    mhmFileRestartOut, mrmFileRestartOut, &
                                    fileLatLon, history, iFlag_cordinate_sys, mHM_details, domainMeta, nLcoverScene, &
                                    nProcesses, nuniqueL0Domains, processMatrix, project_details, resolutionHydrology, &
                                    setup_description, simulation_type, write_restart
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    !> name of file
    character(*), intent(in) :: file_namelist

    !> id of file
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


    ! define namelists
    ! namelist directories
    namelist /project_description/ project_details, setup_description, simulation_type, &
            Conventions, contact, mHM_details, history
    namelist /directories_general/ dirConfigOut, dirCommonFiles, &
            dir_Morpho, dir_LCover, &
            dir_Out, mhm_file_RestartOut, mrm_file_RestartOut, &
            file_LatLon
    ! namelist spatial & temporal resolution, optimization information
    namelist /mainconfig/ iFlag_cordinate_sys, resolution_Hydrology, nDomains, L0Domain, write_restart, &
            read_opt_domain_data
    ! namelist process selection
    namelist /processSelection/ processCase

    ! namelist for land cover scenes
    namelist/LCover/nLcoverScene, LCoverYearStart, LCoverYearEnd, LCoverfName

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
      call error_message('***ERROR: Number of domains is resticted to ', trim(num2str(maxNoDomains)), '!')
    end if

    call check_L0Domain(L0Domain, nDomains)

    ! allocate patharray sizes
    allocate(resolutionHydrology(domainMeta%nDomains))
    allocate(dirMorpho(domainMeta%nDomains))
    allocate(mhmFileRestartOut(domainMeta%nDomains))
    allocate(mrmFileRestartOut(domainMeta%nDomains))
    allocate(dirLCover(domainMeta%nDomains))
    allocate(dirOut(domainMeta%nDomains))
    allocate(fileLatLon(domainMeta%nDomains))
    allocate(domainMeta%L0DataFrom(domainMeta%nDomains))
    allocate(domainMeta%optidata(domainMeta%nDomains))
    allocate(domainMeta%doRouting(domainMeta%nDomains))

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
    if(.NOT. (iFlag_cordinate_sys == 0 .OR. iFlag_cordinate_sys == 1)) then
      call error_message('***ERROR: coordinate system for the model run should be 0 or 1')
    end if

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
      dirMorpho(iDomain)           = dir_Morpho(domainID)
      mhmFileRestartOut(iDomain)   = mhm_file_RestartOut(domainID)
      mrmFileRestartOut(iDomain)   = mrm_file_RestartOut(domainID)
      dirLCover(iDomain)           = dir_LCover(domainID)
      dirOut(iDomain)              = dir_Out(domainID)
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

    call close_nml(unamelist)

  end subroutine common_read_config


  !> \brief Set land cover scenes IDs
  !> \changelog
  !! - Robert Schweppe Dec  2018
  !!   - refactoring and restructuring
  !> \authors Matthias Zink
  !> \date Dec 2012
  subroutine set_land_cover_scenes_id(sim_Per, LCyear_Id)

    use mo_common_constants, only : nodata_i4
    use mo_common_types, only: period
    use mo_common_variables, only : LC_year_end, LC_year_start, domainMeta, nLcoverScene
    use mo_string_utils, only : num2str

    implicit none

    type(period), dimension(:), intent(in) :: sim_Per !< simulation period
    integer(i4), dimension(:, :), allocatable, intent(inout) :: LCyear_Id !< land cover year ID

    integer(i4) :: ii, iDomain


    ! countercheck if land cover covers simulation period
    if (LC_year_start(1) .GT. minval(sim_Per(1 : domainMeta%nDomains)%yStart)) then
      call error_message('***ERROR: Land cover for warming period is missing!', raise=.false.)
      call error_message('   SimStart   : ', trim(num2str(minval(sim_Per(1 : domainMeta%nDomains)%yStart))), raise=.false.)
      call error_message('   LCoverStart: ', trim(num2str(LC_year_start(1))))
    end if
    if (LC_year_end(nLCoverScene) .LT. maxval(sim_Per(1 : domainMeta%nDomains)%yEnd)) then
      call error_message('***ERROR: Land cover period shorter than modelling period!', raise=.false.)
      call error_message('   SimEnd   : ', trim(num2str(maxval(sim_Per(1 : domainMeta%nDomains)%yEnd))), raise=.false.)
      call error_message('   LCoverEnd: ', trim(num2str(LC_year_end(nLCoverScene))))
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


  !> \brief Initialization of the domain variables
  !> \details Initialization of the domain variable for all domain loops and if activated for parallelization
  !! In case of MPI parallelization domainMeta%overAllNumberOfDomains is a
  !! variable where the number of domains from the namelist is stored. By this
  !! every process knows the total number of domains. Then, in a loop the
  !! domains are distributed onto the processes. There is a master process
  !! and several subprocesses. The master process only reads the confings in the
  !! mHM driver.
  !!
  !! The subprocesses get a number of domains. domainMeta%nDomain refers
  !! to the number of domains assigned to a specific process. It is a local
  !! variable and therefore has a different value for each process.
  !!
  !! In case more domains are there than processes, currently the domains
  !! are distributed round robin, i.e. like cards in a card game.
  !!
  !! In case less domains than processes exist, all remaining processes
  !! are assigned to the routing domains round robin. In that case the
  !! local communicator is of interest: It is a group of processes assigned
  !! to a routing domain again with a master process
  !! (domainMeta%isMasterInComLocal) and subprocesses. This communicator can
  !! in future be passed to the routing parallelization.
  !> \author Maren Kaluza
  !> \date Sep 2019
  subroutine init_domain_variable(nDomains, optiData, domainMeta)
    use mo_common_types, only: domain_meta
#ifdef MPI
    use mo_common_variables, only: comm
    use mpi_f08
#endif
    integer(i4),       intent(in)    :: nDomains !< number of domains
    integer(i4), dimension(:), intent(in) :: optiData !< optimization data
    type(domain_meta), intent(inout) :: domainMeta !< domain meta info

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
      call error_message('at least 2 processes are required')
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
    use mo_common_types, only: domain_meta
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
    use mo_common_types, only: domain_meta
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
    use mo_common_types, only: domain_meta
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

  ! \brief check the L0Domain variable from the namelist
  !> \authors Sebastian Müller
  !> \date Mar 2023
  subroutine check_L0Domain(L0Domain, nDomains)
    use mo_common_constants, only : maxNoDomains
    use mo_string_utils, only : num2str

    integer(i4), dimension(maxNoDomains), intent(in) ::L0Domain !< given L0Domain variable
    integer(i4), intent(in) :: nDomains !< number of domains

    integer(i4) :: i

    do i = 1, nDomains
      if (L0Domain(i) < 0) call error_message( &
        "L0Domain values need to be positive: ", &
        "L0Domain(", trim(adjustl(num2str(i))), ") = ", trim(adjustl(num2str(L0Domain(i)))))
      if (L0Domain(i) > i) call error_message( &
        "L0Domain values need to be less or equal to the domain index: ", &
        "L0Domain(", trim(adjustl(num2str(i))), ") = ", trim(adjustl(num2str(L0Domain(i)))))
      ! check for increasing values
      if (i > 1) then
        if (L0Domain(i) < L0Domain(i-1)) call error_message( &
          "L0Domain values need to be increasing: ", &
          "L0Domain(", trim(adjustl(num2str(i-1))), ") = ", trim(adjustl(num2str(L0Domain(i-1)))), &
          ", L0Domain(", trim(adjustl(num2str(i))), ") = ", trim(adjustl(num2str(L0Domain(i)))))
      end if
      ! if lower, check that the reference domain uses its own L0Data
      if (L0Domain(i) < i) then
        if (L0Domain(L0Domain(i)) /= L0Domain(i)) call error_message( &
          "L0Domain values should be taken from a domain with its own L0 data: ", &
          "L0Domain(", trim(adjustl(num2str(i))), ") = ", trim(adjustl(num2str(L0Domain(i)))), &
          ", L0Domain(", trim(adjustl(num2str(L0Domain(i)))), ") = ", trim(adjustl(num2str(L0Domain(L0Domain(i))))))
      end if
    end do

  end subroutine check_L0Domain

END MODULE mo_common_read_config
