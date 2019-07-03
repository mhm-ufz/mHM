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

  PUBLIC :: common_read_config, set_land_cover_scenes_id

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        common_read_config

  !    PURPOSE
  !>       \brief Read main configurations commonly used by mHM, mRM and MPR

  !>       \details Read the main configurations commonly used by mHM, mRM and MPR, namely:
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
    use mo_common_variables, only : Conventions, L0_Domain, LC_year_end, LC_year_start, LCfilename, contact, &
                                    dirCommonFiles, dirConfigOut, dirLCover, dirMorpho, dirOut, dirRestartOut, &
                                    fileLatLon, history, iFlag_cordinate_sys, mHM_details, domainMeta, nLcoverScene, &
                                    nProcesses, nuniqueL0Domains, processMatrix, project_details, resolutionHydrology, &
                                    setup_description, simulation_type, write_restart
    use mo_message, only : message
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    ! name of file
    character(*), intent(in) :: file_namelist

    ! id of file
    integer, intent(in) :: unamelist

    ! Choosen process description number
    integer(i4), dimension(nProcesses) :: processCase

    character(256), dimension(maxNoDomains) :: dir_Morpho

    character(256), dimension(maxNoDomains) :: dir_RestartOut

    character(256), dimension(maxNoDomains) :: dir_LCover

    character(256), dimension(maxNoDomains) :: dir_Out

    character(256), dimension(maxNoDomains) :: file_LatLon

    real(dp), dimension(maxNoDomains) :: resolution_Hydrology

    integer(i4), dimension(maxNoDomains) :: L0Domain

    ! starting year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart

    ! ending year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd

    ! filename of Lcover file
    character(256), dimension(maxNLCovers) :: LCoverfName

    integer(i4) :: i, newDomainID, domainID, iDomain, nDomains


    ! define namelists
    ! namelist directories
    namelist /project_description/ project_details, setup_description, simulation_type, &
            Conventions, contact, mHM_details, history
    namelist /directories_general/ dirConfigOut, dirCommonFiles, &
            dir_Morpho, dir_LCover, &
            dir_Out, dir_RestartOut, &
            file_LatLon
    ! namelist spatial & temporal resolution, optimization information
    namelist /mainconfig/ iFlag_cordinate_sys, resolution_Hydrology, nDomains, L0Domain, write_restart
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

    call init_domain_variable(nDomains, domainMeta)

    if (nDomains .GT. maxNoDomains) then
      call message()
      call message('***ERROR: Number of domains is resticted to ', trim(num2str(maxNoDomains)), '!')
      stop 1
    end if

    ! allocate patharray sizes
    allocate(resolutionHydrology(domainMeta%nDomains))
    allocate(dirMorpho(domainMeta%nDomains))
    allocate(dirRestartOut(domainMeta%nDomains))
    allocate(dirLCover(domainMeta%nDomains))
    allocate(dirOut(domainMeta%nDomains))
    allocate(fileLatLon(domainMeta%nDomains))
    allocate(L0_Domain(domainMeta%nDomains))

    nuniqueL0Domains = 0_i4
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      resolutionHydrology(iDomain) = resolution_Hydrology(domainID)
      ! if a domain uses the same L0 data as a previous one, write
      ! the index into L0_Domain
      ! ToDo: switch L0_Domain with L0_data_from as type part of domainMeta
      newDomainID = L0Domain(domainID)
      L0_Domain(iDomain) = iDomain
      do i = 1, iDomain - 1
        if (newDomainID == domainMeta%indices(i)) then
          L0_Domain(iDomain) = i
          cycle 
        end if
      end do
      nuniqueL0Domains = nuniqueL0Domains + 1_i4
     ! L0_Domain(iDomain) = L0Domain(domainID)
    end do

    ! check for possible options
    if(.NOT. (iFlag_cordinate_sys == 0 .OR. iFlag_cordinate_sys == 1)) then
      call message()
      call message('***ERROR: coordinate system for the model run should be 0 or 1')
      stop 1
    end if

  !  nuniqueL0Domains = 0_i4
  !  do iDomain = 1, domainMeta%nDomains
  !    if (iDomain .gt. 1) then
  !      if (L0_Domain(iDomain) .eq. L0_Domain(iDomain - 1)) then
  !        cycle
  !      end if
  !    end if
  !    nuniqueL0Domains = nuniqueL0Domains + 1_i4
  !  end do


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
    !  Read namelist for mainpaths
    !===============================================================
    call position_nml('directories_general', unamelist)
    read(unamelist, nml = directories_general)

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      dirMorpho(iDomain)     = dir_Morpho(domainID)
      dirRestartOut(iDomain) = dir_RestartOut(domainID)
      dirLCover(iDomain)     = dir_LCover(domainID)
      dirOut(iDomain)        = dir_Out(domainID)
      fileLatLon(iDomain)    = file_LatLon(domainID)
    end do

    !===============================================================
    ! Read process selection list
    !===============================================================
    call position_nml('processselection', unamelist)
    read(unamelist, nml = processSelection)

    processMatrix = 0_i4
    processMatrix(:, 1) = processCase

    call close_nml(unamelist)

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

  subroutine set_land_cover_scenes_id(sim_Per, LCyear_Id, LCfilename)

    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : LC_year_end, LC_year_start, domainMeta, nLcoverScene, period
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    type(period), dimension(:), intent(in) :: sim_Per

    integer(i4), dimension(:, :), allocatable, intent(inout) :: LCyear_Id

    character(256), dimension(:), allocatable, intent(inout) :: LCfilename

    integer(i4) :: ii, iDomain, max_lcs, min_lcs, jj

    character(256), dimension(:), allocatable :: dummy_LCfilenames

    integer(i4), dimension(:,:), allocatable :: dummy_LCyears


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

    ! correct number of input land cover scenes to number of needed scenes
    max_lcs = maxval(LCyear_Id, mask = (LCyear_Id .gt. nodata_i4))
    min_lcs = minval(LCyear_Id, mask = (LCyear_Id .gt. nodata_i4))
    nLCoverScene = max_lcs - min_lcs + 1

    ! select the LC_years for only the needed scenes
    allocate(dummy_LCyears(2, nLCoverScene))
    jj = 1
    do ii = 1, size(LC_year_start)
      if ((LC_year_start(ii) .lt. LC_year_end(max_lcs)) .and. (LC_year_end(ii) .gt. LC_year_start(min_lcs))) then
        dummy_LCyears(1, jj) = LC_year_start(ii)
        dummy_LCyears(2, jj) = LC_year_end(ii)
        jj = jj + 1
      end if
    end do

    ! put land cover scenes to corresponding file name and LuT
    ! this was allocated for MPR before, now update using only needed scenes
    allocate(dummy_LCfilenames(nLCoverScene))
    dummy_LCfilenames(:) = LCfilename(minval(LCyear_Id, mask = (LCyear_Id .gt. nodata_i4)) : &
            maxval(LCyear_Id, mask = (LCyear_Id .gt. nodata_i4)))
    deallocate(LCfilename, LC_year_start, LC_year_end)
    allocate(LCfilename(nLCoverScene))
    allocate(LC_year_start(nLCoverScene))
    allocate(LC_year_end(nLCoverScene))
    LCfilename(:) = dummy_LCfilenames(:)
    LC_year_start(:) = dummy_LCyears(1, :)
    LC_year_end(:) = dummy_LCyears(2, :)

    ! update the ID's
    if (maxval(sim_Per(1 : domainMeta%nDomains)%julStart) .eq. minval(sim_Per(1 : domainMeta%nDomains)%julStart) .and. &
            maxval(sim_Per(1 : domainMeta%nDomains)%julEnd) .eq. minval(sim_Per(1 : domainMeta%nDomains)%julEnd)) then
      if (any(LCyear_Id .EQ. nodata_i4)) then
        call message()
        call message('***ERROR: Intermediate land cover period is missing!')
        stop 1
      end if
    end if

  end subroutine set_land_cover_scenes_id

  subroutine init_domain_variable(nDomains, domainMeta)
    use mo_common_variables, only: domain_meta
    use mo_common_mHM_mRM_variables, only: optimize
#ifdef MPI
    use mo_common_variables, only: comm
    use mpi_f08
#endif
    integer(i4),       intent(in)    :: nDomains
    type(domain_meta), intent(inout) :: domainMeta

    integer             :: ierror
    integer(i4)         :: nproc
    integer(i4)         :: rank
    integer(i4)         :: iDomain

    domainMeta%overAllNumberOfDomains = nDomains
#ifdef MPI
    ! find number of processes nproc
    call MPI_Comm_size(comm, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(comm, rank, ierror)
    if (optimize) then
      if (nproc < domainMeta%overAllNumberOfDomains + 1) then
        !ToDo: message
        write(*,*) "Warning: not all domains will be simulated"
      end if
      if (rank == 0) then
        domainMeta%nDomains = 1
        allocate(domainMeta%indices(domainMeta%nDomains))
        domainMeta%indices(1) = 1
      end if
      do iDomain = 1 , domainMeta%overallNumberOfDomains
        if (rank == iDomain) then
          domainMeta%nDomains = 1
          allocate(domainMeta%indices(domainMeta%nDomains))
          domainMeta%indices(1) = iDomain
        end if
      end do
      if (rank > domainMeta%overallNumberOfDomains) then
        domainMeta%nDomains = 1
        allocate(domainMeta%indices(domainMeta%nDomains))
        domainMeta%indices(1) = 1
      end if
    else
      if (nproc < domainMeta%overAllNumberOfDomains) then
        !ToDo: message
        write(*,*) "Warning: not all domains will be simulated"
      end if
      do iDomain = 1 , domainMeta%overallNumberOfDomains
        if (rank == iDomain - 1) then
          domainMeta%nDomains = 1
          allocate(domainMeta%indices(domainMeta%nDomains))
          domainMeta%indices(1) = iDomain
        end if
      end do
      if (rank > domainMeta%overallNumberOfDomains - 1) then
        domainMeta%nDomains = 1
        allocate(domainMeta%indices(domainMeta%nDomains))
        domainMeta%indices(1) = 1
      end if
    end if

#else
    domainMeta%nDomains = nDomains
    allocate(domainMeta%indices(domainMeta%nDomains))
    do iDomain = 1, domainMeta%nDomains
      domainMeta%indices(iDomain) = iDomain
    end do
#endif

  end subroutine init_domain_variable

END MODULE mo_common_read_config
