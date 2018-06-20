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

    use mo_common_constants, only : maxNLcovers, maxNoBasins
    use mo_common_variables, only : Conventions, L0_Basin, LC_year_end, LC_year_start, LCfilename, contact, &
                                    dirCommonFiles, dirConfigOut, dirLCover, dirMorpho, dirOut, dirRestartOut, &
                                    fileLatLon, history, iFlag_cordinate_sys, mHM_details, nBasins, nLcoverScene, &
                                    nProcesses, nuniquel0Basins, processMatrix, project_details, resolutionHydrology, &
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

    character(256), dimension(maxNoBasins) :: dir_Morpho

    character(256), dimension(maxNoBasins) :: dir_RestartOut

    character(256), dimension(maxNoBasins) :: dir_LCover

    character(256), dimension(maxNoBasins) :: dir_Out

    character(256), dimension(maxNoBasins) :: file_LatLon

    real(dp), dimension(maxNoBasins) :: resolution_Hydrology

    integer(i4), dimension(maxNoBasins) :: L0Basin

    ! starting year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart

    ! ending year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd

    ! filename of Lcover file
    character(256), dimension(maxNLCovers) :: LCoverfName

    integer(i4) :: iBasin


    ! define namelists
    ! namelist directories
    namelist /project_description/ project_details, setup_description, simulation_type, &
            Conventions, contact, mHM_details, history
    namelist /directories_general/ dirConfigOut, dirCommonFiles, &
            dir_Morpho, dir_LCover, &
            dir_Out, dir_RestartOut, &
            file_LatLon
    ! namelist spatial & temporal resolution, optimization information
    namelist /mainconfig/ iFlag_cordinate_sys, resolution_Hydrology, nBasins, L0Basin, write_restart
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

    if (nBasins .GT. maxNoBasins) then
      call message()
      call message('***ERROR: Number of basins is resticted to ', trim(num2str(maxNoBasins)), '!')
      stop 1
    end if

    ! allocate patharray sizes
    allocate(resolutionHydrology(nBasins))
    allocate(dirMorpho(nBasins))
    allocate(dirRestartOut(nBasins))
    allocate(dirLCover(nBasins))
    allocate(dirOut(nBasins))
    allocate(fileLatLon(nBasins))
    allocate(L0_Basin(nBasins))

    resolutionHydrology = resolution_Hydrology(1 : nBasins)
    L0_Basin = L0Basin(1 : nBasins)

    ! check for possible options
    if(.NOT. (iFlag_cordinate_sys == 0 .OR. iFlag_cordinate_sys == 1)) then
      call message()
      call message('***ERROR: coordinate system for the model run should be 0 or 1')
      stop 1
    end if

    nuniquel0Basins = 0_i4
    do iBasin = 1, nBasins
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          cycle
        end if
      end if
      nuniquel0Basins = nuniquel0Basins + 1_i4
    end do


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

    dirMorpho = dir_Morpho(1 : nBasins)
    dirRestartOut = dir_RestartOut(1 : nBasins)
    dirLCover = dir_LCover(1 : nBasins)
    dirOut = dir_Out(1 : nBasins)
    fileLatLon = file_LatLon(1 : nBasins)

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
    use mo_common_variables, only : LC_year_end, LC_year_start, nBasins, nLcoverScene, period
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    type(period), dimension(:), intent(in) :: sim_Per

    integer(i4), dimension(:, :), allocatable, intent(inout) :: LCyear_Id

    character(256), dimension(:), allocatable, intent(inout) :: LCfilename

    integer(i4) :: ii, iBasin, max_lcs, min_lcs, jj

    character(256), dimension(:), allocatable :: dummy_LCfilenames

    integer(i4), dimension(:,:), allocatable :: dummy_LCyears


    ! countercheck if land cover covers simulation period
    if (LC_year_start(1) .GT. minval(sim_Per(1 : nBasins)%yStart)) then
      call message()
      call message('***ERROR: Land cover for warming period is missing!')
      call message('   SimStart   : ', trim(num2str(minval(sim_Per(1 : nBasins)%yStart))))
      call message('   LCoverStart: ', trim(num2str(LC_year_start(1))))
      stop 1
    end if
    if (LC_year_end(nLCoverScene) .LT. maxval(sim_Per(1 : nBasins)%yEnd)) then
      call message()
      call message('***ERROR: Land cover period shorter than modelling period!')
      call message('   SimEnd   : ', trim(num2str(maxval(sim_Per(1 : nBasins)%yEnd))))
      call message('   LCoverEnd: ', trim(num2str(LC_year_end(nLCoverScene))))
      stop 1
    end if
    !
    allocate(LCyear_Id(minval(sim_Per(1 : nBasins)%yStart) : maxval(sim_Per(1 : nBasins)%yEnd), nBasins))
    LCyear_Id = nodata_i4
    do iBasin = 1, nBasins
      do ii = 1, nLCoverScene
        ! land cover before model period or land cover after model period
        if ((LC_year_end(ii) .LT. sim_Per(iBasin)%yStart) .OR. &
                (LC_year_start(ii) .GT. sim_Per(iBasin)%yEnd)) then
          cycle
          ! land cover period fully covers model period
        else if ((LC_year_start(ii) .LE. sim_Per(iBasin)%yStart) .AND. &
                (LC_year_end(ii) .GE. sim_Per(iBasin)%yEnd)) then
          LCyear_Id(sim_Per(iBasin)%yStart : sim_Per(iBasin)%yEnd, iBasin) = ii
          exit
          ! land cover period covers beginning of model period
        else if ((LC_year_start(ii) .LE. sim_Per(iBasin)%yStart) .AND. &
                (LC_year_end(ii) .LT. sim_Per(iBasin)%yEnd)) then
          LCyear_Id(sim_Per(iBasin)%yStart : LC_year_end(ii), iBasin) = ii
          ! land cover period covers end of model period
        else if ((LC_year_start(ii) .GT. sim_Per(iBasin)%yStart) .AND. &
                (LC_year_end(ii) .GE. sim_Per(iBasin)%yEnd)) then
          LCyear_Id(LC_year_start(ii) : sim_Per(iBasin)%yEnd, iBasin) = ii
          ! land cover period covers part of model_period
        else
          LCyear_Id(LC_year_start(ii) : LC_year_end(ii), iBasin) = ii
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
    if (maxval(sim_Per(1 : nBasins)%julStart) .eq. minval(sim_Per(1 : nBasins)%julStart) .and. &
            maxval(sim_Per(1 : nBasins)%julEnd) .eq. minval(sim_Per(1 : nBasins)%julEnd)) then
      if (any(LCyear_Id .EQ. nodata_i4)) then
        call message()
        call message('***ERROR: Intermediate land cover period is missing!')
        stop 1
      end if
    end if

  end subroutine set_land_cover_scenes_id


END MODULE mo_common_read_config
