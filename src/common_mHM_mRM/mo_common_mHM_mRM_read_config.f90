!>       \file mo_common_mHM_mRM_read_config.f90

!>       \brief Reading of main model configurations.

!>       \details This routine reads the configurations of common program parts

!>       \authors Matthias Zink

!>       \date Dec 2012

! Modifications:

MODULE mo_common_mHM_mRM_read_config

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: common_mHM_mRM_read_config, common_check_resolution, check_optimization_settings

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        common_mHM_mRM_read_config

  !    PURPOSE
  !>       \brief Read main configurations for common parts

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist"
  !>       \param[in] "integer :: unamelist"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Dec  2017 - based on mhm_read_config

  subroutine common_mHM_mRM_read_config(file_namelist, unamelist)

    use mo_common_constants, only : maxNoBasins
    use mo_common_mHM_mRM_variables, only : LCyearId, dds_r, dirRestartIn, evalPer, mcmc_error_params, mcmc_opti, nIterations, &
                                            nTStepDay, opti_function, opti_method, optimize, optimize_restart, &
                                            read_restart, resolutionRouting, sa_temp, sce_ngs, sce_npg, sce_nps, seed, &
                                            simPer, timestep, warmPer, warmingDays
    use mo_common_read_config, only : set_land_cover_scenes_id
    use mo_common_variables, only : LCfilename, nBasins, period
    use mo_julian, only : caldat, julday
    use mo_message, only : message
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    character(*), intent(in) :: file_namelist

    integer, intent(in) :: unamelist

    integer(i4) :: jday

    integer(i4) :: iBasin

    integer(i4), dimension(maxNoBasins) :: warming_Days

    type(period), dimension(maxNoBasins) :: eval_Per

    real(dp), dimension(maxNoBasins) :: resolution_Routing

    character(256), dimension(maxNoBasins) :: dir_RestartIn


    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig_mhm_mrm/ timestep, resolution_Routing, optimize, &
            optimize_restart, opti_method, opti_function, &
            read_restart, dir_RestartIn
    ! namelist for optimization settings
    namelist /Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, &
            sce_npg, sce_nps, mcmc_opti, mcmc_error_params
    ! namelist for time settings
    namelist /time_periods/ warming_Days, eval_Per

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist specifying the model configuration
    !===============================================================
    call position_nml('mainconfig_mhm_mrm', unamelist)
    read(unamelist, nml = mainconfig_mhm_mrm)

    allocate(resolutionRouting(nBasins))
    allocate(dirRestartIn(nBasins))
    dirRestartIn = dir_RestartIn(1 : nBasins)

    resolutionRouting = resolution_Routing(1 : nBasins)

    ! check for optimize and read restart
    if ((read_restart) .and. (optimize)) then
      call message()
      call message('***ERROR: cannot read states from restart file when optimizing')
      stop 1
    end if

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
    allocate(simPer(nBasins))
    allocate(evalPer(nBasins))
    allocate(warmingDays(nBasins))
    allocate(warmPer(nBasins))

    !===============================================================
    !  read simulation time periods incl. warming days
    !===============================================================
    call position_nml('time_periods', unamelist)
    read(unamelist, nml = time_periods)
    warmingDays = warming_Days(1 : nBasins)
    evalPer = eval_Per(1 : nBasins)

    !===============================================================
    !  determine simulation time period incl. warming days for each
    !  basin
    !===============================================================
    do iBasin = 1, nBasins
      ! julian days for evaluation period
      jday = julday(dd = evalPer(iBasin)%dStart, mm = evalPer(iBasin)%mStart, yy = evalPer(iBasin)%yStart)
      evalPer(iBasin)%julStart = jday

      jday = julday(dd = evalPer(iBasin)%dEnd, mm = evalPer(iBasin)%mEnd, yy = evalPer(iBasin)%yEnd)
      evalPer(iBasin)%julEnd = jday

      ! determine warming period
      warmPer(iBasin)%julStart = evalPer(iBasin)%julStart - warmingDays(iBasin)
      warmPer(iBasin)%julEnd = evalPer(iBasin)%julStart - 1

      call caldat(warmPer(iBasin)%julStart, dd = warmPer(iBasin)%dStart, mm = warmPer(iBasin)%mStart, yy = warmPer(iBasin)%yStart)
      call caldat(warmPer(iBasin)%julEnd, dd = warmPer(iBasin)%dEnd, mm = warmPer(iBasin)%mEnd, yy = warmPer(iBasin)%yEnd)

      ! simulation Period = warming Period + evaluation Period
      simPer(iBasin)%dStart = warmPer(iBasin)%dStart
      simPer(iBasin)%mStart = warmPer(iBasin)%mStart
      simPer(iBasin)%yStart = warmPer(iBasin)%yStart
      simPer(iBasin)%julStart = warmPer(iBasin)%julStart
      simPer(iBasin)%dEnd = evalPer(iBasin)%dEnd
      simPer(iBasin)%mEnd = evalPer(iBasin)%mEnd
      simPer(iBasin)%yEnd = evalPer(iBasin)%yEnd
      simPer(iBasin)%julEnd = evalPer(iBasin)%julEnd
    end do

    call set_land_cover_scenes_id(simPer, LCyearId, LCfilename)

    !===============================================================
    ! Settings for Optimization
    !===============================================================
    ! namelist for Optimization settings
    call position_nml('Optimization', unamelist)
    read(unamelist, nml = Optimization)
    ! checking of settings and default value initialization moved to new subroutine
    ! because global_parameters need to be set, which is not the case right now
    call close_nml(unamelist)

  end subroutine common_mHM_mRM_read_config

  !    NAME
  !        check_optimization_settings

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine check_optimization_settings

    use mo_common_mHM_mRM_variables, only : dds_r, nIterations, sce_ngs, sce_npg, sce_nps
    use mo_common_variables, only : global_parameters
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

  !    NAME
  !        common_check_resolution

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "logical :: do_message"
  !>       \param[in] "logical :: allow_subgrid_routing"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine common_check_resolution(do_message, allow_subgrid_routing)

    use mo_common_mHM_mRM_variables, only : resolutionRouting
    use mo_common_variables, only : nBasins, resolutionHydrology
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    logical, intent(in) :: do_message

    logical, intent(in) :: allow_subgrid_routing

    integer(i4) :: ii

    ! conversion factor L11 to L1
    real(dp) :: cellFactorRbyH


    !===============================================================
    ! check matching of resolutions: hydrology, forcing and routing
    !===============================================================
    do ii = 1, nBasins
      cellFactorRbyH = resolutionRouting(ii) / resolutionHydrology(ii)
      if (do_message) then
        call message()
        call message('Basin ', trim(adjustl(num2str(ii))), ': ')
        call message('resolution Hydrology (basin ', trim(adjustl(num2str(ii))), ')     = ', &
                trim(adjustl(num2str(resolutionHydrology(ii)))))
        call message('resolution Routing (basin ', trim(adjustl(num2str(ii))), ')       = ', &
                trim(adjustl(num2str(resolutionRouting(ii)))))
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

END MODULE mo_common_mHM_mRM_read_config
