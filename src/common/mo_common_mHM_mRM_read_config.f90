!> \file mo_common_mHM_mRM_read_config.f90
!> \brief \copybrief mo_common_mhm_mrm_read_config
!> \details \copydetails mo_common_mhm_mrm_read_config

!> \brief Reading of main model configurations.
!> \details This routine reads the configurations of common program parts
!> \authors Matthias Zink
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_common_mHM_mRM_read_config

  use mo_kind, only : i4, dp
  use mo_message, only : message, error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: common_mHM_mRM_read_config, common_check_resolution, check_optimization_settings

CONTAINS


  !> \brief Read main configurations for common parts
  !> \changelog
  !! - Robert Schweppe Dec  2017
  !!   - based on mhm_read_config
  !! - Stephan Thober Jan 2022
  !!   - added nTStepForcingDay
  !> \authors Matthias Zink
  !> \date Dec 2012
  subroutine common_mHM_mRM_read_config(file_namelist, unamelist)

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_mHM_mRM_variables, only : LCyearId, dds_r, mhmFileRestartIn, mrmFileRestartIn, evalPer,&
                                            mcmc_error_params, mcmc_opti, nIterations, &
                                            nTStepDay, opti_function, opti_method, optimize, optimize_restart, &
                                            read_restart, mrm_read_river_network, resolutionRouting, sa_temp, &
                                            sce_ngs, sce_npg, sce_nps, seed, &
                                            simPer, timestep, warmPer, warmingDays, read_old_style_restart_bounds, &
                                            restart_reset_fluxes_states
    use mo_common_read_config, only : set_land_cover_scenes_id
    use mo_common_types, only: period
    use mo_common_variables, only : LCfilename, domainMeta, processMatrix
    use mo_julian, only : caldat, julday
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    character(*), intent(in) :: file_namelist !< namelist file name
    integer, intent(in) :: unamelist !< unit to open namelist file

    integer(i4) :: jday

    integer(i4) :: domainID, iDomain

    integer(i4), dimension(maxNoDomains) :: warming_Days

    type(period), dimension(maxNoDomains) :: eval_Per

    real(dp), dimension(maxNoDomains) :: resolution_Routing

    character(256), dimension(maxNoDomains) :: mhm_file_RestartIn
    character(256), dimension(maxNoDomains) :: mrm_file_RestartIn


    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig_mhm_mrm/ timestep, resolution_Routing, optimize, &
            optimize_restart, opti_method, opti_function, &
            read_restart, mrm_read_river_network, read_old_style_restart_bounds, restart_reset_fluxes_states, &
            mhm_file_RestartIn, mrm_file_RestartIn
    ! namelist for optimization settings
    namelist /Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, &
            sce_npg, sce_nps, mcmc_opti, mcmc_error_params
    ! namelist for time settings
    namelist /time_periods/ warming_Days, eval_Per

    ! set default values for optional arguments
    mrm_read_river_network = .false.
    read_old_style_restart_bounds = .false.
    restart_reset_fluxes_states = .false.

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist specifying the model configuration
    !===============================================================
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
      call error_message('***ERROR: cannot read states from restart file when optimizing')
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
      call error_message('mo_startup: timeStep must be a divisor of 24: ', num2str(timeStep))
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
      ! this will be a procedure subroutine
      ! therefore inout first, in second
      call period_copy_period_data(evalPer(iDomain), eval_Per(domainID))
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

  end subroutine common_mHM_mRM_read_config


  !> \brief check optimization settings
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine check_optimization_settings

    use mo_common_mHM_mRM_variables, only : dds_r, nIterations, sce_ngs, sce_npg, sce_nps
    use mo_common_variables, only : global_parameters

    implicit none

    integer(i4) :: n_true_pars


    ! check and set default values
    if (nIterations .le. 0_i4) then
      call error_message('Number of iterations for Optimization (nIterations) must be greater than zero')
    end if
    if (dds_r .lt. 0.0_dp .or. dds_r .gt. 1.0_dp) then
      call error_message('dds_r must be between 0.0 and 1.0')
    end if
    if (sce_ngs .lt. 1_i4) then
      call error_message('number of complexes in SCE (sce_ngs) must be at least 1')
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
      call error_message('number of points per complex (sce_npg) must be greater or', raise=.false.)
      call error_message('equal number of points per sub-complex (sce_nps)')
    end if

  end subroutine check_optimization_settings


  !> \brief check resolution
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine common_check_resolution(do_message, allow_subgrid_routing)

    use mo_common_mHM_mRM_variables, only : resolutionRouting
    use mo_common_variables, only : domainMeta, resolutionHydrology
    use mo_string_utils, only : num2str

    implicit none

    logical, intent(in) :: do_message !< flag to print messages
    logical, intent(in) :: allow_subgrid_routing !< flag to allow subgrid routing

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
          call error_message('***ERROR: Resolution of routing is not a multiple of hydrological model resolution!', raise=.false.)
          call error_message('   FILE: mhm.nml, namelist: mainconfig, variable: resolutionRouting')
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


  ! ToDo: make this a procedure of period
  !> \brief copy period data
  subroutine period_copy_period_data(toPeriod, fromPeriod)
    use mo_common_types, only: period
    type(period), intent(inout) :: toPeriod !< copy to this period
    type(period), intent(in)    :: fromPeriod !< copy from this period

    toPeriod%dStart   = fromPeriod%dStart    ! first day
    toPeriod%mStart   = fromPeriod%mStart    ! first month
    toPeriod%yStart   = fromPeriod%yStart    ! first year
    toPeriod%dEnd     = fromPeriod%dEnd      ! last  day
    toPeriod%mEnd     = fromPeriod%mEnd      ! last  month
    toPeriod%yEnd     = fromPeriod%yEnd      ! last  year
    toPeriod%julStart = 0 ! first julian day
    toPeriod%julEnd   = 0 ! last  julian day
    toPeriod%nObs     = 0 ! total number of observations

  end subroutine period_copy_period_data

END MODULE mo_common_mHM_mRM_read_config
