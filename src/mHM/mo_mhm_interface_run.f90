!> \file    mo_mhm_interface_run.f90
!> \copydoc mo_mhm_interface_run

!> \brief   Module providing interfaces for running preconfigured mHM.
!> \version 0.1
!> \authors Sebastian Mueller, Matthias Kelbling
!> \date    Jan 2022
!> \details Interfaces to control the mHM run from outside (prepare domain, do timestep, ...).
module mo_mhm_interface_run

  ! forces
  use mo_kind, only: i4, dp
  use mo_message, only: message, error_message
  use mo_string_utils, only : num2str
  ! mhm
  use mo_common_run_variables, only : run_cfg
  use mo_optimization_types, only : optidata_sim
  use mo_common_datetime_type, only : datetimeinfo
  use mo_common_mHM_mRM_variables, only : &
    resolutionRouting, &
    LCyearId, &
    mhmFileRestartIn, &
    mrmFileRestartIn, &
    nTstepDay, &
    nTstepForcingDay, &
    optimize, &
    readPer, &
    read_restart, &
    simPer, &
    timeStep, &
    warmingDays, &
    c2TSTu
  use mo_common_variables, only : &
    global_parameters, &
    level1, &
    domainMeta, &
    processMatrix
  use mo_global_variables, only : &
    L1_Throughfall, &
    L1_aETCanopy, &
    L1_aETSealed, &
    L1_aETSoil, &
    L1_absvappress, &
    L1_baseflow, &
    L1_fastRunoff, &
    L1_infilSoil, &
    L1_inter, &
    L1_melt, &
    L1_netrad, &
    L1_neutrons, &
    L1_percol, &
    L1_pet, &
    L1_pet_calc, &
    L1_temp_calc, &
    L1_prec_calc, &
    L1_pet_weights, &
    L1_pre, &
    L1_preEffect, &
    L1_pre_weights, &
    L1_rain, &
    L1_runoffSeal, &
    L1_satSTW, &
    L1_sealSTW, &
    L1_slowRunoff, &
    L1_snow, &
    L1_snowPack, &
    L1_soilMoist, &
    L1_temp, &
    L1_temp_weights, &
    L1_tmax, &
    L1_tmin, &
    L1_total_runoff, &
    L1_unsatSTW, &
    L1_windspeed, &
    L1_twsaObs, &
    L1_etObs, &
    L1_smObs, &
    L1_neutronsObs, &
    L1_tann, &
    L1_ssrd, &
    L1_strd, &
    evap_coeff, &
    fday_pet, &
    fday_prec, &
    fday_temp, &
    fnight_pet, &
    fnight_prec, &
    fnight_temp, &
    nSoilHorizons_sm_input, &
    neutron_integral_AFast, &
    outputFlxState, &
    read_meteo_weights, &
    timeStep_model_inputs, &
    timeStep_model_outputs, &
    fday_ssrd, &
    fnight_ssrd, &
    fday_strd, &
    fnight_strd, &
    BFI_qBF_sum, &
    BFI_qT_sum
  use mo_init_states, only : variables_default_init
  use mo_julian, only : caldat, julday
  use mo_message, only : error_message
  use mo_string_utils, only : num2str
  use mo_meteo_forcings, only : prepare_meteo_forcings_data
  use mo_mhm, only : mhm
  use mo_restart, only : read_restart_states
  use mo_write_fluxes_states, only : OutputDataset
  use mo_constants, only : HourSecs
  use mo_common_variables, only : resolutionHydrology
  use mo_mrm_global_variables, only : &
    InflowGauge, &
    L11_C1, &
    L11_C2, &
    L11_L1_Id, &
    L11_TSrout, &
    L11_fromN, &
    L11_length, &
    L11_nLinkFracFPimp, &
    L11_nOutlets, &
    L11_netPerm, &
    L11_qMod, &
    L11_qOUT, &
    L11_qTIN, &
    L11_qTR, &
    L11_slope, &
    L11_toN, &
    L1_L11_Id, &
    domain_mrm, &
    level11, &
    mRM_runoff, &
    outputFlxState_mrm, &
    timeStep_model_outputs_mrm, &
    gw_coupling, &
    L0_river_head_mon_sum, &
    riv_temp_pcs
  use mo_mpr_global_variables, only : &
    L1_HarSamCoeff, &
    L1_PrieTayAlpha, &
    L1_aeroResist, &
    L1_alpha, &
    L1_degDay, &
    L1_degDayInc, &
    L1_degDayMax, &
    L1_degDayNoPre, &
    L1_fAsp, &
    L1_fRoots, &
    L1_fSealed, &
    L1_jarvis_thresh_c1, &
    L1_kBaseFlow, &
    L1_kPerco, &
    L1_kSlowFlow, &
    L1_karstLoss, &
    L1_kfastFlow, &
    L1_maxInter, &
    L1_petLAIcorFactor, &
    L1_sealedThresh, &
    L1_soilMoistExp, &
    L1_soilMoistFC, &
    L1_soilMoistSat, &
    L1_surfResist, &
    L1_tempThresh, &
    L1_unsatThresh, &
    L1_wiltingPoint, &
    L1_No_Count, &
    L1_bulkDens, &
    L1_latticeWater, &
    L1_COSMICL3, &
    HorizonDepth_mHM, &
    nSoilHorizons_mHM
  use mo_mrm_init, only : variables_default_init_routing
  use mo_mrm_mpr, only : mrm_update_param
  use mo_mrm_restart, only : mrm_read_restart_states
  use mo_mrm_routing, only : mrm_routing
  use mo_mrm_write, only : mrm_write_output_fluxes
  use mo_utils, only : ge
  use mo_mrm_river_head, only: calc_river_head, avg_and_write_timestep
  use mo_mpr_eval, only : mpr_eval

contains

  !> \brief prepare single run of mHM
  subroutine mhm_interface_run_prepare(parameterset, opti_domain_indices, runoff_present, BFI_present)
    implicit none
    !> a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), optional, intent(in) :: parameterset
    !> selected domains for optimization
    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices
    !> whether runoff is present
    logical, optional, intent(in) :: runoff_present
    !> whether BFI is present
    logical, optional, intent(in) :: BFI_present

    integer(i4) :: i, iDomain, domainID

    ! run_cfg%output_runoff = .false. by default
    if (present(runoff_present)) run_cfg%output_runoff = runoff_present
    ! run_cfg%output_BFI = .false. by default
    if (present(BFI_present)) run_cfg%output_BFI = BFI_present

    ! store current parameter set
    allocate(run_cfg%parameterset(size(global_parameters, dim=1)))
    if (.not. present(parameterset) .and. optimize) then
      call error_message("mhm_interface_run_prepare: Can't optimize without parameter!")
    else if (.not. present(parameterset)) then
      run_cfg%parameterset = global_parameters(:, 3)
    else
      run_cfg%parameterset = parameterset
    end if

    ! store domain indices for domain loop
    if (optimize .and. present(opti_domain_indices)) then
      run_cfg%nDomains = size(opti_domain_indices)
      allocate(run_cfg%domain_indices(run_cfg%nDomains))
      run_cfg%domain_indices = opti_domain_indices
    else
      run_cfg%nDomains = domainMeta%nDomains
      allocate(run_cfg%domain_indices(run_cfg%nDomains))
      run_cfg%domain_indices = [(i, i=1, run_cfg%nDomains)]
    end if

    run_cfg%is_hourly_forcing = (nTstepForcingDay .eq. 24_i4)

    !----------------------------------------------------------
    ! Check optionals and initialize
    !----------------------------------------------------------
    if (run_cfg%output_runoff) then
      do i = 1, run_cfg%nDomains
        iDomain = run_cfg%get_domain_index(i)
        domainID = domainMeta%indices(iDomain)
        if (.not. domainMeta%doRouting(iDomain)) then
          call error_message("***ERROR: runoff for domain", trim(num2str(domainID)),&
                        "can not be produced, since routing process is off in Process Matrix")
        end if
      end do
    end if

    ! prepare BFI calculation
    if (run_cfg%output_BFI) then
      allocate(BFI_qBF_sum(run_cfg%nDomains))
      allocate(BFI_qT_sum(run_cfg%nDomains))
      BFI_qBF_sum = 0.0_dp
      BFI_qT_sum = 0.0_dp
    end if

    if (read_restart) then
      do i = 1, run_cfg%nDomains
        iDomain = run_cfg%get_domain_index(i)
        domainID = domainMeta%indices(iDomain)
        ! this reads the eff. parameters and optionally the states and fluxes
        call read_restart_states(iDomain, domainID, mhmFileRestartIn(iDomain))
      end do
    else
      call variables_default_init()
      call mpr_eval(run_cfg%parameterset)
      if (processMatrix(8, 1) > 0) call variables_default_init_routing()
    end if

    allocate(run_cfg%L1_fNotSealed(size(L1_fSealed, 1), size(L1_fSealed, 2), size(L1_fSealed, 3)))
    run_cfg%L1_fNotSealed = 1.0_dp - L1_fSealed

  end subroutine mhm_interface_run_prepare

  !> \brief get number of domains for looping
  subroutine mhm_interface_run_get_ndomains(ndomains)
    implicit none
    integer(i4), intent(inout) :: ndomains
    ndomains = run_cfg%nDomains
  end subroutine mhm_interface_run_get_ndomains

  !> \brief prepare single domain to run mHM on
  subroutine mhm_interface_run_prepare_domain(domain, etOptiSim, twsOptiSim, neutronsOptiSim, smOptiSim)
    implicit none
    !> domain loop counter
    integer(i4), intent(in), optional :: domain
    !> returns soil moisture time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim
    !> dim1=ncells, dim2=time
    type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim
    !> returns evapotranspiration time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim
    !> returns tws time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim

    integer(i4) :: iDomain, domainID

    ! get domain index
    run_cfg%selected_domain = 1_i4
    if (present(domain)) run_cfg%selected_domain = domain
    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)
    domainID = domainMeta%indices(iDomain)

    ! evapotranspiration optimization
    if (present(etOptiSim)) call etOptiSim(iDomain)%init(L1_etObs(iDomain))
    ! total water storage optimization
    if (present(twsOptiSim)) call twsOptiSim(iDomain)%init(L1_twsaObs(iDomain))
    ! neutrons optimization
    if (present(neutronsOptiSim)) call neutronsOptiSim(iDomain)%init(L1_neutronsObs(iDomain))
    ! sm optimization
    if (present(smOptiSim)) call smOptiSim(iDomain)%init(L1_smObs(iDomain))

    ! get Domain information
    run_cfg%nCells = level1(iDomain)%nCells
    run_cfg%mask1 => level1(iDomain)%mask
    run_cfg%s1 = level1(iDomain)%iStart
    run_cfg%e1 = level1(iDomain)%iEnd

    if (domainMeta%doRouting(iDomain)) then
      ! ----------------------------------------
      ! initialize factor between routing resolution and hydrologic model resolution
      ! ----------------------------------------
      run_cfg%tsRoutFactor = 1_i4
      allocate(run_cfg%InflowDischarge(size(InflowGauge%Q, dim = 2)))
      run_cfg%InflowDischarge = 0._dp

      ! read states from restart
      if (read_restart) call mrm_read_restart_states(iDomain, domainID, mrmFileRestartIn(iDomain))

      ! get Domain information at L11 if routing is activated
      run_cfg%s11 = level11(iDomain)%iStart
      run_cfg%e11 = level11(iDomain)%iEnd
      run_cfg%mask11 => level11(iDomain)%mask

      ! initialize routing parameters (has to be called for routing option 2)
      if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) &
          call mrm_update_param( &
            iDomain, run_cfg%parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)))
      ! initialize variable for runoff for routing
      allocate(run_cfg%RunToRout(run_cfg%e1 - run_cfg%s1 + 1))
      run_cfg%RunToRout = 0._dp

      if ( riv_temp_pcs%active ) then
        ! set indices for current L11 domain
        riv_temp_pcs%s11 = run_cfg%s11
        riv_temp_pcs%e11 = run_cfg%e11
        ! allocate current L1 lateral components
        call riv_temp_pcs%alloc_lateral(run_cfg%nCells)
      end if
    end if

    ! init datetime variable
    call run_cfg%domainDateTime%init(iDomain)

    ! init time-loop
    run_cfg%time_step = 0_i4

  end subroutine mhm_interface_run_prepare_domain

  !> \brief check if current time loop is finished
  subroutine mhm_interface_run_finished(time_loop_finished)
    implicit none
    logical, intent(inout) :: time_loop_finished
    time_loop_finished = run_cfg%time_step == run_cfg%domainDateTime%nTimeSteps
  end subroutine mhm_interface_run_finished

  !> \brief do one time-step on current domain
  subroutine mhm_interface_run_do_time_step()
    implicit none

    integer(i4) :: iDomain, domainID, tt, jj

    ! increment time step count (first input is 0)
    run_cfg%time_step = run_cfg%time_step + 1_i4
    tt = run_cfg%time_step

    ! get domain index
    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)
    domainID = domainMeta%indices(iDomain)

    ! time increment is done right after call to mrm (and initially before looping)
    if (timeStep_model_inputs(iDomain) .eq. 0_i4) then
      ! whole meteorology is already read

      ! set start and end of meteo position
      run_cfg%s_meteo = run_cfg%s1
      run_cfg%e_meteo = run_cfg%e1
      ! time step for meteorological variable (daily values)
      ! iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
      run_cfg%iMeteoTS = ceiling(real(tt, dp) / real(nint( 24._dp / real(nTstepForcingDay, dp)), dp))
    else
      ! read chunk of meteorological forcings data (reading, upscaling/downscaling)
      call prepare_meteo_forcings_data(iDomain, domainID, tt)
      ! set start and end of meteo position
      run_cfg%s_meteo = 1
      run_cfg%e_meteo = run_cfg%e1 - run_cfg%s1 + 1
      ! time step for meteorological variable (daily values)
      run_cfg%iMeteoTS = &
        ceiling(real(tt, dp) / real(nint( 24._dp / real(nTstepForcingDay, dp)), dp)) &
        - (readPer%julStart - simPer(iDomain)%julStart)
    end if

    ! preapare vector length specifications depending on the process case
    ! process 5 - PET
    select case (processMatrix(5, 1))
      !      [pet,        tmax,    tmin,  netrad, absVapP,windspeed]
      case(-1 : 0) ! PET is input
        run_cfg%s_p5 = [run_cfg%s_meteo, 1, 1, 1, 1, 1]
        run_cfg%e_p5 = [run_cfg%e_meteo, 1, 1, 1, 1, 1]
      case(1) ! Hargreaves-Samani
        run_cfg%s_p5 = [run_cfg%s_meteo, run_cfg%s_meteo, run_cfg%s_meteo, 1, 1, 1]
        run_cfg%e_p5 = [run_cfg%e_meteo, run_cfg%e_meteo, run_cfg%e_meteo, 1, 1, 1]
      case(2) ! Priestely-Taylor
        run_cfg%s_p5 = [run_cfg%s_meteo, 1, 1, run_cfg%s_meteo, 1, 1]
        run_cfg%e_p5 = [run_cfg%e_meteo, 1, 1, run_cfg%e_meteo, 1, 1]
      case(3) ! Penman-Monteith
        run_cfg%s_p5 = [run_cfg%s_meteo, 1, 1, run_cfg%s_meteo, run_cfg%s_meteo, run_cfg%s_meteo]
        run_cfg%e_p5 = [run_cfg%e_meteo, 1, 1, run_cfg%e_meteo, run_cfg%e_meteo, run_cfg%e_meteo]
    end select

    ! customize iMeteoTS for process 5 - PET
    select case (processMatrix(5, 1))
      !              [     pet,     tmin,     tmax,   netrad,  absVapP,windspeed ]
      case(-1 : 0) ! PET is input
        run_cfg%iMeteo_p5 = [run_cfg%iMeteoTS, 1, 1, 1, 1, 1 ]
      case(1) ! Hargreaves-Samani
        run_cfg%iMeteo_p5 = [run_cfg%iMeteoTS, run_cfg%iMeteoTS, run_cfg%iMeteoTS, 1, 1, 1 ]
      case(2) ! Priestely-Taylor
        run_cfg%iMeteo_p5 = [run_cfg%iMeteoTS, 1, 1, run_cfg%iMeteoTS, 1, 1 ]
      case(3) ! Penman-Monteith
        run_cfg%iMeteo_p5 = [run_cfg%iMeteoTS, 1, 1, run_cfg%iMeteoTS, run_cfg%iMeteoTS, run_cfg%iMeteoTS ]
    end select

    call run_cfg%domainDateTime%update_LAI_timestep()

    ! -------------------------------------------------------------------------
    ! ARGUMENT LIST KEY FOR mHM
    ! -------------------------------------------------------------------------
    !  C    CONFIGURATION
    !  F    FORCING DATA L2
    !  Q    INFLOW FROM UPSTREAM AREAS
    !  L0   MORPHOLOGIC DATA L0
    !  L1   MORPHOLOGIC DATA L1
    !  L11  MORPHOLOGIC DATA L11
    !  P    GLOBAL PARAMETERS
    !  E    EFFECTIVE PARAMETER FIELDS (L1, L11 levels)
    !  S    STATE VARIABLES L1
    !  X    FLUXES (L1, L11 levels)
    ! --------------------------------------------------------------------------
    call mhm( &
      read_restart, run_cfg%is_hourly_forcing, & ! IN C
      tt, run_cfg%domainDateTime%newTime - 0.5_dp, processMatrix, &
      HorizonDepth_mHM, & ! IN C
      run_cfg%nCells, nSoilHorizons_mHM, real(nTstepDay, dp), c2TSTu,  & ! IN C
      neutron_integral_AFast, & ! IN C
      pack(level1(iDomain)%y, level1(iDomain)%mask), & ! IN L1
      evap_coeff, fday_prec, fnight_prec, fday_pet, fnight_pet, & ! IN F
      fday_temp, fnight_temp, & ! IN F
      L1_temp_weights(run_cfg%s1 : run_cfg%e1, :, :), & ! IN F
      L1_pet_weights(run_cfg%s1 : run_cfg%e1, :, :), & ! IN F
      L1_pre_weights(run_cfg%s1 : run_cfg%e1, :, :), & ! IN F
      read_meteo_weights, & ! IN F
      L1_pet(run_cfg%s_p5(1) : run_cfg%e_p5(1), run_cfg%iMeteo_p5(1)), & ! INOUT F:PET
      L1_tmin(run_cfg%s_p5(2) : run_cfg%e_p5(2), run_cfg%iMeteo_p5(2)), & ! IN F:PET
      L1_tmax(run_cfg%s_p5(3) : run_cfg%e_p5(3), run_cfg%iMeteo_p5(3)), & ! IN F:PET
      L1_netrad(run_cfg%s_p5(4) : run_cfg%e_p5(4), run_cfg%iMeteo_p5(4)), & ! IN F:PET
      L1_absvappress(run_cfg%s_p5(5) : run_cfg%e_p5(5), run_cfg%iMeteo_p5(5)), & ! IN F:PET
      L1_windspeed(run_cfg%s_p5(6) : run_cfg%e_p5(6), run_cfg%iMeteo_p5(6)), & ! IN F:PET
      L1_pre(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), & ! IN F:Pre
      L1_temp(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), & ! IN F:Temp
      L1_fSealed(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), & ! INOUT L1
      L1_inter(run_cfg%s1 : run_cfg%e1), &
      L1_snowPack(run_cfg%s1 : run_cfg%e1), &
      L1_sealSTW(run_cfg%s1 : run_cfg%e1), & ! INOUT S
      L1_soilMoist(run_cfg%s1 : run_cfg%e1, :), &
      L1_unsatSTW(run_cfg%s1 : run_cfg%e1), &
      L1_satSTW(run_cfg%s1 : run_cfg%e1), & ! INOUT S
      L1_neutrons(run_cfg%s1 : run_cfg%e1), & ! INOUT S
      L1_pet_calc(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_temp_calc(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_prec_calc(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_aETSoil(run_cfg%s1 : run_cfg%e1, :), &
      L1_aETCanopy(run_cfg%s1 : run_cfg%e1), &
      L1_aETSealed(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_baseflow(run_cfg%s1 : run_cfg%e1), &
      L1_infilSoil(run_cfg%s1 : run_cfg%e1, :), &
      L1_fastRunoff(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_melt(run_cfg%s1 : run_cfg%e1), &
      L1_percol(run_cfg%s1 : run_cfg%e1), &
      L1_preEffect(run_cfg%s1 : run_cfg%e1), &
      L1_rain(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_runoffSeal(run_cfg%s1 : run_cfg%e1), &
      L1_slowRunoff(run_cfg%s1 : run_cfg%e1), &
      L1_snow(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      L1_Throughfall(run_cfg%s1 : run_cfg%e1), &
      L1_total_runoff(run_cfg%s1 : run_cfg%e1), & ! INOUT X
      ! MPR
      L1_alpha(run_cfg%s1 : run_cfg%e1, 1, 1), &
      L1_degDayInc(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), &
      L1_degDayMax(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_degDayNoPre(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), &
      L1_degDay(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_fAsp(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_petLAIcorFactor(run_cfg%s1 : run_cfg%e1, run_cfg%domainDateTime%iLAI, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_HarSamCoeff(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_PrieTayAlpha(run_cfg%s1 : run_cfg%e1, run_cfg%domainDateTime%iLAI, 1), & ! INOUT E1
      L1_aeroResist(run_cfg%s1 : run_cfg%e1, run_cfg%domainDateTime%iLAI, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_surfResist(run_cfg%s1 : run_cfg%e1, run_cfg%domainDateTime%iLAI, 1), &
      L1_fRoots(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_maxInter(run_cfg%s1 : run_cfg%e1, run_cfg%domainDateTime%iLAI, 1), &
      L1_karstLoss(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_kFastFlow(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), &
      L1_kSlowFlow(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_kBaseFlow(run_cfg%s1 : run_cfg%e1, 1, 1), &
      L1_kPerco(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_soilMoistFC(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_soilMoistSat(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_soilMoistExp(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), &
      L1_jarvis_thresh_c1(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_tempThresh(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), &
      L1_unsatThresh(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_sealedThresh(run_cfg%s1 : run_cfg%e1, 1, 1), & ! INOUT E1
      L1_wiltingPoint(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), & ! INOUT E1
      !>> neutron count
      L1_No_Count(run_cfg%s1:run_cfg%e1, 1, 1),  &                     ! INOUT E1
      L1_bulkDens(run_cfg%s1:run_cfg%e1,     :, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_latticeWater(run_cfg%s1:run_cfg%e1, :, run_cfg%domainDateTime%yId), & ! INOUT E1
      L1_COSMICL3(run_cfg%s1:run_cfg%e1,     :, run_cfg%domainDateTime%yId)  & ! INOUT E1
    )

    ! call mRM routing
    run_cfg%doRoute = .false.
    if (domainMeta%doRouting(iDomain)) then
      ! set discharge timestep
      run_cfg%iDischargeTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
      ! set input variables for routing
      if (processMatrix(8, 1) .eq. 1) then
        ! >>>
        ! >>> original Muskingum routing, executed every time
        ! >>>
        run_cfg%doRoute = .True.
        run_cfg%tsRoutFactorIn = 1._dp
        run_cfg%timestep_rout = timestep
        run_cfg%RunToRout = L1_total_runoff(run_cfg%s1 : run_cfg%e1) ! runoff [mm TS-1] mm per timestep
        run_cfg%InflowDischarge = InflowGauge%Q(run_cfg%iDischargeTS, :) ! inflow discharge in [m3 s-1]
        !
      else if ((processMatrix(8, 1) .eq. 2) .or. &
               (processMatrix(8, 1) .eq. 3)) then
        ! >>>
        ! >>> adaptive timestep
        ! >>>
        run_cfg%doRoute = .False.
        ! calculate factor
        run_cfg%tsRoutFactor = L11_tsRout(iDomain) / (timestep * HourSecs)
        ! print *, 'routing factor: ', tsRoutFactor
        ! prepare routing call
        if (run_cfg%tsRoutFactor .lt. 1._dp) then
          ! ----------------------------------------------------------------
          ! routing timesteps are shorter than hydrologic time steps
          ! ----------------------------------------------------------------
          ! set all input variables
          run_cfg%tsRoutFactorIn = run_cfg%tsRoutFactor
          run_cfg%RunToRout = L1_total_runoff(run_cfg%s1 : run_cfg%e1) ! runoff [mm TS-1] mm per timestep
          run_cfg%InflowDischarge = InflowGauge%Q(run_cfg%iDischargeTS, :) ! inflow discharge in [m3 s-1]
          run_cfg%timestep_rout = timestep
          run_cfg%doRoute = .True.
        else
          ! ----------------------------------------------------------------
          ! routing timesteps are longer than hydrologic time steps
          ! ----------------------------------------------------------------
          ! set all input variables
          run_cfg%tsRoutFactorIn = run_cfg%tsRoutFactor
          ! Runoff is accumulated in [mm]
          run_cfg%RunToRout = run_cfg%RunToRout + L1_total_runoff(run_cfg%s1 : run_cfg%e1)
          run_cfg%InflowDischarge = run_cfg%InflowDischarge + InflowGauge%Q(run_cfg%iDischargeTS, :)
          ! reset tsRoutFactorIn if last period did not cover full period
          if ((tt == run_cfg%domainDateTime%nTimeSteps) .and. (mod(tt, nint(run_cfg%tsRoutFactorIn)) /= 0_i4)) &
            run_cfg%tsRoutFactorIn = mod(tt, nint(run_cfg%tsRoutFactorIn))
          if ((mod(tt, nint(run_cfg%tsRoutFactorIn)) .eq. 0_i4) .or. (tt .eq. run_cfg%domainDateTime%nTimeSteps)) then
            ! Inflow discharge is given as flow-rate and has to be converted to [m3]
            run_cfg%InflowDischarge = run_cfg%InflowDischarge / run_cfg%tsRoutFactorIn
            run_cfg%timestep_rout = timestep * nint(run_cfg%tsRoutFactorIn, i4)
            run_cfg%doRoute = .True.
          end if
        end if
      end if
      ! prepare temperature routing
      if ( riv_temp_pcs%active ) then
        ! init riv-temp from current air temp
        if ( tt .eq. 1_i4 ) call riv_temp_pcs%init_riv_temp( &
          run_cfg%domainDateTime%newTime - 0.5_dp, &
          real(nTstepDay, dp), &
          L1_temp(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), &
          read_meteo_weights, &
          L1_temp_weights(run_cfg%s1 : run_cfg%e1, :, :), &
          fday_temp, fnight_temp, &
          ! mapping info
          level1(iDomain)%CellArea * 1.E-6_dp, &
          L1_L11_Id(run_cfg%s1 : run_cfg%e1), &
          level11(iDomain)%CellArea * 1.E-6_dp, &
          L11_L1_Id(run_cfg%s11 : run_cfg%e11), &
          ! map_flag
          ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)) &
        )
        ! accumulate source Energy at L1 level
        call riv_temp_pcs%acc_source_E( &
          run_cfg%domainDateTime%newTime - 0.5_dp, &
          real(nTstepDay, dp), &
          L1_fSealed(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId), &
          L1_fastRunoff(run_cfg%s1 : run_cfg%e1), &
          L1_slowRunoff(run_cfg%s1 : run_cfg%e1), &
          L1_baseflow(run_cfg%s1 : run_cfg%e1), &
          L1_runoffSeal(run_cfg%s1 : run_cfg%e1), &
          L1_temp(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), &
          L1_tann(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), &
          L1_ssrd(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), &
          L1_strd(run_cfg%s_meteo : run_cfg%e_meteo, run_cfg%iMeteoTS), &
          read_meteo_weights, &
          L1_temp_weights(run_cfg%s1 : run_cfg%e1, :, :), &
          fday_temp, fnight_temp, &
          fday_ssrd, fnight_ssrd, &
          fday_strd, fnight_strd  &
        )
        ! if routing should be performed, scale source energy to L11 level
        if ( run_cfg%doRoute ) call riv_temp_pcs%finalize_source_E( &
          level1(iDomain)%CellArea * 1.E-6_dp, &
          L1_L11_Id(run_cfg%s1 : run_cfg%e1), &
          level11(iDomain)%CellArea * 1.E-6_dp, &
          L11_L1_Id(run_cfg%s11 : run_cfg%e11), &
          run_cfg%timestep_rout, &
          ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)) &
        )
      end if
      ! -------------------------------------------------------------------
      ! execute routing
      ! -------------------------------------------------------------------
      if (run_cfg%doRoute) call mRM_routing(&
        ! general INPUT variables
        read_restart, &
        processMatrix(8, 1), & ! parse process Case to be used
        run_cfg%parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)), & ! routing par.
        run_cfg%RunToRout, & ! runoff [mm TS-1] mm per timestep old: L1_total_runoff_in(run_cfg%s1:run_cfg%e1, tt), &
        level1(iDomain)%CellArea * 1.E-6_dp, &
        L1_L11_Id(run_cfg%s1 : run_cfg%e1), &
        level11(iDomain)%CellArea * 1.E-6_dp, &
        L11_L1_Id(run_cfg%s11 : run_cfg%e11), &
        L11_netPerm(run_cfg%s11 : run_cfg%e11), & ! routing order at L11
        L11_fromN(run_cfg%s11 : run_cfg%e11), & ! link source at L11
        L11_toN(run_cfg%s11 : run_cfg%e11), & ! link target at L11
        L11_nOutlets(iDomain), & ! number of outlets
        run_cfg%timestep_rout, & ! timestep of runoff to rout [h]
        run_cfg%tsRoutFactorIn, & ! simulate timestep in [h]
        level11(iDomain)%nCells, & ! number of Nodes
        domain_mrm(iDomain)%nInflowGauges, &
        domain_mrm(iDomain)%InflowGaugeIndexList(:), &
        domain_mrm(iDomain)%InflowGaugeHeadwater(:), &
        domain_mrm(iDomain)%InflowGaugeNodeList(:), &
        run_cfg%InflowDischarge, &
        domain_mrm(iDomain)%nGauges, &
        domain_mrm(iDomain)%gaugeIndexList(:), &
        domain_mrm(iDomain)%gaugeNodeList(:), &
        ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)), &
        ! original routing specific input variables
        L11_length(run_cfg%s11 : run_cfg%e11 - 1), & ! link length
        L11_slope(run_cfg%s11 : run_cfg%e11 - 1), &
        L11_nLinkFracFPimp(run_cfg%s11 : run_cfg%e11, run_cfg%domainDateTime%yId), & ! fraction of impervious layer at L11 scale
        ! general INPUT/OUTPUT variables
        L11_C1(run_cfg%s11 : run_cfg%e11), & ! first muskingum parameter
        L11_C2(run_cfg%s11 : run_cfg%e11), & ! second muskigum parameter
        L11_qOUT(run_cfg%s11 : run_cfg%e11), & ! routed runoff flowing out of L11 cell
        L11_qTIN(run_cfg%s11 : run_cfg%e11, :), & ! inflow water into the reach at L11
        L11_qTR(run_cfg%s11 : run_cfg%e11, :), & !
        L11_qMod(run_cfg%s11 : run_cfg%e11), &
        mRM_runoff(tt, :) &
      )
      ! -------------------------------------------------------------------
      ! groundwater coupling
      ! -------------------------------------------------------------------
      if (gw_coupling) then
          call calc_river_head(iDomain, L11_Qmod, L0_river_head_mon_sum)
          if (run_cfg%domainDateTime%is_new_month .and. tt > 1) then
              call avg_and_write_timestep(iDomain, tt, L0_river_head_mon_sum)
          end if
      end if
      ! -------------------------------------------------------------------
      ! reset variables
      ! -------------------------------------------------------------------
      if (processMatrix(8, 1) .eq. 1) then
        ! reset Input variables
        run_cfg%InflowDischarge = 0._dp
        run_cfg%RunToRout = 0._dp
      else if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) then
        if ((.not. (run_cfg%tsRoutFactorIn .lt. 1._dp)) .and. run_cfg%doRoute) then
          do jj = 1, nint(run_cfg%tsRoutFactorIn) ! BUG: this should start at 2
            mRM_runoff(tt - jj + 1, :) = mRM_runoff(tt, :)
          end do
          ! reset Input variables
          run_cfg%InflowDischarge = 0._dp
          run_cfg%RunToRout = 0._dp
          ! reset lateral fluxes and time-step counter if routing was done
          if ( riv_temp_pcs%active ) call riv_temp_pcs%reset_timestep()
        end if
        ! if routing is done every time-step, reset river-temp time step
        if (run_cfg%tsRoutFactor .lt. 1._dp .and. riv_temp_pcs%active ) call riv_temp_pcs%reset_timestep()
      end if
    end if

    ! output only for evaluation period
    run_cfg%domainDateTime%tIndex_out = (tt - warmingDays(iDomain) * nTstepDay) ! tt if write out of warming period

    call run_cfg%domainDateTime%increment()

    ! update the year-dependent domainDateTime%yId (land cover id)
    if (run_cfg%domainDateTime%is_new_year .and. tt < run_cfg%domainDateTime%nTimeSteps) then
      run_cfg%domainDateTime%yId = LCyearId(run_cfg%domainDateTime%year, iDomain)
    end if

    ! calculate BFI releated after warming days if wanted
    if ( run_cfg%output_BFI .and. (run_cfg%domainDateTime%tIndex_out > 0_i4) ) then
      BFI_qBF_sum(iDomain) = BFI_qBF_sum(iDomain) &
        + sum(L1_baseflow(run_cfg%s1 : run_cfg%e1) * level1(iDomain)%CellArea) / level1(iDomain)%nCells
      BFI_qT_sum(iDomain) = BFI_qT_sum(iDomain) &
        + sum(L1_total_runoff(run_cfg%s1 : run_cfg%e1) * level1(iDomain)%CellArea) / level1(iDomain)%nCells
    end if

  end subroutine mhm_interface_run_do_time_step

  !> \brief write output after current time-step
  subroutine mhm_interface_run_write_output()
    implicit none
    integer(i4) :: iDomain, tt

    ! get time step
    tt = run_cfg%time_step

    ! get domain index
    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)

    if ( .not. optimize ) then
      if (any(outputFlxState_mrm) .AND. (domainMeta%doRouting(iDomain))) then
        call mrm_write_output_fluxes( &
          iDomain, & ! Domain id
          level11(iDomain)%nCells, & ! nCells in Domain
          timeStep_model_outputs_mrm, & ! output specification
          run_cfg%domainDateTime, tt, timestep, & ! time specification
          run_cfg%mask11, & ! mask specification
          L11_qmod(run_cfg%s11 : run_cfg%e11) & ! output variables
        )
      end if

      if ((any(outputFlxState)) .and. (run_cfg%domainDateTime%tIndex_out > 0_i4)) then

        if (run_cfg%domainDateTime%tIndex_out == 1) then
          run_cfg%nc = OutputDataset(iDomain, run_cfg%mask1, level1(iDomain)%nCells)
        end if

        call run_cfg%nc%updateDataset(&
          run_cfg%s1, run_cfg%e1, &
          L1_fSealed(:, 1, run_cfg%domainDateTime%yId), &
          run_cfg%L1_fNotSealed(:, 1, run_cfg%domainDateTime%yId), &
          L1_inter, &
          L1_snowPack, &
          L1_soilMoist, &
          L1_soilMoistSat(:, :, run_cfg%domainDateTime%yId), &
          L1_sealSTW, &
          L1_unsatSTW, &
          L1_satSTW, &
          L1_neutrons, &
          L1_pet_calc, &
          L1_aETSoil, &
          L1_aETCanopy, &
          L1_aETSealed, &
          L1_total_runoff, &
          L1_runoffSeal, &
          L1_fastRunoff, &
          L1_slowRunoff, &
          L1_baseflow, &
          L1_percol, &
          L1_infilSoil, &
          L1_preEffect, &
          L1_melt &
        )

        ! write data
        if (run_cfg%domainDateTime%writeout(timeStep_model_outputs, tt)) then
          call run_cfg%nc%writeTimestep(run_cfg%domainDateTime%tIndex_out * timestep - 1)
        end if

        if(tt == run_cfg%domainDateTime%nTimeSteps) then
          call run_cfg%nc%close()
        end if

      end if
    end if ! <-- if (.not. optimize)

  end subroutine mhm_interface_run_write_output

  !> \brief add simulation data to optimization data types
  subroutine mhm_interface_run_update_optisim(etOptiSim, twsOptiSim, neutronsOptiSim, smOptiSim)
    implicit none
    !> returns soil moisture time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim
    !> dim1=ncells, dim2=time
    type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim
    !> returns evapotranspiration time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim
    !> returns tws time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim

    integer(i4) :: gg, s1, e1, lcId

    integer(i4) :: iDomain, tt

    ! get time step
    tt = run_cfg%time_step

    ! get domain index
    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)

    ! slice settings
    s1 = run_cfg%s1
    e1 = run_cfg%e1
    lcId = run_cfg%domainDateTime%yId

    !----------------------------------------------------------------------
    ! FOR SOIL MOISTURE
    ! NOTE:: modeled soil moisture is averaged according to input time step
    !        soil moisture (timeStep_sm_input)
    !----------------------------------------------------------------------
    if (present(smOptiSim)) then
      ! only for evaluation period - ignore warming days
      if ((tt - warmingDays(iDomain) * nTstepDay) .GT. 0) then
        ! decide for daily, monthly or yearly aggregation
        call smOptiSim(iDomain)%average_per_timestep(L1_smObs(iDomain)%timeStepInput, &
                     run_cfg%domainDateTime%is_new_day, run_cfg%domainDateTime%is_new_month, run_cfg%domainDateTime%is_new_year)
        ! last timestep is already done - write_counter exceeds size(smOptiSim(iDomain)%dataSim, dim=2)
        if (tt /= run_cfg%domainDateTime%nTimeSteps) then
          ! aggregate soil moisture to needed time step for optimization
          call smOptiSim(iDomain)%average_add(sum(L1_soilMoist(:, 1 : nSoilHorizons_sm_input), dim = 2) / &
                          sum(L1_soilMoistSat(:, 1 : nSoilHorizons_sm_input, lcId), dim = 2))
        end if
      end if
    end if

    !----------------------------------------------------------------------
    ! FOR NEUTRONS
    ! NOTE:: modeled neutrons are averaged daily
    !----------------------------------------------------------------------
    if (present(neutronsOptiSim)) then
      ! only for evaluation period - ignore warming days
      if ((tt - warmingDays(iDomain) * nTstepDay) .GT. 0) then
        ! decide for daily, monthly or yearly aggregation
        ! daily
        if (run_cfg%domainDateTime%is_new_day)   then
          call neutronsOptiSim(iDomain)%average()
        end if

        ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
        if (tt /= run_cfg%domainDateTime%nTimeSteps) then
          ! aggregate neutrons to needed time step for optimization
          call neutronsOptiSim(iDomain)%average_add(L1_neutrons(s1 : e1))
        end if
      end if
    end if

    !----------------------------------------------------------------------
    ! FOR EVAPOTRANSPIRATION
    ! NOTE:: modeled evapotranspiration is averaged according to input time step
    !        evapotranspiration (timeStep_et_input)
    !----------------------------------------------------------------------
    if (present(etOptiSim)) then
      ! only for evaluation period - ignore warming days
      if ((tt - warmingDays(iDomain) * nTstepDay) .GT. 0) then
        ! decide for daily, monthly or yearly aggregation
        call etOptiSim(iDomain)%increment_counter(L1_etObs(iDomain)%timeStepInput, &
          run_cfg%domainDateTime%is_new_day, run_cfg%domainDateTime%is_new_month, run_cfg%domainDateTime%is_new_year)

        ! last timestep is already done - write_counter exceeds size(etOptiSim(iDomain)%dataSim, dim=2)
        if (tt /= run_cfg%domainDateTime%nTimeSteps) then
          ! aggregate evapotranspiration to needed time step for optimization
          call etOptiSim(iDomain)%add(sum(L1_aETSoil(s1 : e1, :), dim = 2) * &
            run_cfg%L1_fNotSealed(s1 : e1, 1, lcId) + &
            L1_aETCanopy(s1 : e1) + &
            L1_aETSealed(s1 : e1) * L1_fSealed(s1 : e1, 1, lcId))
        end if
      end if
    end if

    !----------------------------------------------------------------------
    ! FOR TWS
    ! NOTE:: modeled tws is averaged according to input time step
    !        (timeStepInput)
    !----------------------------------------------------------------------
    if (present(twsOptiSim)) then
      ! only for evaluation period - ignore warming days
      if ((tt - warmingDays(iDomain) * nTstepDay) > 0) then
        ! decide for daily, monthly or yearly aggregation
        call twsOptiSim(iDomain)%average_per_timestep(L1_twsaObs(iDomain)%timeStepInput, &
          run_cfg%domainDateTime%is_new_day, run_cfg%domainDateTime%is_new_month, run_cfg%domainDateTime%is_new_year)

        ! last timestep is already done - write_counter exceeds size(twsOptiSim(iDomain)%dataSim, dim=2)
        if (tt /= run_cfg%domainDateTime%nTimeSteps) then
          ! aggregate evapotranspiration to needed time step for optimization
          call twsOptiSim(iDomain)%average_add(L1_inter(s1 : e1) + L1_snowPack(s1 : e1) + L1_sealSTW(s1 : e1) + &
            L1_unsatSTW(s1 : e1) + L1_satSTW(s1 : e1))
          do gg = 1, nSoilHorizons_mHM
            call twsOptiSim(iDomain)%add(L1_soilMoist (s1 : e1, gg))
          end do
        end if
      end if
    end if

    ! TODO-RIV-TEMP: add OptiSim for river temperature

  end subroutine mhm_interface_run_update_optisim

  !> \brief finalize current domain after running
  subroutine mhm_interface_run_finalize_domain()
    implicit none

    integer(i4) :: iDomain

    ! get domain index
    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)

    if (allocated(run_cfg%InflowDischarge)) deallocate(run_cfg%InflowDischarge)
      if (domainMeta%doRouting(iDomain)) then
      ! clean runoff variable
      deallocate(run_cfg%RunToRout)
      if ( riv_temp_pcs%active ) call riv_temp_pcs%dealloc_lateral()
    end if

  end subroutine mhm_interface_run_finalize_domain

  !> \brief finalize run
  subroutine mhm_interface_run_finalize(runoff, BFI)
    implicit none
    !> returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: runoff
    !> baseflow index, dim1=domainID
    real(dp), dimension(:), allocatable, optional, intent(out) :: BFI

    if (present(runoff) .and. (processMatrix(8, 1) > 0)) runoff = mRM_runoff
    if (present(BFI)) then
      BFI = BFI_qBF_sum / BFI_qT_sum
      deallocate(BFI_qBF_sum)
      deallocate(BFI_qT_sum)
    end if

    if (allocated(run_cfg%parameterset)) deallocate(run_cfg%parameterset)
    if (allocated(run_cfg%L1_fNotSealed)) deallocate(run_cfg%L1_fNotSealed)
    if (allocated(run_cfg%domain_indices)) deallocate(run_cfg%domain_indices)

  end subroutine mhm_interface_run_finalize

end module mo_mhm_interface_run
