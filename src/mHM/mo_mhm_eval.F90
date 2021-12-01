!>       \file mo_mhm_eval.f90

!>       \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

!>       \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

!>       \authors Juliane Mai, Rohini Kumar

!>       \date Feb 2013

! Modifications:

MODULE mo_mhm_eval

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mhm_eval

  !    PURPOSE
  !>       \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !>       \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset" a set of global parameter (gamma) to run mHM, DIMENSION
  !>       [no. of global_Parameters]

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), dimension(:, :), optional :: runoff"        returns runoff time series, DIMENSION
  !>       [nTimeSteps, nGaugesTotal]
  !>       \param[out] "real(dp), dimension(:, :), optional :: sm_opti"       returns soil moisture time series for all
  !>       grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
  !>       time series, DIMENSION [nTimeSteps, nDomains]
  !>       \param[out] "real(dp), dimension(:, :), optional :: neutrons_opti" dim1=ncells, dim2=time
  !>       \param[out] "real(dp), dimension(:, :), optional :: et_opti"       returns evapotranspiration time series for
  !>       \param[out] "real(dp), dimension(:, :), optional :: tws_opti"      returns tws time series
  !>       all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]

  !    HISTORY
  !>       \authors Juliane Mai, Rohini Kumar

  !>       \date Feb 2013

  ! Modifications:
  ! R. Kumar             Jun 2013 - restart_flag_states_read is passed to mhm call for the soil moisture initalisation
  ! R. Kumar             Jun 2013 - frac_sealed_city_area is added
  ! R. Kumar & S. Thober Aug 2013 - code change to incorporate output timestep during writing of the netcdf file
  ! R. Kumar             Aug 2013 - added iFlag_LAI_data_format to handle LAI options, and changed within the code made accordingly
  ! R. Kumar, J. Mai     Sep 2013 - Splitting allocation and initialization of arrays
  ! R. Kumar             Nov 2013 - update intent variables in documentation
  ! L. Samaniego         Nov 2013 - relational statements == to .eq., etc.
  ! M. Zink              Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  ! M. Zink              Mar 2014 - added inflow from upstream areas
  ! Stephan Thober       Jun 2014 - added chunk read for meteorological input
  ! Stephan Thober       Jun 2014 - updated flag for read_restart
  ! M. Cuntz & J. Mai    Nov 2014 - LAI input from daily, monthly or yearly files
  ! Matthias Zink        Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober       Aug 2015 - moved writing of daily discharge to mo_write_routing, included routing related variables from mRM
  ! David Schaefer       Aug 2015 - changed to new netcdf-writing scheme
  ! Stephan Thober       Sep 2015 - updated mrm_routing call
  ! O. Rakovec, R. Kumar Oct 2015 - added optional output for Domain averaged TWS
  ! Rohini Kumar         Mar 2016 - changes for handling multiple soil database options
  ! Stephan Thober       Nov 2016 - added two options for routing
  ! Rohini Kuamr         Dec 2016 - option to handle monthly mean gridded fields of LAI
  ! Stephan Thober       Jan 2017 - added prescribed weights for tavg and pet
  ! Zink M. Demirel C.   Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  ! Robert Schweppe      Dec 2017 - extracted call to mpr from inside mhm
  ! Robert Schweppe      Jun 2018 - refactoring and reformatting

  SUBROUTINE mhm_eval(parameterset, opti_domain_indices, runoff, smOptiSim, neutronsOptiSim, etOptiSim, twsOptiSim)

    use mo_optimization_types, only : optidata_sim
    use mo_common_datetime_type, only : datetimeinfo, nTstepDay, simPer, timeStep, landCoverPeriods, laiPeriods
    use mo_common_variables, only : mhmFileRestartIn, mrmFileRestartIn, global_parameters_name, &
                                            optimize, readPer, read_restart, &
                                            warmingDays, c2TSTu, level1, domainMeta, processMatrix
    use mo_global_variables, only : L1_Throughfall, L1_aETCanopy, L1_aETSealed, L1_aETSoil, &
                                    L1_absvappress, L1_baseflow, L1_fastRunoff, L1_infilSoil, L1_inter, L1_melt, &
                                    L1_netrad, L1_neutrons, L1_percol, L1_pet, L1_pet_calc, L1_pet_weights, L1_pre, &
                                    L1_preEffect, L1_pre_weights, L1_rain, L1_runoffSeal, L1_satSTW, L1_sealSTW, &
                                    L1_slowRunoff, L1_snow, L1_snowPack, L1_soilMoist, L1_temp, L1_temp_weights, L1_tmax, &
                                    L1_tmin, L1_total_runoff, L1_unsatSTW, L1_windspeed, evap_coeff, &
                                    fday_pet, fday_prec, fday_temp, fnight_pet, fnight_prec, fnight_temp, &
                                    nSoilHorizons_sm_input, &
                                    neutron_integral_AFast, outputFlxState, read_meteo_weights, &
                                    timeStep_model_inputs, timeStep_model_outputs, &
                                    L1_twsaObs, L1_etObs, L1_smObs, L1_neutronsObs, &
                                    L1_tann, L1_ssrd, L1_strd, fday_ssrd, fnight_ssrd, fday_strd, fnight_strd, & ! meteo for riv-temp
                                    L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_alpha, L1_degDay, &
                                    nSoilHorizons, soilHorizonBoundaries, lowestDepth, &
                                        L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, &
                                        L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, &
                                        L1_kfastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                        L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, &
                                        L1_wiltingPoint, are_parameter_initialized, L1_latitude
    use mo_init_states, only : variables_default_init
    use mo_julian, only : caldat, julday
    use mo_message, only : error_message
    use mo_string_utils, only : num2str
    use mo_meteo_forcings, only : prepare_meteo_forcings_data
    use mo_mhm, only : mhm
    use mo_restart, only : read_restart_states
    use mo_write_fluxes_states, only : OutputDataset
    use mo_constants, only : HourSecs
    use mo_common_variables, only : resolutionRouting
    use mo_common_variables, only : resolutionHydrology
    use mo_mrm_global_variables, only : InflowGauge, L11_C1, L11_C2, &
                                        L11_L1_Id, L11_TSrout, L11_fromN, L11_length, L11_nLinkFracFPimp, L11_nOutlets, &
                                        L11_netPerm, L11_qMod, L11_qOUT, L11_qTIN, L11_qTR, L11_slope, L11_toN, &
                                        L1_L11_Id, domain_mrm, level11, mRM_runoff, outputFlxState_mrm, &
                                        timeStep_model_outputs_mrm, gw_coupling, L0_river_head_mon_sum, &
                                        riv_temp_pcs
    use mo_mrm_mpr, only : mrm_update_param
    use mo_mrm_restart, only : mrm_read_restart_states
    use mo_mrm_routing, only : mrm_routing
    use mo_mrm_write, only : mrm_write_output_fluxes
    use mo_utils, only : ge
    use mo_mrm_river_head, only: calc_river_head, avg_and_write_timestep
    use mo_mhm_mpr_interface, only: call_mpr
#ifdef pgiFortran154
    use mo_write_fluxes_states, only : newOutputDataset
#endif

    implicit none

    ! a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), intent(in) :: parameterset

    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices
    ! returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: runoff

    ! returns soil moisture time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells,
    ! nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim

    ! dim1=ncells, dim2=time
    type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim

    ! returns evapotranspiration time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells,
    ! nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim

    ! returns tws time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells,
    ! nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim

    real(dp), dimension(:, :), allocatable :: L1_fNotSealed

    type(OutputDataset) :: nc

    ! Counters
    integer(i4) :: domainID, iDomain, tt, lcId, laiId

    ! No. of cells at level 1 for current Domain
    integer(i4) :: nCells

    ! start and end index at level 1 for current Domain
    integer(i4) :: s1, e1
    ! start and end index at level 1 for current Domain
    integer(i4) :: s1_param, e1_param

    ! meteorological time step for process 5 (PET)
    integer(i4), dimension(6) :: iMeteo_p5

    ! process 5: start and end index of vectors
    ! index 1: pet
    ! index 2: tmin
    ! index 3: tmax
    ! index 4: netrad
    ! index 5: absolute vapour pressure
    ! index 6: windspeed
    integer(i4), dimension(6) :: s_p5, e_p5

    integer(i4) :: s_meteo, e_meteo

    logical, dimension(:, :), pointer :: mask1

    integer(i4) :: iMeteoTS

    ! datetimeinfo variable for everything that has to do with time dependend
    ! calculations
    type(datetimeinfo) :: domainDateTime

    integer(i4) :: jj

    ! discharge timestep
    integer(i4) :: iDischargeTS

    ! start and end index at L11
    integer(i4) :: s11, e11

    ! factor between routing and hydrological modelling resolution
    real(dp) :: tsRoutFactor

    ! factor between routing and hydrological modelling resolution (dummy)
    real(dp) :: tsRoutFactorIn

    ! timestep of runoff to rout [h]
    ! - identical to timestep of input if
    ! tsRoutFactor is less than 1
    ! - tsRoutFactor * timestep if
    ! tsRoutFactor is greater than 1
    integer(i4) :: timestep_rout

    ! Runoff that is input for routing
    real(dp), allocatable, dimension(:) :: RunToRout

    real(dp), allocatable, dimension(:) :: mhmHorizons

    ! inflowing discharge
    real(dp), allocatable, dimension(:) :: InflowDischarge

    logical, pointer, dimension(:, :) :: mask11

    ! flag for performing routing
    logical :: do_rout

    integer(i4) :: gg

    ! number of domains simulated in this mhm_eval run. Depends on opti_function
    integer(i4) :: nDomains, ii


    if (optimize .and. present(opti_domain_indices)) then
      nDomains = size(opti_domain_indices)
    else
      nDomains = domainMeta%nDomains
    end if
    !----------------------------------------------------------
    ! Check optionals and initialize
    !----------------------------------------------------------
    if (present(runoff)) then
      do ii = 1, nDomains
        if (present(opti_domain_indices)) then
          iDomain = opti_domain_indices(ii)
        else
          iDomain = ii
        end if
        domainID = domainMeta%indices(iDomain)
        if (.not. domainMeta%doRouting(iDomain)) then
          call error_message("***ERROR: runoff for domain", trim(num2str(domainID)),&
                        "can not be produced, since routing process is off in Process Matrix")
        end if
      end do
    end if


    if (read_restart) then
      do ii = 1, nDomains
        if (optimize .and. present(opti_domain_indices)) then
          iDomain = opti_domain_indices(ii)
        else
          iDomain = ii
        end if
        ! this reads the eff. parameters and optionally the states and fluxes
        call read_restart_states(iDomain, mhmFileRestartIn(iDomain), do_read_dims_arg=.false.)
      end do
    end if
    if (.not. are_parameter_initialized) then
      !-------------------------------------------------------------------
      ! All variables had been allocated to the required
      ! space before this point (see, mo_startup: initialise) and initialised
      !-------------------------------------------------------------------
      ! as default values,
      ! all cells for all modeled Domains are simultenously initalized ONLY ONCE
      call variables_default_init()

      call call_mpr(parameterset, global_parameters_name, level1, .false., opti_domain_indices)
    end if

    allocate(L1_fNotSealed(size(L1_fSealed, 1), size(L1_fSealed, 2)))
    L1_fNotSealed = 1.0_dp - L1_fSealed
    !----------------------------------------
    ! loop over Domains
    !----------------------------------------
    DomainLoop: do ii = 1, nDomains
      if (optimize .and. present(opti_domain_indices)) then
        iDomain = opti_domain_indices(ii)
      else
        iDomain = ii
      end if

      ! evapotranspiration optimization
      if (present(etOptiSim)) call etOptiSim(iDomain)%init(L1_etObs(iDomain))
      ! total water storage optimization
      if (present(twsOptiSim)) call twsOptiSim(iDomain)%init(L1_twsaObs(iDomain))
      ! neutrons optimization
      if (present(neutronsOptiSim)) call neutronsOptiSim(iDomain)%init(L1_neutronsObs(iDomain))
      ! sm optimization
      if (present(smOptiSim)) call smOptiSim(iDomain)%init(L1_smObs(iDomain))

      ! get Domain information
      nCells = level1(iDomain)%nCells
      mask1 => level1(iDomain)%mask
      s1 = level1(iDomain)%iStart
      e1 = level1(iDomain)%iEnd
      ! TODO: at some point, parameters that share L0 might be put in a seperate structure
      ! for this MPR would be only called for each unique domainMeta%L0DataFrom(iDomain) and then also indexed
      s1_param = s1
      e1_param = e1

      ! this is done for correct handling of soil horizons in mHM subroutine, it is rather hacky, so it gets a TODO
      allocate(mhmHorizons(nSoilHorizons))
      mhmHorizons(:) = soilHorizonBoundaries(2: nSoilHorizons+1)
      mhmHorizons(nSoilHorizons) = lowestDepth

      if (domainMeta%doRouting(iDomain)) then
        ! ----------------------------------------
        ! initialize factor between routing resolution and hydrologic model resolution
        ! ----------------------------------------
        tsRoutFactor = 1_i4
        allocate(InflowDischarge(size(InflowGauge%Q, dim = 2)))
        InflowDischarge = 0._dp

        ! read states from restart
        if (read_restart) call mrm_read_restart_states(iDomain, mrmFileRestartIn(iDomain))
        !
        ! get Domain information at L11 if routing is activated
        s11 = level11(iDomain)%iStart
        e11 = level11(iDomain)%iEnd
        mask11 => level11(iDomain)%mask

        ! initialize routing parameters (has to be called for routing option 2)
        if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) &
            call mrm_update_param(iDomain, parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)))
        ! initialize variable for runoff for routing
        allocate(RunToRout(e1 - s1 + 1))
        RunToRout = 0._dp

        if ( riv_temp_pcs%active ) then
          ! set indices for current L11 domain
          riv_temp_pcs%s11 = s11
          riv_temp_pcs%e11 = e11
          ! allocate current L1 lateral components
          call riv_temp_pcs%alloc_lateral(nCells)
        end if
      end if

      ! init datetime variable
      call domainDateTime%init(iDomain)

      ! Loop over time
      TimeLoop: do tt = 1, domainDateTime%nTimeSteps
        ! time increment is done right after call to mrm (and initially before looping)
        if (timeStep_model_inputs(iDomain) .eq. 0_i4) then
          ! whole meteorology is already read

          ! set start and end of meteo position
          s_meteo = s1
          e_meteo = e1
          ! time step for meteorological variable (daily values)
          iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
        else
          ! read chunk of meteorological forcings data (reading, upscaling/downscaling)
          call prepare_meteo_forcings_data(iDomain, tt)
          ! set start and end of meteo position
          s_meteo = 1
          e_meteo = e1 - s1 + 1
          ! time step for meteorological variable (daily values)
          iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp)) &
                  - (readPer%julStart - simPer(iDomain)%julStart)
        end if

        ! preapare vector length specifications depending on the process case
        ! process 5 - PET
        select case (processMatrix(5, 1))
          !      [pet,        tmax,    tmin,  netrad, absVapP,windspeed]
        case(-1 : 0) ! PET is input
          s_p5 = [s_meteo, 1, 1, 1, 1, 1]
          e_p5 = [e_meteo, 1, 1, 1, 1, 1]
        case(1) ! Hargreaves-Samani
          s_p5 = [s_meteo, s_meteo, s_meteo, 1, 1, 1]
          e_p5 = [e_meteo, e_meteo, e_meteo, 1, 1, 1]
        case(2) ! Priestely-Taylor
          s_p5 = [s_meteo, 1, 1, s_meteo, 1, 1]
          e_p5 = [e_meteo, 1, 1, e_meteo, 1, 1]
        case(3) ! Penman-Monteith
          s_p5 = [s_meteo, 1, 1, s_meteo, s_meteo, s_meteo]
          e_p5 = [e_meteo, 1, 1, e_meteo, e_meteo, e_meteo]
        end select

        ! customize iMeteoTS for process 5 - PET
        select case (processMatrix(5, 1))
          !              [     pet,     tmin,     tmax,   netrad,  absVapP,windspeed ]
          case(-1 : 0) ! PET is input
          iMeteo_p5 = [iMeteoTS, 1, 1, 1, 1, 1 ]
          case(1) ! Hargreaves-Samani
          iMeteo_p5 = [iMeteoTS, iMeteoTS, iMeteoTS, 1, 1, 1 ]
          case(2) ! Priestely-Taylor
          iMeteo_p5 = [iMeteoTS, 1, 1, iMeteoTS, 1, 1 ]
          case(3) ! Penman-Monteith
            iMeteo_p5 = [iMeteoTS, 1, 1, iMeteoTS, iMeteoTS, iMeteoTS ]
        end select

        if (tt < domainDateTime%nTimeSteps) then
          call laiPeriods(iDomain)%increment(domainDateTime)
        end if
        lcId = landCoverPeriods(iDomain)%i
        laiId = laiPeriods(iDomain)%i

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
        call mhm(read_restart, & ! IN C
                tt, domainDateTime%newTime - 0.5_dp, processMatrix, &
                mhmHorizons, & ! IN C
                nCells, nSoilHorizons, real(nTstepDay, dp), c2TSTu,  & ! IN C
                neutron_integral_AFast, & ! IN C
                parameterset, & ! IN
                ! TODO MPR: reactivate this in version 6: pack(level1(iDomain)%y, level1(iDomain)%mask), & ! IN L1
                ! for backward compat and check cases, use the latitude from the extra parameter
                L1_latitude(s1 : e1), & ! IN L1
                evap_coeff, fday_prec, fnight_prec, fday_pet, fnight_pet, & ! IN F
                fday_temp, fnight_temp, & ! IN F
                L1_temp_weights(s1 : e1, :, :), & ! IN F
                L1_pet_weights(s1 : e1, :, :), & ! IN F
                L1_pre_weights(s1 : e1, :, :), & ! IN F
                read_meteo_weights, & ! IN F
                L1_pet(s_p5(1) : e_p5(1), iMeteo_p5(1)), & ! INOUT F:PET
                L1_tmin(s_p5(2) : e_p5(2), iMeteo_p5(2)), & ! IN F:PET
                L1_tmax(s_p5(3) : e_p5(3), iMeteo_p5(3)), & ! IN F:PET
                L1_netrad(s_p5(4) : e_p5(4), iMeteo_p5(4)), & ! IN F:PET
                L1_absvappress(s_p5(5) : e_p5(5), iMeteo_p5(5)), & ! IN F:PET
                L1_windspeed(s_p5(6) : e_p5(6), iMeteo_p5(6)), & ! IN F:PET
                L1_pre(s_meteo : e_meteo, iMeteoTS), & ! IN F:Pre
                L1_temp(s_meteo : e_meteo, iMeteoTS), & ! IN F:Temp
                L1_fSealed(s1_param : e1_param, lcId), & ! INOUT L1
                L1_inter(s1 : e1), L1_snowPack(s1 : e1), L1_sealSTW(s1 : e1), & ! INOUT S
                L1_soilMoist(s1 : e1, :), L1_unsatSTW(s1 : e1), L1_satSTW(s1 : e1), & ! INOUT S
                L1_neutrons(s1 : e1), & ! INOUT S
                L1_pet_calc(s1 : e1), & ! INOUT X
                L1_aETSoil(s1 : e1, :), L1_aETCanopy(s1 : e1), L1_aETSealed(s1 : e1), & ! INOUT X
                L1_baseflow(s1 : e1), L1_infilSoil(s1 : e1, :), L1_fastRunoff(s1 : e1), & ! INOUT X
                L1_melt(s1 : e1), L1_percol(s1 : e1), L1_preEffect(s1 : e1), L1_rain(s1 : e1), & ! INOUT X
                L1_runoffSeal(s1 : e1), L1_slowRunoff(s1 : e1), L1_snow(s1 : e1), & ! INOUT X
                L1_Throughfall(s1 : e1), L1_total_runoff(s1 : e1), & ! INOUT X
                ! TODO: MPR comment more distributed parameters
                ! L1_alpha(s1_param : e1_param, lcId), L1_degDayInc(s1_param : e1_param, lcId), &
                L1_alpha(s1_param : e1_param, 1), L1_degDayInc(s1_param : e1_param, lcId), &
                L1_degDayMax(s1_param : e1_param, lcId), & ! INOUT E1
                L1_degDayNoPre(s1_param : e1_param, lcId), L1_degDay(s1 : e1), & ! INOUT E1
                L1_fAsp(s1_param : e1_param), & ! INOUT E1
                L1_petLAIcorFactor(s1_param : e1_param, laiId, lcId), & ! INOUT E1
                L1_HarSamCoeff(s1_param : e1_param), & ! INOUT E1
                L1_PrieTayAlpha(s1_param : e1_param, laiId), & ! INOUT E1
                L1_aeroResist(s1_param : e1_param, laiId, lcId), & ! INOUT E1
                L1_surfResist(s1_param : e1_param, laiId), L1_fRoots(s1_param : e1_param, :, lcId), & ! INOUT E1
                L1_maxInter(s1_param : e1_param, laiId), L1_karstLoss(s1_param : e1_param), & ! INOUT E1
                ! L1_kFastFlow(s1_param : e1_param, lcId), L1_kSlowFlow(s1_param : e1_param, lcId), & ! INOUT E1
                L1_kFastFlow(s1_param : e1_param, lcId), L1_kSlowFlow(s1_param : e1_param, 1), & ! INOUT E1
                ! L1_kBaseFlow(s1_param : e1_param, lcId), L1_kPerco(s1_param : e1_param, lcId), & ! INOUT E1
                L1_kBaseFlow(s1_param : e1_param, 1), L1_kPerco(s1_param : e1_param, 1), & ! INOUT E1
                L1_soilMoistFC(s1_param : e1_param, :, lcId), & ! INOUT E1
                L1_soilMoistSat(s1_param : e1_param, :, lcId), & ! INOUT E1
                L1_soilMoistExp(s1_param : e1_param, :, lcId), L1_jarvis_thresh_c1(s1_param : e1_param), & ! INOUT E1
                ! L1_tempThresh(s1_param : e1_param, lcId), L1_unsatThresh(s1_param : e1_param, lcId), & ! INOUT E1
                L1_tempThresh(s1_param : e1_param, lcId), L1_unsatThresh(s1_param : e1_param, 1), & ! INOUT E1
                L1_sealedThresh(s1_param : e1_param), & ! INOUT E1
                L1_wiltingPoint(s1_param : e1_param, :, lcId)) ! INOUT E1

        ! call mRM routing
        if (domainMeta%doRouting(iDomain)) then
          ! set discharge timestep
          iDischargeTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
          ! set input variables for routing
          if (processMatrix(8, 1) .eq. 1) then
            ! >>>
            ! >>> original Muskingum routing, executed every time
            ! >>>
            do_rout = .True.
            tsRoutFactorIn = 1._dp
            timestep_rout = timestep
            RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TS-1] mm per timestep
            InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
            !
          else if ((processMatrix(8, 1) .eq. 2) .or. &
                   (processMatrix(8, 1) .eq. 3)) then
            ! >>>
            ! >>> adaptive timestep
            ! >>>
            do_rout = .False.
            ! calculate factor
            tsRoutFactor = L11_tsRout(iDomain) / (timestep * HourSecs)
            ! print *, 'routing factor: ', tsRoutFactor
            ! prepare routing call
            if (tsRoutFactor .lt. 1._dp) then
              ! ----------------------------------------------------------------
              ! routing timesteps are shorter than hydrologic time steps
              ! ----------------------------------------------------------------
              ! set all input variables
              tsRoutFactorIn = tsRoutFactor
              RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TS-1] mm per timestep
              InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
              timestep_rout = timestep
              do_rout = .True.
            else
              ! ----------------------------------------------------------------
              ! routing timesteps are longer than hydrologic time steps
              ! ----------------------------------------------------------------
              ! set all input variables
              tsRoutFactorIn = tsRoutFactor
              ! Runoff is accumulated in [mm]
              RunToRout = RunToRout + L1_total_runoff(s1 : e1)
              InflowDischarge = InflowDischarge + InflowGauge%Q(iDischargeTS, :)
              ! reset tsRoutFactorIn if last period did not cover full period
              if ((tt == domainDateTime%nTimeSteps) .and. (mod(tt, nint(tsRoutFactorIn)) /= 0_i4)) &
                      tsRoutFactorIn = mod(tt, nint(tsRoutFactorIn))
              if ((mod(tt, nint(tsRoutFactorIn)) .eq. 0_i4) .or. (tt .eq. domainDateTime%nTimeSteps)) then
                ! Inflow discharge is given as flow-rate and has to be converted to [m3]
                InflowDischarge = InflowDischarge / tsRoutFactorIn
                timestep_rout = timestep * nint(tsRoutFactorIn, i4)
                do_rout = .True.
              end if
            end if
          end if
          ! prepare temperature routing
          if ( riv_temp_pcs%active ) then
            ! init riv-temp from current air temp
            if ( tt .eq. 1_i4 ) call riv_temp_pcs%init_riv_temp( &
              domainDateTime%newTime - 0.5_dp, &
              real(nTstepDay, dp), &
              L1_temp(s_meteo : e_meteo, iMeteoTS), &
              read_meteo_weights, &
              L1_temp_weights(s1 : e1, :, :), &
              fday_temp, fnight_temp, &
              ! mapping info
              level1(iDomain)%CellArea * 1.E-6_dp, &
              L1_L11_Id(s1 : e1), &
              level11(iDomain)%CellArea * 1.E-6_dp, &
              L11_L1_Id(s11 : e11), &
              ! map_flag
              ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)) &
            )
            ! accumulate source Energy at L1 level
            call riv_temp_pcs%acc_source_E( &
              domainDateTime%newTime - 0.5_dp, &
              real(nTstepDay, dp), &
              L1_fSealed(s1_param : e1_param, lcId), &
              L1_fastRunoff(s1 : e1), &
              L1_slowRunoff(s1 : e1), &
              L1_baseflow(s1 : e1), &
              L1_runoffSeal(s1 : e1), &
              L1_temp(s_meteo : e_meteo, iMeteoTS), &
              L1_tann(s_meteo : e_meteo, iMeteoTS), &
              L1_ssrd(s_meteo : e_meteo, iMeteoTS), &
              L1_strd(s_meteo : e_meteo, iMeteoTS), &
              read_meteo_weights, &
              L1_temp_weights(s1 : e1, :, :), &
              fday_temp, fnight_temp, &
              fday_ssrd, fnight_ssrd, &
              fday_strd, fnight_strd  &
            )
            ! if routing should be performed, scale source energy to L11 level
            if ( do_rout ) call riv_temp_pcs%finalize_source_E( &
              level1(iDomain)%CellArea * 1.E-6_dp, &
              L1_L11_Id(s1 : e1), &
              level11(iDomain)%CellArea * 1.E-6_dp, &
              L11_L1_Id(s11 : e11), &
              timestep_rout, &
              ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)) &
            )
          end if
          ! -------------------------------------------------------------------
          ! execute routing
          ! -------------------------------------------------------------------
          if (do_rout) call mRM_routing(&
            ! general INPUT variables
            read_restart, &
            processMatrix(8, 1), & ! parse process Case to be used
            parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)), & ! routing par.
            RunToRout, & ! runoff [mm TS-1] mm per timestep old: L1_total_runoff_in(s1:e1, tt), &
            level1(iDomain)%CellArea * 1.E-6_dp, &
            L1_L11_Id(s1 : e1), &
            level11(iDomain)%CellArea * 1.E-6_dp, &
            L11_L1_Id(s11 : e11), &
            L11_netPerm(s11 : e11), & ! routing order at L11
            L11_fromN(s11 : e11), & ! link source at L11
            L11_toN(s11 : e11), & ! link target at L11
            L11_nOutlets(iDomain), & ! number of outlets
            timestep_rout, & ! timestep of runoff to rout [h]
            tsRoutFactorIn, & ! simulate timestep in [h]
            level11(iDomain)%nCells, & ! number of Nodes
            domain_mrm(iDomain)%nInflowGauges, &
            domain_mrm(iDomain)%InflowGaugeIndexList(:), &
            domain_mrm(iDomain)%InflowGaugeHeadwater(:), &
            domain_mrm(iDomain)%InflowGaugeNodeList(:), &
            InflowDischarge, &
            domain_mrm(iDomain)%nGauges, &
            domain_mrm(iDomain)%gaugeIndexList(:), &
            domain_mrm(iDomain)%gaugeNodeList(:), &
            ge(resolutionRouting(iDomain), resolutionHydrology(iDomain)), &
            ! original routing specific input variables
            L11_length(s11 : e11 - 1), & ! link length
            L11_slope(s11 : e11 - 1), &
            L11_nLinkFracFPimp(s11 : e11, lcId), & ! fraction of impervious layer at L11 scale
            ! general INPUT/OUTPUT variables
            L11_C1(s11 : e11), & ! first muskingum parameter
            L11_C2(s11 : e11), & ! second muskigum parameter
            L11_qOUT(s11 : e11), & ! routed runoff flowing out of L11 cell
            L11_qTIN(s11 : e11, :), & ! inflow water into the reach at L11
            L11_qTR(s11 : e11, :), & !
            L11_qMod(s11 : e11), &
            mRM_runoff(tt, :) &
          )
          ! -------------------------------------------------------------------
          ! groundwater coupling
          ! -------------------------------------------------------------------
          if (gw_coupling) then
              call calc_river_head(iDomain, L11_Qmod, L0_river_head_mon_sum)
              if (domainDateTime%is_new_month .and. tt > 1) then
                  call avg_and_write_timestep(iDomain, tt, L0_river_head_mon_sum)
              end if
          end if
          ! -------------------------------------------------------------------
          ! reset variables
          ! -------------------------------------------------------------------
          if (processMatrix(8, 1) .eq. 1) then
            ! reset Input variables
            InflowDischarge = 0._dp
            RunToRout = 0._dp
          else if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) then
            if ((.not. (tsRoutFactorIn .lt. 1._dp)) .and. do_rout) then
              do jj = 1, nint(tsRoutFactorIn) ! BUG: this should start at 2
                mRM_runoff(tt - jj + 1, :) = mRM_runoff(tt, :)
              end do
              ! reset Input variables
              InflowDischarge = 0._dp
              RunToRout = 0._dp
              ! reset lateral fluxes and time-step counter if routing was done
              if ( riv_temp_pcs%active ) call riv_temp_pcs%reset_timestep()
            end if
            ! if routing is done every time-step, reset river-temp time step
            if (tsRoutFactor .lt. 1._dp .and. riv_temp_pcs%active ) call riv_temp_pcs%reset_timestep()
          end if
        end if

        ! output only for evaluation period
        domainDateTime%tIndex_out = (tt - warmingDays(iDomain) * nTstepDay) ! tt if write out of warming period

        call domainDateTime%increment()

        if ( .not. optimize ) then
          if (any(outputFlxState_mrm) .AND. (processMatrix(8, 1) .NE. 0) ) then
            call mrm_write_output_fluxes( &
              iDomain, & ! Domain id
              level11(iDomain)%nCells, & ! nCells in Domain
              timeStep_model_outputs_mrm, & ! output specification
              domainDateTime, tt, timestep, & ! time specification
              mask11, & ! mask specification
              L11_qmod(s11 : e11) & ! output variables
            )
          end if

        if ((any(outputFlxState)) .and. (domainDateTime%tIndex_out > 0_i4)) then

          if (domainDateTime%tIndex_out == 1) then
#ifdef pgiFortran154
            nc = newOutputDataset(iDomain, mask1, level1(iDomain)%nCells)
#else
            nc = OutputDataset(iDomain, mask1, level1(iDomain)%nCells)
#endif
          end if

          call nc%updateDataset(&
            s1, e1, &
            s1_param, e1_param, &
            L1_fSealed(:, lcId), &
            L1_fNotSealed(:, lcId), &
            L1_inter, &
            L1_snowPack, &
            L1_soilMoist, &
            L1_soilMoistSat(:, :, lcId), &
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
            L1_preEffect &
          )

          ! write data
          if (domainDateTime%writeout(timeStep_model_outputs, tt)) then
            call nc%writeTimestep(domainDateTime%tIndex_out * timestep - 1)
          end if

          if(tt == domainDateTime%nTimeSteps) then
            call nc%close()
          end if

        end if
        end if ! <-- if (.not. optimize)

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
                         domainDateTime%is_new_day, domainDateTime%is_new_month, domainDateTime%is_new_year)
            ! last timestep is already done - write_counter exceeds size(smOptiSim(iDomain)%dataSim, dim=2)
            if (tt /= domainDateTime%nTimeSteps) then
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
            if (domainDateTime%is_new_day)   then
              call neutronsOptiSim(iDomain)%average()
            end if

            ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
            if (tt /= domainDateTime%nTimeSteps) then
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
                   domainDateTime%is_new_day, domainDateTime%is_new_month, domainDateTime%is_new_year)

            ! last timestep is already done - write_counter exceeds size(etOptiSim(iDomain)%dataSim, dim=2)
            if (tt /= domainDateTime%nTimeSteps) then
              ! aggregate evapotranspiration to needed time step for optimization
              call etOptiSim(iDomain)%add(sum(L1_aETSoil(s1 : e1, :), dim = 2) * &
                      L1_fNotSealed(s1_param : e1_param, lcId) + &
                      L1_aETCanopy(s1 : e1) + &
                      L1_aETSealed(s1 : e1) * L1_fSealed(s1_param : e1_param, lcId))
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
                               domainDateTime%is_new_day, domainDateTime%is_new_month, domainDateTime%is_new_year)

            ! last timestep is already done - write_counter exceeds size(twsOptiSim(iDomain)%dataSim, dim=2)
            if (tt /= domainDateTime%nTimeSteps) then
              ! aggregate evapotranspiration to needed time step for optimization
              call twsOptiSim(iDomain)%average_add(L1_inter(s1 : e1) + L1_snowPack(s1 : e1) + L1_sealSTW(s1 : e1) + &
                   L1_unsatSTW(s1 : e1) + L1_satSTW(s1 : e1))
              do gg = 1, nSoilHorizons
                call twsOptiSim(iDomain)%add(L1_soilMoist (s1 : e1, gg))
              end do
            end if
          end if
        end if

        ! TODO-RIV-TEMP: add OptiSim for river temperature

        if (tt < domainDateTime%nTimeSteps) then
          call landCoverPeriods(iDomain)%increment(domainDateTime)
        end if

      end do TimeLoop !<< TIME STEPS LOOP

      if (allocated(InflowDischarge)) deallocate(InflowDischarge)
       if (domainMeta%doRouting(iDomain)) then
        ! clean runoff variable
        deallocate(RunToRout)
        if ( riv_temp_pcs%active ) call riv_temp_pcs%dealloc_lateral()
      end if
      deallocate(mhmHorizons)

    end do DomainLoop !<< Domain LOOP

    ! =========================================================================
    ! SET RUNOFF OUTPUT VARIABLE
    ! =========================================================================
    if (present(runoff) .and. (processMatrix(8, 1) > 0)) runoff = mRM_runoff
    ! reset to false
    are_parameter_initialized = .false.

  end SUBROUTINE mhm_eval


END MODULE mo_mhm_eval
