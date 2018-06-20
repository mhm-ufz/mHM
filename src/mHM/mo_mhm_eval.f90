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
  !>       grid cells (of multiple basins concatenated),DIMENSION [nCells, nTimeSteps]
  !>       \param[out] "real(dp), dimension(:, :), optional :: basin_avg_tws" returns basin averaged total water storage
  !>       time series, DIMENSION [nTimeSteps, nBasins]
  !>       \param[out] "real(dp), dimension(:, :), optional :: neutrons_opti" dim1=ncells, dim2=time
  !>       \param[out] "real(dp), dimension(:, :), optional :: et_opti"       returns evapotranspiration time series for
  !>       all grid cells (of multiple basins concatenated),DIMENSION [nCells, nTimeSteps]

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
  ! O. Rakovec, R. Kumar Oct 2015 - added optional output for basin averaged TWS
  ! Rohini Kumar         Mar 2016 - changes for handling multiple soil database options
  ! Stephan Thober       Nov 2016 - added two options for routing
  ! Rohini Kuamr         Dec  2016 - option to handle monthly mean gridded fields of LAI
  ! Stephan Thober       Jan 2017 - added prescribed weights for tavg and pet
  ! Zink M. Demirel C.   Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  ! Robert Schweppe      Dec 2017 - extracted call to mpr from inside mhm
  ! Robert Schweppe      Jun 2018 - refactoring and reformatting

  SUBROUTINE mhm_eval(parameterset, runoff, sm_opti, basin_avg_tws, neutrons_opti, et_opti)

    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : LCyearId, dirRestartIn, nTstepDay, optimize, readPer, read_restart, simPer, timeStep, &
                                            warmingDays
    use mo_common_variables, only : level1, nBasins, processMatrix
    use mo_global_variables, only : L1_Throughfall, L1_aETCanopy, L1_aETSealed, L1_aETSoil, &
                                    L1_absvappress, L1_baseflow, L1_fastRunoff, L1_infilSoil, L1_inter, L1_melt, &
                                    L1_netrad, L1_neutrons, L1_percol, L1_pet, L1_pet_calc, L1_pet_weights, L1_pre, &
                                    L1_preEffect, L1_pre_weights, L1_rain, L1_runoffSeal, L1_satSTW, L1_sealSTW, &
                                    L1_slowRunoff, L1_snow, L1_snowPack, L1_soilMoist, L1_temp, L1_temp_weights, L1_tmax, &
                                    L1_tmin, L1_total_runoff, L1_unsatSTW, L1_windspeed, basin_avg_TWS_sim, evap_coeff, &
                                    fday_pet, fday_prec, fday_temp, fnight_pet, fnight_prec, fnight_temp, &
                                    nSoilHorizons_sm_input, nTimeSteps_L1_et, nTimeSteps_L1_neutrons, nTimeSteps_L1_sm, &
                                    neutron_integral_AFast, outputFlxState, read_meteo_weights, timeStep_et_input, &
                                    timeStep_model_inputs, timeStep_model_outputs, timeStep_sm_input
    use mo_init_states, only : variables_default_init
    use mo_julian, only : caldat, julday
    use mo_message, only : message
    use mo_meteo_forcings, only : prepare_meteo_forcings_data
    use mo_mhm, only : mhm
    use mo_mpr_eval, only : mpr_eval
    use mo_mpr_global_variables, only : HorizonDepth_mHM, &
                                        L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_alpha, L1_degDay, &
                                        L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, &
                                        L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, &
                                        L1_kfastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                        L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, &
                                        L1_wiltingPoint, nSoilHorizons_mHM, timeStep_LAI_input
    use mo_restart, only : read_restart_states
    use mo_write_fluxes_states, only : OutputDataset
#ifdef MRM2MHM
    use mo_common_constants, only : HourSecs
    use mo_common_mHM_mRM_variables, only : resolutionRouting
    use mo_common_variables, only : resolutionHydrology
    use mo_mrm_global_variables, only : InflowGauge, L11_C1, L11_C2, &
                                        L11_L1_Id, L11_TSrout, L11_fromN, L11_length, L11_nLinkFracFPimp, L11_nOutlets, &
                                        L11_netPerm, L11_qMod, L11_qOUT, L11_qTIN, L11_qTR, L11_slope, L11_toN, &
                                        L1_L11_Id, basin_mrm, level11, mRM_runoff, outputFlxState_mrm, &
                                        timeStep_model_outputs_mrm
    use mo_mrm_init, only : mrm_update_param, variables_default_init_routing
    use mo_mrm_restart, only : mrm_read_restart_states
    use mo_mrm_routing, only : mrm_routing
    use mo_mrm_write, only : mrm_write_output_fluxes
    use mo_utils, only : ge
#endif
#ifdef pgiFortran154
    use mo_write_fluxes_states, only : newOutputDataset
#endif

    implicit none

    ! a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), intent(in) :: parameterset

    ! returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: runoff

    ! returns soil moisture time series for all grid cells (of multiple basins concatenated),DIMENSION [nCells,
    ! nTimeSteps]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: sm_opti

    ! returns basin averaged total water storage time series, DIMENSION [nTimeSteps, nBasins]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: basin_avg_tws

    ! dim1=ncells, dim2=time
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: neutrons_opti

    ! returns evapotranspiration time series for all grid cells (of multiple basins concatenated),DIMENSION [nCells,
    ! nTimeSteps]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: et_opti

    ! for writing netcdf file
    integer(i4) :: tIndex_out

    real(dp), dimension(size(L1_fSealed, 1), size(L1_fSealed, 2), size(L1_fSealed, 3)) :: L1_fNotSealed

    type(OutputDataset) :: nc

    integer(i4) :: nTimeSteps

    ! Counters
    integer(i4) :: iBasin, tt

    ! No. of cells at level 1 for current basin
    integer(i4) :: nCells

    ! start and end index at level 1 for current basin
    integer(i4) :: s1, e1

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

    logical, dimension(:, :), allocatable :: mask1

    integer(i4) :: day, month, year, hour, prev_day, prev_month, prev_year

    integer(i4) :: iMeteoTS

    integer(i4) :: iLAI

    integer(i4) :: yId

    real(dp) :: newTime

    ! flags for stepping into new period
    logical :: is_new_day, is_new_month, is_new_year

    ! for averaging output
    integer(i4) :: average_counter

    ! if true write out netcdf files
    logical :: writeout

    ! write out time step
    integer(i4) :: writeout_counter

#ifdef MRM2MHM
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

    ! inflowing discharge
    real(dp), allocatable, dimension(:) :: InflowDischarge

    logical, allocatable, dimension(:, :) :: mask11

    ! flag for performing routing
    logical :: do_rout

#endif
    integer(i4) :: gg

    ! field of TWS
    real(dp), dimension(:), allocatable :: TWS_field

    real(dp) :: area_basin


    !----------------------------------------------------------
    ! Check optionals and initialize
    !----------------------------------------------------------
    if (present(runoff)) then
      if (processMatrix(8, 1) .eq. 0) then
        call message("***ERROR: runoff can not be produced, since routing process is off in Process Matrix")
        stop
      end if
    end if
    ! soil moisture optimization
    !--------------------------
    if (present(sm_opti)) then
      !                ! total No of cells, No of timesteps
      !                ! of all basins    , in soil moist input
      allocate(sm_opti(size(L1_pre, dim = 1), nTimeSteps_L1_sm))
      sm_opti(:, :) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! neutrons optimization
    !--------------------------
    if (present(neutrons_opti)) then
      !                ! total No of cells, No of timesteps
      !                ! of all basins    , in neutrons input
      allocate(neutrons_opti(size(L1_pre, dim = 1), nTimeSteps_L1_neutrons))
      neutrons_opti(:, :) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! evapotranspiration optimization
    !--------------------------
    if (present(et_opti)) then
      !                ! total No of cells, No of timesteps
      !                ! of all basins    , in evapotranspiration input
      allocate(et_opti(size(L1_pre, dim = 1), nTimeSteps_L1_et))
      et_opti(:, :) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! add other optionals...

    !-------------------------------------------------------------------
    ! All variables had been allocated to the required
    ! space before this point (see, mo_startup: initialise) and initialised
    !-------------------------------------------------------------------
    if (.NOT. read_restart) then
      ! as default values,
      ! all cells for all modeled basins are simultenously initalized ONLY ONCE
      call variables_default_init()
      call mpr_eval(parameterset)

#ifdef MRM2MHM
       if (processMatrix(8, 1) .gt. 0) then
        !-------------------------------------------
        ! L11 ROUTING STATE VARIABLES, FLUXES AND
        !             PARAMETERS
        !-------------------------------------------
        call variables_default_init_routing()
      end if
#endif
    else
      do iBasin = 1, nBasins
        ! this reads the eff. parameters and optionally the states and fluxes
        call read_restart_states(iBasin, dirRestartIn(iBasin))
      end do
    end if

#ifdef MRM2MHM
    if (processMatrix(8, 1) .gt. 0) then
      ! ----------------------------------------
      ! initialize factor between routing resolution and hydrologic model resolution
      ! ----------------------------------------
      tsRoutFactor = 1_i4
      allocate(InflowDischarge(size(InflowGauge%Q, dim = 2)))
      InflowDischarge = 0._dp
    end if
#endif

    L1_fNotSealed = 1.0_dp - L1_fSealed
    !----------------------------------------
    ! loop over basins
    !----------------------------------------
    do iBasin = 1, nBasins

      ! get basin information
      nCells = level1(iBasin)%nCells
      mask1 = level1(iBasin)%mask
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

#ifdef MRM2MHM
       if (processMatrix(8, 1) .gt. 0) then
        ! read states from restart
        if (read_restart) call mrm_read_restart_states(iBasin, dirRestartIn(iBasin))
        !
        ! get basin information at L11 and L110 if routing is activated
        s11 = level11(iBasin)%iStart
        e11 = level11(iBasin)%iEnd
        mask11 = level11(iBasin)%mask

        ! initialize routing parameters (has to be called for routing option 2)
        if (processMatrix(8, 1) .eq. 2) call mrm_update_param(iBasin, &
                parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)))
        ! initialize variable for runoff for routing
        allocate(RunToRout(e1 - s1 + 1))
        RunToRout = 0._dp
      end if
#endif

      ! allocate space for local tws field
      if (present(basin_avg_tws)) then
        allocate(TWS_field(s1 : e1))
        TWS_field(s1 : e1) = nodata_dp
      end if

      ! calculate NtimeSteps for this basin
      nTimeSteps = (simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1) * nTstepDay

      ! reinitialize time counter for LCover and MPR
      ! -0.5 is due to the fact that dec2date routine
      !   changes the day at 12:00 in NOON
      ! Whereas mHM needs day change at 00:00 h
      ! initialize the julian day as real
      newTime = real(simPer(iBasin)%julStart, dp)
      ! initialize the date
      call caldat(int(newTime), yy = year, mm = month, dd = day)
      ! initialize flags for period changes, they are true for first time step
      is_new_day = .true.
      is_new_month = .true.
      is_new_year = .true.

      ! initialize arrays and counters
      yId = LCyearId(year, iBasin)
      average_counter = 0
      writeout_counter = 0
      hour = -timestep
      iLAI = 0

      ! Loop over time
      do tt = 1, nTimeSteps
        ! time increment is done right after call to mrm (and initially before looping)
        if (timeStep_model_inputs(iBasin) .eq. 0_i4) then
          ! whole meteorology is already read

          ! set start and end of meteo position
          s_meteo = s1
          e_meteo = e1
          ! time step for meteorological variable (daily values)
          iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
        else
          ! read chunk of meteorological forcings data (reading, upscaling/downscaling)
          call prepare_meteo_forcings_data(iBasin, tt)
          ! set start and end of meteo position
          s_meteo = 1
          e_meteo = e1 - s1 + 1
          ! time step for meteorological variable (daily values)
          iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp)) &
                  - (readPer%julStart - simPer(iBasin)%julStart)
        end if

        ! preapare vector length specifications depending on the process case
        ! process 5 - PET
        select case (processMatrix(5, 1))
          !      (/pet,        tmax,    tmin,  netrad, absVapP,windspeed/)
        case(-1 : 0) ! PET is input
          s_p5 = (/s_meteo, 1, 1, 1, 1, 1/)
          e_p5 = (/e_meteo, 1, 1, 1, 1, 1/)
        case(1) ! Hargreaves-Samani
          s_p5 = (/s_meteo, s_meteo, s_meteo, 1, 1, 1/)
          e_p5 = (/e_meteo, e_meteo, e_meteo, 1, 1, 1/)
        case(2) ! Priestely-Taylor
          s_p5 = (/s_meteo, 1, 1, s_meteo, 1, 1/)
          e_p5 = (/e_meteo, 1, 1, e_meteo, 1, 1/)
        case(3) ! Penman-Monteith
          s_p5 = (/s_meteo, 1, 1, s_meteo, s_meteo, s_meteo/)
          e_p5 = (/e_meteo, 1, 1, e_meteo, e_meteo, e_meteo/)
        end select

        ! customize iMeteoTS for process 5 - PET
        select case (processMatrix(5, 1))
          !              (/     pet,     tmin,     tmax,   netrad,  absVapP,windspeed /)
        case(-1 : 0) ! PET is input
          iMeteo_p5 = (/iMeteoTS, 1, 1, 1, 1, 1 /)
        case(1) ! Hargreaves-Samani
          iMeteo_p5 = (/iMeteoTS, iMeteoTS, iMeteoTS, 1, 1, 1 /)
        case(2) ! Priestely-Taylor
          iMeteo_p5 = (/iMeteoTS, 1, 1, iMeteoTS, 1, 1 /)
        case(3) ! Penman-Monteith
          iMeteo_p5 = (/iMeteoTS, 1, 1, iMeteoTS, iMeteoTS, iMeteoTS /)
        end select

        select case (timeStep_LAI_input)
        case(0 : 1) ! long term mean monthly gridded fields or LUT-based values
          iLAI = month
        case(-1) ! daily timestep
          if (is_new_day) then
            iLAI = iLAI + 1
          end if
        case(-2) ! monthly timestep
          if (is_new_month) then
            iLAI = iLAI + 1
          end if
        case(-3) ! yearly timestep
          if (is_new_year) then
            iLAI = iLAI + 1
          end if
        end select

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
                tt, newTime - 0.5_dp, processMatrix, HorizonDepth_mHM, & ! IN C
                nCells, nSoilHorizons_mHM, real(nTstepDay, dp), & ! IN C
                neutron_integral_AFast, & ! IN C
                parameterset, & ! IN
                pack(level1(iBasin)%y, level1(iBasin)%mask), & ! IN L1
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
                L1_fSealed(s1 : e1, 1, yId), & ! INOUT L1
                L1_inter(s1 : e1), L1_snowPack(s1 : e1), L1_sealSTW(s1 : e1), & ! INOUT S
                L1_soilMoist(s1 : e1, :), L1_unsatSTW(s1 : e1), L1_satSTW(s1 : e1), & ! INOUT S
                L1_neutrons(s1 : e1), & ! INOUT S
                L1_pet_calc(s1 : e1), & ! INOUT X
                L1_aETSoil(s1 : e1, :), L1_aETCanopy(s1 : e1), L1_aETSealed(s1 : e1), & ! INOUT X
                L1_baseflow(s1 : e1), L1_infilSoil(s1 : e1, :), L1_fastRunoff(s1 : e1), & ! INOUT X
                L1_melt(s1 : e1), L1_percol(s1 : e1), L1_preEffect(s1 : e1), L1_rain(s1 : e1), & ! INOUT X
                L1_runoffSeal(s1 : e1), L1_slowRunoff(s1 : e1), L1_snow(s1 : e1), & ! INOUT X
                L1_Throughfall(s1 : e1), L1_total_runoff(s1 : e1), & ! INOUT X
                L1_alpha(s1 : e1, 1, 1), L1_degDayInc(s1 : e1, 1, yId), L1_degDayMax(s1 : e1, 1, yId), & ! INOUT E1
                L1_degDayNoPre(s1 : e1, 1, yId), L1_degDay(s1 : e1, 1, 1), L1_fAsp(s1 : e1, 1, 1), & ! INOUT E1
                L1_petLAIcorFactor(s1 : e1, iLAI, yId), L1_HarSamCoeff(s1 : e1, 1, 1), & ! INOUT E1
                L1_PrieTayAlpha(s1 : e1, iLAI, 1), L1_aeroResist(s1 : e1, iLAI, yId), & ! INOUT E1
                L1_surfResist(s1 : e1, iLAI, 1), L1_fRoots(s1 : e1, :, yId), & ! INOUT E1
                L1_maxInter(s1 : e1, iLAI, 1), L1_karstLoss(s1 : e1, 1, 1), & ! INOUT E1
                L1_kFastFlow(s1 : e1, 1, yId), L1_kSlowFlow(s1 : e1, 1, 1), & ! INOUT E1
                L1_kBaseFlow(s1 : e1, 1, 1), L1_kPerco(s1 : e1, 1, 1), & ! INOUT E1
                L1_soilMoistFC(s1 : e1, :, yId), L1_soilMoistSat(s1 : e1, :, yId), & ! INOUT E1
                L1_soilMoistExp(s1 : e1, :, yId), L1_jarvis_thresh_c1(s1 : e1, 1, 1), & ! INOUT E1
                L1_tempThresh(s1 : e1, 1, yId), L1_unsatThresh(s1 : e1, 1, 1), & ! INOUT E1
                L1_sealedThresh(s1 : e1, 1, 1), & ! INOUT E1
                L1_wiltingPoint(s1 : e1, :, yId)) ! INOUT E1

        ! call mRM routing
#ifdef MRM2MHM
        if (processMatrix(8, 1) .gt. 0) then
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
            RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TST-1] mm per timestep
            InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
            timestep_rout = timestep
            !
          else if (processMatrix(8, 1) .eq. 2) then
            ! >>>
            ! >>> adaptive timestep
            ! >>>
            do_rout = .False.
            ! calculate factor
            tsRoutFactor = L11_tsRout(iBasin) / (timestep * HourSecs)
            ! print *, 'routing factor: ', tsRoutFactor
            ! prepare routing call
            if (tsRoutFactor .lt. 1._dp) then
              ! ----------------------------------------------------------------
              ! routing timesteps are shorter than hydrologic time steps
              ! ----------------------------------------------------------------
              ! set all input variables
              tsRoutFactorIn = tsRoutFactor
              RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TST-1] mm per timestep
              InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
              timestep_rout = timestep
              do_rout = .True.
            else
              ! ----------------------------------------------------------------
              ! routing timesteps are longer than hydrologic time steps
              ! ----------------------------------------------------------------
              ! set all input variables
              tsRoutFactorIn = tsRoutFactor
              RunToRout = RunToRout + L1_total_runoff(s1 : e1)
              InflowDischarge = InflowDischarge + InflowGauge%Q(iDischargeTS, :)
              ! reset tsRoutFactorIn if last period did not cover full period
              if ((tt .eq. nTimeSteps) .and. (mod(tt, nint(tsRoutFactorIn)) .ne. 0_i4)) &
                      tsRoutFactorIn = mod(tt, nint(tsRoutFactorIn))
              if ((mod(tt, nint(tsRoutFactorIn)) .eq. 0_i4) .or. (tt .eq. nTimeSteps)) then
                InflowDischarge = InflowDischarge / tsRoutFactorIn
                timestep_rout = timestep * nint(tsRoutFactor, i4)
                do_rout = .True.
              end if
            end if
          end if
          ! -------------------------------------------------------------------
          ! execute routing
          ! -------------------------------------------------------------------
          if (do_rout) call mRM_routing(&
                  ! general INPUT variables
                  read_restart, &
                  processMatrix(8, 1), & ! parse process Case to be used
                  parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)), & ! routing par.
                  RunToRout, & ! runoff [mm TST-1] mm per timestep old: L1_total_runoff_in(s1:e1, tt), &
                  level1(iBasin)%CellArea * 1.E-6_dp, &
                  L1_L11_Id(s1 : e1), &
                  level11(iBasin)%CellArea * 1.E-6_dp, &
                  L11_L1_Id(s11 : e11), &
                  L11_netPerm(s11 : e11), & ! routing order at L11
                  L11_fromN(s11 : e11), & ! link source at L11
                  L11_toN(s11 : e11), & ! link target at L11
                  L11_nOutlets(iBasin), & ! number of outlets
                  timestep_rout, & ! timestep of runoff to rout [h]
                  tsRoutFactorIn, & ! simulate timestep in [h]
                  level11(iBasin)%nCells, & ! number of Nodes
                  basin_mrm(iBasin)%nInflowGauges, &
                  basin_mrm(iBasin)%InflowGaugeIndexList(:), &
                  basin_mrm(iBasin)%InflowGaugeHeadwater(:), &
                  basin_mrm(iBasin)%InflowGaugeNodeList(:), &
                  InflowDischarge, &
                  basin_mrm(iBasin)%nGauges, &
                  basin_mrm(iBasin)%gaugeIndexList(:), &
                  basin_mrm(iBasin)%gaugeNodeList(:), &
                  ge(resolutionRouting(iBasin), resolutionHydrology(iBasin)), &
                  ! original routing specific input variables
                  L11_length(s11 : e11 - 1), & ! link length
                  L11_slope(s11 : e11 - 1), &
                  L11_nLinkFracFPimp(s11 : e11, yID), & ! fraction of impervious layer at L11 scale
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
            ! reset variables
            ! -------------------------------------------------------------------
            if (processMatrix(8, 1) .eq. 1) then
              ! reset Input variables
              InflowDischarge = 0._dp
              RunToRout = 0._dp
            else if (processMatrix(8, 1) .eq. 2) then
              if ((.not. (tsRoutFactorIn .lt. 1._dp)) .and. do_rout) then
                do jj = 1, nint(tsRoutFactorIn)
                  mRM_runoff(tt - jj + 1, :) = mRM_runoff(tt, :)
                end do
                ! reset Input variables
                InflowDischarge = 0._dp
                RunToRout = 0._dp
              end if
            end if
          end if
#endif

        ! prepare the date and time information for next iteration step...
        ! TODO: turn datetime information into type with procedure "increment"
        ! set the current year as previous
        prev_day = day
        prev_month = month
        prev_year = year
        ! set the flags to false
        is_new_day = .false.
        is_new_month = .false.
        is_new_year = .false.

        ! increment of timestep
        hour = mod(hour + timestep, 24)
        newTime = julday(day, month, year) + real(hour + timestep, dp) / 24._dp
        ! calculate new year, month and day
        call caldat(int(newTime), yy = year, mm = month, dd = day)
        ! update the flags
        if (prev_day .ne. day) is_new_day = .true.
        if (prev_month .ne. month) is_new_month = .true.
        if (prev_year .ne. year) is_new_year = .true.

        if (.not. optimize) then
#ifdef MRM2MHM
          if (any(outputFlxState_mrm)) then
            call mrm_write_output_fluxes(&
                  ! basin id
                  iBasin, &
                  ! nCells in basin
                  level11(iBasin)%nCells, &
                  ! output specification
                  timeStep_model_outputs_mrm, &
                  ! time specification
                  warmingDays(iBasin), newTime, nTimeSteps, nTstepDay, &
                  tt, &
                  ! parse previous date to mRM writer
                  prev_day, prev_month, prev_year, &
                  timestep, &
                  ! mask specification
                  mask11, &
                  ! output variables
                  L11_qmod(s11 : e11))
        end if
#endif

        ! output only for evaluation period
        tIndex_out = (tt - warmingDays(iBasin) * nTstepDay) ! tt if write out of warming period

        if ((any(outputFlxState)) .and. (tIndex_out .gt. 0_i4)) then

          average_counter = average_counter + 1

          if (tIndex_out .EQ. 1) then
#ifdef pgiFortran154
            nc = newOutputDataset(iBasin, mask1, level1(iBasin)%nCells)
#else
            nc = OutputDataset(iBasin, mask1, level1(iBasin)%nCells)
#endif
          end if

          call nc%updateDataset(&
                  s1, &
                  e1, &
                  L1_fSealed(:, 1, yId), &
                  L1_fNotSealed(:, 1, yId), &
                  L1_inter, &
                  L1_snowPack, &
                  L1_soilMoist, &
                  L1_soilMoistSat(:, :, yId), &
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
                  L1_preEffect      &
                  )

          ! write data
          writeout = .false.
          if (timeStep_model_outputs .gt. 0) then
            if ((mod(tIndex_out, timeStep_model_outputs) .eq. 0) .or. (tt .eq. nTimeSteps)) writeout = .true.
          else
            select case(timeStep_model_outputs)
            case(0) ! only at last time step
              if (tt .eq. nTimeSteps) writeout = .true.
            case(-1) ! daily
              if (((tIndex_out .gt. 1) .and. is_new_day) .or. (tt .eq. nTimeSteps))     writeout = .true.
            case(-2) ! monthly
              if (((tIndex_out .gt. 1) .and. is_new_month) .or. (tt .eq. nTimeSteps)) writeout = .true.
            case(-3) ! yearly
              if (((tIndex_out .gt. 1) .and. is_new_year) .or. (tt .eq. nTimeSteps))   writeout = .true.
            case default ! no output at all

            end select
          end if

          if (writeout) then
            call nc%writeTimestep(tIndex_out * timestep - 1)
          end if

          if(tt .eq. nTimeSteps) then
            call nc%close()
          end if

        end if
        end if ! <-- if (.not. optimize)

        !----------------------------------------------------------------------
        ! FOR SOIL MOISTURE
        ! NOTE:: modeled soil moisture is averaged according to input time step
        !        soil moisture (timeStep_sm_input)
        !----------------------------------------------------------------------
        if (present(sm_opti)) then
          if (tt .EQ. 1) writeout_counter = 1
          ! only for evaluation period - ignore warming days
          if ((tt - warmingDays(iBasin) * nTstepDay) .GT. 0) then
            ! decide for daily, monthly or yearly aggregation
            select case(timeStep_sm_input)
            case(-1) ! daily
              if (is_new_day)   then
                sm_opti(s1 : e1, writeout_counter) = sm_opti(s1 : e1, writeout_counter) / real(average_counter, dp)
                writeout_counter = writeout_counter + 1
                average_counter = 0
              end if
            case(-2) ! monthly
              if (is_new_month) then
                sm_opti(s1 : e1, writeout_counter) = sm_opti(s1 : e1, writeout_counter) / real(average_counter, dp)
                writeout_counter = writeout_counter + 1
                average_counter = 0
              end if
            case(-3) ! yearly
              if (is_new_year)  then
                sm_opti(s1 : e1, writeout_counter) = sm_opti(s1 : e1, writeout_counter) / real(average_counter, dp)
                writeout_counter = writeout_counter + 1
                average_counter = 0
              end if
            end select

            ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
            if (.not. (tt .eq. nTimeSteps)) then
              ! aggregate soil moisture to needed time step for optimization
              sm_opti(s1 : e1, writeout_counter) = sm_opti(s1 : e1, writeout_counter) + &
                      sum(L1_soilMoist   (s1 : e1, 1 : nSoilHorizons_sm_input), dim = 2) / &
                              sum(L1_soilMoistSat(s1 : e1, 1 : nSoilHorizons_sm_input, yId), dim = 2)
            end if

            ! increase average counter by one
            average_counter = average_counter + 1
          end if
        end if

        !----------------------------------------------------------------------
        ! FOR TOTAL WATER STORAGE
        if(present(basin_avg_tws)) then
          area_basin = sum(level1(iBasin)%CellArea) * 1.E-6_dp
          TWS_field(s1 : e1) = L1_inter(s1 : e1) + L1_snowPack(s1 : e1) + L1_sealSTW(s1 : e1) + &
                  L1_unsatSTW(s1 : e1) + L1_satSTW(s1 : e1)
          do gg = 1, nSoilHorizons_mHM
            TWS_field(s1 : e1) = TWS_field(s1 : e1) + L1_soilMoist (s1 : e1, gg)
          end do
          basin_avg_TWS_sim(tt, iBasin) = (dot_product(TWS_field (s1 : e1), level1(iBasin)%CellArea * 1.E-6_dp) / area_basin)
        end if
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        ! FOR NEUTRONS
        ! NOTE:: modeled neutrons are averaged daily
        !----------------------------------------------------------------------
        if (present(neutrons_opti)) then
          if (tt .EQ. 1) writeout_counter = 1
          ! only for evaluation period - ignore warming days
          if ((tt - warmingDays(iBasin) * nTstepDay) .GT. 0) then
            ! decide for daily, monthly or yearly aggregation
            ! daily
            if (is_new_day)   then
              neutrons_opti(s1 : e1, writeout_counter) = neutrons_opti(s1 : e1, writeout_counter) / real(average_counter, dp)
              writeout_counter = writeout_counter + 1
              average_counter = 0
            end if

            ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
            if (.not. (tt .eq. nTimeSteps)) then
              ! aggregate neutrons to needed time step for optimization
              neutrons_opti(s1 : e1, writeout_counter) = neutrons_opti(s1 : e1, writeout_counter) + L1_neutrons(s1 : e1)
            end if

            average_counter = average_counter + 1
          end if
        end if

        !----------------------------------------------------------------------
        ! FOR EVAPOTRANSPIRATION
        ! NOTE:: modeled evapotranspiration is averaged according to input time step
        !        evapotranspiration (timeStep_sm_input)
        !----------------------------------------------------------------------
        if (present(et_opti)) then
          if (tt .EQ. 1) then
            writeout_counter = 1
          end if

          ! only for evaluation period - ignore warming days
          if ((tt - warmingDays(iBasin) * nTstepDay) .GT. 0) then
            ! decide for daily, monthly or yearly aggregation
            select case(timeStep_et_input)
            case(-1) ! daily
              if (is_new_day)   then
                writeout_counter = writeout_counter + 1
              end if
            case(-2) ! monthly
              if (is_new_month) then
                writeout_counter = writeout_counter + 1
              end if
            case(-3) ! yearly
              if (is_new_year)  then
                writeout_counter = writeout_counter + 1
              end if
            end select

            ! last timestep is already done - write_counter exceeds size(et_opti, dim=2)
            if (.not. (tt .eq. nTimeSteps)) then
              ! aggregate evapotranspiration to needed time step for optimization
              et_opti(s1 : e1, writeout_counter) = et_opti(s1 : e1, writeout_counter) + &
                      sum(L1_aETSoil(s1 : e1, :), dim = 2) * L1_fNotSealed(s1 : e1, 1, yId) + &
                      L1_aETCanopy(s1 : e1) + &
                      L1_aETSealed(s1 : e1) * L1_fSealed(s1 : e1, 1, yId)
            end if
          end if
        end if

        ! update the year-dependent yID (land cover id)
        if (is_new_year .and. tt .lt. nTimeSteps) then
          yId = LCyearId(year, iBasin)
        end if

      end do !<< TIME STEPS LOOP

      ! deallocate TWS field temporal variable
      if (allocated(TWS_field)) deallocate(TWS_field)
#ifdef MRM2MHM
       if (processMatrix(8, 1) .ne. 0) then
        ! clean runoff variable
        deallocate(RunToRout)
      end if
#endif

    end do !<< BASIN LOOP
#ifdef MRM2MHM
    ! =========================================================================
    ! SET RUNOFF OUTPUT VARIABLE
    ! =========================================================================
    if (present(runoff) .and. (processMatrix(8, 1) .gt. 0)) runoff = mRM_runoff
#endif

    ! =========================================================================
    ! SET TWS OUTPUT VARIABLE
    ! =========================================================================
    if(present(basin_avg_tws)) basin_avg_tws = basin_avg_TWS_sim

  end SUBROUTINE mhm_eval


END MODULE mo_mhm_eval
