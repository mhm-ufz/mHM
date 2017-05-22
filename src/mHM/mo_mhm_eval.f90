!> \file mo_mhm_eval.f90

!> \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

!> \details  Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

!> \authors Juliane Mai, Rohini Kumar
!> \date Feb 2013

MODULE mo_mhm_eval

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          mhm_eval

  !>        \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !>        \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !     INTENT(IN)
  !>       \param[in] "real(dp), dimension(:)    :: parameterset"
  !>          a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]

  !     INTENT(INOUT)
  !         None
  

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), dimension(:,:), optional  :: runoff"
  !>           returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: sm_opti"
  !>           returns soil moisture time series for all grid cells (of multiple basins concatenated), DIMENSION [nCells, nTimeSteps]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: basin_avg_tws"
  !>           returns basin averaged total water storage time series, DIMENSION [nTimeSteps, nBasins]

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai, Rohini Kumar
  !>        \date Feb 2013
  !         Modified, R. Kumar,             Jun 2013 - restart_flag_states_read is passed to mhm call
  !                                                    for the soil moisture initalisation
  !                   R. Kumar,             Jun 2013 - frac_sealed_city_area is added
  !                   R. Kumar & S. Thober, Aug 2013 - code change to incorporate output timestep
  !                                                    during writing of the netcdf file
  !                   R. Kumar,             Aug 2013 - added iFlag_LAI_data_format to handle LAI options,
  !                                                    and changed within the code made accordingly
  !                   R. Kumar, J. Mai,     Sep 2013 - Splitting allocation and initialization of arrays
  !                   R. Kumar              Nov 2013 - update intent variables in documentation
  !                   L. Samaniego,         Nov 2013 - relational statements == to .eq., etc.
  !                   M. Zink,              Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  !                   M. Zink,              Mar 2014 - added inflow from upstream areas
  !                   Stephan Thober,       Jun 2014 - added chunk read for meteorological input
  !                   Stephan Thober,       Jun 2014 - updated flag for read_restart
  !                   M. Cuntz & J. Mai,    Nov 2014 - LAI input from daily, monthly or yearly files
  !                   Matthias Zink,        Dec 2014 - adopted inflow gauges to ignore headwater cells
  !                   Stephan Thober,       Aug 2015 - moved writing of daily discharge to mo_write_routing,
  !                                                    included routing related variables from mRM
  !                   David Schaefer,       Aug 2015 - changed to new netcdf-writing scheme
  !                   Stephan Thober,       Sep 2015 - updated mrm_routing call
  !                   O. Rakovec, R. Kumar, Oct 2015 - added optional output for basin averaged TWS
  !                   Rohini Kumar,         Mar 2016 - changes for handling multiple soil database options
  !                   Stephan Thober,       Nov 2016 - added two options for routing
  !                   Rohini Kuamr,         Dec  2016 - option to handle monthly mean gridded fields of LAI
  !                   Stephan Thober,       Jan 2017 - added prescribed weights for tavg and pet
  !                   Zink M. Demirel C.,   Mar 2017 - Added Jarvis soil water stress function at SM process(3)  

  
  SUBROUTINE mhm_eval(parameterset, runoff, sm_opti, basin_avg_tws, neutrons_opti, et_opti)

    use mo_init_states,         only : get_basin_info
    use mo_init_states,         only : variables_default_init   ! default initalization of variables
    use mo_julian,              only : caldat, julday
    use mo_message,             only : message
    use mo_mhm,                 only : mhm

    use mo_mhm_constants,       only : nodata_dp
    use mo_restart,             only : read_restart_states      ! read initial values of variables

    use mo_meteo_forcings,      only : prepare_meteo_forcings_data
    use mo_write_fluxes_states, only : OutputDataset
#ifdef pgiFortran154
    use mo_write_fluxes_states, only : newOutputDataset
#endif
    use mo_global_variables,    only : &
         nTstepDay,                                          &
         timeStep_model_outputs, outputFlxState,             &  ! definition which output to write
         read_restart, perform_mpr, fracSealed_CityArea,     &
         timeStep_model_inputs,                              &
         timeStep, nBasins, simPer, readPer,                 & ! [h] simulation time step, No. of basins
         c2TSTu, HorizonDepth_mHM,            &
         nSoilHorizons_mHM, NTSTEPDAY, timeStep,             &
         LCyearId, LAIUnitList, LAILUT,                      &
         GeoUnitList, GeoUnitKar, soilDB,                    &
         iFlag_soilDB,                                       & ! options to handle different types of soil databases
         L0_Id, L0_soilId,                                   &
         L0_LCover, L0_asp, L0_LCover_LAI, L0_geoUnit,       &
         soilDB, L1_nTCells_L0,                              &
         L0_slope_emp,                                       &
         L1_upBound_L0, L1_downBound_L0, L1_leftBound_L0,    &
         L1_rightBound_L0, L0_latitude, L1_latitude,         &
         evap_coeff, fday_prec,                              &
         fnight_prec, fday_pet, fnight_pet, fday_temp,       &
         fnight_temp, L1_pet, L1_tmin, L1_tmax, L1_netrad,   &
         L1_temp_weights, L1_pet_weights, L1_pre_weights,    &
         read_meteo_weights,                                 &
         L1_absvappress, L1_windspeed,                       &
         L1_pre, L1_temp , L1_fForest,                       &
         L1_fPerm, L1_fSealed, L1_inter,                     &
         L1_snowPack, L1_sealSTW, L1_soilMoist, L1_unsatSTW, &
         L1_satSTW, L1_pet_calc,                             &
         L1_aETSoil, L1_aETCanopy, L1_aETSealed,             &
         L1_baseflow, L1_infilSoil, L1_fastRunoff, L1_melt,  &
         L1_percol, L1_preEffect, L1_rain, L1_runoffSeal,    &
         L1_slowRunoff, L1_snow, L1_Throughfall,             &
         L1_total_runoff, L1_alpha, L1_degDayInc,            &
         L1_degDayMax,                                       &
         L1_degDayNoPre, L1_degDay, L1_fAsp, L1_HarSamCoeff, &
         L1_PrieTayAlpha, L1_aeroResist, L1_surfResist,      &
         L1_fRoots, L1_maxInter, L1_karstLoss, L1_kfastFlow, &
         L1_kSlowFlow, L1_kBaseFlow, L1_kPerco,              &
         L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp,   &
         L1_jarvis_thresh_c1,                                &
         L1_tempThresh, L1_unsatThresh, L1_sealedThresh,     &
         L1_wiltingPoint, L1_neutrons,                       &
         basin_avg_TWS_sim,                                  &
         warmingDays,                                        &
         timeStep_LAI_input,                                 & ! flag on how LAI data has to be read
         L0_gridded_LAI, dirRestartIn,                       & ! restart directory location
         timeStep_sm_input,                                  & ! time step of soil moisture input (day, month, year)
         timeStep_et_input,                                  & ! time step of soil moisture input (day, month, year)
         nSoilHorizons_sm_input,                             & ! no. of mhm soil horizons equivalent to sm input
         nTimeSteps_L1_sm,                                   & ! total number of timesteps in soil moisture input
         nTimeSteps_L1_neutrons,                             & ! total number of timesteps in neutrons input
         nTimeSteps_L1_et                                      ! total number of timesteps in evapotranspiration input
    use mo_common_variables, only: &
         optimize, &
         processMatrix
#ifdef MRM2MHM
    use mo_utils, only: ge
    use mo_mrm_global_variables, only: &
         ! INPUT variables for mRM routing ====================================
         resolutionRouting,          &
         resolutionHydrology,        &
         L0_LCover_mRM,              & ! L0 land cover
         L0_floodPlain,              & ! flood plains at L0 level
         L0_areaCell,                &
         L1_areaCell,                &
         L1_L11_Id,                  &
         L11_areaCell,               &
         L11_L1_Id,                  &
         L11_aFloodPlain,            & ! flood plains at L11 level
         L11_length,                 & ! link length
         L11_slope,                  &
         L11_netPerm,                & ! routing order at L11
         L11_fromN,                  & ! link source at L11
         L11_toN,                    & ! link target at L11
         L11_nOutlets,               & ! number of outlets at the L11
         basin_mrm,                  & ! basin_mrm structure
         InflowGauge,                &
         outputFlxState_mrm,         &
         ! INPUT variables for writing output =================================
         warmingDays_mrm,            &
         timeStep_model_outputs_mrm, &
         ! INPUT/OUTPUT variables for mRM routing =============================
         L11_TSrout,                 & ! routing timestep in seconds
         L11_C1,                     & ! first muskingum parameter
         L11_C2,                     & ! second muskigum parameter
         L11_qOUT,                   & ! routed runoff flowing out of L11 cell
         L11_qTIN,                   & ! inflow water into the reach at L11
         L11_qTR,                    & !
         L11_FracFPimp,              & ! fraction of impervious layer at L11 scale
         L11_qMod,                   &
         mRM_runoff ! global variable containing runoff for every gauge
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_restart, only: mrm_read_restart_states
    use mo_mrm_routing, only: mrm_routing
    use mo_mrm_write, only: mrm_write_output_fluxes
    use mo_mrm_init, only: variables_default_init_routing, mrm_update_param
    use mo_mrm_constants, only: hourSecs
#endif

    implicit none

    real(dp), dimension(:),                          intent(in)  :: parameterset
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: runoff        ! dim1=time dim2=gauge
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: sm_opti       ! dim1=ncells, dim2=time
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: basin_avg_tws ! dim1=time dim2=nBasins
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: neutrons_opti ! dim1=ncells, dim2=time
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: et_opti       ! dim1=ncells, dim2=time

    ! -------------------------------------
    ! local variables
    !
    ! FOR WRITING GRIDDED STATES AND FLUXES
    integer(i4)                               :: tIndex_out          ! for writing netcdf file
    real(dp), dimension(size(L1_fSealed))     :: L1_fNotSealed
    type(OutputDataset)                       :: nc
    !
    ! counters and indexes
    integer(i4)                               :: nTimeSteps
    integer(i4)                               :: ii, tt, ll       ! Counters
    integer(i4)                               :: nCells           ! No. of cells at level 1 for current basin
    integer(i4)                               :: Lcover_yID       ! Land cover year ID
    integer(i4)                               :: s0, e0           ! start and end index at level 0 for current basin
    integer(i4)                               :: s1, e1           ! start and end index at level 1 for current basin
    !
    ! process case dependent length specifiers of vectors to pass to mHM
    integer(i4), dimension(6)                 :: iMeteo_p5        ! meteorological time step for process 5 (PET)
    integer(i4), dimension(6)                 :: s_p5, e_p5       ! process 5: start and end index of vectors
    !                                                             ! index 1: pet
    !                                                             ! index 2: tmin
    !                                                             ! index 3: tmax
    !                                                             ! index 4: netrad
    !                                                             ! index 5: absolute vapour pressure
    !                                                             ! index 6: windspeed
    integer(i4)                               :: s_meteo, e_meteo
    logical, dimension(:,:), allocatable      :: mask0, mask1
    integer(i4)                               :: nrows, ncols
    integer(i4)                               :: day, month, year, hour
    integer(i4)                               :: iMeteoTS
    integer(i4)                               :: iGridLAI_TS
    integer(i4)                               :: yId
    real(dp)                                  :: newTime
    integer(i4)                               :: year_counter     ! for yearly output
    integer(i4)                               :: average_counter  ! for averaging output
    logical                                   :: writeout         ! if true write out netcdf files
    integer(i4)                               :: writeout_counter ! write out time step
    !
#ifdef MRM2MHM
    ! for routing
    logical               :: do_mpr
    integer(i4)           :: jj
    integer(i4)           :: iDischargeTS ! discharge timestep
    integer(i4)           :: s11, e11 ! start and end index at L11
    integer(i4)           :: s110, e110 ! start and end index of L11 at L0
    real(dp)              :: tsRoutFactor ! factor between routing and hydrological modelling resolution
    real(dp)              :: tsRoutFactorIn ! factor between routing and hydrological modelling resolution (dummy)
    integer(i4) :: timestep_rout            ! timestep of runoff to rout [h]
    !                                       ! - identical to timestep of input if
    !                                       !   tsRoutFactor is less than 1
    !                                       ! - tsRoutFactor * timestep if
    !                                       !   tsRoutFactor is greater than 1
    !

    real(dp), allocatable :: RunToRout(:) ! Runoff that is input for routing
    real(dp), allocatable :: InflowDischarge(:) ! inflowing discharge
    logical,  allocatable :: mask11(:,:)
    logical               :: do_rout ! flag for performing routing
#endif
    !

    ! for basin average tws timeseries
    integer(i4)                               :: gg
    real(dp), dimension(:), allocatable       :: TWS_field      ! field of TWS
    real(dp)                                  :: area_basin

    ! LAI options
    integer(i4)                               :: day_counter
    integer(i4)                               :: month_counter
    real(dp), dimension(:), allocatable       :: LAI            ! local variable for leaf area index

    !----------------------------------------------------------
    ! Check optionals and initialize
    !----------------------------------------------------------
    if ( present(runoff) ) then
       if ( processMatrix(8, 1) .eq. 0 ) then
          call message("***ERROR: runoff can not be produced, since routing process is off in Process Matrix")
          stop
       end if
    end if
    ! soil mosiure optimization
    !--------------------------
    if ( present(sm_opti) ) then
       !                ! total No of cells, No of timesteps
       !                ! of all basins    , in soil moist input
       allocate(sm_opti(size(L1_pre, dim=1), nTimeSteps_L1_sm))
       sm_opti(:,:) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! neutrons optimization
    !--------------------------
    if ( present(neutrons_opti) ) then
       !                ! total No of cells, No of timesteps
       !                ! of all basins    , in neutrons input
       allocate(neutrons_opti(size(L1_pre, dim=1), nTimeSteps_L1_neutrons))
       neutrons_opti(:,:) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! evapotranspiration optimization
    !--------------------------
    if ( present(et_opti) ) then
       !                ! total No of cells, No of timesteps
       !                ! of all basins    , in evapotranspiration input
       allocate(et_opti(size(L1_pre, dim=1), nTimeSteps_L1_et))
       et_opti(:,:) = 0.0_dp ! has to be intialized with zero because later summation
    end if
    ! add other optionals...

    !-------------------------------------------------------------------
    ! Initalize State variables either to the default value or
    ! from the restart_files.
    ! All variables to be initalized had been allocated to the required
    ! space before this point (see, mo_startup: initialise)
    !-------------------------------------------------------------------
    if (.NOT. read_restart ) then
       ! as default values,
       ! all cells for all modeled basins are simultenously initalized ONLY ONCE
       call variables_default_init()
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
       ! read from restart files, basin wise ...
       do ii = 1, nBasins
          call read_restart_states(ii, dirRestartIn(ii) )
       end do
    end if

#ifdef MRM2MHM
    if (processMatrix(8, 1) .gt. 0) then
       ! ----------------------------------------
       ! initialize factor between routing resolution and hydrologic model resolution
       ! ----------------------------------------
       tsRoutFactor = 1_i4
       allocate(InflowDischarge(size(InflowGauge%Q, dim=2)))
       InflowDischarge = 0._dp
    end if
#endif

    !----------------------------------------
    ! loop over basins
    !----------------------------------------
    do ii = 1, nBasins

       ! calculate NtimeSteps for this basin
       nTimeSteps = ( simPer(ii)%julEnd - simPer(ii)%julStart + 1 ) * NTSTEPDAY

       ! reinitialize time counter for LCover and MPR
       ! -0.5 is due to the fact that dec2date routine
       !   changes the day at 12:00 in NOON
       ! Whereas mHM needs day change at 00:00 h
       newTime = real(simPer(ii)%julStart,dp)
       yId     = 0

       ! get basin information
       call get_basin_info ( ii,  0, nrows, ncols,                iStart=s0,  iEnd=e0, mask=mask0 )
       call get_basin_info ( ii,  1, nrows, ncols, ncells=nCells, iStart=s1,  iEnd=e1, mask=mask1 )

#ifdef MRM2MHM
       if (processMatrix(8, 1) .gt. 0) then
          ! read states from restart
          if (read_restart) call mrm_read_restart_states(ii, dirRestartIn(ii))
          !
          ! get basin information at L11 and L110 if routing is activated
          call get_basin_info_mrm ( ii,  11, nrows, ncols,  iStart=s11,  iEnd=e11, mask=mask11  )
          call get_basin_info_mrm ( ii, 110, nrows, ncols, iStart=s110,  iEnd=e110 )
          ! initialize routing parameters (has to be called for routing option 2)
          if (processMatrix(8, 1) .eq. 2) call mrm_update_param(ii, &
               parameterset(processMatrix(8,3) - processMatrix(8,2) + 1:processMatrix(8,3)))
          ! initialize variable for runoff for routing
          allocate(RunToRout(e1 - s1 + 1))
          RunToRout = 0._dp
       end if
#endif
       ! allocate space for local LAI grid
       allocate( LAI(s0:e0) )
       LAI(:) = nodata_dp

       ! allocate space for local tws field
       if (present(basin_avg_tws)) then
          allocate(TWS_field(s1:e1))
          TWS_field(s1:e1) = nodata_dp
       end if

       ! Loop over time
       average_counter  = 0
       writeout_counter = 0
       hour = -timestep
       iGridLAI_TS = 0
       do tt = 1, nTimeSteps
          if ( timeStep_model_inputs(ii) .eq. 0_i4 ) then
             ! whole meteorology is already read

             ! set start and end of meteo position
             s_meteo = s1
             e_meteo = e1
             ! time step for meteorological variable (daily values)
             iMeteoTS = ceiling( real(tt,dp) / real(NTSTEPDAY,dp) )
          else
             ! read chunk of meteorological forcings data (reading, upscaling/downscaling)
             call prepare_meteo_forcings_data(ii, tt)
             ! set start and end of meteo position
             s_meteo = 1
             e_meteo = e1 - s1 + 1
             ! time step for meteorological variable (daily values)
             iMeteoTS = ceiling( real(tt,dp) / real(NTSTEPDAY,dp) ) &
                        - ( readPer%julStart - simPer(ii)%julStart )
          end if

          hour = mod(hour+timestep, 24)

          ! year needed to be passed to mHM call, e.g., yearly LC scene
          ! month needed for LAI process
          call caldat(int(newTime), yy=year, mm=month, dd=day)

          ! preapare vector length specifications depending on the process case
          ! process 5 - PET
          select case (processMatrix(5,1))
             !      (/pet,        tmax,    tmin,  netrad, absVapP,windspeed/)
          case(0) ! PET is input
             s_p5 = (/s_meteo,       1,       1,       1,       1,       1/)
             e_p5 = (/e_meteo,       1,       1,       1,       1,       1/)
          case(1) ! Hargreaves-Samani
             s_p5 = (/s_meteo, s_meteo, s_meteo,       1,       1,       1/)
             e_p5 = (/e_meteo, e_meteo, e_meteo,       1,       1,       1/)
          case(2) ! Priestely-Taylor
             s_p5 = (/s_meteo,       1,       1, s_meteo,       1,       1/)
             e_p5 = (/e_meteo,       1,       1, e_meteo,       1,       1/)
          case(3) ! Penman-Monteith
             s_p5 = (/s_meteo,       1,       1, s_meteo, s_meteo, s_meteo/)
             e_p5 = (/e_meteo,       1,       1, e_meteo, e_meteo, e_meteo/)
          end select

          ! customize iMeteoTS for process 5 - PET
          select case (processMatrix(5,1))
          !              (/     pet,     tmin,     tmax,   netrad,  absVapP,windspeed /)
          case(0) ! PET is input
             iMeteo_p5 = (/iMeteoTS,        1,        1,        1,        1,        1 /)
          case(1) ! Hargreaves-Samani
             iMeteo_p5 = (/iMeteoTS, iMeteoTS, iMeteoTS,        1,        1,        1 /)
          case(2) ! Priestely-Taylor
             iMeteo_p5 = (/iMeteoTS,        1,        1, iMeteoTS,        1,        1 /)
          case(3) ! Penman-Monteith
             iMeteo_p5 = (/iMeteoTS,        1,        1, iMeteoTS, iMeteoTS, iMeteoTS /)
          end select

          !--------------------------------------------------------------------
          ! call LAI function to get LAI fields for this timestep and basin
          !--------------------------------------------------------------------
          ! initalise counters
          if ( tt .EQ. 1 ) then
             day_counter   = day
             month_counter = month
             year_counter  = year
          end if
          !
          select case(timeStep_LAI_input)
          case(0)
             ! create gridded fields of LAI using long term monthly mean LAI values
             ! and the corresponding LC file
             ! update LAI --> for 1st timestep and when month changes
             if( (tt .EQ. 1) .OR. (month_counter .NE. month) ) then
                do ll = 1, size(LAIUnitList)
                   where( L0_LCover_LAI(s0:e0) .EQ. LAIUnitList(ll) ) LAI(:) = LAILUT(ll, month)
                end do
             end if
             
          case(1) ! long term mean monthly gridded lai 
             LAI(:) = L0_gridded_LAI(s0:e0, month)
             
          case(-1) ! daily
             if ( (tt .EQ. 1) .OR. (day .NE. day_counter) ) then
                iGridLAI_TS = iGridLAI_TS + 1_i4
                LAI(:) = L0_gridded_LAI(s0:e0, iGridLAI_TS)
             endif
          case(-2) ! monthly
             if ( (tt .EQ. 1) .OR. (month .NE. month_counter) ) then
                iGridLAI_TS = iGridLAI_TS + 1_i4
                LAI(:) = L0_gridded_LAI(s0:e0, iGridLAI_TS)
             endif
          case(-3) ! yearly
             if ( (tt .EQ. 1) .OR. (year .NE. year_counter) ) then
                iGridLAI_TS = iGridLAI_TS + 1_i4
                LAI(:) = L0_gridded_LAI(s0:e0, iGridLAI_TS)
             endif
          case default ! no output at all
             continue
          end select
          !
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
          call mhm(perform_mpr, read_restart, fracSealed_cityArea,                          & ! IN C
               timeStep_LAI_input, year_counter, month_counter, day_counter,                & ! IN C
               tt, newTime-0.5_dp, processMatrix, c2TSTu, HorizonDepth_mHM,                 & ! IN C
               nCells, nSoilHorizons_mHM, real(NTSTEPDAY,dp), mask0,                        & ! IN C
               iFlag_soilDB,                                                                & ! IN C
               parameterset,                                                                & ! IN P
               LCyearId(year,ii), GeoUnitList, GeoUnitKar, LAIUnitList, LAILUT,             & ! IN L0
               L0_slope_emp(s0:e0), L0_Latitude(s0:e0),                                     & ! IN L0
               L0_Id(s0:e0), L0_soilId(s0:e0,:), L0_LCover_LAI(s0:e0),                      & ! IN L0
               L0_LCover(s0:e0, LCyearId(year,ii)), L0_asp(s0:e0), LAI(s0:e0),              & ! IN L0
               L0_geoUnit(s0:e0),                                                           & ! IN L0
               soilDB%is_present, soilDB%nHorizons, soilDB%nTillHorizons,                   & ! IN L0
               soilDB%sand, soilDB%clay, soilDB%DbM, soilDB%Wd, soilDB%RZdepth,             & ! IN L0
               L1_nTCells_L0(s1:e1),                                                        & ! IN L1
               L1_upBound_L0(s1:e1), L1_downBound_L0(s1:e1),                                & ! IN L1
               L1_leftBound_L0(s1:e1), L1_rightBound_L0(s1:e1),                             & ! IN L1
               L1_latitude(s_p5(1):e_p5(1)),                                                & ! IN L1
               evap_coeff, fday_prec, fnight_prec, fday_pet, fnight_pet,                    & ! IN F
               fday_temp, fnight_temp,                                                      & ! IN F
               L1_temp_weights(s1:e1,:,:),                                                  & ! IN F
               L1_pet_weights(s1:e1,:,:),                                                   & ! IN F
               L1_pre_weights(s1:e1,:,:),                                                   & ! IN F
               read_meteo_weights,                                                          & ! IN F
               L1_pet(s_p5(1):e_p5(1), iMeteo_p5(1)),                                       & ! INOUT F:PET
               L1_tmin(s_p5(2):e_p5(2), iMeteo_p5(2)),                                      & ! IN F:PET
               L1_tmax(s_p5(3):e_p5(3), iMeteo_p5(3)),                                      & ! IN F:PET
               L1_netrad(s_p5(4):e_p5(4), iMeteo_p5(4)),                                    & ! IN F:PET
               L1_absvappress(s_p5(5):e_p5(5), iMeteo_p5(5)),                               & ! IN F:PET
               L1_windspeed(s_p5(6):e_p5(6), iMeteo_p5(6)),                                 & ! IN F:PET
               L1_pre(s_meteo:e_meteo,iMeteoTS),                                            & ! IN F:Pre
               L1_temp(s_meteo:e_meteo,iMeteoTS),                                           & ! IN F:Temp
               yId,                                                                         & ! INOUT C
               L1_fForest(s1:e1), L1_fPerm(s1:e1),  L1_fSealed(s1:e1),                      & ! INOUT L1
               L1_inter(s1:e1), L1_snowPack(s1:e1), L1_sealSTW(s1:e1),                      & ! INOUT S
               L1_soilMoist(s1:e1,:), L1_unsatSTW(s1:e1), L1_satSTW(s1:e1),                 & ! INOUT S
               L1_neutrons(s1:e1),                                                          & ! INOUT S
               L1_pet_calc(s1:e1),                                                          & ! INOUT X
               L1_aETSoil(s1:e1,:), L1_aETCanopy(s1:e1), L1_aETSealed(s1:e1),               & ! INOUT X
               L1_baseflow(s1:e1), L1_infilSoil(s1:e1,:), L1_fastRunoff(s1:e1),             & ! INOUT X
               L1_melt(s1:e1), L1_percol(s1:e1), L1_preEffect(s1:e1), L1_rain(s1:e1),       & ! INOUT X
               L1_runoffSeal(s1:e1), L1_slowRunoff(s1:e1), L1_snow(s1:e1),                  & ! INOUT X
               L1_Throughfall(s1:e1), L1_total_runoff(s1:e1),                               & ! INOUT X
               L1_alpha(s1:e1), L1_degDayInc(s1:e1), L1_degDayMax(s1:e1),                   & ! INOUT E1
               L1_degDayNoPre(s1:e1), L1_degDay(s1:e1), L1_fAsp(s1:e1),                     & ! INOUT E1
               L1_HarSamCoeff(s1:e1), L1_PrieTayAlpha(s1:e1,:), L1_aeroResist(s1:e1,:),     & ! INOUT E1
               L1_surfResist(s1:e1,:), L1_fRoots(s1:e1,:),                                  & ! INOUT E1
               L1_maxInter(s1:e1), L1_karstLoss(s1:e1),  L1_kFastFlow(s1:e1),               & ! INOUT E1
               L1_kSlowFlow(s1:e1), L1_kBaseFlow(s1:e1), L1_kPerco(s1:e1),                  & ! INOUT E1
               L1_soilMoistFC(s1:e1,:), L1_soilMoistSat(s1:e1,:), L1_soilMoistExp(s1:e1,:), & ! INOUT E1
               L1_jarvis_thresh_c1(s1:e1),                                                  & ! INOUT E1
               L1_tempThresh(s1:e1), L1_unsatThresh(s1:e1), L1_sealedThresh(s1:e1),         & ! INOUT E1
               L1_wiltingPoint(s1:e1,:)                                                     ) ! INOUT E1

          ! call mRM routing
#ifdef MRM2MHM
          if (processMatrix(8, 1) .gt. 0) then
             ! set discharge timestep
             iDischargeTS = ceiling(real(tt,dp) / real(NTSTEPDAY,dp))
             ! determine whether mpr is to be executed
             if( ( LCyearId(year,ii) .NE. yId) .or. (tt .EQ. 1) ) then
                do_mpr = perform_mpr
                Lcover_yID = LCyearId(year, ii)
             else
                do_mpr = .false.
             end if
             !
             ! set input variables for routing
             if (processMatrix(8, 1) .eq. 1) then
                ! >>>
                ! >>> original Muskingum routing, executed every time
                ! >>>
                do_rout         = .True.
                tsRoutFactorIn  = 1._dp
                timestep_rout   = timestep
                RunToRout       = L1_total_runoff(s1:e1) ! runoff [mm TST-1] mm per timestep
                InflowDischarge = InflowGauge%Q(iDischargeTS,:) ! inflow discharge in [m3 s-1]
                timestep_rout   = timestep
                !
             else if (processMatrix(8, 1) .eq. 2) then
                ! >>>
                ! >>> adaptive timestep
                ! >>>
                do_rout = .False.
                ! set dummy lcover_yID
                Lcover_yID = 1_i4
                ! calculate factor
                tsRoutFactor = L11_tsRout(ii) / (timestep * HourSecs)
                ! print *, 'routing factor: ', tsRoutFactor
                ! prepare routing call
                if (tsRoutFactor .lt. 1._dp) then
                   ! ----------------------------------------------------------------
                   ! routing timesteps are shorter than hydrologic time steps
                   ! ----------------------------------------------------------------
                   ! set all input variables
                   tsRoutFactorIn = tsRoutFactor
                   RunToRout = L1_total_runoff(s1:e1) ! runoff [mm TST-1] mm per timestep
                   InflowDischarge = InflowGauge%Q(iDischargeTS,:) ! inflow discharge in [m3 s-1]
                   timestep_rout = timestep
                   do_rout = .True.
                else
                   ! ----------------------------------------------------------------
                   ! routing timesteps are longer than hydrologic time steps
                   ! ----------------------------------------------------------------
                   ! set all input variables
                   tsRoutFactorIn = tsRoutFactor
                   RunToRout = RunToRout + L1_total_runoff(s1:e1)
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
             if (do_rout) call mRM_routing( &
                  ! general INPUT variables
                  processMatrix(8, 1), & ! parse process Case to be used
                  parameterset(processMatrix(8, 3) - processMatrix(8, 2) + 1 : processMatrix(8, 3)), & ! routing par.
                  RunToRout, & ! runoff [mm TST-1] mm per timestep old: L1_total_runoff_in(s1:e1, tt), &
                  L1_areaCell(s1:e1), &
                  L1_L11_Id(s1:e1), &
                  L11_areaCell(s11:e11), &
                  L11_L1_Id(s11:e11), &
                  L11_netPerm(s11:e11), & ! routing order at L11
                  L11_fromN(s11:e11), & ! link source at L11
                  L11_toN(s11:e11), & ! link target at L11
                  L11_nOutlets(ii), & ! number of outlets
                  timestep_rout, & ! timestep of runoff to rout [h]
                  tsRoutFactorIn, & ! simulate timestep in [h]
                  basin_mrm%L11_iEnd(ii) - basin_mrm%L11_iStart(ii) + 1, & ! number of Nodes
                  basin_mrm%nInflowGauges(ii), &
                  basin_mrm%InflowGaugeIndexList(ii,:), &
                  basin_mrm%InflowGaugeHeadwater(ii,:), &
                  basin_mrm%InflowGaugeNodeList(ii,:), &
                  InflowDischarge, &
                  basin_mrm%nGauges(ii), &
                  basin_mrm%gaugeIndexList(ii,:), &
                  basin_mrm%gaugeNodeList(ii,:), &
                  ge(resolutionRouting(ii), resolutionHydrology(ii)), &
                  ! original routing specific input variables
                  L0_LCover_mRM(s0:e0, Lcover_yID), & ! L0 land cover
                  L0_floodPlain(s110:e110), & ! flood plains at L0 level
                  L0_areaCell(s0:e0), &
                  L11_aFloodPlain(s11:e11), & ! flood plains at L11 level
                  L11_length(s11:e11 - 1), & ! link length
                  L11_slope(s11:e11 - 1), &
                  ! general INPUT/OUTPUT variables
                  L11_C1(s11:e11), & ! first muskingum parameter
                  L11_C2(s11:e11), & ! second muskigum parameter
                  L11_qOUT(s11:e11), & ! routed runoff flowing out of L11 cell
                  L11_qTIN(s11:e11,:), & ! inflow water into the reach at L11
                  L11_qTR(s11:e11,:), & !
                  L11_qMod(s11:e11), &
                  mRM_runoff(tt, :), &
                  ! original routing specific input/output variables
                  L11_FracFPimp(s11:e11), & ! fraction of impervious layer at L11 scale
                  ! OPTIONAL INPUT variables
                  do_mpr)
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

          ! update the counters
          if (day_counter   .NE. day  ) day_counter   = day
          if (month_counter .NE. month) month_counter = month
          if (year_counter  .NE. year)  year_counter  = year

          ! increment of timestep
          newTime = julday(day,month,year) + real(hour+timestep,dp)/24._dp
          ! calculate new year, month and day
          call caldat(int(newTime), yy=year, mm=month, dd=day)

          if (.not. optimize) then
#ifdef MRM2MHM
             if (any(outputFlxState_mrm)) then
                call mrm_write_output_fluxes( &
                     ! basin id
                     ii, &
                     ! output specification
                     timeStep_model_outputs_mrm, &
                     ! time specification
                     warmingDays_mrm(ii), newTime, nTimeSteps, nTStepDay, &
                     tt, &
                     ! parse previous date to mRM writer
                     day_counter, month_counter, year_counter, &
                     timestep, &
                     ! mask specification
                     mask11, &
                     ! output variables
                     L11_qmod(s11:e11))
             end if
#endif

             ! output only for evaluation period
             tIndex_out = (tt-warmingDays(ii)*NTSTEPDAY) ! tt if write out of warming period

             if ((any(outputFlxState)) .and. (tIndex_out .gt. 0_i4)) then

                average_counter = average_counter + 1

                if ( tIndex_out .EQ. 1 ) then
                   L1_fNotSealed = 1.0_dp - L1_fSealed
#ifdef pgiFortran154
                   nc = newOutputDataset(ii, mask1)
#else
                   nc = OutputDataset(ii, mask1)
#endif
                end if

                call nc%updateDataset( &
                     s1              , &
                     e1              , &
                     L1_fSealed      , &
                     L1_fNotSealed   , &
                     L1_inter        , &
                     L1_snowPack     , &
                     L1_soilMoist    , &
                     L1_soilMoistSat , &
                     L1_sealSTW      , &
                     L1_unsatSTW     , &
                     L1_satSTW       , &
                     L1_neutrons     , &
                     L1_pet_calc     , &
                     L1_aETSoil      , &
                     L1_aETCanopy    , &
                     L1_aETSealed    , &
                     L1_total_runoff , &
                     L1_runoffSeal   , &
                     L1_fastRunoff   , &
                     L1_slowRunoff   , &
                     L1_baseflow     , &
                     L1_percol       , &
                     L1_infilSoil    , &
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
                      if (((tIndex_out .gt. 1) .and. (day_counter .ne. day)) .or. (tt .eq. nTimeSteps))     writeout = .true.
                   case(-2) ! monthly
                      if (((tIndex_out .gt. 1) .and. (month_counter .ne. month)) .or. (tt .eq. nTimeSteps)) writeout = .true.
                   case(-3) ! yearly
                      if (((tIndex_out .gt. 1) .and. (year_counter .ne. year)) .or. (tt .eq. nTimeSteps))   writeout = .true.
                   case default ! no output at all
                      continue
                   end select
                endif

                if (writeout) then
                   call nc%writeTimestep(tIndex_out*timestep-1)
                end if

                if( tt .eq. nTimeSteps ) then
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
             if ( tt .EQ. 1 ) writeout_counter = 1
             ! only for evaluation period - ignore warming days
             if ( (tt-warmingDays(ii)*NTSTEPDAY) .GT. 0 ) then
                ! decide for daily, monthly or yearly aggregation
                select case(timeStep_sm_input)
                case(-1) ! daily
                   if (day   .NE. day_counter)   then
                      sm_opti(s1:e1,writeout_counter) = sm_opti(s1:e1,writeout_counter) / real(average_counter,dp)
                      writeout_counter = writeout_counter + 1
                      average_counter = 0
                   end if
                case(-2) ! monthly
                   if (month .NE. month_counter) then
                      sm_opti(s1:e1,writeout_counter) = sm_opti(s1:e1,writeout_counter) / real(average_counter,dp)
                      writeout_counter = writeout_counter + 1
                      average_counter = 0
                   end if
                case(-3) ! yearly
                   if (year  .NE. year_counter)  then
                      sm_opti(s1:e1,writeout_counter) = sm_opti(s1:e1,writeout_counter) / real(average_counter,dp)
                      writeout_counter = writeout_counter + 1
                      average_counter = 0
                   end if
                end select

                ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
                if (.not. (tt .eq. nTimeSteps) ) then
                   ! aggregate soil moisture to needed time step for optimization
                   sm_opti(s1:e1,writeout_counter) = sm_opti(s1:e1,writeout_counter) + &
                        sum(L1_soilMoist   (s1:e1, 1:nSoilHorizons_sm_input), dim=2) / &
                        sum(L1_soilMoistSat(s1:e1, 1:nSoilHorizons_sm_input), dim=2)
                end if

                ! increase average counter by one
                average_counter = average_counter + 1
             end if
          end if

          !----------------------------------------------------------------------
          ! FOR TOTAL WATER STORAGE
          if( present(basin_avg_tws) ) then
               area_basin = sum(L1_areaCell(s1:e1) )
               TWS_field(s1:e1) = L1_inter(s1:e1) + L1_snowPack(s1:e1) + L1_sealSTW(s1:e1) + &
                                        L1_unsatSTW(s1:e1) + L1_satSTW(s1:e1)
               do gg = 1, nSoilHorizons_mHM
                  TWS_field(s1:e1) =   TWS_field(s1:e1) + L1_soilMoist (s1:e1,gg)
                            end do
             basin_avg_TWS_sim(tt,ii) = ( dot_product( TWS_field (s1:e1), L1_areaCell(s1:e1) ) / area_basin )
          end if
          !----------------------------------------------------------------------
          
          !----------------------------------------------------------------------
          ! FOR NEUTRONS
          ! NOTE:: modeled neutrons are averaged daily
          !----------------------------------------------------------------------
          if (present(neutrons_opti)) then
             if ( tt .EQ. 1 ) writeout_counter = 1
             ! only for evaluation period - ignore warming days
             if ( (tt-warmingDays(ii)*NTSTEPDAY) .GT. 0 ) then
                ! decide for daily, monthly or yearly aggregation
                ! daily
                if (day   .NE. day_counter)   then
                   neutrons_opti(s1:e1,writeout_counter) = neutrons_opti(s1:e1,writeout_counter) / real(average_counter,dp)
                   writeout_counter = writeout_counter + 1
                   average_counter = 0
                end if

                ! last timestep is already done - write_counter exceeds size(sm_opti, dim=2)
                if (.not. (tt .eq. nTimeSteps) ) then
                   ! aggregate neutrons to needed time step for optimization
                   neutrons_opti(s1:e1,writeout_counter) = neutrons_opti(s1:e1,writeout_counter) + L1_neutrons(s1:e1)
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
             if ( tt .EQ. 1 ) then
                writeout_counter = 1
                L1_fNotSealed = 1.0_dp - L1_fSealed
             end if
             
             ! only for evaluation period - ignore warming days
             if ( (tt-warmingDays(ii)*NTSTEPDAY) .GT. 0 ) then
                ! decide for daily, monthly or yearly aggregation
                select case(timeStep_et_input)
                case(-1) ! daily
                   if (day   .NE. day_counter)   then
                      writeout_counter = writeout_counter + 1
                   end if
                case(-2) ! monthly
                   if (month .NE. month_counter) then
                      writeout_counter = writeout_counter + 1
                   end if
                case(-3) ! yearly
                   if (year  .NE. year_counter)  then
                      writeout_counter = writeout_counter + 1
                   end if
                end select

                ! last timestep is already done - write_counter exceeds size(et_opti, dim=2)
                if (.not. (tt .eq. nTimeSteps) ) then
                   ! aggregate evapotranspiration to needed time step for optimization
                   et_opti(s1:e1,writeout_counter) = et_opti(s1:e1,writeout_counter) + &
                        sum(L1_aETSoil(s1:e1,:), dim=2) * L1_fNotSealed(s1:e1) + &
                        L1_aETCanopy(s1:e1) + &
                        L1_aETSealed(s1:e1) * L1_fSealed(s1:e1)
                end if
             end if
          end if

       end do !<< TIME STEPS LOOP

       ! deallocate space for temprory LAI fields
       deallocate(LAI)
       ! deallocate TWS field temporal variable
       if (allocated(TWS_field) ) deallocate(TWS_field)

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
    if( present(basin_avg_tws) ) basin_avg_tws = basin_avg_TWS_sim

  end SUBROUTINE mhm_eval


END MODULE mo_mhm_eval
