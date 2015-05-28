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
  !                   Matthias Cuntz & Juliane Mai, Nov 2014 - LAI input from daily, monthly or yearly files
  !                   Matthias Zink,        Dec 2014 - adopted inflow gauges to ignore headwater cells

  SUBROUTINE mhm_eval(parameterset, runoff, sm_opti)

    use mo_init_states,         only : get_basin_info
    use mo_init_states,         only : variables_default_init   ! default initalization of variables
    use mo_julian,              only : caldat, julday
    use mo_message,             only : message
    use mo_mhm,                 only : mhm
    use mo_mhm_constants,       only : nodata_dp
    use mo_restart,             only : read_restart_states      ! read initial values of variables
    use mo_meteo_forcings,      only : prepare_meteo_forcings_data
    use mo_write_ascii,         only : write_daily_obs_sim_discharge
    use mo_write_fluxes_states, only : CloseFluxState_file
    use mo_write_fluxes_states, only : WriteFluxState
    use mo_write_fluxes_states, only : WriteFluxStateInit
    use mo_global_variables,    only : &
         timeStep_model_outputs, outputFlxState,             &  ! definition which output to write
         read_restart, perform_mpr, fracSealed_CityArea,     &
         timeStep_model_inputs,                              &
         timeStep, nBasins, basin, simPer, readPer,          & ! [h] simulation time step, No. of basins
         nGaugesTotal,                                       &
         processMatrix, c2TSTu, HorizonDepth_mHM,            & 
         nSoilHorizons_mHM, NTSTEPDAY, timeStep,             & 
         LCyearId, LAIUnitList, LAILUT,                      & 
         GeoUnitList, GeoUnitKar, soilDB,                    &
         L0_Id, L0_soilId,                                   & 
         L0_LCover, L0_asp, L0_LCover_LAI, L0_geoUnit,       &
         L0_areaCell,L0_floodPlain,                          &        
         soilDB, L1_areaCell, L1_nTCells_L0, L1_L11_Id,      & 
         L0_slope_emp,                                       &
         L1_upBound_L0, L1_downBound_L0, L1_leftBound_L0,    & 
         L1_rightBound_L0, latitude,                         &
         L11_netPerm, L11_fromN, L11_toN,                    & 
         L11_length, L11_slope, evap_coeff, fday_prec,       & 
         fnight_prec, fday_pet, fnight_pet, fday_temp,       & 
         fnight_temp, L1_pet, L1_tmin, L1_tmax, L1_netrad,   &
         L1_absvappress, L1_windspeed,                       &
         L1_pre, L1_temp , L1_fForest,                       & 
         L1_fPerm, L1_fSealed, L11_FracFPimp,                & 
         L11_aFloodPlain, L1_inter,                          & 
         L1_snowPack, L1_sealSTW, L1_soilMoist, L1_unsatSTW, & 
         L1_satSTW, L1_pet_calc,                             &
         L1_aETSoil, L1_aETCanopy, L1_aETSealed,             &
         L1_baseflow, L1_infilSoil, L1_fastRunoff, L1_melt,  & 
         L1_percol, L1_preEffect, L1_rain, L1_runoffSeal,    & 
         L1_slowRunoff, L1_snow, L1_Throughfall,             & 
         L1_total_runoff, L11_Qmod, L11_qOUT, L11_qTIN,      & 
         L11_qTR, L1_alpha, L1_degDayInc, L1_degDayMax,      & 
         L1_degDayNoPre, L1_degDay, L1_fAsp, L1_HarSamCoeff, & 
         L1_PrieTayAlpha, L1_aeroResist, L1_surfResist,      &
         L1_fRoots, L1_maxInter, L1_karstLoss, L1_kfastFlow, & 
         L1_kSlowFlow, L1_kBaseFlow, L1_kPerco,              & 
         L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp,   & 
         L1_tempThresh, L1_unsatThresh, L1_sealedThresh,     & 
         L1_wiltingPoint, L11_C1, L11_C2, L1_neutrons,       &
         warmingDays, evalPer, gauge, InflowGauge,           &  
         optimize,  nMeasPerDay,                             &
         timeStep_LAI_input,                                 & ! flag on how LAI data has to be read
         L0_gridded_LAI, dirRestartIn,                       & ! restart directory location
         timeStep_sm_input,                                  & ! time step of soil moisture input (day, month, year)
         nSoilHorizons_sm_input,                             & ! no. of mhm soil horizons equivalent to sm input 
         nTimeSteps_L1_sm                                      ! total number of timesteps in soil moisture input
    
    implicit none

    real(dp), dimension(:),                          intent(in)  :: parameterset
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: runoff       ! dim1=time dim2=gauge
    real(dp), dimension(:,:), allocatable, optional, intent(out) :: sm_opti      ! dim1=ncells, dim2=time

    ! -------------------------------------
    ! local variables
    !
    ! FOR WRITING GRIDDED STATES AND FLUXES 
    integer(i4)                           :: hh                  ! Counter
    integer(i4)                           :: ncid                ! netcdf fileID
    integer(i4)                           :: tIndex_out          ! for writing netcdf file
    ! States L1
    real(dp), dimension(:),   allocatable :: L1_inter_out        ! Interception
    real(dp), dimension(:),   allocatable :: L1_snowPack_out     ! Snowpack
    real(dp), dimension(:,:), allocatable :: L1_soilMoist_out    ! Soil moisture of each horizon
    real(dp), dimension(:),   allocatable :: L1_sealSTW_out      ! Retention storage of impervious areas
    real(dp), dimension(:),   allocatable :: L1_unsatSTW_out     ! Upper soil storage
    real(dp), dimension(:),   allocatable :: L1_satSTW_out       ! Groundwater storage
    real(dp), dimension(:),   allocatable :: L1_neutrons_out     ! Ground albedo neutrons
    ! Fluxes L1
    real(dp), dimension(:),   allocatable :: L1_pet_out          ! potential evapotranpiration (PET)
    real(dp), dimension(:,:), allocatable :: L1_aETSoil_out      ! actual ET of each horizon
    real(dp), dimension(:),   allocatable :: L1_aETCanopy_out    ! Real evaporation intensity from canopy
    real(dp), dimension(:),   allocatable :: L1_aETSealed_out    ! Actual ET from free-water surfaces
    real(dp), dimension(:),   allocatable :: L1_total_runoff_out ! Generated runoff
    real(dp), dimension(:),   allocatable :: L1_runoffSeal_out   ! Direct runoff from impervious areas
    real(dp), dimension(:),   allocatable :: L1_fastRunoff_out   ! Fast runoff component
    real(dp), dimension(:),   allocatable :: L1_slowRunoff_out   ! Slow runoff component
    real(dp), dimension(:),   allocatable :: L1_baseflow_out     ! Baseflow
    real(dp), dimension(:),   allocatable :: L1_percol_out       ! Percolation
    real(dp), dimension(:,:), allocatable :: L1_infilSoil_out    ! Infiltration 
    !
    ! counters and indexes
    integer(i4)                               :: nTimeSteps
    integer(i4)                               :: maxTimeSteps
    integer(i4)                               :: ii, tt, gg, ll   ! Counters
    integer(i4)                               :: nCells           ! No. of cells at level 1 for current basin
    integer(i4)                               :: nNodes           !
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
    integer(i4)                               :: s11, e11         ! process 8: start and end index of vectors (on or off)
    integer(i4)                               :: s110, e110
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
    real(dp)                                  :: multiplier       ! for averaging output
    logical                                   :: writeout         ! if true write out netcdf files
    integer(i4)                               :: writeout_counter ! write out time step
    !
    ! for discharge timeseries
    integer(i4)                               :: iday, iS, iE
    real(dp), dimension(:,:), allocatable     :: d_Qmod
    !
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
       else 
          !----------------------------------------------------------
          ! estimate maximum modeling timesteps including warming days
          !----------------------------------------------------------
          maxTimeSteps = maxval( simPer(1:nBasins)%julEnd - simPer(1:nBasins)%julStart + 1 ) * NTSTEPDAY
          allocate( runoff(maxTimeSteps, nGaugesTotal) )
          runoff = nodata_dp
       end if
    else 
       if ( (processMatrix(8,1) .gt. 0) .AND. (.NOT. optimize)) then
          call message("***ERROR: runoff can not be produced, since runoff variable is not present")
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
    ! add other optionals...


    !-------------------------------------------------------------------
    ! Initalize State variables either to the default value or 
    ! from the restrat_files.
    ! All variables to be initalized had been allocated to the required 
    ! space before this point (see, mo_startup: initialise)
    !-------------------------------------------------------------------
    if (.NOT. read_restart ) then
       ! as default values, 
       ! all cells for all modeled basins are simultenously initalized ONLY ONCE
       call variables_default_init()
    else 
       ! read from restart files, basin wise ...
       do ii = 1, nBasins
          call read_restart_states(ii, dirRestartIn(ii) )
       end do
    end if


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
       call get_basin_info ( ii,110, nrows, ncols,                iStart=s110,iEnd=e110 ) 
       call get_basin_info ( ii,  1, nrows, ncols, ncells=nCells, iStart=s1,  iEnd=e1, mask=mask1 ) 

       ! process 8 - routing process (on or off)
       if( processMatrix(8, 1) .eq. 0 ) then
          s11 = 1
          e11 = 1
          nNodes = 1
       else
          call get_basin_info ( ii, 11, nrows, ncols, ncells=nNodes, iStart=s11, iEnd=e11 ) 
       end if

       ! allocate space for local LAI grid
       allocate( LAI(s0:e0) )
       LAI(:) = nodata_dp

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
          case(1) ! HarSam
             s_p5 = (/s_meteo, s_meteo, s_meteo,       1,       1,       1/)
             e_p5 = (/e_meteo, e_meteo, e_meteo,       1,       1,       1/)
          case(2) ! PrieTay
             s_p5 = (/s_meteo,       1,       1, s_meteo,       1,       1/)
             e_p5 = (/e_meteo,       1,       1, e_meteo,       1,       1/)
          case(3) ! PenMon
             s_p5 = (/s_meteo,       1,       1, s_meteo, s_meteo, s_meteo/)
             e_p5 = (/e_meteo,       1,       1, e_meteo, e_meteo, e_meteo/)
          end select

          ! customize iMeteoTS for process 5 - PET
          select case (processMatrix(5,1))
          !              (/     pet,     tmin,     tmax,   netrad,  absVapP,windspeed /)  
          case(0) ! PET is input
             iMeteo_p5 = (/iMeteoTS,        1,        1,        1,        1,        1 /)
          case(1) ! HarSam
             iMeteo_p5 = (/iMeteoTS, iMeteoTS, iMeteoTS,        1,        1,        1 /)
          case(2) ! PrieTay
             iMeteo_p5 = (/iMeteoTS,        1,        1, iMeteoTS,        1,        1 /)
          case(3) ! PenMon
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
             !              
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
               nCells, nNodes, nSoilHorizons_mHM, real(NTSTEPDAY,dp), timeStep, mask0,      & ! IN C 
               basin%nInflowGauges(ii), basin%InflowGaugeIndexList(ii,:),                   & ! IN C
               basin%InflowGaugeHeadwater(ii,:), basin%InflowGaugeNodeList(ii,:),           & ! IN C
               parameterset,                                                                & ! IN P
               LCyearId(year,ii), GeoUnitList, GeoUnitKar, LAIUnitList, LAILUT,             & ! IN L0
               L0_slope_emp(s0:e0), L0_Id(s0:e0), L0_soilId(s0:e0), L0_LCover_LAI(s0:e0),   & ! IN L0
               L0_LCover(s0:e0, LCyearId(year,ii)), L0_asp(s0:e0), LAI(s0:e0),              & ! IN L0
               L0_geoUnit(s0:e0), L0_areaCell(s0:e0),L0_floodPlain(s110:e110),              & ! IN L0
               soilDB%is_present, soilDB%nHorizons, soilDB%nTillHorizons,                   & ! IN L0
               soilDB%sand, soilDB%clay, soilDB%DbM, soilDB%Wd, soilDB%RZdepth,             & ! IN L0
               L1_areaCell(s1:e1), L1_nTCells_L0(s1:e1),  L1_L11_Id(s1:e1),                 & ! IN L1
               L1_upBound_L0(s1:e1), L1_downBound_L0(s1:e1),                                & ! IN L1
               L1_leftBound_L0(s1:e1), L1_rightBound_L0(s1:e1),                             & ! IN L1
               latitude(s_p5(1):e_p5(1)),                                                   & ! IN L1
               L11_netPerm(s11:e11), L11_fromN(s11:e11), L11_toN(s11:e11),                  & ! IN L11
               L11_length(s11:e11), L11_slope(s11:e11),                                     & ! IN L11
               evap_coeff, fday_prec, fnight_prec, fday_pet, fnight_pet,                    & ! IN F
               fday_temp, fnight_temp,                                                      & ! IN F
               L1_pet(s_p5(1):e_p5(1), iMeteo_p5(1)),                                       & ! INOUT F:PET
               L1_tmin(s_p5(2):e_p5(2), iMeteo_p5(2)),                                      & ! IN F:PET
               L1_tmax(s_p5(3):e_p5(3), iMeteo_p5(3)),                                      & ! IN F:PET
               L1_netrad(s_p5(4):e_p5(4), iMeteo_p5(4)),                                    & ! IN F:PET
               L1_absvappress(s_p5(5):e_p5(5), iMeteo_p5(5)),                               & ! IN F:PET
               L1_windspeed(s_p5(6):e_p5(6), iMeteo_p5(6)),                                 & ! IN F:PET
               L1_pre(s_meteo:e_meteo,iMeteoTS),                                            & ! IN F:Pre 
               L1_temp(s_meteo:e_meteo,iMeteoTS),                                           & ! IN F:Temp
               InflowGauge%Q(iMeteoTS,:),                                                   & ! IN Q
               yId,                                                                         & ! INOUT C
               L1_fForest(s1:e1), L1_fPerm(s1:e1),  L1_fSealed(s1:e1),                      & ! INOUT L1 
               L11_FracFPimp(s11:e11), L11_aFloodPlain(s11:e11),                            & ! INOUT L11
               L1_inter(s1:e1), L1_snowPack(s1:e1), L1_sealSTW(s1:e1),                      & ! INOUT S 
               L1_soilMoist(s1:e1,:), L1_unsatSTW(s1:e1), L1_satSTW(s1:e1), L1_neutrons,    & ! INOUT S 
               L1_pet_calc(s1:e1),                                                          & ! INOUT X
               L1_aETSoil(s1:e1,:), L1_aETCanopy(s1:e1), L1_aETSealed(s1:e1),               & ! INOUT X
               L1_baseflow(s1:e1), L1_infilSoil(s1:e1,:), L1_fastRunoff(s1:e1),             & ! INOUT X
               L1_melt(s1:e1), L1_percol(s1:e1), L1_preEffect(s1:e1), L1_rain(s1:e1),       & ! INOUT X
               L1_runoffSeal(s1:e1), L1_slowRunoff(s1:e1), L1_snow(s1:e1),                  & ! INOUT X
               L1_Throughfall(s1:e1), L1_total_runoff(s1:e1),                               & ! INOUT X
               L11_Qmod(s11:e11), L11_qOUT(s11:e11),L11_qTIN(s11:e11,:),L11_qTR(s11:e11,:), & ! INOUT X11
               L1_alpha(s1:e1), L1_degDayInc(s1:e1), L1_degDayMax(s1:e1),                   & ! INOUT E1
               L1_degDayNoPre(s1:e1), L1_degDay(s1:e1), L1_fAsp(s1:e1),                     & ! INOUT E1
               L1_HarSamCoeff(s1:e1), L1_PrieTayAlpha(s1:e1,:), L1_aeroResist(s1:e1,:),     & ! INOUT E1
               L1_surfResist(s1:e1,:), L1_fRoots(s1:e1,:),                                  & ! INOUT E1
               L1_maxInter(s1:e1), L1_karstLoss(s1:e1),  L1_kFastFlow(s1:e1),               & ! INOUT E1
               L1_kSlowFlow(s1:e1), L1_kBaseFlow(s1:e1), L1_kPerco(s1:e1),                  & ! INOUT E1
               L1_soilMoistFC(s1:e1,:), L1_soilMoistSat(s1:e1,:), L1_soilMoistExp(s1:e1,:), & ! INOUT E1
               L1_tempThresh(s1:e1), L1_unsatThresh(s1:e1), L1_sealedThresh(s1:e1),         & ! INOUT E1
               L1_wiltingPoint(s1:e1,:),                                                    & ! INOUT E1
               L11_C1(s11:e11), L11_C2(s11:e11)                                             ) ! INOUT E11


          ! update the counters
          if (day_counter   .NE. day  ) day_counter   = day
          if (month_counter .NE. month) month_counter = month
          if (year_counter  .NE. year)  year_counter  = year

          ! increment of timestep
          newTime = julday(day,month,year) + real(hour+timestep,dp)/24._dp
          call caldat(int(newTime), yy=year, mm=month, dd=day)

          if (.not. optimize) then

             ! output only for evaluation period
             tIndex_out = (tt-warmingDays(ii)*NTSTEPDAY) ! tt if write out of warming period

             if ((any(outputFlxState)) .and. (tIndex_out .gt. 0_i4)) then

                average_counter = average_counter + 1

                ! check the compatibility of timeStep_model_outputs and also initialize NetCDF variables
                if ( tIndex_out .EQ. 1 ) then

                   ! initalize
                   call WriteFluxStateInit( &
                        ! Input
                        ii                     , &
                        timeStep_model_outputs , &
                        ! Inout: States L1
                        L1_inter_out           , & ! Interception
                        L1_snowPack_out        , & ! Snowpack
                        L1_soilMoist_out       , & ! Soil moisture of each horizon
                        L1_sealSTW_out         , & ! Retention storage of impervious areas
                        L1_unsatSTW_out        , & ! Upper soil storage
                        L1_satSTW_out          , & ! Groundwater storage
                        L1_neutrons_out        , & ! ground albedo neutrons
                        ! Inout: Fluxes L1
                        L1_pet_out             , & ! potential evapotranspiration (PET)
                        L1_aETSoil_out         , & ! actual ET
                        L1_aETCanopy_out       , & ! Real evaporation intensity from canopy
                        L1_aETSealed_out       , & ! Actual ET from free-water surfaces
                        L1_total_runoff_out    , & ! Generated runoff
                        L1_runoffSeal_out      , & ! Direct runoff from impervious areas
                        L1_fastRunoff_out      , & ! Fast runoff component
                        L1_slowRunoff_out      , & ! Slow runoff component
                        L1_baseflow_out        , & ! Baseflow
                        L1_percol_out          , & ! Percolation 
                        L1_infilSoil_out       , & ! Infiltration 
                        ! Output
                        ncid )


                end if ! <- tIndex_out CONDITION

                ! PREPARE FIELDS ON REQUIRED OUTPUT TIME STEP

                ! States L1 --> AVERAGE
                if (outputFlxState(1) ) L1_inter_out    (:)   = L1_inter_out    (:)   + L1_inter    (s1:e1)
                if (outputFlxState(2) ) L1_snowPack_out (:)   = L1_snowPack_out (:)   + L1_snowPack (s1:e1)
                if (outputFlxState(3) .OR. &
                     outputFlxState(4) .OR. &
                     outputFlxState(5) ) L1_soilMoist_out(:,:) = L1_soilMoist_out(:,:) + L1_soilMoist(s1:e1,:)
                if (outputFlxState(6) ) L1_sealSTW_out  (:)   = L1_sealSTW_out  (:)   + L1_sealSTW  (s1:e1)
                if (outputFlxState(7) ) L1_unsatSTW_out (:)   = L1_unsatSTW_out (:)   + L1_unsatSTW (s1:e1)
                if (outputFlxState(8) ) L1_satSTW_out   (:)   = L1_satSTW_out   (:)   + L1_satSTW   (s1:e1)
                if (outputFlxState(18) ) L1_neutrons_out(:)   = L1_neutrons_out (:)   + L1_neutrons (s1:e1)

                ! Fluxes L1  --> AGGREGATED
                if (outputFlxState(9) ) &
                     L1_pet_out(:)          = L1_pet_out(:)        + L1_pet_calc(s1:e1)
                if (outputFlxState(10)      ) then
                   do hh = 1, nSoilHorizons_mHM
                      L1_aETSoil_out(:,hh)  = L1_aETSoil_out(:,hh) + L1_aETSoil(s1:e1,hh)*(1.0_dp - L1_fSealed(s1:e1))
                   end do
                end if
                if (outputFlxState(10) ) &
                     L1_aETCanopy_out(:)    = L1_aETCanopy_out(:)    + L1_aETCanopy(s1:e1)
                if (outputFlxState(10) ) &
                     L1_aETSealed_out(:)    = L1_aETSealed_out(:)    + L1_aETSealed(s1:e1)*L1_fSealed(s1:e1)
                if (outputFlxState(11)) &
                     L1_total_runoff_out(:) = L1_total_runoff_out(:) + L1_total_runoff(s1:e1)
                if (outputFlxState(12)) &
                     L1_runoffSeal_out(:)   = L1_runoffSeal_out(:)   + L1_runoffSeal(s1:e1)*L1_fSealed(s1:e1)
                if (outputFlxState(13)) &
                     L1_fastRunoff_out(:)   = L1_fastRunoff_out(:)   + L1_fastRunoff(s1:e1)*(1.0_dp - L1_fSealed(s1:e1))
                if (outputFlxState(14)) &
                     L1_slowRunoff_out(:)   = L1_slowRunoff_out(:)   + L1_slowRunoff(s1:e1)*(1.0_dp - L1_fSealed(s1:e1))
                if (outputFlxState(15)) &
                     L1_baseflow_out(:)     = L1_baseflow_out(:)     + L1_baseflow(s1:e1)*(1.0_dp - L1_fSealed(s1:e1))
                if (outputFlxState(16)) &
                     L1_percol_out(:)       = L1_percol_out(:)       + L1_percol(s1:e1)*(1.0_dp - L1_fSealed(s1:e1))
                if (outputFlxState(17)      ) then 
                   do hh = 1, nSoilHorizons_mHM 
                      L1_infilSoil_out(:,hh) = L1_infilSoil_out(:,hh) + L1_infilSoil(s1:e1,hh)*(1.0_dp - L1_fSealed(s1:e1)) 
                   end do
                end if

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
                   ! Average States
                   multiplier = 1.0_dp/real(average_counter,dp)
                   if (outputFlxState(1)) L1_inter_out(:)    = L1_inter_out(:)    * multiplier
                   if (outputFlxState(2)) L1_snowPack_out(:) = L1_snowPack_out(:) * multiplier
                   if (outputFlxState(3) .OR. outputFlxState(4) .OR. outputFlxState(5)) &
                        L1_soilMoist_out(:,:) = L1_soilMoist_out(:,:) * multiplier
                   if (outputFlxState(6)) L1_sealSTW_out(:)  = L1_sealSTW_out(:)  * multiplier
                   if (outputFlxState(7)) L1_unsatSTW_out(:) = L1_unsatSTW_out(:) * multiplier
                   if (outputFlxState(8)) L1_satSTW_out(:)   = L1_satSTW_out(:)   * multiplier
                   if (outputFlxState(18)) L1_neutrons_out(:)= L1_neutrons_out(:) * multiplier

                   average_counter = 0

                   ! Write netcdf file
                   writeout_counter = writeout_counter + 1
                   call WriteFluxState(tIndex_out*timestep-1, writeout_counter, ncid, ii, mask1, &
                        ! States L1
                        L1_inter_out             , & ! Interception
                        L1_snowPack_out          , & ! Snowpack
                        L1_soilMoist_out         , & ! Soil moisture of each horizon
                        L1_sealSTW_out           , & ! Retention storage of impervious areas
                        L1_unsatSTW_out          , & ! Upper soil storage
                        L1_satSTW_out            , & ! Groundwater storage
                        L1_neutrons_out          , & ! ground albedo neutrons
                        ! Fluxes L1
                        L1_pet_out               , & ! potential evapotranspiration (PET)
                        L1_aETSoil_out           , & ! actual ET
                        L1_aETCanopy_out         , & ! Real evaporation intensity from canopy
                        L1_aETSealed_out         , & ! Actual ET from free-water surfaces
                        L1_total_runoff_out      , & ! Generated runoff
                        L1_runoffSeal_out        , & ! Direct runoff from impervious areas
                        L1_fastRunoff_out        , & ! Fast runoff component
                        L1_slowRunoff_out        , & ! Slow runoff component
                        L1_baseflow_out          , & ! Baseflow
                        L1_percol_out            , & ! Percolation
                        L1_soilMoistSat(s1:e1 ,:), & ! Saturation soil moisture for each horizon [mm] 
                        L1_infilSoil_out)            ! Infiltration

                   ! set variables to zero
                   ! States L1
                   if (outputFlxState(1)  ) L1_inter_out(:)        = 0.0_dp       
                   if (outputFlxState(2)  ) L1_snowPack_out(:)     = 0.0_dp      
                   if( outputFlxState(3) .OR. &
                        outputFlxState(4) .OR. &
                        outputFlxState(5)  ) L1_soilMoist_out(:,:) = 0.0_dp       
                   if (outputFlxState(6)  ) L1_sealSTW_out(:)      = 0.0_dp      
                   if (outputFlxState(7)  ) L1_unsatSTW_out(:)     = 0.0_dp      
                   if (outputFlxState(8)  ) L1_satSTW_out(:)       = 0.0_dp      
                   if (outputFlxState(18)  ) L1_neutrons_out(:)    = 0.0_dp      
                   ! Fluxes L1
                   if (outputFlxState(9)  ) L1_pet_out(:)          = 0.0_dp     
                   if (outputFlxState(10) ) L1_aETSoil_out(:,:)    = 0.0_dp     
                   if (outputFlxState(10) ) L1_aETCanopy_out(:)    = 0.0_dp    
                   if (outputFlxState(10) ) L1_aETSealed_out(:)    = 0.0_dp    
                   if (outputFlxState(11) ) L1_total_runoff_out(:) = 0.0_dp   
                   if (outputFlxState(12) ) L1_runoffSeal_out(:)   = 0.0_dp   
                   if (outputFlxState(13) ) L1_fastRunoff_out(:)   = 0.0_dp   
                   if (outputFlxState(14) ) L1_slowRunoff_out(:)   = 0.0_dp   
                   if (outputFlxState(15) ) L1_baseflow_out(:)     = 0.0_dp   
                   if (outputFlxState(16) ) L1_percol_out(:)       = 0.0_dp
                   if (outputFlxState(17) ) L1_infilSoil_out(:,:)  = 0.0_dp   
                end if
                !    
                ! close file and deallocate variables
                if( tt .eq. nTimeSteps ) then
                   call CloseFluxState_file(ii, ncid)
                   ! States L1
                   if (outputFlxState(1)  ) deallocate( L1_inter_out       )         
                   if (outputFlxState(2)  ) deallocate( L1_snowPack_out    )        
                   if( outputFlxState(3) .OR. &
                       outputFlxState(4) .OR. &
                       outputFlxState(5)  ) deallocate( L1_soilMoist_out   )  
                   if (outputFlxState(6)  ) deallocate( L1_sealSTW_out     )        
                   if (outputFlxState(7)  ) deallocate( L1_unsatSTW_out    )        
                   if (outputFlxState(8)  ) deallocate( L1_satSTW_out      )        
                   if (outputFlxState(18)  ) deallocate( L1_neutrons_out   )        
                   ! Fluxes L1
                   if (outputFlxState(9)   ) deallocate( L1_pet_out        )    
                   if (outputFlxState(10)  ) deallocate( L1_aETSoil_out    )    
                   if (outputFlxState(10)  ) deallocate( L1_aETCanopy_out  )   
                   if (outputFlxState(10)  ) deallocate( L1_aETSealed_out  )   
                   if (outputFlxState(11) ) deallocate( L1_total_runoff_out)  
                   if (outputFlxState(12) ) deallocate( L1_runoffSeal_out  )  
                   if (outputFlxState(13) ) deallocate( L1_fastRunoff_out  )  
                   if (outputFlxState(14) ) deallocate( L1_slowRunoff_out  )  
                   if (outputFlxState(15) ) deallocate( L1_baseflow_out    )  
                   if (outputFlxState(16) ) deallocate( L1_percol_out      )  
                   if (outputFlxState(17) ) deallocate( L1_infilSoil_out   ) 
                   !
                end if
                !
             end if
          end if ! <-- if (.not. optimize)

          !----------------------------------------------------------------------
          ! FOR STORING the optional arguments
          ! 
          ! FOR RUNOFF
          ! NOTE:: Node ID for a given gauging station is stored at gaugeindex's
          !        index in runoff. In consequence the gauges in runoff are 
          !        ordered corresponing to gauge%Q(:,:)
          !----------------------------------------------------------------------
          if( present(runoff) ) then
             do gg = 1, basin%nGauges(ii)
                runoff(tt,basin%gaugeIndexList(ii,gg)) = L11_Qmod( basin%gaugeNodeList(ii,gg) + s11 - 1 )
             end do
          end if

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

       end do !<< TIME STEPS LOOP

       ! deallocate space for temprory LAI fields
       deallocate(LAI)

    end do !<< BASIN LOOP

    ! --------------------------------------------------------------------------
    ! STORE DAILY DISCHARGE TIMESERIES OF EACH GAUGING STATION 
    ! FOR SIMULATIONS DURING THE EVALUATION PERIOD
    !
    !  **** AT DAILY TIME STEPS ****
    ! Note:: Observed Q are stored only for the evaluation period and not for
    !        the warming days
    ! --------------------------------------------------------------------------
    if( (.not. optimize) .AND. present(runoff) .AND. (nMeasPerDay .eq. 1) ) then
       !
       ii = maxval( evalPer(1:nBasins)%julEnd - evalPer(1:nBasins)%julStart + 1 )
       allocate( d_Qmod(ii, nGaugesTotal) ) 
       d_Qmod = 0.0_dp

       ! loop over basins
       do ii = 1, nBasins
          iDay = 0
          ! loop over timesteps
          do tt = warmingDays(ii)*NTSTEPDAY+1, nTimeSteps, NTSTEPDAY
             iS = tt
             iE = tt + NTSTEPDAY - 1
             iDay = iDay + 1
             ! over gauges
             do gg = 1, basin%nGauges(ii)
                d_Qmod(iDay, basin%gaugeIndexList(ii,gg) ) = &
                     sum( runoff(iS:iE, basin%gaugeIndexList(ii,gg)) )/ real(NTSTEPDAY,dp)
             end do
             !
          end do
       end do
       ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
       call write_daily_obs_sim_discharge( gauge%Q(:,:), d_Qmod(:,:) )
       ! free space
       deallocate(d_Qmod)        
       !
    end if


  end SUBROUTINE mhm_eval


END MODULE mo_mhm_eval
