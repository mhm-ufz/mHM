!> \file mo_read_config.f90

!> \brief Reading of main model configurations.

!> \details This routine reads the configurations of mHM including, input and
!>          output directories, module usage specification, simulation time periods,
!>          global parameters, ...

!> \authors Matthias Zink
!> \date Dec 2012

MODULE mo_read_config

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_config ! read main directories

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_config

  !     PURPOSE
  !>        \brief Read main configurations for mHM

  !>        \details The main configurations in mHM are read from three files:
  !>                 <ol>
  !>                   <li> mhm.nml
  !>                   <li> mhm_parameters.nml
  !>                   <li> mhm_outputs.nml
  !>                 </ol>
  !>                 For details please refer to the above mentioned namelist files.

  !     CALLING SEQUENCE
  !         None

  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Dec 2012
  !         Modified Luis Samaniego, Jan 2013 - messages
  !                  Rohini Kumar
  !                  Matthias Cuntz, Jan  2013 - namelist consolidation and positioning
  !                  Matthias Zink,  Jan  2013 - bug fix, added gaugeinfo reading
  !                  Rohini Kumar,   Jun  2013 - added restart flags
  !                  R. Kumar &
  !                  S. Thober,      Aug. 2013 - code change to incorporate output timestep
  !                                              during writing of the netcdf file
  !                  Rohini Kumar,   Aug  2013 - name changed from "inputFormat" to inputFormat_meteo_forcings
  !                  Rohini Kumar,   Aug  2013 - added dirSoil_LUT and dirGeology_LUT, and changed in
  !                                              namelist made accordingly
  !                  Rohini Kumar,   Aug  2013 - added new namelist for LAI related datasets, and changed in
  !                                              within the code made accordingly
  !                  Matthias Zink,  Aug  2013 - changed read in for land cover period
  !                  Juliane Mai,    Oct  2013 - adding global_parameters_name
  !                  Matthias Zink,  Nov  2013 - edited documentation and included DEFAULT cases for ptocess Matrix
  !                  Stephan Thober, Nov  2013 - added read of directories where latitude longitude fields are located
  !                  Matthias Zink,  Feb  2014 - added multiple options for PET process
  !                  Matthias Zink,  Mar  2014 - added inflow from upstream areas and gauge information as namelist
  !                  Rohini Kumar,   May  2014 - added options for the model run coordinate system
  !                  Stephan Thober, May  2014 - added switch for chunk read in
  !                  Stephan Thober, Jun  2014 - added option for switching off mpr
  !                  Matthias Cuntz & Juliane Mai Nov 2014 - LAI input from daily, monthly or yearly files
  !                  Matthias Zink,  Dec  2014 - adopted inflow gauges to ignore headwater cells
  !                  Matthias Zink,  Mar  2015 - added optional soil moisture read in for calibration
  !                  Matthias Cuntz, Jul  2015 - removed adjustl from trim(adjustl()) of Geoparams for compilation with PGI
  !                  Stephan Thober, Aug  2015 - added read_config_routing and read_routing_params from mRM
  !                  Oldrich Rakovec,Oct  2015 - added reading of the basin average TWS data
  !                     Rohini Kumar, Mar 2016 - options to handle different soil databases

  subroutine read_config()

    use mo_julian,           only: date2dec, dec2date
    use mo_append,           only: append
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_nml,              only: open_nml, close_nml, position_nml
    use mo_constants,        only: eps_dp                     ! epsilon(1.0) in double precision 
    use mo_mhm_constants,    only:                          &
         nodata_i4, nodata_dp,                              & ! nodata values
         nColPars,                                          & ! number of properties of the global variables
         maxNoSoilHorizons,                                 & ! maximum number of allowed soil layers
         maxNoBasins,                                       & ! maximum number of allowed basins
         maxNLcovers,                                       & ! maximum number of allowed LCover scenes
         maxGeoUnit                                           ! maximum number of allowed geological classes
    use mo_file,             only:                          &
         file_namelist, unamelist,                          & ! file containing main configurations
         file_namelist_param, unamelist_param,              & ! file containing parameter values
         file_defOutput, udefOutput,                        & ! file specifying which output to write
         file_geolut, ugeolut                                 ! file specifying geological formations
    use mo_global_variables, only:                          &
         timestep,                                          & ! model time step
         period,                                            & ! data structure for period
         timestep_model_inputs,                             & ! read input frequency
         resolutionHydrology,                               & ! resolutions of hydrology
         resolutionRouting,                                 & ! resolution of routing
         L0_Basin,                                          & ! L0_Basin ID
         dirMorpho, dirLCover,                              & ! input directory of morphological
         dirPrecipitation, dirTemperature,                  & ! directory of meteo input
         dirReferenceET,                                    & ! PET input path  if process 5 is 'PET is input' (case 0)
         dirMinTemperature, dirMaxTemperature,              & ! PET input paths if process 5 is Hargreaves-Samani (case 1)
         dirNetRadiation,                                   & ! PET input paths if process 5 is Priestely-Taylor (case 2)
         dirabsVapPressure, dirwindspeed,                   & ! PET input paths if process 5 is Penman-Monteith (case 3)
         inputFormat_meteo_forcings,                        & ! input format either bin or nc
         fileLatLon,                                        & ! directory of latitude and longitude files
         dirConfigOut,                                      & ! configuration run output directory
         dirCommonFiles,                                    & ! directory where common files are located
         dirOut,                                            & ! output directory basin wise
         dirRestartOut,                                     & ! output directory of restart file basin wise
         dirRestartIn,                                      & ! input directory of restart file basin wise
         dirgridded_LAI,                                    & ! directory where gridded LAI is located
         dirSoil_moisture, timeStep_sm_input,               & ! directory and time stepping of soil moisture data
         nSoilHorizons_sm_input,                            & ! No. of mhm soil horizons equivalent to soil moisture input
         basin_avg_TWS_obs,                                 & ! basin avg TWS data
         fileTWS,                                           & ! directory with basin average tws data
         dirNeutrons, timeStep_neutrons_input,              & ! directory where neutron data is located
         iFlag_soilDB,                                      & ! options to handle different types of soil databases
         HorizonDepth_mHM, nSoilHorizons_mHM, tillageDepth, & ! soil horizons info for mHM
         fracSealed_cityArea, nLcoverScene,                 & ! land cover information
         LCfilename, LCyearId,                              & !
         nBasins,                                           & ! number of basins
         read_restart,                                      & ! flag reading restart
         write_restart,                                     & ! flag writing restart
         perform_mpr,                                       & ! switch for performing mpr
         warmingDays, warmPer,                              & ! warming days and warming period
         evalPer, simPer,                                   & ! model eval. & sim. periods
                                !                                                    ! (sim. = wrm. + eval.)
         evap_coeff,                                        & ! pan evaporation
         fday_prec, fnight_prec, fday_pet,                  & ! day-night fraction
         fnight_pet, fday_temp, fnight_temp,                & ! day-night fraction
         nProcesses, processMatrix,                         & ! process configuration
         nGeoUnits,                                         & ! number of geological classes
                                !                                                    ! for parameter read-in
         timeStep_model_outputs,                            & ! timestep for writing model outputs
         outputFlxState,                                    & ! definition which output to write
         inputFormat_gridded_LAI,                           & ! format of gridded LAI data(bin or nc)
         timeStep_LAI_input,                                & ! time step of gridded LAI input
         iFlag_cordinate_sys                                  ! model run cordinate system

    use mo_common_variables, only: &
         global_parameters,                                 & ! global parameters
         global_parameters_name,                            & ! clear names of global parameters
         optimize,                                          & ! if mhm runs in optimization mode or not
         optimize_restart,                                  & ! Optimization will be restarted from
         !                                                    ! mo_<opti_method>.restart file (.true.)
         opti_method,                                       & ! optimization algorithm used
         opti_function,                                     & ! objective function to be optimized
         nIterations,                                       & ! number of iterations in optimization
         seed,                                              & ! seed used for optimization
         dds_r,                                             & ! DDS: perturbation rate
         sa_temp,                                           & ! SA: initial temperature
         sce_ngs, sce_npg, sce_nps,                         & ! SCE: # complexes, # points per complex,
         !                                                    !      # points per subcomplex
         mcmc_opti,                                         & ! MCMC: if optimization mode of MCMC or only uncertainty estimation
         mcmc_error_params                                    !       parameters of error model used in likelihood

    implicit none

    ! LOCAL variables
    ! PARAMETERS
    integer(i4), dimension(nProcesses)              :: processCase             ! Choosen process description number
    real(dp), dimension(5, nColPars)                :: dummy_2d_dp = nodata_dp ! space holder for routing parameters
    ! interception
    real(dp), dimension(nColPars)                   :: canopyInterceptionFactor
    ! snow
    real(dp), dimension(nColPars)                   :: snowTreshholdTemperature
    real(dp), dimension(nColPars)                   :: degreeDayFactor_forest
    real(dp), dimension(nColPars)                   :: degreeDayFactor_impervious
    real(dp), dimension(nColPars)                   :: degreeDayFactor_pervious
    real(dp), dimension(nColPars)                   :: increaseDegreeDayFactorByPrecip
    real(dp), dimension(nColPars)                   :: maxDegreeDayFactor_forest
    real(dp), dimension(nColPars)                   :: maxDegreeDayFactor_impervious
    real(dp), dimension(nColPars)                   :: maxDegreeDayFactor_pervious
    ! soilmoisture
    real(dp), dimension(nColPars)                   :: orgMatterContent_forest
    real(dp), dimension(nColPars)                   :: orgMatterContent_impervious
    real(dp), dimension(nColPars)                   :: orgMatterContent_pervious
    real(dp), dimension(nColPars)                   :: PTF_lower66_5_constant
    real(dp), dimension(nColPars)                   :: PTF_lower66_5_clay
    real(dp), dimension(nColPars)                   :: PTF_lower66_5_Db
    real(dp), dimension(nColPars)                   :: PTF_higher66_5_constant
    real(dp), dimension(nColPars)                   :: PTF_higher66_5_clay
    real(dp), dimension(nColPars)                   :: PTF_higher66_5_Db
    real(dp), dimension(nColPars)                   :: infiltrationShapeFactor
    real(dp), dimension(nColPars)                   :: PTF_Ks_constant
    real(dp), dimension(nColPars)                   :: PTF_Ks_sand
    real(dp), dimension(nColPars)                   :: PTF_Ks_clay
    real(dp), dimension(nColPars)                   :: PTF_Ks_curveSlope
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_forest
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_impervious
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_pervious
    ! directRunoff
    real(dp), dimension(nColPars)                   :: imperviousStorageCapacity
    ! PET0
    real(dp), dimension(nColPars)                   :: minCorrectionFactorPET
    real(dp), dimension(nColPars)                   :: maxCorrectionFactorPET
    real(dp), dimension(nColPars)                   :: aspectTresholdPET
    real(dp), dimension(nColPars)                   :: HargreavesSamaniCoeff
    real(dp), dimension(nColPars)                   :: PriestleyTaylorCoeff
    real(dp), dimension(nColPars)                   :: PriestleyTaylorLAIcorr
    real(dp), dimension(nColPars)                   :: canopyheigth_forest
    real(dp), dimension(nColPars)                   :: canopyheigth_impervious
    real(dp), dimension(nColPars)                   :: canopyheigth_pervious
    real(dp), dimension(nColPars)                   :: displacementheight_coeff
    real(dp), dimension(nColPars)                   :: roughnesslength_momentum_coeff
    real(dp), dimension(nColPars)                   :: roughnesslength_heat_coeff
    real(dp), dimension(nColPars)                   :: stomatal_resistance
    ! interflow
    real(dp), dimension(nColPars)                   :: interflowStorageCapacityFactor
    real(dp), dimension(nColPars)                   :: interflowRecession_slope
    real(dp), dimension(nColPars)                   :: fastInterflowRecession_forest
    real(dp), dimension(nColPars)                   :: slowInterflowRecession_Ks
    real(dp), dimension(nColPars)                   :: exponentSlowInterflow
    ! percolation
    real(dp), dimension(nColPars)                   :: rechargeCoefficient
    real(dp), dimension(nColPars)                   :: rechargeFactor_karstic
    real(dp), dimension(nColPars)                   :: gain_loss_GWreservoir_karstic
    ! routing moved to mRM
    ! neutrons
    real(dp), dimension(nColPars)                   :: Desilets_N0
    real(dp), dimension(nColPars)                   :: COSMIC_N0
    real(dp), dimension(nColPars)                   :: COSMIC_N1
    real(dp), dimension(nColPars)                   :: COSMIC_N2
    real(dp), dimension(nColPars)                   :: COSMIC_alpha0
    real(dp), dimension(nColPars)                   :: COSMIC_alpha1
    real(dp), dimension(nColPars)                   :: COSMIC_L30
    real(dp), dimension(nColPars)                   :: COSMIC_L31
    !
    integer(i4)                                     :: ii, iBasin, n_true_pars
    real(dp)                                        :: cellFactorRbyH            ! conversion factor L11 to L1
    !
    ! some dummy arrays for namelist read in (allocatables not allowed in namelists)
    character(256)                                  :: dummy
    character(256)                                  :: fname

    integer(i4),dimension(maxNoSoilHorizons)        :: soil_Depth           ! depth of the single horizons
    character(256), dimension(maxNoBasins)          :: dir_Morpho
    character(256), dimension(maxNoBasins)          :: dir_LCover
    character(256), dimension(maxNoBasins)          :: dir_Precipitation
    character(256), dimension(maxNoBasins)          :: dir_Temperature
    character(256), dimension(maxNoBasins)          :: dir_MinTemperature
    character(256), dimension(maxNoBasins)          :: dir_MaxTemperature
    character(256), dimension(maxNoBasins)          :: dir_NetRadiation
    character(256), dimension(maxNoBasins)          :: dir_windspeed
    character(256), dimension(maxNoBasins)          :: dir_absVapPressure
    character(256), dimension(maxNoBasins)          :: dir_ReferenceET
    character(256), dimension(maxNoBasins)          :: dir_Out
    character(256), dimension(maxNoBasins)          :: dir_RestartOut
    character(256), dimension(maxNoBasins)          :: dir_RestartIn
    character(256), dimension(maxNoBasins)          :: file_LatLon
    character(256), dimension(maxNoBasins)          :: dir_gridded_LAI           ! directory of gridded LAI data
    !                                                                            ! used when timeStep_LAI_input<0
    character(256), dimension(maxNoBasins)          :: dir_soil_moisture         ! soil moisture input
    character(256), dimension(maxNoBasins)          :: file_TWS                  ! total water storage input file
    character(256), dimension(maxNoBasins)          :: dir_neutrons              ! ground albedo neutron input

    !
  integer(i4)                                       :: nLCover_scene   ! given number of land cover scenes
  integer(i4),    dimension(maxNLCovers)            :: LCoverYearStart ! starting year LCover
  integer(i4),    dimension(maxNLCovers)            :: LCoverYearEnd   ! ending year LCover
  character(256), dimension(maxNLCovers)            :: LCoverfName     ! filename of Lcover file
    real(dp),       dimension(maxGeoUnit, nColPars) :: GeoParam                  ! geological parameters
    !
    real(dp)                                        :: jday_frac
    integer(i4),    dimension(maxNoBasins)          :: warming_Days
    type(period),   dimension(maxNoBasins)          :: eval_Per
    integer(i4),    dimension(maxNoBasins)          :: time_step_model_inputs
    !
    real(dp), dimension(maxNoBasins)                :: resolution_Hydrology
    real(dp), dimension(maxNoBasins)                :: resolution_Routing
    integer(i4), dimension(maxNoBasins)             :: L0Basin

    ! define namelists
    ! namelist directories
    namelist /directories_general/ dirConfigOut, dirCommonFiles, &
         dir_Morpho, dir_LCover,                                 &
         dir_Out, dir_RestartOut,                                &
         dir_RestartIn, file_LatLon
    namelist /directories_mHM/ inputFormat_meteo_forcings, &
         dir_Precipitation, dir_Temperature, dir_ReferenceET, dir_MinTemperature,    &
         dir_MaxTemperature, dir_absVapPressure, dir_windspeed,                      &
         dir_NetRadiation, dir_gridded_LAI
    ! optional data used for optimization
    namelist /optional_data/ &
        dir_soil_moisture, &
        nSoilHorizons_sm_input, &
        timeStep_sm_input, &
        file_TWS, &
        dir_neutrons
    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig/ timestep, iFlag_cordinate_sys, resolution_Hydrology, resolution_Routing, &
         L0Basin, optimize, optimize_restart, opti_method, opti_function, nBasins, read_restart,   &
         write_restart, perform_mpr
    ! namelist for time settings
    namelist /time_periods/ warming_Days, eval_Per, time_step_model_inputs
    ! namelist soil layering
    namelist /soilLayer/ iFlag_soilDB, tillageDepth, nSoilHorizons_mHM, soil_Depth
    ! namelist for land cover scenes
    namelist/LCover/fracSealed_cityArea,nLcover_scene,LCoverYearStart,LCoverYearEnd,LCoverfName
    ! namelist for LAI related data
    namelist /LAI_data_information/ inputFormat_gridded_LAI, timeStep_LAI_input
    ! namelist for pan evaporation
    namelist/panEvapo/evap_coeff
    ! namelist for night-day ratio of precipitation, referenceET and temperature
    namelist/nightDayRatio/fnight_prec,fnight_pet,fnight_temp
    ! namelsit process selection
    namelist /processSelection/ processCase
    ! namelist parameters
    namelist /interception1/ canopyInterceptionFactor
    namelist /snow1/snowTreshholdTemperature, degreeDayFactor_forest, degreeDayFactor_impervious, &
         degreeDayFactor_pervious, increaseDegreeDayFactorByPrecip, maxDegreeDayFactor_forest,    &
         maxDegreeDayFactor_impervious, maxDegreeDayFactor_pervious
    namelist/soilmoisture1/ orgMatterContent_forest, orgMatterContent_impervious, orgMatterContent_pervious, &
         PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant,              &
         PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant,                                            &
         PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope,                                                        &
         rootFractionCoefficient_forest, rootFractionCoefficient_impervious,                                 &
         rootFractionCoefficient_pervious, infiltrationShapeFactor
    namelist /directRunoff1/ imperviousStorageCapacity
    namelist /PET0/  minCorrectionFactorPET, maxCorrectionFactorPET, aspectTresholdPET
    ! Hargreaves-Samani
    namelist /PET1/  minCorrectionFactorPET, maxCorrectionFactorPET, aspectTresholdPET, HargreavesSamaniCoeff
    ! Priestely-Taylor
    namelist /PET2/  PriestleyTaylorCoeff, PriestleyTaylorLAIcorr
    ! Penman-Monteith
    namelist /PET3/  canopyheigth_forest, canopyheigth_impervious, canopyheigth_pervious, displacementheight_coeff, &
         roughnesslength_momentum_coeff, roughnesslength_heat_coeff, stomatal_resistance
    namelist /interflow1/ interflowStorageCapacityFactor, interflowRecession_slope, fastInterflowRecession_forest, &
         slowInterflowRecession_Ks, exponentSlowInterflow
    namelist /percolation1/ rechargeCoefficient, rechargeFactor_karstic, gain_loss_GWreservoir_karstic
    namelist /neutrons1/ Desilets_N0, COSMIC_N0, COSMIC_N1, COSMIC_N2, COSMIC_alpha0, COSMIC_alpha1, COSMIC_L30, COSMIC_L31
    !
    namelist /geoparameter/ GeoParam
    ! name list regarding output
    namelist/NLoutputResults/timeStep_model_outputs, outputFlxState
    ! namelist for optimization settings
    namelist/Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, sce_npg, sce_nps, mcmc_opti, mcmc_error_params

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet=.true.)

    !===============================================================
    !  Read namelist specifying the model configuration
    !===============================================================
    call position_nml('mainconfig', unamelist)
    read(unamelist, nml=mainconfig)

    if (nBasins .GT. maxNoBasins) then
       call message()
       call message('***ERROR: Number of basins is resticted to ', trim(num2str(maxNoBasins)),'!')
       stop
    end if
    ! allocate patharray sizes
    allocate(resolutionHydrology(nBasins))
    allocate(resolutionRouting(nBasins))
    allocate(L0_Basin(nBasins))
    allocate(dirMorpho(nBasins))
    allocate(dirLCover(nBasins))
    allocate(dirPrecipitation(nBasins))
    allocate(dirTemperature(nBasins))
    allocate(dirwindspeed(nBasins))
    allocate(dirabsVapPressure(nBasins))
    allocate(dirReferenceET(nBasins))
    allocate(dirMinTemperature(nBasins))
    allocate(dirMaxTemperature(nBasins))
    allocate(dirNetRadiation(nBasins))
    allocate(dirOut(nBasins))
    allocate(dirRestartOut(nBasins))
    allocate(dirRestartIn(nBasins))
    allocate(fileLatLon(nBasins))
    allocate(dirgridded_LAI(nBasins))
    allocate(dirSoil_Moisture(nBasins))
    allocate(dirNeutrons(nBasins))
    allocate(fileTWS(nBasins))
    !
    resolutionHydrology = resolution_Hydrology(1:nBasins)
    resolutionRouting   = resolution_Routing(1:nBasins)
    L0_Basin            = L0Basin(1:nBasins)
    !
    ! check for possible options
    if( .NOT. (iFlag_cordinate_sys == 0 .OR. iFlag_cordinate_sys == 1) ) then
       call message()
       call message('***ERROR: coordinate system for the model run should be 0 or 1')
       stop
    end if
    ! check for optimize and read restart
    if ( ( read_restart ) .and. ( optimize ) ) then
       call message()
       call message('***ERROR: cannot read states from restart file when optimizing')
       stop
    end if
    ! check for perform_mpr and read restart
    if ( ( read_restart ) .and. ( perform_mpr ) ) then
       call message()
       call message('***WARNING: MPR has not effect when reading states from restart file')
       call message('            MPR should only be activated if a land cover change accurs.')
    end if
    if ( ( .not. read_restart ) .and. ( .not. perform_mpr ) ) then
       call message()
       call message('***ERROR: cannot omit mpr when read_restart is set to .false.')
       stop
    end if

    ! allocate time periods
    allocate(simPer(nBasins))
    allocate(evalPer(nBasins))
    allocate(warmingDays(nBasins))
    allocate(warmPer(nBasins))
    allocate(timestep_model_inputs(nBasins))

    !===============================================================
    !  read simulation time periods incl. warming days
    !===============================================================
    call position_nml('time_periods', unamelist)
    read(unamelist, nml=time_periods)
    warmingDays = warming_Days(1:nBasins)
    evalPer     = eval_Per(1:nBasins)
    timestep_model_inputs = time_step_model_inputs(1:nBasins)

    ! consistency check for timestep_model_inputs
    if (       any( timestep_model_inputs .ne. 0 ) .and. &
         .not. all( timestep_model_inputs .ne. 0 ) ) then
       call message()
       call message('***ERROR: timestep_model_inputs either have to be all zero or all non-zero')
       stop
    end if
    ! check for optimzation and timestep_model_inputs options
    if ( optimize .and. ( any(timestep_model_inputs .ne. 0) ) ) then
       call message()
       call message('***ERROR: optimize and chunk read is switched on! (set timestep_model_inputs to zero)')
       stop
    end if

    !===============================================================
    !  determine simulation time period incl. warming days for each
    !  basin
    !===============================================================
    do ii = 1, nBasins
       ! julain days for evaluation period
       jday_frac = date2dec(dd=evalPer(ii)%dStart, mm=evalPer(ii)%mStart, yy=evalPer(ii)%yStart)
       evalPer(ii)%julStart = nint(jday_frac)

       jday_frac = date2dec(dd=evalPer(ii)%dEnd, mm=evalPer(ii)%mEnd, yy=evalPer(ii)%yEnd)
       evalPer(ii)%julEnd  = nint(jday_frac, i4 )

       ! determine warming period
       warmPer(ii)%julStart = evalPer(ii)%julStart - warmingDays(ii)
       warmPer(ii)%julEnd   = evalPer(ii)%julStart - 1

       jday_frac = real(warmPer(ii)%julStart,dp)
       call dec2date(jday_frac, dd=warmPer(ii)%dStart, mm=warmPer(ii)%mStart, yy=warmPer(ii)%yStart)

       jday_frac = real(warmPer(ii)%julEnd,dp)
       call dec2date(jday_frac, dd=warmPer(ii)%dEnd,   mm=warmPer(ii)%mEnd,   yy=warmPer(ii)%yEnd  )

       ! sumulation Period = warming Period + evaluation Period
       simPer(ii)%dStart   = warmPer(ii)%dStart
       simPer(ii)%mStart   = warmPer(ii)%mStart
       simPer(ii)%yStart   = warmPer(ii)%yStart
       simPer(ii)%julStart = warmPer(ii)%julStart
       simPer(ii)%dEnd     = evalPer(ii)%dEnd
       simPer(ii)%mEnd     = evalPer(ii)%mEnd
       simPer(ii)%yEnd     = evalPer(ii)%yEnd
       simPer(ii)%julEnd   = evalPer(ii)%julEnd
    end do

    !===============================================================
    !  Read namelist for mainpaths
    !===============================================================
    call position_nml('directories_general', unamelist)
    read(unamelist, nml=directories_general)
    call position_nml('directories_mHM', unamelist)
    read(unamelist, nml=directories_mHM)

    dirMorpho         = dir_Morpho(1:nBasins)
    dirLCover         = dir_LCover(1:nBasins)
    dirPrecipitation  = dir_Precipitation(1:nBasins)
    dirTemperature    = dir_Temperature(1:nBasins)
    dirReferenceET    = dir_ReferenceET(1:nBasins)
    dirMinTemperature = dir_MinTemperature(1:nBasins)
    dirMaxTemperature = dir_MaxTemperature(1:nBasins)
    dirNetRadiation   = dir_NetRadiation(1:nBasins)
    dirwindspeed      = dir_windspeed(1:nBasins)
    dirabsVapPressure = dir_absVapPressure(1:nBasins)
    dirOut            = dir_Out(1:nBasins)
    dirRestartOut     = dir_RestartOut(1:nBasins)
    dirRestartIn      = dir_RestartIn(1:nBasins)
    fileLatLon        = file_LatLon(1:nBasins)
    dirgridded_LAI    = dir_gridded_LAI(1:nBasins)
    
    ! counter checks -- soil horizons
    if (nSoilHorizons_mHM .GT. maxNoSoilHorizons) then
       call message()
       call message('***ERROR: Number of soil horizons is resticted to ', trim(num2str(maxNoSoilHorizons)),'!')
       stop
    end if

    !===============================================================
    ! Read soil layering information
    !===============================================================
    call position_nml('soillayer', unamelist)
    read(unamelist, nml=soillayer)

    allocate( HorizonDepth_mHM(nSoilHorizons_mHM) )
    HorizonDepth_mHM(:) = 0.0_dp
   
    if( iFlag_soilDB .eq. 0 ) then
       ! classical mhm soil database
       HorizonDepth_mHM(1:nSoilHorizons_mHM-1) = soil_Depth(1:nSoilHorizons_mHM-1)
    else if( iFlag_soilDB .eq. 1 ) then
       ! for option-1 where horizon specific information are taken into consideration
       HorizonDepth_mHM(1:nSoilHorizons_mHM) = soil_Depth(1:nSoilHorizons_mHM)
    else
       call message()
       call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
       stop
    end if

    ! some consistency checks for the specification of the tillage depth
    if(iFlag_soilDB .eq. 1) then
       if( count( abs(HorizonDepth_mHM(:) - tillageDepth) .lt. eps_dp )  .eq. 0 ) then
          call message()
          call message('***ERROR: Soil tillage depth must conform with one of the specified horizon (lower) depth.')
          stop
       end if
    end if
   
    !===============================================================
    !  Read namelist of optional input data
    !===============================================================
    call position_nml('optional_data', unamelist)
    read(unamelist, nml=optional_data)
    dirSoil_moisture = dir_Soil_moisture(1:nBasins)
    if ( nSoilHorizons_sm_input .GT. nSoilHorizons_mHM ) then
       call message()
       call message('***ERROR: Number of soil horizons representative for input soil moisture exceeded')
       call message('          defined number of soil horizions: ', trim(num2str(maxNoSoilHorizons)),'!')
       stop
    end if
    dirNeutrons = dir_neutrons(1:nBasins)
    timeStep_neutrons_input = -1 ! daily, hard-coded, to be flexibilized

    !===============================================================
    ! Read evaluation basin average TWS data
    !===============================================================
    fileTWS = file_TWS (1:nBasins)

    ! if (opti_function .EQ. 15) then !OROROR
    allocate(basin_avg_TWS_obs%basinId(nBasins)); basin_avg_TWS_obs%basinId = nodata_i4
    allocate(basin_avg_TWS_obs%fName  (nBasins)); basin_avg_TWS_obs%fName(:)= num2str(nodata_i4)
    
    do iBasin = 1, nBasins
       if (trim(fileTWS(iBasin)) .EQ. trim(num2str(nodata_i4))) then
          call message()
          call message('***ERROR: mhm.nml: Filename of evaluation TWS data ', &
               ' for subbasin ', trim(adjustl(num2str(iBasin))), &
               ' is not defined!')
          call message('          Error occured in namelist: evaluation_tws')
          stop
       end if
       
       basin_avg_TWS_obs%basinId(iBasin) = iBasin
       basin_avg_TWS_obs%fname(iBasin)   = trim(file_TWS(iBasin))
    end do

    !===============================================================
    ! Read process selection list
    !===============================================================
    call position_nml('processselection', unamelist)
    read(unamelist, nml=processSelection)

    !===============================================================
    ! Read land cover information
    !===============================================================
    call position_nml('LCover', unamelist)
    read(unamelist, nml=LCover)

    !===============================================================
    ! Read LAI related information
    !===============================================================
    call position_nml('LAI_data_information', unamelist)
    read(unamelist, nml=LAI_data_information)

    if (timeStep_LAI_input .ne. 0) then
       if ( (timeStep_LAI_input .ne. -1) .and. (trim(inputFormat_gridded_LAI) .eq. 'bin') ) then
          call message()
          call message('***ERROR: Gridded LAI input in bin format must be daily.')
          stop
       end if
       if (timeStep_LAI_input > 0) then
          call message()
          call message('***ERROR: timeStep_LAI_input must be <= 0.')
          stop
       end if
    end if

    !===============================================================
    ! Read night-day ratios and pan evaporation
    !===============================================================
    ! Evap. coef. for free-water surfaces
    call position_nml('panEvapo', unamelist)
    read(unamelist, nml=panEvapo)
    ! namelist for night-day ratio of precipitation, referenceET and temperature
    call position_nml('nightDayRatio', unamelist)
    read(unamelist, nml=nightDayRatio)
    !
    fday_prec =  1.0_dp - fnight_prec
    fday_pet  =  1.0_dp - fnight_pet
    fday_temp = -1.0_dp * fnight_temp

    !===============================================================
    !  determine land cover periods
    !===============================================================
    ! countercheck if land cover covers simulation period
    if (LCoverYearStart(1) .GT. minval(evalPer(1:nBasins)%yStart) ) then
       call message()
       call message('***ERROR: Land cover for warming period is missing!')
       call message('   FILE: mhm.nml, namelist: LCover')
       call message('   SimStart   : ', trim(num2str(minval(evalPer(1:nBasins)%yStart))))
       call message('   LCoverStart: ', trim(num2str(LCoverYearStart(1))))
       stop
    end if
    if (LCoverYearEnd(nLcover_scene) .LT. maxval(evalPer(1:nBasins)%yEnd) ) then
       call message()
       call message('***ERROR: Land cover period shorter than modelling period!')
       call message('   FILE: mhm.nml, namelist: LCover')
       call message('   SimEnd   : ', trim(num2str(maxval(evalPer(1:nBasins)%yEnd))))
       call message('   LCoverEnd: ', trim(num2str(LCoverYearEnd(nLcover_scene))))
       stop
    end if
    !
    allocate(LCYearId(minval(simPer(1:nBasins)%yStart):maxval(simPer(1:nBasins)%yEnd),nBasins))
    LCYearId = nodata_i4
    do iBasin = 1, nBasins
       do ii = 1, nLcover_scene
          ! land cover before model period                        ! land cover after model period
          if ((LCoverYearEnd(ii)        .LT. evalPer(iBasin)%yStart) .OR. &
               (LCoverYearStart(ii)      .GT. evalPer(iBasin)%yEnd)) then
             cycle
          else if ((LCoverYearStart(ii) .LE. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii)   .GE. evalPer(iBasin)%yEnd)) then
             LCyearId(simPer(iBasin)%yStart:simPer(iBasin)%yEnd, iBasin) = ii
             exit
          else if ((LCoverYearStart(ii) .LE. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii)   .LT. evalPer(iBasin)%yEnd)) then
             LCyearId(simPer(iBasin)%yStart:LCoverYearEnd(ii), iBasin) = ii
          else if ((LCoverYearStart(ii) .GT. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii)   .GE. evalPer(iBasin)%yEnd)) then
             LCyearId(LCoverYearStart(ii):simPer(iBasin)%yEnd, iBasin) = ii
          else
             LCyearId(LCoverYearStart(ii):LCoverYearEnd(ii), iBasin) = ii
          end if
       end do
    end do
    !
    ! correct number of input land cover scenes to number of needed scenes
    nLCoverScene = maxval(LCyearId, mask = (LCyearId .gt. nodata_i4) ) - &
         minval(LCyearId, mask = (LCyearId .gt. nodata_i4) ) + 1
    ! put land cover scenes to corresponding file name and LuT
    allocate(LCfilename(nLCoverScene))
    LCfilename(:) = LCoverfName( minval(LCyearId, mask = ( LCyearId .gt. nodata_i4 ) ) : &
         maxval(LCyearId))
    ! update the ID's
    ! use next line because of Intel11 bug: LCyearId = LCyearId - minval(LCyearId) + 1
    LCyearId(:,:) = LCyearId(:,:) - minval(LCyearId, mask = ( LCyearId .gt. nodata_i4 ) ) + 1
    !
    if ( maxval( simPer(1:nBasins)%julStart ) .eq. minval( simPer(1:nBasins)%julStart) .and. &
         maxval( simPer(1:nBasins)%julEnd   ) .eq. minval( simPer(1:nBasins)%julEnd  ) ) then
       if (any(LCyearId .EQ. nodata_i4)) then
          call message()
          call message('***ERROR: Intermidiate land cover period is missing!')
          call message('   FILE: mhm.nml, namelist: LCover')
          stop
       end if
    else
       call message()
       call message('***WARNING: No check on missing land cover period is performed!')
    end if
    !
    !===============================================================
    ! check matching of resolutions: hydrology, forcing and routing
    !===============================================================
    !
    do ii = 1, nBasins
       cellFactorRbyH = resolutionRouting(ii) / resolutionHydrology(ii)
       call message()
       call message('Basin ', trim(adjustl(num2str(ii))), ': ')
       call message('resolution Hydrology (basin ', trim(adjustl(num2str(ii))), ')     = ', &
            trim(adjustl(num2str(resolutionHydrology(ii)))))
       call message('resolution Routing (basin ', trim(adjustl(num2str(ii))), ')       = ', &
            trim(adjustl(num2str(resolutionRouting(ii)))))
       !
       if(       nint(cellFactorRbyH * 100.0_dp) .eq. 100) then
          call message()
          call message('Resolution of routing and hydrological modeling are equal!')

       else if ( nint(cellFactorRbyH * 100.0_dp) .lt. 100) then
          call message()
          call message('***ERROR: Resolution of routing is smaller than hydrological model resolution!')
          call message('   FILE: mhm.nml, namelist: mainconfig, variable: resolutionRouting')
          STOP

       else if ( nint(cellFactorRbyH * 100.0_dp) .gt. 100) then
          if( nint(mod(cellFactorRbyH, 2.0_dp) * 100.0_dp) .ne. 0) then
             call message()
             call message('***ERROR: Resolution of routing is not a multiple of hydrological model resolution!')
             call message('   FILE: mhm.nml, namelist: mainconfig, variable: resolutionRouting')
             STOP
          end if
          !
          call message()
          call message('Resolution of routing is bigger than hydrological model resolution by ', &
               trim(adjustl(num2str(nint(cellFactorRbyH)))), ' times !')
       end if
       !
    end do
    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    call open_nml(file_namelist_param, unamelist_param, quiet=.true.)
    processMatrix = 0_i4
    ! decide which parameters to read depending on specified processes

    ! Process 1 - interception
    select case (processCase(1))
       ! 1 - maximum Interception
    case(1)
       call position_nml('interception1', unamelist_param)
       read(unamelist_param, nml=interception1)

       processMatrix(1, 1) = processCase(1)
       processMatrix(1, 2) = 1_i4
       processMatrix(1, 3) = 1_i4
       call append(global_parameters, reshape(canopyInterceptionFactor,(/1, nColPars/)))

       call append(global_parameters_name, (/  &
            'canopyInterceptionFactor'/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "interception1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "interception" does not exist!')
       stop
    end select

    ! Process 2 - snow
    select case (processCase(2))
       ! 1 - degree-day approach
    case(1)
       call position_nml('snow1', unamelist_param)
       read(unamelist_param, nml=snow1)

       processMatrix(2, 1) = processCase(2)
       processMatrix(2, 2) = 8_i4
       processMatrix(2, 3) = sum(processMatrix(1:2, 2))
       call append(global_parameters, reshape(snowTreshholdTemperature,        (/1, nColPars/)))
       call append(global_parameters, reshape(degreeDayFactor_forest,          (/1, nColPars/)))
       call append(global_parameters, reshape(degreeDayFactor_impervious,      (/1, nColPars/)))
       call append(global_parameters, reshape(degreeDayFactor_pervious,        (/1, nColPars/)))
       call append(global_parameters, reshape(increaseDegreeDayFactorByPrecip, (/1, nColPars/)))
       call append(global_parameters, reshape(maxDegreeDayFactor_forest,       (/1, nColPars/)))
       call append(global_parameters, reshape(maxDegreeDayFactor_impervious,   (/1, nColPars/)))
       call append(global_parameters, reshape(maxDegreeDayFactor_pervious,     (/1, nColPars/)))

       call append(global_parameters_name, (/  &
            'snowTreshholdTemperature       ', &
            'degreeDayFactor_forest         ', &
            'degreeDayFactor_impervious     ', &
            'degreeDayFactor_pervious       ', &
            'increaseDegreeDayFactorByPrecip', &
            'maxDegreeDayFactor_forest      ', &
            'maxDegreeDayFactor_impervious  ', &
            'maxDegreeDayFactor_pervious    '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "snow1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "snow" does not exist!')
       stop
    end select

    ! Process 3 - soilmoisture
    select case (processCase(3))
       ! 1 - bucket approach, Brooks-Corey like
    case(1)
       call position_nml('soilmoisture1', unamelist_param)
       read(unamelist_param, nml=soilmoisture1)
       processMatrix(3, 1) = processCase(3)
       processMatrix(3, 2) = 17_i4
       processMatrix(3, 3) = sum(processMatrix(1:3, 2))
       call append(global_parameters, reshape(orgMatterContent_forest,     (/1, nColPars/)))
       call append(global_parameters, reshape(orgMatterContent_impervious, (/1, nColPars/)))
       call append(global_parameters, reshape(orgMatterContent_pervious,   (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_lower66_5_constant,      (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_lower66_5_clay,          (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_lower66_5_Db,            (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_higher66_5_constant,     (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_higher66_5_clay,         (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_higher66_5_Db,           (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_Ks_constant,             (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_Ks_sand,                 (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_Ks_clay,                 (/1, nColPars/)))
       call append(global_parameters, reshape(PTF_Ks_curveSlope,           (/1, nColPars/)))
       call append(global_parameters, reshape(rootFractionCoefficient_forest,     (/1, nColPars/)))
       call append(global_parameters, reshape(rootFractionCoefficient_impervious, (/1, nColPars/)))
       call append(global_parameters, reshape(rootFractionCoefficient_pervious,   (/1, nColPars/)))
       call append(global_parameters, reshape(infiltrationShapeFactor,     (/1, nColPars/)))

       call append(global_parameters_name, (/     &
            'orgMatterContent_forest           ', &
            'orgMatterContent_impervious       ', &
            'orgMatterContent_pervious         ', &
            'PTF_lower66_5_constant            ', &
            'PTF_lower66_5_clay                ', &
            'PTF_lower66_5_Db                  ', &
            'PTF_higher66_5_constant           ', &
            'PTF_higher66_5_clay               ', &
            'PTF_higher66_5_Db                 ', &
            'PTF_Ks_constant                   ', &
            'PTF_Ks_sand                       ', &
            'PTF_Ks_clay                       ', &
            'PTF_Ks_curveSlope                 ', &
            'rootFractionCoefficient_forest    ', &
            'rootFractionCoefficient_impervious', &
            'rootFractionCoefficient_pervious  ', &
            'infiltrationShapeFactor           '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "soilmoisture1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "soilmoisture" does not exist!')
       stop
    end select

    ! Process 4 - sealed area directRunoff
    select case (processCase(4))
       ! 1 - bucket exceedance approach
    case(1)
       call position_nml('directRunoff1', unamelist_param)
       read(unamelist_param, nml=directRunoff1)
       processMatrix(4, 1) = processCase(4)
       processMatrix(4, 2) = 1_i4
       processMatrix(4, 3) = sum(processMatrix(1:4, 2))
       call append(global_parameters, reshape(imperviousStorageCapacity, (/1, nColPars/)))

       call append(global_parameters_name, (/'imperviousStorageCapacity'/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "directRunoff1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "directRunoff" does not exist!')
       stop
    end select

    ! Process 5 - potential evapotranspiration (PET)
    select case (processCase(5))
    case(0) ! 0 - PET is input, correct PET by aspect
       call position_nml('PET0', unamelist_param)
       read(unamelist_param, nml=PET0)
       processMatrix(5, 1) = processCase(5)
       processMatrix(5, 2) = 3_i4
       processMatrix(5, 3) = sum(processMatrix(1:5, 2))
       call append(global_parameters, reshape(minCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(maxCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(aspectTresholdPET,                  (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'minCorrectionFactorPET ', &
            'maxCorrectionFactorPET ', &
            'aspectTresholdPET      '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "PET0" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case(1) ! 1 - Hargreaves-Samani method (HarSam) - additional input needed: Tmin, Tmax
       call position_nml('PET1', unamelist_param)
       read(unamelist_param, nml=PET1)
       processMatrix(5, 1) = processCase(5)
       processMatrix(5, 2) = 4_i4
       processMatrix(5, 3) = sum(processMatrix(1:5, 2))
       call append(global_parameters, reshape(minCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(maxCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(aspectTresholdPET,                  (/1, nColPars/)))
       call append(global_parameters, reshape(HargreavesSamaniCoeff,              (/1, nColPars/)))
       call append(global_parameters_name, (/ &
            'minCorrectionFactorPET', &
            'maxCorrectionFactorPET', &
            'aspectTresholdPET     ', &
            'HargreavesSamaniCoeff '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "PET1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case(2) ! 2 - Priestley-Taylor method (PrieTay) - additional input needed: net_rad
       ! check which LAI input is specified
       if (timeStep_LAI_input .NE. 0) then
          call message('***ERROR: The specified option of process 5 does only work with LAI from LUT.')
          call message('          For process 5 the options 0 and 1 work with timeStep_LAI_input unequal to 0.')
          stop
       end if

       call position_nml('PET2', unamelist_param)
       read(unamelist_param, nml=PET2)
       processMatrix(5, 1) = processCase(5)
       processMatrix(5, 2) = 2_i4
       processMatrix(5, 3) = sum(processMatrix(1:5, 2))
       call append(global_parameters, reshape(PriestleyTaylorCoeff,               (/1, nColPars/)))
       call append(global_parameters, reshape(PriestleyTaylorLAIcorr,             (/1, nColPars/)))
       call append(global_parameters_name, (/ &
            'PriestleyTaylorCoeff  ', &
            'PriestleyTaylorLAIcorr'/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "PET2" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case(3) ! 3 - Penman-Monteith method - additional input needed: net_rad, abs. vapour pressue, windspeed
       ! check which LAI input is specified
       if (timeStep_LAI_input .NE. 0) then
          call message('***ERROR: The specified option of process 5 does only work with LAI from LUT.')
          call message('          For process 5 the options 0 and 1 work with timeStep_LAI_input unequal to 0.')
          stop
       end if

       call position_nml('PET3', unamelist_param)
       read(unamelist_param, nml=PET3)
       processMatrix(5, 1) = processCase(5)
       processMatrix(5, 2) = 7_i4
       processMatrix(5, 3) = sum(processMatrix(1:5, 2))

       call append(global_parameters, reshape(canopyheigth_forest,                (/1, nColPars/)))
       call append(global_parameters, reshape(canopyheigth_impervious,            (/1, nColPars/)))
       call append(global_parameters, reshape(canopyheigth_pervious,              (/1, nColPars/)))
       call append(global_parameters, reshape(displacementheight_coeff,           (/1, nColPars/)))
       call append(global_parameters, reshape(roughnesslength_momentum_coeff,     (/1, nColPars/)))
       call append(global_parameters, reshape(roughnesslength_heat_coeff,         (/1, nColPars/)))
       call append(global_parameters, reshape(stomatal_resistance,                (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'canopyheigth_forest           ', &
            'canopyheigth_impervious       ', &
            'canopyheigth_pervious         ', &
            'displacementheight_coeff      ', &
            'roughnesslength_momentum_coeff', &
            'roughnesslength_heat_coeff    ', &
            'stomatal_resistance           '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "PET3" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "actualET" does not exist!')
       stop
    end select


    ! Process 6 - interflow
    select case (processCase(6))
       ! 1 - parallel soil reservoir approach
    case(1)
       call position_nml('interflow1', unamelist_param)
       read(unamelist_param, nml=interflow1)
       processMatrix(6, 1) = processCase(6)
       processMatrix(6, 2) = 5_i4
       processMatrix(6, 3) = sum(processMatrix(1:6, 2))
       call append(global_parameters, reshape(interflowStorageCapacityFactor, (/1, nColPars/)))
       call append(global_parameters, reshape(interflowRecession_slope,       (/1, nColPars/)))
       call append(global_parameters, reshape(fastInterflowRecession_forest,  (/1, nColPars/)))
       call append(global_parameters, reshape(slowInterflowRecession_Ks,      (/1, nColPars/)))
       call append(global_parameters, reshape(exponentSlowInterflow,          (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'interflowStorageCapacityFactor', &
            'interflowRecession_slope      ', &
            'fastInterflowRecession_forest ', &
            'slowInterflowRecession_Ks     ', &
            'exponentSlowInterflow         '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "interflow1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "interflow" does not exist!')
       stop
    end select

    ! Process 7 - percolation
    select case (processCase(7))
       ! 1 - GW layer is assumed as bucket
    case(1)
       call position_nml('percolation1', unamelist_param)
       read(unamelist_param, nml=percolation1)
       processMatrix(7, 1) = processCase(7)
       processMatrix(7, 2) = 3_i4
       processMatrix(7, 3) = sum(processMatrix(1:7, 2))
       call append(global_parameters, reshape(rechargeCoefficient,           (/1, nColPars/)))
       call append(global_parameters, reshape(rechargeFactor_karstic,        (/1, nColPars/)))
       call append(global_parameters, reshape(gain_loss_GWreservoir_karstic, (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'rechargeCoefficient          ', &
            'rechargeFactor_karstic       ', &
            'gain_loss_GWreservoir_karstic'/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "percolation1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "percolation" does not exist!')
       stop
    end select

    ! Process 8 - routing
    select case (processCase(8))
    case(0)
       ! 0 - deactivated
       call message()
       call message('***CAUTION: Routing is deativated! ')

       processMatrix(8, 1) = processCase(8)
       processMatrix(8, 2) = 0_i4
       processMatrix(8, 3) = sum(processMatrix(1:8, 2))
    case(1)
       ! parameter values and names are set in mRM
       ! 1 - Muskingum approach
#ifndef mrm2mhm
       call message('***ERROR processCase(8) equals 1, but mrm2mhm preprocessor flag is not given in Makefile')
       stop
#endif
       processMatrix(8, 1) = processCase(8)
       processMatrix(8, 2) = 5_i4
       processMatrix(8, 3) = sum(processMatrix(1:8, 2))
       call append(global_parameters, dummy_2d_dp)
       call append(global_parameters_name, (/'dummy', 'dummy', 'dummy', 'dummy', 'dummy'/))

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "routing" does not exist!')
       stop
    end select

    !===============================================================
    ! Geological formations
    !===============================================================

    ! read in of number of geological formations
    fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_geolut))
    open( unit=ugeolut, file=fname, action='read', status='old')
    read(ugeolut, *) dummy, nGeoUnits
    close(ugeolut)
    dummy = dummy//''   ! only to avoid warning

    ! Process 9 - geoparameter
    select case (processCase(9))
    case(1)
       ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
       call position_nml('geoparameter', unamelist_param)
       GeoParam = nodata_dp
       read(unamelist_param, nml=geoparameter)

       ! for geology parameters
       processMatrix(9,1) = processCase(9)
       processMatrix(9,2) = nGeoUnits
       processMatrix(9,3) = sum(processMatrix(1:9, 2))

       call append(global_parameters, GeoParam(1:nGeoUnits,:))

       ! create names
       do ii=1, nGeoUnits
          dummy = 'GeoParam('//trim(adjustl(num2str(ii)))//',:)'
          call append(global_parameters_name, (/ trim(dummy) /))
       end do

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "geoparameter" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "geoparameter" does not exist!')
       stop
    end select

    ! Process 10 - neutrons
    !   0 - deactivated
    !   1 - inverse N0 based on Desilets et al. 2010
    !   2 - COSMIC forward operator by Shuttlworth et al. 2013
    if (processCase(10) .gt. 0) then

       call position_nml('neutrons1', unamelist_param)
       read(unamelist_param, nml=neutrons1)

       processMatrix(10, 1) = processCase(10)
       processMatrix(10, 2) = 8_i4
       processMatrix(10, 3) = sum(processMatrix(1:10, 2))
       call append(global_parameters, reshape(Desilets_N0,   (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_N0,     (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_N1,     (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_N2,     (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_alpha0, (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_alpha1, (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_L30,    (/1, nColPars/)))
       call append(global_parameters, reshape(COSMIC_L31,    (/1, nColPars/)))

       call append(global_parameters_name, (/  &
            'Desilets_N0   ', &
            'COSMIC_N0     ', &
            'COSMIC_N1     ', &
            'COSMIC_N2     ', &
            'COSMIC_alpha0 ', &
            'COSMIC_alpha1 ', &
            'COSMIC_L30    ', &
            'COSMIC_L31    '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "neutrons1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if
    else
       call message(' INFO: Process (10, neutrons) is deactivated, so output will be suppressed.')
       ! this is done below, where nml_output is read
       processMatrix(10, 1) = processCase(10)
       processMatrix(10, 2) = 0_i4
       processMatrix(10, 3) = sum(processMatrix(1:10, 2))
    end if

    call close_nml(unamelist_param)

    !===============================================================
    ! Settings for Optimization
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet=.true.)
    ! namelist for Optimization settings
    call position_nml('Optimization', unamelist)
    read(unamelist, nml=Optimization)
    call close_nml(unamelist)
    ! check and set default values
    if (nIterations .le. 0_i4) then
       call message('Number of iterations for Optimization (nIterations) must be greater than zero')
       stop
    end if
    if (dds_r .lt. 0.0_dp .or. dds_r .gt. 1.0_dp) then
       call message('dds_r must be between 0.0 and 1.0')
       stop
    end if
    if (sce_ngs .lt. 1_i4) then
       call message ('number of complexes in SCE (sce_ngs) must be at least 1')
       stop
    end if
    ! number of points in each complex: default = 2n+1
    if (sce_npg .lt. 0_i4) then
       n_true_pars = count(nint(global_parameters(:,4)) .eq. 1)
       sce_npg = 2 * n_true_pars + 1_i4
    end if
    ! number of points in each sub-complex: default = n+1
    if (sce_nps .lt. 0_i4) then
       n_true_pars = count(nint(global_parameters(:,4)) .eq. 1)
       sce_nps = n_true_pars + 1_i4
    end if
    if (sce_npg .lt. sce_nps) then
       call message ('number of points per complex (sce_npg) must be greater or')
       call message ('equal number of points per sub-complex (sce_nps)')
       stop
    end if

    call close_nml(unamelist)

    !===============================================================
    ! Read output specifications for mHM
    !===============================================================
    call open_nml(file_defOutput, udefOutput, quiet=.true.)
    outputFlxState = .FALSE.
    call position_nml('NLoutputResults', udefOutput)
    read(udefOutput, nml=NLoutputResults)
    call close_nml(udefOutput)

    call message( '' )
    call message( 'Following output will be written:' )
    call message( '  STATES:' )
    if (outputFlxState(1)) then
       call message( '    interceptional storage                      (L1_inter)     [mm]')
    end if
    if (outputFlxState(2)) then
       call message( '    height of snowpack                          (L1_snowpack)  [mm]')
    end if
    if (outputFlxState(3)) then
       call message( '    soil water content in the single layers     (L1_soilMoist) [mm]')
    end if
    if (outputFlxState(4)) then
       call message( '    volumetric soil moisture in the single layers              [mm/mm]')
    end if
    if (outputFlxState(5)) then
       call message( '    mean volum. soil moisture averaged over all soil layers    [mm/mm]')
    end if
    if (outputFlxState(6)) then
       call message( '    waterdepth in reservoir of sealed areas     (L1_sealSTW)   [mm]')
    end if
    if (outputFlxState(7)) then
       call message( '    waterdepth in reservoir of unsat. soil zone (L1_unsatSTW)  [mm]')
    end if
    if (outputFlxState(8)) then
       call message( '    waterdepth in reservoir of sat. soil zone   (L1_satSTW)    [mm]')
    end if
    if (processCase(10) .eq. 0) outputFlxState(18) = .false. ! suppress output if process is off
    if (outputFlxState(18)) then
       call message( '    ground albedo neutrons                      (L1_neutrons)  [cph]')
    end if

    call message( '  FLUXES:' )
    if (outputFlxState(9)) then
       call message( '    actual evapotranspiration aET      (L1_pet)                [mm/T]')
    end if
    if (outputFlxState(10)) then
       call message( '    total discharge generated per cell (L1_total_runoff)       [mm/T]')
    end if
    if (outputFlxState(11)) then
       call message( '    direct runoff generated per cell   (L1_runoffSeal)         [mm/T]')
    end if
    if (outputFlxState(12)) then
       call message( '    fast interflow generated per cell  (L1_fastRunoff)         [mm/T]')
    end if
    if (outputFlxState(13)) then
       call message( '    slow interflow generated per cell  (L1_slowRunoff)         [mm/T]')
    end if
    if (outputFlxState(14)) then
       call message( '    baseflow generated per cell        (L1_baseflow)           [mm/T]')
    end if
    if (outputFlxState(15)) then
       call message( '    groundwater recharge               (L1_percol)             [mm/T]')
    end if
    if (outputFlxState(16)) then
       call message( '    infiltration                       (L1_infilSoil)          [mm/T]')
    end if
    call message( '' )
    call message( 'FINISHED reading config' )

    ! warning message
    if ( any(outputFlxState) .and. optimize ) then
       call message( 'WARNING: FLUXES and STATES netCDF will be not written since optimization flag is TRUE ' )
    end if

  end subroutine read_config

  ! --------------------------------------------------------------------------------
  ! private funtions and subroutines
  ! --------------------------------------------------------------------------------

  function in_bound(params)
    real(dp), dimension(:,:), intent(in) :: params ! parameter:
    !                                              !   col_1=Lower bound,
    !                                              !   col_2=Upper bound
    !                                              !   col_3=initial
    logical :: in_bound

    if ( any(params(:,3) .lt. params(:,1)) .or. any(params(:,3) .gt. params(:,2)) ) then
       in_bound=.false.
    else
       in_bound=.true.
    end if

  end function in_bound


END MODULE mo_read_config
