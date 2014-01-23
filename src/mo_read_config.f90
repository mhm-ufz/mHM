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

  subroutine read_config()

    use mo_julian,           only: date2dec, dec2date
    use mo_append,           only: append
    use mo_message,          only: message 
    use mo_string_utils,     only: num2str
    use mo_nml,              only: open_nml, close_nml, position_nml
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
         resolutionHydrology, resolutionRouting,            & ! resolutions of hydrology and routing 
         dirMorpho, dirLCover,                              & ! input directory of morphological
         dirGauges,                                         & ! and discharge files
         dirPrecipitation, dirTemperature, dirReferenceET,  & ! directory of meteo input
         inputFormat_meteo_forcings,                        & ! input format either bin or nc
         dirLatLon,                                         & ! directory of latitude and longitude files
         dirConfigOut,                                      & ! configuration run output directory
         dirCommonFiles_In,                                 & ! directory where common files are located
         dirOut,                                            & ! output directory basin wise 
         dirRestartOut,                                     & ! output directory of restart file basin wise
         dirRestartIn,                                      & ! input directory of restart file basin wise
         optimize,                                          & ! if mhm runs in optimization mode or not
         opti_method,                                       & ! optimization algorithm used    
         opti_function,                                     & ! objective function to be optimized
         nIterations,                                       & ! number of iterations in optimization
         seed,                                              & ! seed used for optimization
         dds_r,                                             & ! DDS: perturbation rate
         sa_temp,                                           & ! SA: initial temperature
         sce_ngs, sce_npg, sce_nps,                         & ! SCE: # complexes, # points per complex,
         !                                                    !      # points per subcomplex
         HorizonDepth_mHM, nSoilHorizons_mHM, tillageDepth, & ! soil horizons info for mHM
         fracSealed_cityArea, nLcover_scene,                & ! land cover information
         LCfilename, LCyearId,                              & ! 
         nBasins,                                           & ! number of basins
         restart_flag_states_read,                          & ! flag reading restart (state variables)
         restart_flag_states_write,                         & ! flag writing restart (state variables)
         restart_flag_config_read,                          & ! flag reading restart (config variables)
         restart_flag_config_write,                         & ! flag writing restart (config variables)
         warmingDays, warmPer,                              & ! warming days and warming period
         evalPer, simPer,                                   & ! model eval. & sim. periods  
         !                                                    ! (sim. = wrm. + eval.)
         evap_coeff,                                        & ! pan evaporation
         fday_prec, fnight_prec, fday_pet,                  & ! day-night fraction
         fnight_pet, fday_temp, fnight_temp,                & ! day-night fraction
         nProcesses, processMatrix,                         & ! process configuration
         nGeoUnits,                                         & ! number of geological classes
         !                                                    ! for parameter read-in
         global_parameters,                                 & ! global parameters
         global_parameters_name,                            & ! clear names of global parameters
         timeStep_model_outputs,                            & ! timestep for writing model outputs
         outputFlxState,                                    & ! definition which output to write
         iFlag_LAI_data_format,                             & ! flag on how LAI data has to be read
         !                                                    ! used when iFlag_LAI_data_format = 1
         inputFormat_gridded_LAI,                           & ! format of gridded LAI data(bin or nc)
         dir_gridded_LAI                                      ! Directory where gridded LAI is located

    implicit none

    ! LOCAL variables
    ! PARAMETERS
    integer(i4), dimension(nProcesses)              :: processCase               ! Choosen process description number

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
    ! directRunoff
    real(dp), dimension(nColPars)                   :: imperviousStorageCapacity         
    ! actualET
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_forest    
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_impervious
    real(dp), dimension(nColPars)                   :: rootFractionCoefficient_pervious  
    real(dp), dimension(nColPars)                   :: minCorrectionFactorPET            
    real(dp), dimension(nColPars)                   :: maxCorrectionFactorPET            
    real(dp), dimension(nColPars)                   :: aspectTresholdPET                 
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
    ! routing
    real(dp), dimension(nColPars)                   :: muskingumTravelTime_constant      
    real(dp), dimension(nColPars)                   :: muskingumTravelTime_riverLength   
    real(dp), dimension(nColPars)                   :: muskingumTravelTime_riverSlope    
    real(dp), dimension(nColPars)                   :: muskingumTravelTime_impervious    
    real(dp), dimension(nColPars)                   :: muskingumAttenuation_riverSlope    
    !
    integer(i4)                                     :: ii
    real(dp)                                        :: cellFactorRbyH            ! conversion factor L11 to L1
    !
    ! some dummy arrays for namelist read in (allocatables not allowed in namelists)
    character(256)                                  :: dummy 
    character(256)                                  :: fname
    integer(i4),dimension(maxNoSoilHorizons)        :: soilDepth_dummy           ! depth of the single horizons
    character(256), dimension(maxNoBasins)          :: dirMorpho_dummy
    character(256), dimension(maxNoBasins)          :: dirLCover_dummy
    character(256), dimension(maxNoBasins)          :: dirGauges_dummy
    character(256), dimension(maxNoBasins)          :: dirPrecipitation_dummy
    character(256), dimension(maxNoBasins)          :: dirTemperature_dummy
    character(256), dimension(maxNoBasins)          :: dirReferenceET_dummy
    character(256), dimension(maxNoBasins)          :: dirOut_dummy
    character(256), dimension(maxNoBasins)          :: dirRestartOut_dummy
    character(256), dimension(maxNoBasins)          :: dirRestartIn_dummy
    character(256), dimension(maxNoBasins)          :: dirLatLon_dummy
    integer(i4),    dimension(maxNLCovers)          :: LCoverYearStart           ! starting year of LCover
    integer(i4),    dimension(maxNLCovers)          :: LCoverYearEnd             ! ending year  of LCover
    character(256), dimension(maxNLCovers)          :: LCoverfName_dummy         ! filename of Lcover file
    real(dp),       dimension(maxGeoUnit, nColPars) :: GeoParam                  ! geological parameters
    !
    character(256), dimension(maxNoBasins)          :: dir_gridded_LAI_dummy     ! directory of gridded LAI data 
    !                                                                            ! used when iFlag_LAI_data_format = 1
    real(dp)                                        :: jday_frac
    ! define namelists
    ! namelist directories
    namelist /directories/ dirConfigOut, dirCommonFiles_In, inputFormat_meteo_forcings,            &
         dirMorpho_dummy,dirLCover_dummy,dirGauges_dummy,dirPrecipitation_dummy, &
         dirTemperature_dummy, dirReferenceET_dummy, dirOut_dummy, dirRestartOut_dummy,&
         dirRestartIn_dummy, dirLatLon_dummy
    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig/ timestep, resolutionHydrology, resolutionRouting, optimize, opti_method,  &
         opti_function, nBasins, restart_flag_states_read, restart_flag_states_write, &
         restart_flag_config_read, restart_flag_config_write, warmingDays, evalPer
    ! namelsit soil layering
    namelist /soilLayer/ tillageDepth, nSoilHorizons_mHM, soilDepth_dummy
    ! namelist for land cover scenes
    namelist/LCover/fracSealed_cityArea,nLcover_scene,LCoverYearStart,LCoverYearEnd,LCoverfName_dummy
    ! namelist for LAI related data
    namelist/LAI_data_information/iFlag_LAI_data_format,inputFormat_gridded_LAI,dir_gridded_LAI_dummy
    ! namelist for pan evaporation
    namelist/panEvapo/evap_coeff
    ! namelist for night-day ratio of precipitation, referenceET and temperature
    namelist/nightDayRatio/fnight_prec,fnight_pet,fnight_temp
    ! namelsit process selection
    namelist /processSelection/ processCase
    ! namelist parameters
    namelist /interception1/ canopyInterceptionFactor
    namelist /snow1/snowTreshholdTemperature, degreeDayFactor_forest, degreeDayFactor_impervious,                      &
         degreeDayFactor_pervious, increaseDegreeDayFactorByPrecip, maxDegreeDayFactor_forest,               &
         maxDegreeDayFactor_impervious, maxDegreeDayFactor_pervious       
    namelist/soilmoisture1/ orgMatterContent_forest, orgMatterContent_impervious, orgMatterContent_pervious,           &         
         PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant,      &           
         PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant,                                    &
         PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope,                                                &
         rootFractionCoefficient_forest, rootFractionCoefficient_impervious,                         &
         rootFractionCoefficient_pervious, infiltrationShapeFactor
    namelist /directRunoff1/ imperviousStorageCapacity
    namelist /actualET1/  minCorrectionFactorPET, maxCorrectionFactorPET,                                              &
         aspectTresholdPET 
    namelist /interflow1/ interflowStorageCapacityFactor, interflowRecession_slope, fastInterflowRecession_forest,     &     
         slowInterflowRecession_Ks, exponentSlowInterflow    
    namelist /percolation1/ rechargeCoefficient, rechargeFactor_karstic, gain_loss_GWreservoir_karstic     
    namelist /routing1/ muskingumTravelTime_constant, muskingumTravelTime_riverLength, muskingumTravelTime_riverSlope, &
         muskingumTravelTime_impervious, muskingumAttenuation_riverSlope
    namelist /geoparameter/ GeoParam
    ! name list regarding output
    namelist/NLoutputResults/timeStep_model_outputs, outputFlxState
    ! namelist for optimization settings
    namelist/Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, sce_npg, sce_nps

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
    allocate(dirMorpho       (nBasins))
    allocate(dirLCover       (nBasins))
    allocate(dirGauges       (nBasins))
    allocate(dirPrecipitation(nBasins))
    allocate(dirTemperature  (nBasins))
    allocate(dirReferenceET  (nBasins))
    allocate(dirOut          (nBasins))
    allocate(dirRestartOut   (nBasins))
    allocate(dirRestartIn    (nBasins))
    allocate(dirLatLon       (nBasins))

    !===============================================================
    !  determine simulation time period incl. warming days
    !===============================================================
    ! julain days for evaluation period
    jday_frac = date2dec(dd=evalPer%dStart, mm=evalPer%mStart, yy=evalPer%yStart)
    evalPer%julStart = nint(jday_frac) 

    jday_frac = date2dec(dd=evalPer%dEnd, mm=evalPer%mEnd, yy=evalPer%yEnd)
    evalPer%julEnd  = nint(jday_frac, i4 ) 

    ! determine warming period
    warmPer%julStart = evalPer%julStart - warmingDays 
    warmPer%julEnd   = evalPer%julStart - 1 

    jday_frac = real(warmPer%julStart,dp)
    call dec2date(jday_frac, dd=warmPer%dStart, mm=warmPer%mStart, yy=warmPer%yStart)

    jday_frac = real(warmPer%julEnd,dp)
    call dec2date(jday_frac, dd=warmPer%dEnd,   mm=warmPer%mEnd,   yy=warmPer%yEnd  )

    ! sumulation Period = warming Period + evaluation Period
    simPer%dStart   = warmPer%dStart
    simPer%mStart   = warmPer%mStart
    simPer%yStart   = warmPer%yStart
    simPer%julStart = warmPer%julStart
    simPer%dEnd     = evalPer%dEnd
    simPer%mEnd     = evalPer%mEnd
    simPer%yEnd     = evalPer%yEnd
    simPer%julEnd   = evalPer%julEnd

    !===============================================================
    !  Read namelist for mainpaths
    !===============================================================
    call position_nml('directories', unamelist)
    read(unamelist, nml=directories)
    dirMorpho       = dirMorpho_dummy       (1:nBasins)
    dirLCover       = dirLCover_dummy       (1:nBasins)
    dirGauges       = dirGauges_dummy       (1:nBasins)      
    dirPrecipitation= dirPrecipitation_dummy(1:nBasins)
    dirTemperature  = dirTemperature_dummy  (1:nBasins)
    dirReferenceET  = dirReferenceET_dummy  (1:nBasins)
    dirOut          = dirOut_dummy          (1:nBasins)
    dirRestartOut   = dirRestartOut_dummy   (1:nBasins)
    dirRestartIn    = dirRestartIn_dummy    (1:nBasins)
    dirLatLon       = dirLatLon_dummy       (1:nBasins)
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

    allocate(HorizonDepth_mHM(nSoilHorizons_mHM))
    HorizonDepth_mHM = 0.0_dp
    HorizonDepth_mHM(1:nSoilHorizons_mHM-1)  = soilDepth_dummy(1:nSoilHorizons_mHM-1)

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

    if(iFlag_LAI_data_format .EQ. 1) then
       allocate( dir_gridded_LAI(nBasins) )
       dir_gridded_LAI(:) = dir_gridded_LAI_dummy(1:nBasins)
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
    if (LCoverYearStart(1) .GT. evalPer%yStart) then
       call message()
       call message('***ERROR: Land cover for warming period is missing!')
       call message('   FILE: mhm.nml, namelist: LCover')
       call message('   SimStart   : ', trim(num2str(simPer%yStart)))
       call message('   LCoverStart: ', trim(num2str(LCoverYearStart(1))))
       stop       
    end if
    if (LCoverYearEnd(nLcover_scene) .LT. evalPer%yEnd) then
       call message()
       call message('***ERROR: Land cover period shorter than modelling period!')
       call message('   FILE: mhm.nml, namelist: LCover')
       call message('   SimEnd   : ', trim(num2str(simPer%yEnd)))
       call message('   LCoverEnd: ', trim(num2str(LCoverYearEnd(nLcover_scene))))
       stop
    end if
    !
    allocate(LcYearId(simPer%yStart:simPer%yEnd))
    do ii = 1, nLcover_scene
       ! land cover before model period                        ! land cover after model period
       if ((LCoverYearEnd(ii) .LT. evalPer%yStart)        .OR.  (LCoverYearStart(ii) .GT. evalPer%yEnd)) then
          cycle
       else if ((LCoverYearStart(ii) .LE. evalPer%yStart) .AND. (LCoverYearEnd(ii) .GE. evalPer%yEnd)) then
          LCyearId(simPer%yStart:simPer%yEnd)          = ii
          exit
       else if ((LCoverYearStart(ii) .LE. evalPer%yStart) .AND. (LCoverYearEnd(ii) .LT. evalPer%yEnd)) then
          LCyearId(simPer%yStart:LCoverYearEnd(ii))      = ii
       else if ((LCoverYearStart(ii) .GT. evalPer%yStart) .AND. (LCoverYearEnd(ii) .GE. evalPer%yEnd)) then
          LCyearId(LCoverYearStart(ii):simPer%yEnd) = ii
       else
          LCyearId(LCoverYearStart(ii):LCoverYearEnd(ii)) = ii
       end if
    end do
    !
    ! correct number of input land cover scenes to number of needed scenes
    nLcover_scene = maxval(LCyearId) - minval(LCyearId) + 1
    ! put land cover scenes to corresponding file name and LuT
    allocate(LCfilename(nLcover_scene))
    LCfilename(:) = LCoverfName_dummy(minval(LCyearId):maxval(LCyearId))
    ! update the ID's
    LCyearId = LCyearId - minval(LCyearId) + 1
    !
    if (any(LCyearId .EQ. nodata_i4)) then 
       call message()
       call message('***ERROR: Intermidiate land cover period is missing!')
       call message('   FILE: mhm.nml, namelist: LCover')
       stop
    end if
    !
    !===============================================================
    ! check matching of resolutions: hydrology, forcing and routing
    !===============================================================
    cellFactorRbyH = resolutionRouting / resolutionHydrology
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

    ! Process 5 - actualET (meteo correction  factor)
    select case (processCase(5))
    ! 1 - root fraction approach
    case(1)
       call position_nml('actualET1', unamelist_param)
       read(unamelist_param, nml=actualET1)
       processMatrix(5, 1) = processCase(5)
       processMatrix(5, 2) = 3_i4
       processMatrix(5, 3) = sum(processMatrix(1:5, 2))
       call append(global_parameters, reshape(minCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(maxCorrectionFactorPET,             (/1, nColPars/)))
       call append(global_parameters, reshape(aspectTresholdPET,                  (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'minCorrectionFactorPET', &
            'maxCorrectionFactorPET', &
            'aspectTresholdPET     '/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "actualET1" out of bound in ', &
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
    case(1)
       ! 1 - Muskingum approach
       call position_nml('routing1', unamelist_param)
       read(unamelist_param, nml=routing1)
       processMatrix(8, 1) = processCase(8)
       processMatrix(8, 2) = 5_i4
       processMatrix(8, 3) = sum(processMatrix(1:8, 2))
       call append(global_parameters, reshape(muskingumTravelTime_constant,    (/1, nColPars/)))
       call append(global_parameters, reshape(muskingumTravelTime_riverLength, (/1, nColPars/)))
       call append(global_parameters, reshape(muskingumTravelTime_riverSlope,  (/1, nColPars/)))
       call append(global_parameters, reshape(muskingumTravelTime_impervious,  (/1, nColPars/)))
       call append(global_parameters, reshape(muskingumAttenuation_riverSlope, (/1, nColPars/)))

       call append(global_parameters_name, (/ &
            'muskingumTravelTime_constant   ', &
            'muskingumTravelTime_riverLength', &
            'muskingumTravelTime_riverSlope ', &
            'muskingumTravelTime_impervious ', &
            'muskingumAttenuation_riverSlope'/))

       ! check if parameter are in range
       if ( .not. in_bound(global_parameters) ) then
          call message('***ERROR: parameter in namelist "routing1" out of bound in ', &
               trim(adjustl(file_namelist_param)))
          stop
       end if

    case DEFAULT
       call message()          
       call message('***ERROR: Process description for process "routing" does not exist!')
       stop
    end select

    !===============================================================
    ! Geological formations
    !===============================================================

    ! read in of number of geological formations
    fName = trim(adjustl(dirMorpho(1))) // trim(adjustl(file_geolut))
    open( unit=ugeolut, file=fname, action='read', status='old')
    read(ugeolut, *) dummy, nGeoUnits
    close(ugeolut)
    dummy = dummy//''   ! only to avoid warning

    ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
    call position_nml('geoparameter', unamelist_param)
    GeoParam = nodata_dp
    read(unamelist_param, nml=geoparameter)

    call append(global_parameters, GeoParam(1:nGeoUnits,:))

    do ii=1, nGeoUnits
       dummy = 'GeoParam('//trim(adjustl(num2str(ii)))//',:)'
       call append(global_parameters_name, (/ trim(adjustl(dummy)) /)) 
    end do

    ! check if parameter are in range
    if ( .not. in_bound(global_parameters) ) then
       call message('***ERROR: parameter in namelist "geoparameter" out of bound in ', &
            trim(adjustl(file_namelist_param)))
       stop
    end if

    ! for baseflow parameters
    processMatrix(9,1) = 1
    processMatrix(9,2) = nGeoUnits
    processMatrix(9,3) = sum(processMatrix(1:9, 2))

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
       sce_npg = 2 * size(global_parameters,1) + 1_i4
    end if
    ! number of points in each sub-complex: default = n+1
    if (sce_nps .lt. 0_i4) then
       sce_nps = size(global_parameters,1) + 1_i4
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

    call message( '  FLUXES:' )
    if (outputFlxState(9)) then
      call message( '    actual evapotranspiration aET                              [mm/T]')
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
    call message( '' )
    call message( 'FINISHED readin config' )

    ! warning message  
    if( any(outputFlxState) .and. optimize ) then
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
