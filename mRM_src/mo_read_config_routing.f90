module mo_read_config_routing
  use mo_kind, only : i4, dp
  implicit none
  public :: read_mrm_config
  public :: read_config_routing
  public :: read_routing_params
contains

  !
  subroutine read_mrm_config
    implicit none

    !===============================================================
    ! Read evaluation gauge information
    !===============================================================
    call read_config_routing()

  end subroutine read_mrm_config

  subroutine read_config_routing()
    use mo_julian, only: dec2date, date2dec
    use mo_message, only: message
    use mo_nml, only: open_nml, position_nml, close_nml
    use mo_mrm_constants, only: nodata_i4, &
         maxNoGauges, & ! maximum number of allowed gauges
         maxNLcovers, & ! maximum number of allowed LCover scenes
         maxNoBasins ! maximum number of allowed basins
    use mo_mrm_file, only: file_namelist, unamelist
    use mo_string_utils, only: num2str
    use mo_global_variables_routing, only : &
         timestep, & ! timestep of routing [h]
         iFlag_cordinate_sys, & ! model run cordinate system
         nGaugesTotal, gauge, & ! number of evaluation gauges and gauge informations 
         nInflowGaugesTotal, InflowGauge, & ! number of inflow gauges and gauge informations
         dirCommonFiles, & ! directory where common input files should be located for all modeled basins
         dirConfigOut, &
         dirMorpho, & ! Directory where morphological files are located
         dirLCover, & ! Directory where land cover files are located
         dirGauges, & ! Directory where discharge files are located
         dirOut, & ! Directory where output is written to
         dirRestartOut, & ! Directory where output of restart is written
         dirRestartIn, & ! Directory where input of restart is read from
         is_start, & ! flag for first timestep
         resolutionRouting, & ! resolution of routing
         resolutionHydrology, & ! resolution of Hydrology
         L0_Basin, & ! L0_Basin ID
         optimize, & ! if mhm runs in optimization mode or not
         opti_method, & ! optimization algorithm used
         opti_function, & ! objective function to be optimized
         read_restart, & ! flag reading restart
         write_restart, & ! flag writing restart
         perform_mpr, & ! switch for performing mpr
         nBasins, & ! number of basins
         simPer, &
         evalPer, &
         warmingDays, &
         warmPer, &
         timestep_model_inputs, &
         fracSealed_cityArea, nLcover_scene, & ! land cover information
         LCYearId, &
         LCfilename, &
         nIterations, & ! number of iterations in optimization
         seed, & ! seed used for optimization
         dds_r, & ! DDS: perturbation rate
         sa_temp, & ! SA: initial temperature
         sce_ngs, sce_npg, sce_nps, & ! SCE: # complexes, # points per complex,
         !                            !      # points per subcomplex
         mrm_global_parameters, &
         basin_mrm, &
         period ! structure for time periods

    implicit none
    ! input variables
    !
    ! local variables
    integer(i4),    dimension(maxNoBasins)             :: NoGauges_basin
    integer(i4),    dimension(maxNoBasins,maxNoGauges) :: Gauge_id
    character(256), dimension(maxNoGauges,maxNoGauges) :: Gauge_filename
    integer(i4),    dimension(maxNoBasins)             :: NoInflowGauges_basin
    integer(i4),    dimension(maxNoBasins,maxNoGauges) :: InflowGauge_id
    character(256), dimension(maxNoGauges,maxNoGauges) :: InflowGauge_filename
    logical,        dimension(maxNoBasins,maxNoGauges) :: InflowGauge_Headwater
    integer(i4)                                        :: iBasin
    integer(i4)                                        :: iGauge
    integer(i4)                                        :: idx
    integer(i4)                                        :: ii
    integer(i4)                                        :: n_true_pars
    real(dp)                                           :: cellFactorRbyH ! conversion factor L11 to L1

    ! namelist variables: direcoteries
    character(256), dimension(maxNoBasins) :: dir_Morpho
    character(256), dimension(maxNoBasins) :: dir_LCover
    character(256), dimension(maxNoBasins) :: dir_Gauges
    character(256), dimension(maxNoBasins) :: dir_Out
    character(256), dimension(maxNoBasins) :: dir_RestartOut
    character(256), dimension(maxNoBasins) :: dir_RestartIn
    ! namelist variables: mainconfig
    real(dp), dimension(maxNoBasins) :: resolution_Routing
    real(dp), dimension(maxNoBasins) :: resolution_Hydrology
    integer(i4), dimension(maxNoBasins) :: L0Basin
    ! namelist variables: LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart ! starting year of LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd ! ending year  of LCover
    character(256), dimension(maxNLCovers) :: LCoverfName ! filename of Lcover file
    ! namelist variables: time_periods
    real(dp) :: jday_frac
    integer(i4), dimension(maxNoBasins) :: warming_Days
    type(period), dimension(maxNoBasins) :: eval_Per
    integer(i4), dimension(maxNoBasins) :: time_step_model_inputs

    !

    ! namelist spatial & temporal resolution, otmization information
    namelist /mainconfig/ timestep, iFlag_cordinate_sys, resolution_Routing, resolution_Hydrology, &
         L0Basin, optimize, opti_method, opti_function, nBasins, read_restart,                     &
         write_restart, perform_mpr
    ! namelist for time settings
    namelist /time_periods/ warming_Days, eval_Per, time_step_model_inputs
    ! namelist for evaluation gauges
    ! define namelists
    ! namelist directories
    namelist /directories/ dirConfigOut, dirCommonFiles, &
         dir_Morpho, dir_LCover, dir_Gauges,             &
         dir_Out, dir_RestartOut,                        &
         dir_RestartIn
    namelist/LCover/ fracSealed_cityArea, nLcover_scene, LCoverYearStart, LCoverYearEnd, LCoverfName
    namelist /evaluation_gauges/ nGaugesTotal, NoGauges_basin, Gauge_id, gauge_filename
    ! namelist for inflow gauges
    namelist /inflow_gauges/ nInflowGaugesTotal, NoInflowGauges_basin, InflowGauge_id, &
         InflowGauge_filename, InflowGauge_Headwater
    ! namelist for optimization settings
    namelist/Optimization/ nIterations, seed, dds_r, sa_temp, sce_ngs, sce_npg, sce_nps

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    is_start = .True.
    nGaugesTotal   = nodata_i4
    NoGauges_basin = nodata_i4
    Gauge_id       = nodata_i4
    gauge_filename = num2str(nodata_i4)

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
    allocate(resolutionRouting(nBasins))
    allocate(resolutionHydrology(nBasins))
    allocate(L0_Basin(nBasins))
    allocate(dirMorpho(nBasins))
    allocate(dirLCover(nBasins))
    allocate(dirGauges(nBasins))
    allocate(dirOut(nBasins))
    allocate(dirRestartOut(nBasins))
    allocate(dirRestartIn(nBasins))
    resolutionRouting = resolution_Routing(1:nBasins)
    resolutionHydrology = resolution_Hydrology(1:nBasins)
    L0_Basin = L0Basin(1:nBasins)

    !
    ! check for possible options
    if( .NOT. (iFlag_cordinate_sys == 0 .OR. iFlag_cordinate_sys == 1) ) then
       call message()
       call message('***ERROR: coordinate system for the model run should be 0 or 1')
       stop
    end if
    ! check for perform_mpr
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
    evalPer = eval_Per(1:nBasins)
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
    call position_nml('directories', unamelist)
    read(unamelist, nml=directories)

    dirMorpho = dir_Morpho(1:nBasins)
    dirLCover = dir_LCover(1:nBasins)
    dirGauges = dir_Gauges(1:nBasins)
    dirOut = dir_Out(1:nBasins)
    dirRestartOut = dir_RestartOut(1:nBasins)
    dirRestartIn = dir_RestartIn(1:nBasins)

    !===============================================================
    ! Read land cover information
    !===============================================================
    call position_nml('LCover', unamelist)
    read(unamelist, nml=LCover)

    !===============================================================
    ! READ EVALUATION GAUGES
    !===============================================================
    call position_nml('evaluation_gauges', unamelist)
    read(unamelist, nml=evaluation_gauges)
    
    if (nGaugesTotal .GT. maxNoGauges) then
       call message()
       call message('***ERROR: mhm.nml: Total number of evaluation gauges is restricted to', num2str(maxNoGauges))
       call message('          Error occured in namlist: evaluation_gauges')
       stop
    end if

    allocate(gauge%gaugeId(nGaugesTotal)) ; gauge%gaugeId  = nodata_i4
    allocate(gauge%basinId(nGaugesTotal)) ; gauge%basinId  = nodata_i4
    allocate(gauge%fName  (nGaugesTotal)) ; gauge%fName(1) = num2str(nodata_i4)
    allocate(basin_mrm%nGauges        (nBasins                           )) ; basin_mrm%nGauges        = nodata_i4
    allocate(basin_mrm%gaugeIdList    (nBasins, maxval(NoGauges_basin(:)))) ; basin_mrm%gaugeIdList    = nodata_i4
    allocate(basin_mrm%gaugeIndexList (nBasins, maxval(NoGauges_basin(:)))) ; basin_mrm%gaugeIndexList = nodata_i4
    allocate(basin_mrm%gaugeNodeList  (nBasins, maxval(NoGauges_basin(:)))) ; basin_mrm%gaugeNodeList  = nodata_i4

    idx = 0
    do iBasin = 1, nBasins
       ! check if NoGauges_basin has a valid value
       if ( NoGauges_basin(iBasin) .EQ. nodata_i4 ) then
          call message()
          call message('***ERROR: mhm.nml: Number of evaluation gauges for subbasin ', &
               trim(adjustl(num2str(iBasin))),' is not defined!')
          call message('          Error occured in namlist: evaluation_gauges')
          stop
       end if

       basin_mrm%nGauges(iBasin) = NoGauges_basin(iBasin)
       
       do iGauge = 1, NoGauges_basin(iBasin)
          ! check if NoGauges_basin has a valid value
          if (Gauge_id(iBasin,iGauge) .EQ. nodata_i4) then
             call message()
             call message('***ERROR: mhm.nml: ID of evaluation gauge ',        &
                  trim(adjustl(num2str(iGauge))),' for subbasin ', &
                  trim(adjustl(num2str(iBasin))),' is not defined!')
             call message('          Error occured in namlist: evaluation_gauges')
             stop
          else if (trim(gauge_filename(iBasin,iGauge)) .EQ. trim(num2str(nodata_i4))) then
             call message()
             call message('***ERROR: mhm.nml: Filename of evaluation gauge ', &
                  trim(adjustl(num2str(iGauge))),' for subbasin ',  &
                  trim(adjustl(num2str(iBasin))),' is not defined!')
             call message('          Error occured in namlist: evaluation_gauges')
             stop
          end if
          !
          idx = idx + 1
          gauge%basinId(idx) = iBasin
          gauge%gaugeId(idx) = Gauge_id(iBasin,iGauge)
          gauge%fname(idx)   = trim(dirGauges(iBasin)) // trim(gauge_filename(iBasin,iGauge))
          basin_mrm%gaugeIdList(iBasin,iGauge)    = Gauge_id(iBasin,iGauge)
          basin_mrm%gaugeIndexList(iBasin,iGauge) = idx
       end do
    end do

    if ( nGaugesTotal .NE. idx) then
       call message()
       call message('***ERROR: mhm.nml: Total number of evaluation gauges (', trim(adjustl(num2str(nGaugesTotal))), &
            ') different from sum of gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
       call message('          Error occured in namlist: evaluation_gauges')
       stop
    end if
    
    !===============================================================
    ! Read inflow gauge information
    !===============================================================
    
    nInflowGaugesTotal   = 0
    NoInflowGauges_basin = 0
    InflowGauge_id       = nodata_i4
    InflowGauge_filename = num2str(nodata_i4)

    call position_nml('inflow_gauges', unamelist)
    read(unamelist, nml=inflow_gauges)

    if (nInflowGaugesTotal .GT. maxNoGauges) then
       call message()
       call message('***ERROR: mhm.nml:read_gauge_lut: Total number of inflow gauges is restricted to', num2str(maxNoGauges))
       call message('          Error occured in namlist: inflow_gauges')
       stop
    end if


    ! allocation - max() to avoid allocation with zero, needed for mhm call
    allocate(InflowGauge%gaugeId (max(1,nInflowGaugesTotal)))
    allocate(InflowGauge%basinId (max(1,nInflowGaugesTotal)))
    allocate(InflowGauge%fName   (max(1,nInflowGaugesTotal)))
    allocate(basin_mrm%nInflowGauges        (nBasins                                 ))
    allocate(basin_mrm%InflowGaugeIdList    (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin_mrm%InflowGaugeHeadwater (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin_mrm%InflowGaugeIndexList (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin_mrm%InflowGaugeNodeList  (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    ! dummy initialization
    InflowGauge%gaugeId = nodata_i4
    InflowGauge%basinId = nodata_i4
    InflowGauge%fName   = num2str(nodata_i4)
    basin_mrm%nInflowGauges        = 0
    basin_mrm%InflowGaugeIdList    = nodata_i4
    basin_mrm%InflowGaugeHeadwater = .FALSE.
    basin_mrm%InflowGaugeIndexList = nodata_i4
    basin_mrm%InflowGaugeNodeList  = nodata_i4

    idx = 0
    do iBasin = 1, nBasins

       ! no inflow gauge for subbasin i
       if (NoInflowGauges_basin(iBasin) .EQ. nodata_i4) then
          NoInflowGauges_basin(iBasin)       = 0
       end if

       basin_mrm%nInflowGauges(iBasin) = NoInflowGauges_basin(iBasin)

       do iGauge = 1, NoInflowGauges_basin(iBasin)
          ! check if NoInflowGauges_basin has a valid value
          if (InflowGauge_id(iBasin,iGauge) .EQ. nodata_i4) then
             call message()
             call message('***ERROR: mhm.nml:ID of inflow gauge ',        &
                  trim(adjustl(num2str(iGauge))),' for subbasin ', &
                  trim(adjustl(num2str(iBasin))),' is not defined!')
             call message('          Error occured in namlist: inflow_gauges')
             stop
          else if (trim(InflowGauge_filename(iBasin,iGauge)) .EQ. trim(num2str(nodata_i4))) then
             call message()
             call message('***ERROR: mhm.nml:Filename of inflow gauge ', &
                  trim(adjustl(num2str(iGauge))),' for subbasin ',  &
                  trim(adjustl(num2str(iBasin))),' is not defined!')
             call message('          Error occured in namlist: inflow_gauges')
             stop
          end if
          !
          idx = idx + 1
          InflowGauge%basinId(idx) = iBasin
          InflowGauge%gaugeId(idx) = InflowGauge_id(iBasin,iGauge)
          InflowGauge%fname(idx)   = trim(dirGauges(iBasin)) // trim(InflowGauge_filename(iBasin,iGauge))
          basin_mrm%InflowGaugeIdList(iBasin,iGauge)    = InflowGauge_id(iBasin,iGauge)
          basin_mrm%InflowGaugeHeadwater(iBasin,iGauge) = InflowGauge_Headwater(iBasin,iGauge)
          basin_mrm%InflowGaugeIndexList(iBasin,iGauge) = idx
       end do
    end do

    if ( nInflowGaugesTotal .NE. idx) then
       call message()
       call message('***ERROR: mhm.nml: Total number of inflow gauges (', trim(adjustl(num2str(nInflowGaugesTotal))), &
            ') different from sum of inflow gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
       call message('          Error occured in namlist: inflow_gauges')
       stop
    end if

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
          if ((LCoverYearEnd(ii) .LT. evalPer(iBasin)%yStart) .OR. &
               (LCoverYearStart(ii) .GT. evalPer(iBasin)%yEnd)) then
             cycle
          else if ((LCoverYearStart(ii) .LE. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii) .GE. evalPer(iBasin)%yEnd)) then
             LCyearId(simPer(iBasin)%yStart:simPer(iBasin)%yEnd, iBasin) = ii
             exit
          else if ((LCoverYearStart(ii) .LE. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii) .LT. evalPer(iBasin)%yEnd)) then
             LCyearId(simPer(iBasin)%yStart:LCoverYearEnd(ii), iBasin) = ii
          else if ((LCoverYearStart(ii) .GT. evalPer(iBasin)%yStart) .AND. &
               (LCoverYearEnd(ii) .GE. evalPer(iBasin)%yEnd)) then
             LCyearId(LCoverYearStart(ii):simPer(iBasin)%yEnd, iBasin) = ii
          else
             LCyearId(LCoverYearStart(ii):LCoverYearEnd(ii), iBasin) = ii
          end if
       end do
    end do
    !
    ! correct number of input land cover scenes to number of needed scenes
    nLcover_scene = maxval(LCyearId, mask = (LCyearId .gt. nodata_i4) ) - &
         minval(LCyearId, mask = (LCyearId .gt. nodata_i4) ) + 1
    ! put land cover scenes to corresponding file name and LuT
    allocate(LCfilename(nLcover_scene))
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

    call close_nml(unamelist)
    
    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    call read_routing_params(1)

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
       n_true_pars = count(nint(mrm_global_parameters(:,4)) .eq. 1)
       sce_npg = 2 * n_true_pars + 1_i4
    end if
    ! number of points in each sub-complex: default = n+1
    if (sce_nps .lt. 0_i4) then
       n_true_pars = count(nint(mrm_global_parameters(:,4)) .eq. 1)
       sce_nps = n_true_pars + 1_i4
    end if
    if (sce_npg .lt. sce_nps) then
       call message ('number of points per complex (sce_npg) must be greater or')
       call message ('equal number of points per sub-complex (sce_nps)')
       stop
    end if

    call close_nml(unamelist)

    
end subroutine read_config_routing

  subroutine read_routing_params(processCase)
    use mo_mrm_file, only: file_namelist_param, unamelist_param ! file containing parameter values
    use mo_append, only: append
    use mo_nml, only: open_nml, position_nml, close_nml
    use mo_message, only: message
    use mo_mrm_constants, only: nColPars ! number of properties of the global variables
    use mo_global_variables_routing, only: &
         mrm_global_parameters, &
         mrm_global_parameters_name
#ifdef mrm2mhm    
    use mo_global_variables, only :  processMatrix, & ! process configuration
         global_parameters, & ! global parameters
         global_parameters_name ! clear names of global parameters
#endif

    implicit none
    ! input variables
    integer(i4), intent(in) :: processCase ! it is the default case should be one
    ! local variables
    integer(i4) :: start_index ! equals sum of previous parameters
    real(dp), dimension(nColPars) :: muskingumTravelTime_constant
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope
    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious
    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope

    namelist /routing1/ muskingumTravelTime_constant, muskingumTravelTime_riverLength, &
         muskingumTravelTime_riverSlope, muskingumTravelTime_impervious, muskingumAttenuation_riverSlope
    !
    call open_nml(file_namelist_param, unamelist_param, quiet=.true.)

    call position_nml('routing1', unamelist_param)
    read(unamelist_param, nml=routing1)

#ifdef mrm2mhm
    ! set variables of mhm
    processMatrix(8, 1) = processCase
    processMatrix(8, 2) = 5_i4
    processMatrix(8, 3) = sum(processMatrix(1:8, 2))
    ! insert parameter values and names at position required by mhm
    start_index         = processMatrix(8, 3) - processMatrix(8, 2)
    global_parameters(start_index + 1, :) = muskingumTravelTime_constant
    global_parameters(start_index + 2, :) = muskingumTravelTime_riverLength
    global_parameters(start_index + 3, :) = muskingumTravelTime_riverSlope
    global_parameters(start_index + 4, :) = muskingumTravelTime_impervious
    global_parameters(start_index + 5, :) = muskingumAttenuation_riverSlope
    
    global_parameters_name(start_index + 1 : start_index + processMatrix(8,2)) = (/ &
         'muskingumTravelTime_constant   ', &
         'muskingumTravelTime_riverLength', &
         'muskingumTravelTime_riverSlope ', &
         'muskingumTravelTime_impervious ', &
         'muskingumAttenuation_riverSlope'/)

    ! check if parameter are in range
    if ( .not. in_bound(global_parameters) ) then
       call message('***ERROR: parameter in namelist "routing1" out of bound in ', &
            trim(adjustl(file_namelist_param)))
       stop
    end if
#endif

    ! set variables of mrm (redundant in case of coupling to mhm)
    call append(mrm_global_parameters, reshape(muskingumTravelTime_constant,    (/1, nColPars/)))
    call append(mrm_global_parameters, reshape(muskingumTravelTime_riverLength, (/1, nColPars/)))
    call append(mrm_global_parameters, reshape(muskingumTravelTime_riverSlope,  (/1, nColPars/)))
    call append(mrm_global_parameters, reshape(muskingumTravelTime_impervious,  (/1, nColPars/)))
    call append(mrm_global_parameters, reshape(muskingumAttenuation_riverSlope, (/1, nColPars/)))

    call append(mrm_global_parameters_name, (/ &
         'muskingumTravelTime_constant   ', &
         'muskingumTravelTime_riverLength', &
         'muskingumTravelTime_riverSlope ', &
         'muskingumTravelTime_impervious ', &
         'muskingumAttenuation_riverSlope'/))

    ! check if parameter are in range
    if ( .not. in_bound(mrm_global_parameters) ) then
       call message('***ERROR: parameter in namelist "routing1" out of bound in ', &
            trim(adjustl(file_namelist_param)))
       stop
    end if

    call close_nml(unamelist_param)

  end subroutine read_routing_params

  ! --------------------------------------------------------------------------------
  ! private funtions and subroutines, DUPLICATED FROM mo_read_config.f90
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

end module mo_read_config_routing
