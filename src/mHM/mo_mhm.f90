!> \file mo_mhm.f90

!> \brief Call all main processes of mHM.

!> \details This module calls all processes of mHM for a given configuration.
!>          The configuration of the model is stored in the a process matrix.
!>          This configuration is specified in the namelist mhm.nml.
!>
!>          The processes are executed in ascending order. At the moment only
!>          process 5 and 8 have options.\n
!>          The MPR technique is only called either if the land cover has been
!>          changed or for very first time step.\n
!>
!>          Currently the following processes are implemented: \n
!>
!>          Process    | Name                      | Flag  | Description
!>          ---------- | ------------------------- | ----- | ------------------------------------------
!>          1          | interception              |   1   | Maximum interception
!>          2          | snow and melting          |   1   | Degree-day
!>          3          | soil moisture             |   1   | Feddes equation for ET reduction, Brooks-Corey like
!>          3          | soil moisture             |   2   | Jarvis equation for ET reduction, Brooks-Corey like
!>          3          | soil moisture             |   3   | Jarvis eq. for ET red. + FC dependency on root frac. coef.
!>          4          | direct runoff             |   1   | Linear reservoir exceedance
!>          5          | PET                       |  -1   | PET is input, LAI based correction, dynamic scaling func.
!>          5          | PET                       |   0   | PET is input, Aspect based correction
!>          5          | PET                       |   1   | Hargreaves-Samani
!>          5          | PET                       |   2   | Priestley-Taylor
!>          5          | PET                       |   3   | Penman-Monteith
!>          6          | interflow                 |   1   | Nonlinear reservoir with saturation excess
!>          7          | percolation and base flow |   1   | GW linear reservoir
!>          8          | routing                   |   0   | no routing
!>          8          | routing                   |   1   | use mRM i.e. Muskingum
!>          8          | routing                   |   2   | use mRM i.e. adaptive timestep

!>

!> \author Luis Samaniego
!> \date Dec 2012

MODULE mo_mHM

  use mo_kind,          only: i4, dp
  use mo_mhm_constants, only: nodata_dp
  use mo_message,       only: message
  !$ USE omp_lib

  IMPLICIT NONE

  PUBLIC :: mHM      ! initialization sequence

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------

  !      NAME
  !         mHM

  !     PURPOSE
  !>        \brief Pure mHM calculations.

  !>        \details Pure mHM calculations. All variables are allocated and initialized. \n
  !>                 They will be local variables within this call. \n
  !>

  !     INTENT(IN)
  !         Has to be updated...

  !     INTENT(INOUT)
  !         Has to be updated...

  !     INTENT(OUT)
  !         Has to be updated...

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note Fields must be consistent to DEM.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Luis Samaniego & Rohini Kumar
  !>        \date    Dec 2012

  !         Modified Luis Samaniego, Rohini Kumar     Dec 2012 - modularization
  !                  Luis Samaniego,                  Feb 2013 - call routine
  !                  Rohini Kumar,                    Feb 2013 - MPR call and other pre-requisite
  !                                                              variables for this call
  !                  Rohini Kumar,                    May 2013 - Error checks
  !                  Rohini Kumar,                    Jun 2013 - sealed area correction in total runoff
  !                                                            - initalization of soil moist. at first timestep
  !                  Rohini Kumar,                    Aug 2013 - dynamic LAI option included, and changed within
  !                                                              the code made accordingly (e.g., canopy intecpt.)
  !                                                            - max. canopy interception is estimated outside of MPR
  !                                                              call
  !                  Matthias Zink,                   Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  !                  Matthias Zink,                   Mar 2014 - added inflow from upstream areas
  !                  Matthias Zink,                   Apr 2014 - added PET calculation: Priestley-Taylor and Penamn-Monteith
  !                                                              and its parameterization (Process 5)
  !                  Rohini Kumar,                    Apr 2014 - mHM run with a single L0 grid cell, also in the routing mode
  !                  Stephan Thober,                  Jun 2014 - added flag for switching of MPR
  !                  Matthias Cuntz & Juliane Mai     Nov 2014 - LAI input from daily, monthly or yearly files
  !                  Matthias Zink,                   Dec 2014 - adopted inflow gauges to ignore headwater cells
  !                  Stephan Thober,                  Aug 2015 - moved routing to mRM
  !                  Rohini Kumar,                    Mar 2016 - changes for handling multiple soil database options
  !                  Rohini Kumar,                    Dec 2016 - changes for reading gridded mean monthly LAI fields
  !                  Stephan Thober,                  Jan 2017 - added prescribed weights for tavg and pet
  !                  Zink M. Demirel C.,              Mar 2017 - added Jarvis soil water stress function at SM process(3)  
  !                  M.Cuneyd Demirel & Simon Stisen  May 2017 - added FC dependency on root fraction coef. at SM process(3)
  !                  M.Cuneyd Demirel & Simon Stisen  Jun 2017 - added PET correction based on LAI at PET process(5)
  !                                                  
  ! ------------------------------------------------------------------

  subroutine mHM(  &
      ! Input -----------------------------------------------------------------
      ! Configuration
      perform_mpr,          & ! flag for reading restart files for state variables
      read_states,          & ! flag indicating whether states have been read
      fSealedInCity       , & ! sealed area fraction within cities
      timeStep_LAI_input  , & ! time step of gridded LAI input
      counter_year        , & ! counter to tackle the change of year
      counter_month       , & ! counter to tackle the change of month
      counter_day         , & ! counter to tackle the change of day
      tt                  , & ! simulation time step
      time                , & ! current decimal Julian day
      processMatrix       , & ! mHM process configuration matrix
      c2TSTu              , & ! unit transformation coefficient
      horizon_depth       , & ! Depth of each horizon in mHM
      nCells1             , & ! number of cells in a given basin at level L1
      nHorizons_mHM       , & ! Number of Horizons in mHM
      ntimesteps_day      , & ! number of time intervals per day, transformed in dp
      mask0               , & ! mask 0 for MPR
      neutron_integral_AFast,&! tabular for neutron flux approximation
      iflag_soil_option   , & ! flags for handling multiple soil databases      
      global_parameters   , & ! global mHM parameters
      ! LUT
      LCyearId            , & ! mapping of landcover scenes
      GeoUnitList         , & ! List of Ids for geological units
      GeoUnitKar          , & ! List of Ids for geological units with Karstic formation
      LAIUnitList         , & ! List of ids of each LAI class in LAILUT
      LAILUT              , & ! List of Ids for LAI
      ! Physiographic L0
      slope_emp0          , &
      L0_Latitude         , & ! latitude on level 0
      cellId0             , & ! cell Ids at level 0
      soilId0             , & ! soil Ids at level 0
      L0_LCover_LAI       , & ! land cover ID for LAI estimation
      LCover0             , & ! land use cover at level 0
      Asp0                , & ! [degree] Aspect at Level 0
      LAI0                , & ! LAI at level 0
      geoUnit0            , & ! geological units at level 0
      SDB_is_present      , & ! indicates whether soiltype exists
      SDB_nHorizons       , & ! Number of Horizons per soiltype
      SDB_nTillHorizons   , & ! Number of Tillage Horizons
      SDB_sand            , & ! sand content from LUT
      SDB_clay            , & ! clay content from LUT
      SDB_DbM             , & ! mineral Bulk density from LUT
      SDB_Wd              , & ! soil weighing vertical column of mHM
      SDB_RZdepth         , & ! soil depth from LUT
      ! Physiographic L1
      nTCells0_inL1       , & ! total number of valid L0 cells in a given L1 cell
      L0upBound_inL1      , & ! upper row of L0 block within L1 cell
      L0downBound_inL1    , & ! lower row of L0 block within L1 cell
      L0leftBound_inL1    , & ! left column of L0 block within L1 cell
      L0rightBound_inL1   , & ! right column of L0 block within L1 cell
      latitude            , & ! latitude on level 1
      ! Forcings
      evap_coeff          , & ! Evaporation coefficent for free-water surface of that current month
      fday_prec           , & ! [-] day ratio precipitation < 1
      fnight_prec         , & ! [-] night ratio precipitation < 1
      fday_pet            , & ! [-] day ratio PET  < 1
      fnight_pet          , & ! [-] night ratio PET  < 1
      fday_temp           , & ! [-] day factor mean temp
      fnight_temp         , & ! [-] night factor mean temp
      temp_weights        , & ! multiplicative weights for temperature (deg K)
      pet_weights         , & ! multiplicative weights for potential evapotranspiration
      pre_weights         , & ! multiplicative weights for precipitation
      read_meteo_weights  , & ! flag whether weights for tavg and pet have read and should be used
      pet_in              , & ! [mm d-1] Daily potential evapotranspiration (input)
      tmin_in             , & ! [degc]   Daily minimum temperature
      tmax_in             , & ! [degc]   Daily maxumum temperature
      netrad_in           , & ! [w m2]   Daily average net radiation
      absvappres_in       , & ! [Pa]     Daily average absolute vapour pressure
      windspeed_in        , & ! [m s-1]  Daily average wind speed
      prec_in             , & ! [mm d-1] Daily mean precipitation
      temp_in             , & ! [degc]   Daily average temperature
      ! In-Out -----------------------------------------------------------------
      ! Configuration
      yId                 , & ! Current Id of the LCover year scene
      ! Land cover L1 and L11
      fForest1            , & ! fraction of forest cover at scale L1
      fPerm1              , & ! fraction of permeable area at scale L1
      fSealed1            , & ! fraction of sealed area at scale L1
      ! States
      interc              , & ! Interception
      snowpack            , & ! Snowpack
      sealedStorage       , & ! Retention storage of impervious areas
      soilMoisture        , & ! Soil moisture of each horizon
      unsatStorage        , & ! Upper soil storage
      satStorage          , & ! Groundwater storage
      neutrons            , & ! Ground albedo neutrons
      ! Fluxes L1
      pet_calc            , & ! [mm TST-1] estimated PET (if PET is input = corrected values (fAsp*PET))
      aet_soil            , & ! actual ET
      aet_canopy          , & ! Real evaporation intensity from canopy
      aet_sealed          , & ! Actual ET from free-water surfaces
      baseflow            , & ! Baseflow
      infiltration        , & ! Recharge, infiltration intensity or effective precipitation of each horizon
      fast_interflow      , & ! Fast runoff component
      melt                , & ! Melting snow depth
      perc                , & ! Percolation
      prec_effect         , & ! Effective precipitation depth (snow melt + rain)
      rain                , & ! Rain precipitation depth
      runoff_sealed       , & ! Direct runoff from impervious areas
      slow_interflow      , & ! Slow runoff component
      snow                , & ! Snow precipitation depth
      throughfall         , & ! Throughfall
      total_runoff        , & ! Generated runoff
      ! Effective Parameters
      alpha               , & ! Exponent for the upper reservoir
      deg_day_incr        , & ! Increase of the Degree-day factor per mm of increase in precipitation
      deg_day_max         , & ! Maximum Degree-day factor
      deg_day_noprec      , & ! Degree-day factor with no precipitation
      deg_day             , & ! Degree-day factor
      fAsp                , & ! [1]     PET correction for Aspect at level 1
      petLAIcorFactorL1   , & ! PET correction factor based on LAI at level 1    
      HarSamCoeff         , & ! [1]     PET Hargreaves Samani coefficient at level 1
      PrieTayAlpha        , & ! [1]     PET Priestley Taylor coefficient at level 1
      aeroResist          , & ! [s m-1] PET aerodynamical resitance at level 1
      surfResist          , & ! [s m-1] PET bulk surface resitance at level 1
      frac_roots          , & ! Fraction of Roots in soil horizon
      interc_max          , & ! Maximum interception
      karst_loss          , & ! Karstic percolation loss
      k0                  , & ! Recession coefficient of the upper reservoir, upper outlet
      k1                  , & ! Recession coefficient of the upper reservoir, lower outlet
      k2                  , & ! Baseflow recession coefficient
      kp                  , & ! Percolation coefficient
      soil_moist_FC       , & ! Soil moisture below which actual ET is reduced
      soil_moist_sat      , & ! Saturation soil moisture for each horizon [mm]
      soil_moist_exponen  , & ! Exponential parameter to how non-linear is the soil water retention
      jarvis_thresh_c1    , & ! jarvis critical value for normalized soil water content
      temp_thresh         , & ! Threshold temperature for snow/rain
      unsat_thresh        , & ! Threshold water depth in upper reservoir
      water_thresh_sealed , & ! Threshold water depth in impervious areas
      wilting_point         ) ! Permanent wilting point for each horizon
    ! subroutines required to estimate variables prior to the MPR call
    use mo_upscaling_operators,     only: L0_fractionalCover_in_Lx         ! land cover fraction
    use mo_multi_param_reg,         only: mpr,canopy_intercept_param       ! reg. and scaling
    use mo_pet,                     only: pet_hargreaves, pet_priestly,  & ! calc. of pot. evapotranspiration
                                          pet_penman
    use mo_Temporal_Disagg_Forcing, only: Temporal_Disagg_Forcing
    use mo_canopy_interc ,          only: canopy_interc
    use mo_snow_accum_melt,         only: snow_accum_melt
    use mo_soil_moisture,           only: soil_moisture
    use mo_neutrons,                only: DesiletsN0, COSMIC
    use mo_runoff,                  only: runoff_unsat_zone
    use mo_runoff,                  only: runoff_sat_zone, L1_total_runoff
    use mo_julian,                  only: dec2date, date2dec
    use mo_string_utils,            only: num2str
    use mo_mhm_constants,           only: HarSamConst ! parameters for Hargreaves-Samani Equation
    use mo_mpr_pet,                 only: pet_correctbyLAI

    implicit none

    ! Intent
    logical,                     intent(in) :: perform_mpr          ! flag for reading restart files for state variables
    logical,                     intent(in) :: read_states          ! indicated whether states have been read from file
    real(dp),                    intent(in) :: fSealedInCity        ! fraction of perfectly sealed area within cities
    integer(i4),                 intent(in) :: timeStep_LAI_input   ! time step of gridded LAI input
    integer(i4),                 intent(in) :: counter_year         ! counter to tackle the change of year
    integer(i4),                 intent(in) :: counter_month        ! counter to tackle the change of month
    integer(i4),                 intent(in) :: counter_day          ! counter to tackle the change of day
    integer(i4),                 intent(in) :: tt
    real(dp),                    intent(in) :: time
    integer(i4), dimension(:,:), intent(in) :: processMatrix
    real(dp),                    intent(in) :: c2TSTu
    real(dp),    dimension(:),   intent(in) :: horizon_depth
    integer(i4),                 intent(in) :: nCells1
    integer(i4),                 intent(in) :: nHorizons_mHM
    real(dp),                    intent(in) :: ntimesteps_day
    logical,     dimension(:,:), intent(in) :: mask0
    real(dp),    dimension(:),   intent(in) :: neutron_integral_AFast
    integer(i4),                 intent(in) :: iflag_soil_option                
    real(dp),    dimension(:),   intent(in) :: global_parameters

    ! LUT
    integer(i4),                   intent(in) :: LCyearId
    integer(i4), dimension(:),     intent(in) :: GeoUnitList
    integer(i4), dimension(:),     intent(in) :: GeoUnitKar
    integer(i4), dimension(:),     intent(in) :: LAIUnitList
    real(dp),    dimension(:,:),   intent(in) :: LAILUT

    ! Physiographic L0
    real(dp),    dimension(:),     intent(in) :: slope_emp0
    real(dp),    dimension(:),     intent(in) :: l0_latitude ! l1 ids of l0 cells
    integer(i4), dimension(:),     intent(in) :: cellId0
    integer(i4), dimension(:,:),   intent(in) :: soilId0
    integer(i4), dimension(:),     intent(in) :: L0_LCover_LAI
    integer(i4), dimension(:),     intent(in) :: LCover0
    real(dp),    dimension(:),     intent(in) :: Asp0
    real(dp),    dimension(:),     intent(in) :: LAI0
    integer(i4), dimension(:),     intent(in) :: geoUnit0

    integer(i4), dimension(:),     intent(in) :: SDB_is_present
    integer(i4), dimension(:),     intent(in) :: SDB_nHorizons
    integer(i4), dimension(:),     intent(in) :: SDB_nTillHorizons
    real(dp),    dimension(:,:),   intent(in) :: SDB_sand
    real(dp),    dimension(:,:),   intent(in) :: SDB_clay
    real(dp),    dimension(:,:),   intent(in) :: SDB_DbM
    real(dp),    dimension(:,:,:), intent(in) :: SDB_Wd
    real(dp),    dimension(:),     intent(in) :: SDB_RZdepth

    ! Physiographic L1
    integer(i4), dimension(:),     intent(in) :: nTCells0_inL1
    integer(i4), dimension(:),     intent(in) :: L0upBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0downBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0leftBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0rightBound_inL1
    real(dp),    dimension(:),     intent(in) :: latitude

    ! Forcings
    real(dp),    dimension(:),     intent(in) :: evap_coeff
    real(dp),    dimension(:),     intent(in) :: fday_prec
    real(dp),    dimension(:),     intent(in) :: fnight_prec
    real(dp),    dimension(:),     intent(in) :: fday_pet
    real(dp),    dimension(:),     intent(in) :: fnight_pet
    real(dp),    dimension(:),     intent(in) :: fday_temp
    real(dp),    dimension(:),     intent(in) :: fnight_temp
    real(dp),    dimension(:,:,:), intent(in) :: temp_weights
    real(dp),    dimension(:,:,:), intent(in) :: pet_weights
    real(dp),    dimension(:,:,:), intent(in) :: pre_weights
    logical,                       intent(in) :: read_meteo_weights
    real(dp),    dimension(:),     intent(in) :: pet_in
    real(dp),    dimension(:),     intent(in) :: tmin_in
    real(dp),    dimension(:),     intent(in) :: tmax_in
    real(dp),    dimension(:),     intent(in) :: netrad_in
    real(dp),    dimension(:),     intent(in) :: absvappres_in
    real(dp),    dimension(:),     intent(in) :: windspeed_in
    real(dp),    dimension(:),     intent(in) :: prec_in
    real(dp),    dimension(:),     intent(in) :: temp_in

    ! Configuration
    integer(i4),                   intent(inout) ::  yId

    ! Land cover L1 and
    real(dp), dimension(:),        intent(inout) :: fForest1
    real(dp), dimension(:),        intent(inout) :: fPerm1
    real(dp), dimension(:),        intent(inout) :: fSealed1

    ! States
    real(dp),  dimension(:),       intent(inout) :: interc
    real(dp),  dimension(:),       intent(inout) :: snowpack
    real(dp),  dimension(:),       intent(inout) :: sealedStorage
    real(dp),  dimension(:,:),     intent(inout) :: soilMoisture
    real(dp),  dimension(:),       intent(inout) :: unsatStorage
    real(dp),  dimension(:),       intent(inout) :: satStorage
    real(dp),  dimension(:),       intent(inout) :: neutrons

    ! Fluxes L1
    real(dp),  dimension(:),       intent(inout) :: pet_calc
    real(dp),  dimension(:,:),     intent(inout) :: aet_soil
    real(dp),  dimension(:),       intent(inout) :: aet_canopy
    real(dp),  dimension(:),       intent(inout) :: aet_sealed
    real(dp),  dimension(:),       intent(inout) :: baseflow
    real(dp),  dimension(:,:),     intent(inout) :: infiltration
    real(dp),  dimension(:),       intent(inout) :: fast_interflow
    real(dp),  dimension(:),       intent(inout) :: melt
    real(dp),  dimension(:),       intent(inout) :: perc
    real(dp),  dimension(:),       intent(inout) :: prec_effect
    real(dp),  dimension(:),       intent(inout) :: rain
    real(dp),  dimension(:),       intent(inout) :: runoff_sealed
    real(dp),  dimension(:),       intent(inout) :: slow_interflow
    real(dp),  dimension(:),       intent(inout) :: snow
    real(dp),  dimension(:),       intent(inout) :: throughfall
    real(dp),  dimension(:),       intent(inout) :: total_runoff

    ! Effective Parameters
    real(dp), dimension(:),        intent(inout) ::  alpha
    real(dp), dimension(:),        intent(inout) ::  deg_day_incr
    real(dp), dimension(:),        intent(inout) ::  deg_day_max
    real(dp), dimension(:),        intent(inout) ::  deg_day_noprec
    real(dp), dimension(:),        intent(inout) ::  deg_day
    real(dp), dimension(:),        intent(inout) ::  fAsp
    real(dp), dimension(:),        intent(inout) ::  petLAIcorFactorL1
    real(dp), dimension(:),        intent(inout) ::  HarSamCoeff
    real(dp), dimension(:,:),      intent(inout) ::  PrieTayAlpha
    real(dp), dimension(:,:),      intent(inout) ::  aeroResist
    real(dp), dimension(:,:),      intent(inout) ::  surfResist
    real(dp), dimension(:,:),      intent(inout) ::  frac_roots
    real(dp), dimension(:),        intent(inout) ::  interc_max
    real(dp), dimension(:),        intent(inout) ::  karst_loss
    real(dp), dimension(:),        intent(inout) ::  k0
    real(dp), dimension(:),        intent(inout) ::  k1
    real(dp), dimension(:),        intent(inout) ::  k2
    real(dp), dimension(:),        intent(inout) ::  kp
    real(dp), dimension(:,:),      intent(inout) ::  soil_moist_FC
    real(dp), dimension(:,:),      intent(inout) ::  soil_moist_sat
    real(dp), dimension(:,:),      intent(inout) ::  soil_moist_exponen
    real(dp), dimension(:),        intent(inout) ::  jarvis_thresh_c1
    real(dp), dimension(:),        intent(inout) ::  temp_thresh
    real(dp), dimension(:),        intent(inout) ::  unsat_thresh
    real(dp), dimension(:),        intent(inout) ::  water_thresh_sealed
    real(dp), dimension(:,:),      intent(inout) ::  wilting_point

    ! local
    logical                :: isday       ! is day or night
    integer(i4)            :: hour        ! current hour of a given day
    integer(i4)            :: day         ! day of the month     [1-28 or 1-29 or 1-30 or 1-31]
    integer(i4)            :: month       ! Month of current day [1-12]
    integer(i4)            :: year        ! year
    integer(i4)            :: doy         ! doy of the year [1-365 or 1-366]
    integer(i4)            :: k           ! cell index
    integer(i4)            :: iStart, iEnd 
    real(dp)               :: pet         !
    real(dp)               :: prec        !
    real(dp)               :: temp        !

    ! temporary arrays so that inout of routines is contiguous array
    real(dp), dimension(size(infiltration,2)) :: tmp_infiltration
    real(dp), dimension(size(soilMoisture,2)) :: tmp_soilMoisture
    real(dp), dimension(size(aet_soil,2))     :: tmp_aet_soil

    !-------------------------------------------------------------------
    ! date and month of this timestep
    !-------------------------------------------------------------------
    call dec2date(time, yy=year, mm=month, dd=day, hh=hour)

   !-------------------------------------------------------------------
    ! MPR CALL
    !   -> call only when LC had changed
    !   -> or for very first time step
    !
    ! Variables required prior to the MPR call:
    ! 1) LC fractions
    !      --> time independent variable: to be initalized every time
    !          with landcover change
    ! 2) Impervious area fraction in flood plains
    !      --> time independent variable: to be initalized every time
    !          with landcover change
    !-------------------------------------------------------------------
    if( (LCyearId .NE. yId) .or. (tt .EQ. 1) ) then

       ! abort if land cover change is there and mpr is switched off
       if ( (tt .ne. 1) .and. (.not. perform_mpr) ) then
          call message()
          call message('***ERROR: land cover change detected and mpr is switched off!')
          stop
       end if

        ! update yId to keep track of LC change
        yId = LCyearId

        ! estimate land cover fractions for dominant landcover class
        ! --> time independent variable: to be initalized every time
        !     with landcover change
        fForest1(:) = L0_fractionalCover_in_Lx( LCover0, 1, mask0, &
                                                L0upBound_inL1,    &
                                                L0downBound_inL1,  &
                                                L0leftBound_inL1,  &
                                                L0rightBound_inL1, &
                                                nTCells0_inL1      )
        fSealed1(:) = L0_fractionalCover_in_Lx( LCover0, 2, mask0, &
                                                L0upBound_inL1,    &
                                                L0downBound_inL1,  &
                                                L0leftBound_inL1,  &
                                                L0rightBound_inL1, &
                                                nTCells0_inL1      )
        fPerm1(:)   = L0_fractionalCover_in_Lx( LCover0, 3, mask0, &
                                                L0upBound_inL1,    &
                                                L0downBound_inL1,  &
                                                L0leftBound_inL1,  &
                                                L0rightBound_inL1, &
                                                nTCells0_inL1      )
        !---------------------------------------------------------
        ! Update fractions of sealed area fractions
        ! based on the sealing fraction[0-1] in cities
        !---------------------------------------------------------
        fSealed1(:) = fSealedInCity*fSealed1(:)
        fPerm1(:)   = fPerm1(:) + (1.0_dp - fSealedInCity)*fSealed1(:)

        ! to make sure everything happens smoothly
        fForest1(:) = fForest1(:) / ( fForest1(:) + fSealed1(:)  + fPerm1(:) )
        fSealed1(:) = fSealed1(:) / ( fForest1(:) + fSealed1(:)  + fPerm1(:) )
        fPerm1(:)   = fPerm1(:)   / ( fForest1(:) + fSealed1(:)  + fPerm1(:) )

        !-------------------------------------------------------------------
        ! NOW call MPR
        !-------------------------------------------------------------------
        if ( perform_mpr ) then
           call mpr( processMatrix, iflag_soil_option, global_parameters(:), nodata_dp,   &
                mask0, geoUnit0, GeoUnitList, GeoUnitKar, LAILUT, LAIUnitList,            &
                SDB_is_present, SDB_nHorizons,                                            &
                SDB_nTillHorizons, SDB_sand, SDB_clay, SDB_DbM, SDB_Wd, SDB_RZdepth,      &
                nHorizons_mHM, horizon_depth, c2TSTu, fForest1, fSealed1, fPerm1,         &
                soilId0, Asp0, L0_LCover_LAI, LCover0,                                    &
                slope_emp0, cellId0,                                                      &
                L0upBound_inL1, L0downBound_inL1, L0leftBound_inL1,                       &
                L0rightBound_inL1, nTCells0_inL1, l0_latitude,                            &
                alpha, deg_day_incr, deg_day_max, deg_day_noprec,                         &
                fAsp, HarSamCoeff(:), PrieTayAlpha(:,:), aeroResist(:,:),                 &
                surfResist(:,:), frac_roots, k0, k1, k2, kp, karst_loss,                  &
                soil_moist_FC, soil_moist_sat, soil_moist_exponen, jarvis_thresh_c1(:),   &
                temp_thresh, unsat_thresh, water_thresh_sealed, wilting_point            )
        end if
    
        !-------------------------------------------------------------------
        ! Update the inital states of soil water content for the first time
        ! step and when perform_mpr = FALSE
        ! based on the half of the derived values of Field capacity
        ! other states are kept at their inital values
        !-------------------------------------------------------------------
        if( (tt .EQ. 1) .AND. ( .not. read_states ) ) then
          soilMoisture(:,:) = 0.5_dp*soil_moist_FC(:,:)
        end if

    end if

    !-------------------------------------------------------------------
    ! CALL regionalization of parameters related to LAI
    ! IT is now outside of mHM since LAI is now dynamic variable
    !-------------------------------------------------------------------
    select case(timeStep_LAI_input)
    case(0:1) ! in both case of option with 0 and 1
       ! Estimate max. intecept. capacity based on long term monthly mean LAI values
       ! Max. interception is updated every month rather than every day
       if( (tt .EQ. 1) .OR. (month .NE. counter_month) ) then
          call canopy_intercept_param( processMatrix, global_parameters(:), &
               LAI0, nTCells0_inL1, L0upBound_inL1, &
               L0downBound_inL1, L0leftBound_inL1,  &
               L0rightBound_inL1, cellId0, mask0,   &
               nodata_dp,  interc_max               )

          if (processMatrix(5,1) == -1) then 
             iStart = processMatrix(5,3) - processMatrix(5,2) + 1
             iEnd   = processMatrix(5,3)  

             call pet_correctbyLAI(global_parameters(iStart:iEnd), nodata_dp, LCOVER0, LAI0, mask0, cellId0, &
                  L0upBound_inL1, L0downBound_inL1, L0leftBound_inL1, L0rightBound_inL1, nTCells0_inL1, petLAIcorFactorL1)
          end if
       end if
       ! Estimate max. inteception based on daily LAI values
    case(-1) ! daily
       if ( (tt .EQ. 1) .OR. (day .NE. counter_day) ) then
          call canopy_intercept_param( processMatrix, global_parameters(:), &
               LAI0, nTCells0_inL1, L0upBound_inL1, &
               L0downBound_inL1, L0leftBound_inL1,  &
               L0rightBound_inL1, cellId0, mask0,   &
               nodata_dp,  interc_max               )

          if (processMatrix(5,1) == -1) then 
             iStart = processMatrix(5,3) - processMatrix(5,2) + 1
             iEnd   = processMatrix(5,3)  

             call pet_correctbyLAI(global_parameters(iStart:iEnd), nodata_dp, LCOVER0, LAI0, mask0, cellId0, &
                  L0upBound_inL1, L0downBound_inL1, L0leftBound_inL1, L0rightBound_inL1, nTCells0_inL1, petLAIcorFactorL1)
          end if
       end if
    case(-2) ! monthly
       if ( (tt .EQ. 1) .OR. (month .NE. counter_month) ) then
          call canopy_intercept_param( processMatrix, global_parameters(:), &
               LAI0, nTCells0_inL1, L0upBound_inL1, &
               L0downBound_inL1, L0leftBound_inL1,  &
               L0rightBound_inL1, cellId0, mask0,   &
               nodata_dp,  interc_max               )

          if (processMatrix(5,1) == -1) then 
             iStart = processMatrix(5,3) - processMatrix(5,2) + 1
             iEnd   = processMatrix(5,3)  

             call pet_correctbyLAI(global_parameters(iStart:iEnd), nodata_dp, LCOVER0, LAI0, mask0, cellId0, &
                  L0upBound_inL1, L0downBound_inL1, L0leftBound_inL1, L0rightBound_inL1, nTCells0_inL1, petLAIcorFactorL1)
          end if
       end if
    case(-3) ! yearly
       if ( (tt .EQ. 1) .OR. (year .NE. counter_year) ) then
          call canopy_intercept_param( processMatrix, global_parameters(:), &
               LAI0, nTCells0_inL1, L0upBound_inL1, &
               L0downBound_inL1, L0leftBound_inL1,  &
               L0rightBound_inL1, cellId0, mask0,   &
               nodata_dp,  interc_max               )

          if (processMatrix(5,1) == -1) then
             iStart = processMatrix(5,3) - processMatrix(5,2) + 1
             iEnd   = processMatrix(5,3)  
             
             call pet_correctbyLAI(global_parameters(iStart:iEnd), nodata_dp, LCOVER0, LAI0, mask0, cellId0, &
                  L0upBound_inL1, L0downBound_inL1, L0leftBound_inL1, L0rightBound_inL1, nTCells0_inL1, petLAIcorFactorL1)
          end if
       end if

    case default ! no output at all
       continue
    end select

    !-------------------------------------------------------------------
    ! flag for day or night depending on hours of the day
    !-------------------------------------------------------------------
    isday = ( hour .gt. 6 ) .AND. ( hour .le. 18 )

    !-------------------------------------------------------------------
    ! HYDROLOGICAL PROCESSES at L1-LEVEL
    !-------------------------------------------------------------------

    !$OMP parallel default(shared) &
    !$OMP private(k, prec, pet, temp, tmp_soilmoisture, tmp_infiltration, tmp_aet_soil)
    !$OMP do SCHEDULE(STATIC)
    do k = 1, nCells1

       ! PET calculation
       select case (processMatrix(5,1))
       case(-1) ! PET is input ! correct pet for every day only once at the first time step
          pet =  petLAIcorFactorL1(k) * pet_in(k)

       case(0) ! PET is input ! correct pet for every day only once at the first time step
          pet =  fAsp(k) * pet_in(k)

       case(1) ! Hargreaves-Samani
          ! estimate day of the year (doy) for approximation of the extraterrestrial radiation
          doy = nint(date2dec(day,month,year,12) - date2dec(1,1,year,12) ) + 1
          
          if (tmax_in(k) .lt. tmin_in(k)) call message('WARNING: tmax smaller than tmin at doy ', &
               num2str(doy), ' in year ', num2str(year),' at cell', num2str(k),'!')

          pet = fAsp(k) * pet_hargreaves(HarSamCoeff(k), HarSamConst,  temp_in(k), tmax_in(k),   &
               tmin_in(k), latitude(k), doy)

       case(2) ! Priestley-Taylor
           ! Priestley Taylor is not defined for values netrad < 0.0_dp
          pet = pet_priestly( PrieTayAlpha(k,month), max(netrad_in(k), 0.0_dp), temp_in(k))

       case(3) ! Penman-Monteith
          pet = pet_penman  (max(netrad_in(k), 0.0_dp), temp_in(k), absvappres_in(k)/1000.0_dp, &
               aeroResist(k,month) / windspeed_in(k), surfResist(k,month), 1.0_dp, 1.0_dp)
  
       end select

       ! temporal disaggreagtion of forcing variables
       call temporal_disagg_forcing( isday, ntimesteps_day, prec_in(k),                        & ! Intent IN
            pet, temp_in(k), fday_prec(month), fday_pet(month),                                & ! Intent IN
            fday_temp(month), fnight_prec(month), fnight_pet(month), fnight_temp(month),       & ! Intent IN
            temp_weights(k, month, hour + 1), pet_weights(k, month, hour + 1),                 & ! Intent IN
            pre_weights(k, month, hour + 1),                                                   & ! Intent IN
            read_meteo_weights,                                                                & ! Intent IN
            prec, pet_calc(k), temp )                                                            ! Intent OUT

       call canopy_interc( pet_calc(k), interc_max(k), prec,                                   & ! Intent IN
            interc(k),                                                                         & ! Intent INOUT
            throughfall(k), aet_canopy(k) )                                                      ! Intent OUT

       call snow_accum_melt( deg_day_incr(k), deg_day_max(k),                                  & ! Intent IN
            deg_day_noprec(k), prec, temp, temp_thresh(k), throughfall(k),                     & ! Intent IN
            snowpack(k),                                                                       & ! Intent INOUT
            deg_day(k),                                                                        & ! Intent OUT
            melt(k), prec_effect(k), rain(k), snow(k) )                                          ! Intent OUT

       tmp_soilMoisture(:) = soilMoisture(k,:)
       tmp_infiltration(:) = infiltration(k,:)
       
       call soil_moisture(processMatrix(3,1),                                                  & ! Intent IN
            fSealed1(k), water_thresh_sealed(k),                                               & ! Intent IN
            pet_calc(k), evap_coeff(month), soil_moist_sat(k,:), frac_roots(k,:),              & ! Intent IN
            soil_moist_FC(k,:), wilting_point(k,:), soil_moist_exponen(k,:),                   & ! Intent IN
            jarvis_thresh_c1(k), aet_canopy(k),                                                & ! Intent IN
            prec_effect(k), runoff_sealed(k), sealedStorage(k),                                & ! Intent INOUT
            tmp_infiltration(:), tmp_soilMoisture(:),                                          & ! Intent INOUT
            tmp_aet_soil(:), aet_sealed(k) )                                                     ! Intent OUT
       infiltration(k,:) = tmp_infiltration(:)
       soilMoisture(k,:) = tmp_soilMoisture(:)
       aet_soil(k,:)     = tmp_aet_soil(:)

       call runoff_unsat_zone( k1(k), kp(k), k0(k), alpha(k), karst_loss(k),                   & ! Intent IN
            infiltration(k, nHorizons_mHM),  unsat_thresh(k),                                  & ! Intent IN
            satStorage(k), unsatStorage(k),                                                    & ! Intent INOUT
            slow_interflow(k), fast_interflow(k), perc(k) )                                      ! Intent OUT

       call runoff_sat_zone( k2(k),                                                            & ! Intent IN
            satStorage(k),                                                                     & ! Intent INOUT
            baseflow(k) )                                                                        ! Intent OUT

       call L1_total_runoff( fSealed1(k), fast_interflow(k), slow_interflow(k), baseflow(k),   & ! Intent IN
            runoff_sealed(k),                                                                  & ! Intent IN
            total_runoff(k) )                                                                    ! Intent OUT

       !-------------------------------------------------------------------
       ! Nested model: Neutrons state variable, related to soil moisture   
       !-------------------------------------------------------------------

       ! based on soilMoisture
       ! TODO they again loop over all cells. Maybe move this to line 680 in the loop used above?
       if ( processMatrix(10, 1) .eq. 1 ) &
           call DesiletsN0( soilMoisture(k,:), horizon_depth(:), &
                           global_parameters(processMatrix(10,3)-processMatrix(10,2)+1), &
                           neutrons(k))
       if ( processMatrix(10, 1) .eq. 2 ) &
           call COSMIC( soilMoisture(k,:), horizon_depth(:), &
                       global_parameters(processMatrix(10,3)-processMatrix(10,2)+2:processMatrix(10,3)), &
                       neutron_integral_AFast(:), &
                       neutrons(k))
    end do
    !$OMP end do
    !$OMP end parallel

  end subroutine mHM

END MODULE mo_mHM
