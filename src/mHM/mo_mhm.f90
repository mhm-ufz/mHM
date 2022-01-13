!>       \file mo_mhm.f90

!>       \brief Call all main processes of mHM.

!>       \details This module calls all processes of mHM for a given configuration.
!>       The configuration of the model is stored in the a process matrix.
!>       This configuration is specified in the namelist mhm.nml.

!>       The processes are executed in ascending order. At the moment only
!>       process 5 and 8 have options.

!>       Currently the following processes are implemented:

!>       Process    | Name                      | Flag  | Description
!>       ---------- | ------------------------- | ----- | ------------------------------------------
!>       1          | interception              |   1   | Maximum interception
!>       2          | snow and melting          |   1   | Degree-day
!>       3          | soil moisture             |   1   | Feddes equation for ET reduction, Brooks-Corey like
!>       3          | soil moisture             |   2   | Jarvis equation for ET reduction, Brooks-Corey like
!>       3          | soil moisture             |   3   | Jarvis eq. for ET red. + FC dependency on root frac. coef.
!>       4          | direct runoff             |   1   | Linear reservoir exceedance
!>       5          | PET                       |  -1   | PET is input, LAI based correction, dynamic scaling func.
!>       5          | PET                       |   0   | PET is input, Aspect based correction
!>       5          | PET                       |   1   | Hargreaves-Samani
!>       5          | PET                       |   2   | Priestley-Taylor
!>       5          | PET                       |   3   | Penman-Monteith
!>       6          | interflow                 |   1   | Nonlinear reservoir with saturation excess
!>       7          | percolation and base flow |   1   | GW linear reservoir
!>       8          | routing                   |   0   | no routing
!>       8          | routing                   |   1   | use mRM i.e. Muskingum
!>       8          | routing                   |   2   | use mRM i.e. adaptive timestep

!>       \authors Luis Samaniego

!>       \date Dec 2012

! Modifications:

MODULE mo_mHM

  use mo_kind, only : i4, dp
  use mo_message, only : message
  !$ USE omp_lib

  IMPLICIT NONE

  PUBLIC :: mHM      ! initialization sequence

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------

  !    NAME
  !        mHM

  !    PURPOSE
  !>       \brief Pure mHM calculations.

  !>       \details Pure mHM calculations. All variables are allocated and initialized.
  !>       They will be local variables within this call.

  !    INTENT(IN)
  !>       \param[in] "logical :: read_states"                           indicated whether states have been read from
  !>       file
  !>       \param[in] "integer(i4) :: tt"                                simulation time step
  !>       \param[in] "real(dp) :: time"                                 current decimal Julian day
  !>       \param[in] "integer(i4), dimension(:, :) :: processMatrix"    mHM process configuration matrix
  !>       \param[in] "real(dp), dimension(:) :: horizon_depth"          Depth of each horizon in mHM
  !>       \param[in] "integer(i4) :: nCells1"                           number of cells in a given domain at level L1
  !>       \param[in] "integer(i4) :: nHorizons_mHM"                     Number of Horizons in mHM
  !>       \param[in] "real(dp) :: ntimesteps_day"                       number of time intervals per day, transformed
  !>       in dp
  !>       \param[in] "real(dp), dimension(:) :: neutron_integral_AFast" tabular for neutron flux approximation
  !>       \param[in] "real(dp), dimension(:) :: global_parameters"      global mHM parameters
  !>       \param[in] "real(dp), dimension(:) :: latitude"               latitude on level 1
  !>       \param[in] "real(dp), dimension(:) :: evap_coeff"             Evaporation coefficent for free-water surface
  !>       of that current month
  !>       \param[in] "real(dp), dimension(:) :: fday_prec"              [-] day ratio precipitation < 1
  !>       \param[in] "real(dp), dimension(:) :: fnight_prec"            [-] night ratio precipitation < 1
  !>       \param[in] "real(dp), dimension(:) :: fday_pet"               [-] day ratio PET  < 1
  !>       \param[in] "real(dp), dimension(:) :: fnight_pet"             [-] night ratio PET  < 1
  !>       \param[in] "real(dp), dimension(:) :: fday_temp"              [-] day factor mean temp
  !>       \param[in] "real(dp), dimension(:) :: fnight_temp"            [-] night factor mean temp
  !>       \param[in] "real(dp), dimension(:, :, :) :: temp_weights"     multiplicative weights for temperature (deg K)
  !>       \param[in] "real(dp), dimension(:, :, :) :: pet_weights"      multiplicative weights for potential
  !>       evapotranspiration
  !>       \param[in] "real(dp), dimension(:, :, :) :: pre_weights"      multiplicative weights for precipitation
  !>       \param[in] "logical :: read_meteo_weights"                    flag whether weights for tavg and pet have read
  !>       and should be used
  !>       \param[in] "real(dp), dimension(:) :: pet_in"                 [mm d-1] Daily potential evapotranspiration
  !>       (input)
  !>       \param[in] "real(dp), dimension(:) :: tmin_in"                [degc]   Daily minimum temperature
  !>       \param[in] "real(dp), dimension(:) :: tmax_in"                [degc]   Daily maxumum temperature
  !>       \param[in] "real(dp), dimension(:) :: netrad_in"              [w m2]   Daily average net radiation
  !>       \param[in] "real(dp), dimension(:) :: absvappres_in"          [Pa]     Daily average absolute vapour pressure
  !>       \param[in] "real(dp), dimension(:) :: windspeed_in"           [m s-1]  Daily average wind speed
  !>       \param[in] "real(dp), dimension(:) :: prec_in"                [mm d-1] Daily mean precipitation
  !>       \param[in] "real(dp), dimension(:) :: temp_in"                [degc]   Daily average temperature

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:) :: fSealed1"              fraction of sealed area at scale L1
  !>       \param[inout] "real(dp), dimension(:) :: interc"                Interception
  !>       \param[inout] "real(dp), dimension(:) :: snowpack"              Snowpack
  !>       \param[inout] "real(dp), dimension(:) :: sealedStorage"         Retention storage of impervious areas
  !>       \param[inout] "real(dp), dimension(:, :) :: soilMoisture"       Soil moisture of each horizon
  !>       \param[inout] "real(dp), dimension(:) :: unsatStorage"          Upper soil storage
  !>       \param[inout] "real(dp), dimension(:) :: satStorage"            Groundwater storage
  !>       \param[inout] "real(dp), dimension(:) :: neutrons"              Ground albedo neutrons
  !>       \param[inout] "real(dp), dimension(:) :: pet_calc"              [mm TS-1] estimated PET (if PET is input =
  !>       corrected values (fAsp*PET))
  !>       \param[inout] "real(dp), dimension(:, :) :: aet_soil"           actual ET
  !>       \param[inout] "real(dp), dimension(:) :: aet_canopy"            Real evaporation intensity from canopy
  !>       \param[inout] "real(dp), dimension(:) :: aet_sealed"            Actual ET from free-water surfaces
  !>       \param[inout] "real(dp), dimension(:) :: baseflow"              Baseflow
  !>       \param[inout] "real(dp), dimension(:, :) :: infiltration"       Recharge, infiltration intensity or effective
  !>       precipitation of each horizon
  !>       \param[inout] "real(dp), dimension(:) :: fast_interflow"        Fast runoff component
  !>       \param[inout] "real(dp), dimension(:) :: melt"                  Melting snow depth
  !>       \param[inout] "real(dp), dimension(:) :: perc"                  Percolation
  !>       \param[inout] "real(dp), dimension(:) :: prec_effect"           Effective precipitation depth (snow melt +
  !>       rain)
  !>       \param[inout] "real(dp), dimension(:) :: rain"                  Rain precipitation depth
  !>       \param[inout] "real(dp), dimension(:) :: runoff_sealed"         Direct runoff from impervious areas
  !>       \param[inout] "real(dp), dimension(:) :: slow_interflow"        Slow runoff component
  !>       \param[inout] "real(dp), dimension(:) :: snow"                  Snow precipitation depth
  !>       \param[inout] "real(dp), dimension(:) :: throughfall"           Throughfall
  !>       \param[inout] "real(dp), dimension(:) :: total_runoff"          Generated runoff
  !>       \param[inout] "real(dp), dimension(:) :: alpha"                 Exponent for the upper reservoir
  !>       \param[inout] "real(dp), dimension(:) :: deg_day_incr"          Increase of the Degree-day factor per mm of
  !>       increase in precipitation
  !>       \param[inout] "real(dp), dimension(:) :: deg_day_max"           Maximum Degree-day factor
  !>       \param[inout] "real(dp), dimension(:) :: deg_day_noprec"        Degree-day factor with no precipitation
  !>       \param[inout] "real(dp), dimension(:) :: deg_day"               Degree-day factor
  !>       \param[inout] "real(dp), dimension(:) :: fAsp"                  [1]     PET correction for Aspect at level 1
  !>       \param[inout] "real(dp), dimension(:) :: petLAIcorFactorL1"     PET correction factor based on LAI at level 1
  !>       \param[inout] "real(dp), dimension(:) :: HarSamCoeff"           [1]     PET Hargreaves Samani coefficient at
  !>       level 1
  !>       \param[inout] "real(dp), dimension(:) :: PrieTayAlpha"          [1]     PET Priestley Taylor coefficient at
  !>       level 1
  !>       \param[inout] "real(dp), dimension(:) :: aeroResist"            [s m-1] PET aerodynamical resitance at level
  !>       1
  !>       \param[inout] "real(dp), dimension(:) :: surfResist"            [s m-1] PET bulk surface resitance at level 1
  !>       \param[inout] "real(dp), dimension(:, :) :: frac_roots"         Fraction of Roots in soil horizon
  !>       \param[inout] "real(dp), dimension(:) :: interc_max"            Maximum interception
  !>       \param[inout] "real(dp), dimension(:) :: karst_loss"            Karstic percolation loss
  !>       \param[inout] "real(dp), dimension(:) :: k0"                    Recession coefficient of the upper reservoir,
  !>       upper outlet
  !>       \param[inout] "real(dp), dimension(:) :: k1"                    Recession coefficient of the upper reservoir,
  !>       lower outlet
  !>       \param[inout] "real(dp), dimension(:) :: k2"                    Baseflow recession coefficient
  !>       \param[inout] "real(dp), dimension(:) :: kp"                    Percolation coefficient
  !>       \param[inout] "real(dp), dimension(:, :) :: soil_moist_FC"      Soil moisture below which actual ET is
  !>       reduced
  !>       \param[inout] "real(dp), dimension(:, :) :: soil_moist_sat"     Saturation soil moisture for each horizon
  !>       [mm]
  !>       \param[inout] "real(dp), dimension(:, :) :: soil_moist_exponen" Exponential parameter to how non-linear is
  !>       the soil water retention
  !>       \param[inout] "real(dp), dimension(:) :: jarvis_thresh_c1"      jarvis critical value for normalized soil
  !>       water content
  !>       \param[inout] "real(dp), dimension(:) :: temp_thresh"           Threshold temperature for snow/rain
  !>       \param[inout] "real(dp), dimension(:) :: unsat_thresh"          Threshold water depth in upper reservoir
  !>       \param[inout] "real(dp), dimension(:) :: water_thresh_sealed"   Threshold water depth in impervious areas
  !>       \param[inout] "real(dp), dimension(:, :) :: wilting_point"      Permanent wilting point for each horizon

  !    HISTORY
  !>       \authors Luis Samaniego & Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Luis Samaniego, Rohini Kumar    Dec 2012 - modularization
  ! Luis Samaniego                  Feb 2013 - call routine
  ! Rohini Kumar                    Feb 2013 - MPR call and other pre-requisite variables for this call
  ! Rohini Kumar                    May 2013 - Error checks
  ! Rohini Kumar                    Jun 2013 - sealed area correction in total runoff
  !                                          - initalization of soil moist. at first timestep
  ! Rohini Kumar                    Aug 2013 - dynamic LAI option included, and changed within the code
  !                                            made accordingly (e.g., canopy intecpt.)
  !                                          - max. canopy interception is estimated outside of MPR call
  ! Matthias Zink                   Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  ! Matthias Zink                   Mar 2014 - added inflow from upstream areas
  ! Matthias Zink                   Apr 2014 - added PET calculation: Priestley-Taylor and Penman-Monteith
  !                                            and its parameterization (Process 5)
  ! Rohini Kumar                    Apr 2014 - mHM run with a single L0 grid cell, also in the routing mode
  ! Stephan Thober                  Jun 2014 - added flag for switching of MPR
  ! Matthias Cuntz & Juliane Mai    Nov 2014 - LAI input from daily, monthly or yearly files
  ! Matthias Zink                   Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober                  Aug 2015 - moved routing to mRM
  ! Rohini Kumar                    Mar 2016 - changes for handling multiple soil database options
  ! Rohini Kumar                    Dec 2016 - changes for reading gridded mean monthly LAI fields
  ! Stephan Thober                  Jan 2017 - added prescribed weights for tavg and pet
  ! Zink M. Demirel C.              Mar 2017 - added Jarvis soil water stress function at SM process(3)
  ! M.Cuneyd Demirel & Simon Stisen May 2017 - added FC dependency on root fraction coef. at SM process(3)
  ! M.Cuneyd Demirel & Simon Stisen Jun 2017 - added PET correction based on LAI at PET process(5)
  ! Robert Schweppe, Stephan Thober Nov 2017 - moved call to MPR to mhm_eval
  ! Robert Schweppe                 Jun 2018 - refactoring and reformatting
  ! Robert Schweppe                 Nov 2018 - added c2TSTu for unit conversion (moved here from MPR)
  ! Stephan Thober                  Jan 2022 - added is_hourly_forcing

  subroutine mHM(read_states, is_hourly_forcing, tt, time, processMatrix, horizon_depth, nCells1, nHorizons_mHM, ntimesteps_day, &
                c2TSTu, neutron_integral_AFast, &
                global_parameters, latitude, evap_coeff, fday_prec, fnight_prec, fday_pet, &
                fnight_pet, fday_temp, fnight_temp, temp_weights, pet_weights, pre_weights, read_meteo_weights, pet_in, &
                tmin_in, tmax_in, netrad_in, absvappres_in, windspeed_in, prec_in, temp_in, fSealed1, interc, snowpack, &
                sealedStorage, soilMoisture, unsatStorage, satStorage, neutrons, pet_calc, aet_soil, aet_canopy, &
                aet_sealed, baseflow, infiltration, fast_interflow, melt, perc, prec_effect, rain, runoff_sealed, &
                slow_interflow, snow, throughfall, total_runoff, alpha, deg_day_incr, deg_day_max, deg_day_noprec, &
                deg_day, fAsp, petLAIcorFactorL1, HarSamCoeff, PrieTayAlpha, aeroResist, surfResist, frac_roots, &
                interc_max, karst_loss, k0, k1, k2, kp, soil_moist_FC, soil_moist_sat, soil_moist_exponen, &
                jarvis_thresh_c1, temp_thresh, unsat_thresh, water_thresh_sealed, wilting_point)

    use mo_Temporal_Disagg_Forcing, only : Temporal_Disagg_Forcing
    use mo_canopy_interc, only : canopy_interc
    use mo_julian, only : date2dec, dec2date
    use mo_mhm_constants, only : HarSamConst
    use mo_neutrons, only : COSMIC, DesiletsN0
    use mo_pet, only : pet_hargreaves, pet_penman, pet_priestly
    use mo_runoff, only : L1_total_runoff, runoff_sat_zone, runoff_unsat_zone
    use mo_snow_accum_melt, only : snow_accum_melt
    use mo_soil_moisture, only : soil_moisture
    use mo_string_utils, only : num2str

    implicit none

    ! indicated whether states have been read from file
    logical, intent(in) :: read_states

    ! indicate whether forcing is hourly timestep
    logical, intent(in) :: is_hourly_forcing

    ! simulation time step
    integer(i4), intent(in) :: tt

    ! current decimal Julian day
    real(dp), intent(in) :: time

    ! mHM process configuration matrix
    integer(i4), dimension(:, :), intent(in) :: processMatrix

    ! Depth of each horizon in mHM
    real(dp), dimension(:), intent(in) :: horizon_depth

    ! number of cells in a given domain at level L1
    integer(i4), intent(in) :: nCells1

    ! Number of Horizons in mHM
    integer(i4), intent(in) :: nHorizons_mHM

    ! number of time intervals per day, transformed in dp
    real(dp), intent(in) :: ntimesteps_day

    ! unit conversion
    real(dp), intent(in) :: c2TSTu

    ! tabular for neutron flux approximation
    real(dp), dimension(:), intent(in) :: neutron_integral_AFast

    ! global mHM parameters
    real(dp), dimension(:), intent(in) :: global_parameters

    ! latitude on level 1
    real(dp), dimension(:), intent(in) :: latitude

    ! Evaporation coefficent for free-water surface of that current month
    real(dp), dimension(:), intent(in) :: evap_coeff

    ! [-] day ratio precipitation < 1
    real(dp), dimension(:), intent(in) :: fday_prec

    ! [-] night ratio precipitation < 1
    real(dp), dimension(:), intent(in) :: fnight_prec

    ! [-] day ratio PET  < 1
    real(dp), dimension(:), intent(in) :: fday_pet

    ! [-] night ratio PET  < 1
    real(dp), dimension(:), intent(in) :: fnight_pet

    ! [-] day factor mean temp
    real(dp), dimension(:), intent(in) :: fday_temp

    ! [-] night factor mean temp
    real(dp), dimension(:), intent(in) :: fnight_temp

    ! multiplicative weights for temperature (deg K)
    real(dp), dimension(:, :, :), intent(in) :: temp_weights

    ! multiplicative weights for potential evapotranspiration
    real(dp), dimension(:, :, :), intent(in) :: pet_weights

    ! multiplicative weights for precipitation
    real(dp), dimension(:, :, :), intent(in) :: pre_weights

    ! flag whether weights for tavg and pet have read and should be used
    logical, intent(in) :: read_meteo_weights

    ! [mm d-1] Daily potential evapotranspiration (input)
    real(dp), dimension(:), intent(in) :: pet_in

    ! [degc]   Daily minimum temperature
    real(dp), dimension(:), intent(in) :: tmin_in

    ! [degc]   Daily maxumum temperature
    real(dp), dimension(:), intent(in) :: tmax_in

    ! [w m2]   Daily average net radiation
    real(dp), dimension(:), intent(in) :: netrad_in

    ! [Pa]     Daily average absolute vapour pressure
    real(dp), dimension(:), intent(in) :: absvappres_in

    ! [m s-1]  Daily average wind speed
    real(dp), dimension(:), intent(in) :: windspeed_in

    ! [mm d-1] Daily mean precipitation
    real(dp), dimension(:), intent(in) :: prec_in

    ! [degc]   Daily average temperature
    real(dp), dimension(:), intent(in) :: temp_in

    ! fraction of sealed area at scale L1
    real(dp), dimension(:), intent(inout) :: fSealed1

    ! Interception
    real(dp), dimension(:), intent(inout) :: interc

    ! Snowpack
    real(dp), dimension(:), intent(inout) :: snowpack

    ! Retention storage of impervious areas
    real(dp), dimension(:), intent(inout) :: sealedStorage

    ! Soil moisture of each horizon
    real(dp), dimension(:, :), intent(inout) :: soilMoisture

    ! Upper soil storage
    real(dp), dimension(:), intent(inout) :: unsatStorage

    ! Groundwater storage
    real(dp), dimension(:), intent(inout) :: satStorage

    ! Ground albedo neutrons
    real(dp), dimension(:), intent(inout) :: neutrons

    ! [mm TS-1] estimated PET (if PET is input = corrected values (fAsp*PET))
    real(dp), dimension(:), intent(inout) :: pet_calc

    ! actual ET
    real(dp), dimension(:, :), intent(inout) :: aet_soil

    ! Real evaporation intensity from canopy
    real(dp), dimension(:), intent(inout) :: aet_canopy

    ! Actual ET from free-water surfaces
    real(dp), dimension(:), intent(inout) :: aet_sealed

    ! Baseflow
    real(dp), dimension(:), intent(inout) :: baseflow

    ! Recharge, infiltration intensity or effective precipitation of each horizon
    real(dp), dimension(:, :), intent(inout) :: infiltration

    ! Fast runoff component
    real(dp), dimension(:), intent(inout) :: fast_interflow

    ! Melting snow depth
    real(dp), dimension(:), intent(inout) :: melt

    ! Percolation
    real(dp), dimension(:), intent(inout) :: perc

    ! Effective precipitation depth (snow melt + rain)
    real(dp), dimension(:), intent(inout) :: prec_effect

    ! Rain precipitation depth
    real(dp), dimension(:), intent(inout) :: rain

    ! Direct runoff from impervious areas
    real(dp), dimension(:), intent(inout) :: runoff_sealed

    ! Slow runoff component
    real(dp), dimension(:), intent(inout) :: slow_interflow

    ! Snow precipitation depth
    real(dp), dimension(:), intent(inout) :: snow

    ! Throughfall
    real(dp), dimension(:), intent(inout) :: throughfall

    ! Generated runoff
    real(dp), dimension(:), intent(inout) :: total_runoff

    ! Exponent for the upper reservoir
    real(dp), dimension(:), intent(inout) :: alpha

    ! Increase of the Degree-day factor per mm of increase in precipitation
    real(dp), dimension(:), intent(inout) :: deg_day_incr

    ! Maximum Degree-day factor
    real(dp), dimension(:), intent(inout) :: deg_day_max

    ! Degree-day factor with no precipitation
    real(dp), dimension(:), intent(inout) :: deg_day_noprec

    ! Degree-day factor
    real(dp), dimension(:), intent(inout) :: deg_day

    ! [1]     PET correction for Aspect at level 1
    real(dp), dimension(:), intent(inout) :: fAsp

    ! PET correction factor based on LAI at level 1
    real(dp), dimension(:), intent(inout) :: petLAIcorFactorL1

    ! [1]     PET Hargreaves Samani coefficient at level 1
    real(dp), dimension(:), intent(inout) :: HarSamCoeff

    ! [1]     PET Priestley Taylor coefficient at level 1
    real(dp), dimension(:), intent(inout) :: PrieTayAlpha

    ! [s m-1] PET aerodynamical resitance at level 1
    real(dp), dimension(:), intent(inout) :: aeroResist

    ! [s m-1] PET bulk surface resitance at level 1
    real(dp), dimension(:), intent(inout) :: surfResist

    ! Fraction of Roots in soil horizon
    real(dp), dimension(:, :), intent(inout) :: frac_roots

    ! Maximum interception
    real(dp), dimension(:), intent(inout) :: interc_max

    ! Karstic percolation loss
    real(dp), dimension(:), intent(inout) :: karst_loss

    ! Recession coefficient of the upper reservoir, upper outlet
    real(dp), dimension(:), intent(inout) :: k0

    ! Recession coefficient of the upper reservoir, lower outlet
    real(dp), dimension(:), intent(inout) :: k1

    ! Baseflow recession coefficient
    real(dp), dimension(:), intent(inout) :: k2

    ! Percolation coefficient
    real(dp), dimension(:), intent(inout) :: kp

    ! Soil moisture below which actual ET is reduced
    real(dp), dimension(:, :), intent(inout) :: soil_moist_FC

    ! Saturation soil moisture for each horizon [mm]
    real(dp), dimension(:, :), intent(inout) :: soil_moist_sat

    ! Exponential parameter to how non-linear is the soil water retention
    real(dp), dimension(:, :), intent(inout) :: soil_moist_exponen

    ! jarvis critical value for normalized soil water content
    real(dp), dimension(:), intent(inout) :: jarvis_thresh_c1

    ! Threshold temperature for snow/rain
    real(dp), dimension(:), intent(inout) :: temp_thresh

    ! Threshold water depth in upper reservoir
    real(dp), dimension(:), intent(inout) :: unsat_thresh

    ! Threshold water depth in impervious areas
    real(dp), dimension(:), intent(inout) :: water_thresh_sealed

    ! Permanent wilting point for each horizon
    real(dp), dimension(:, :), intent(inout) :: wilting_point

    ! is day or night
    logical :: isday

    ! current hour of a given day
    integer(i4) :: hour

    ! day of the month     [1-28 or 1-29 or 1-30 or 1-31]
    integer(i4) :: day

    ! Month of current day [1-12]
    integer(i4) :: month

    ! year
    integer(i4) :: year

    ! doy of the year [1-365 or 1-366]
    integer(i4) :: doy

    ! cell index
    integer(i4) :: k
    ! pet in [mm d-1]
    real(dp) :: pet
    ! precipitation in [mm TS-1]
    real(dp) :: prec

    real(dp) :: temp

    real(dp), dimension(size(infiltration, 2)) :: tmp_infiltration

    real(dp), dimension(size(soilMoisture, 2)) :: tmp_soilMoisture

    real(dp), dimension(size(aet_soil, 2)) :: tmp_aet_soil


    !-------------------------------------------------------------------
    ! date and month of this timestep
    !-------------------------------------------------------------------
    call dec2date(time, yy = year, mm = month, dd = day, hh = hour)

    !-------------------------------------------------------------------
    ! Update the inital states of soil water content for the first time
    ! step and when perform_mpr = FALSE
    ! based on the half of the derived values of Field capacity
    ! other states are kept at their inital values
    !-------------------------------------------------------------------
    if((tt .EQ. 1) .AND. (.not. read_states)) then
      soilMoisture(:, :) = 0.5_dp * soil_moist_FC(:, :)
    end if

    !-------------------------------------------------------------------
    ! flag for day or night depending on hours of the day
    !-------------------------------------------------------------------
    isday = (hour .gt. 6) .AND. (hour .le. 18)

    !-------------------------------------------------------------------
    ! HYDROLOGICAL PROCESSES at L1-LEVEL
    !-------------------------------------------------------------------
    !$OMP parallel default(shared) &
    !$OMP private(k, prec, pet, temp, tmp_soilmoisture, tmp_infiltration, tmp_aet_soil)
    !$OMP do SCHEDULE(STATIC)
    do k = 1, nCells1

      ! PET calculation
      select case (processMatrix(5, 1))
      case(-1) ! PET is input ! correct pet for every day only once at the first time step
        pet = petLAIcorFactorL1(k) * pet_in(k)

      case(0) ! PET is input ! correct pet for every day only once at the first time step
        pet = fAsp(k) * pet_in(k)

      case(1) ! Hargreaves-Samani
        ! estimate day of the year (doy) for approximation of the extraterrestrial radiation
        doy = nint(date2dec(day, month, year, 12) - date2dec(1, 1, year, 12)) + 1

        if (tmax_in(k) .lt. tmin_in(k)) call message('WARNING: tmax smaller than tmin at doy ', &
                num2str(doy), ' in year ', num2str(year), ' at cell', num2str(k), '!')

        pet = fAsp(k) * pet_hargreaves(HarSamCoeff(k), HarSamConst, temp_in(k), tmax_in(k), &
                tmin_in(k), latitude(k), doy)

      case(2) ! Priestley-Taylor

        ! Priestley Taylor is not defined for values netrad < 0.0_dp
        pet = pet_priestly(PrieTayAlpha(k), max(netrad_in(k), 0.0_dp), temp_in(k))

      case(3) ! Penman-Monteith
        pet = pet_penman  (max(netrad_in(k), 0.0_dp), temp_in(k), absvappres_in(k) / 1000.0_dp, &
                aeroResist(k) / windspeed_in(k), surfResist(k), 1.0_dp, 1.0_dp)

      end select
      ! temporal disaggreagtion of forcing variables
      if (is_hourly_forcing) then
         prec = prec_in(k)
         pet_calc(k) = pet
         temp = temp_in(k)
      else
         call temporal_disagg_forcing(isday, ntimesteps_day, prec_in(k), & ! Intent IN
              pet, temp_in(k), fday_prec(month), fday_pet(month), & ! Intent IN
              fday_temp(month), fnight_prec(month), fnight_pet(month), fnight_temp(month), & ! Intent IN
              temp_weights(k, month, hour + 1), pet_weights(k, month, hour + 1), & ! Intent IN
              pre_weights(k, month, hour + 1), & ! Intent IN
              read_meteo_weights, & ! Intent IN
              prec, pet_calc(k), temp)                                                            ! Intent OUT
      end if
      call canopy_interc(pet_calc(k), interc_max(k), prec, & ! Intent IN
              interc(k), & ! Intent INOUT
              throughfall(k), aet_canopy(k))                                                      ! Intent OUT
      call snow_accum_melt(deg_day_incr(k), deg_day_max(k) * c2TSTu, & ! Intent IN
              deg_day_noprec(k) * c2TSTu, prec, temp, temp_thresh(k), throughfall(k), & ! Intent IN
              snowpack(k), & ! Intent INOUT
              deg_day(k), & ! Intent OUT
              melt(k), prec_effect(k), rain(k), snow(k))                                          ! Intent OUT

      tmp_soilMoisture(:) = soilMoisture(k, :)
      tmp_infiltration(:) = infiltration(k, :)

      call soil_moisture(processMatrix(3, 1), & ! Intent IN
              fSealed1(k), water_thresh_sealed(k), & ! Intent IN
              pet_calc(k), evap_coeff(month), soil_moist_sat(k, :), frac_roots(k, :), & ! Intent IN
              soil_moist_FC(k, :), wilting_point(k, :), soil_moist_exponen(k, :), & ! Intent IN
              jarvis_thresh_c1(k), aet_canopy(k), & ! Intent IN
              prec_effect(k), runoff_sealed(k), sealedStorage(k), & ! Intent INOUT
              tmp_infiltration(:), tmp_soilMoisture(:), & ! Intent INOUT
              tmp_aet_soil(:), aet_sealed(k))                                                     ! Intent OUT

      infiltration(k, :) = tmp_infiltration(:)
      soilMoisture(k, :) = tmp_soilMoisture(:)
      aet_soil(k, :) = tmp_aet_soil(:)
      call runoff_unsat_zone(c2TSTu / k1(k), c2TSTu / kp(k), c2TSTu / k0(k), alpha(k), karst_loss(k), & ! Intent IN
              infiltration(k, nHorizons_mHM), unsat_thresh(k), & ! Intent IN
              satStorage(k), unsatStorage(k), & ! Intent INOUT
              slow_interflow(k), fast_interflow(k), perc(k))                                      ! Intent OUT
      call runoff_sat_zone(c2TSTu / k2(k), & ! Intent IN
              satStorage(k), & ! Intent INOUT
              baseflow(k))                                                                        ! Intent OUT
      call L1_total_runoff(fSealed1(k), fast_interflow(k), slow_interflow(k), baseflow(k), & ! Intent IN
              runoff_sealed(k), & ! Intent IN
              total_runoff(k))                                                                    ! Intent OUT

      !-------------------------------------------------------------------
      ! Nested model: Neutrons state variable, related to soil moisture
      !-------------------------------------------------------------------

      ! based on soilMoisture
      if (processMatrix(10, 1) .eq. 1) &
              call DesiletsN0(soilMoisture(k, :), horizon_depth(:), &
                      global_parameters(processMatrix(10, 3) - processMatrix(10, 2) + 1), &
                      neutrons(k))
      if (processMatrix(10, 1) .eq. 2) &
              call COSMIC(soilMoisture(k, :), horizon_depth(:), &
                      global_parameters(processMatrix(10, 3) - processMatrix(10, 2) + 2 : processMatrix(10, 3)), &
                      neutron_integral_AFast(:), &
                      neutrons(k))
    end do
    !$OMP end do
    !$OMP end parallel

  end subroutine mHM

END MODULE mo_mHM
