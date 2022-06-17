!>       \file mo_global_variables.f90

!>       \brief Global variables ONLY used in reading, writing and startup.

!>       \details TODO: add description

!>       \authors Luis Samaniego

!>       \date Dec 2012

! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting


MODULE mo_global_variables

  ! This module provides

  !
  ! Written   Luis Samaniego,     Dec 2005
  ! Modified  Luis Samaniego,     Feb 2013 - new variable names, new modules, units
  !           Rohini Kumar,       Jul 2013 - fraction of perfectly sealed area within city added
  !           Rohini Kumar,       Aug 2013 - name changed "inputFormat" to "inputFormat_meteo_forcings"
  !           Rohini Kumar,       Aug 2013 - name changed from "L0_LAI" to "L0_LCover_LAI"
  !           Rohini Kumar,       Aug 2013 - added dirSoil_LUT and dirGeology_LUT
  !           Luis Samaniego,     Nov 2013 - documentation of dimensions
  !           Matthias Zink,      Nov 2013 - added "InflowGauge" and inflow gauge variabels in Domain
  !           Rohini Kumar,       May 2014 - added options for the model run cordinate system
  !           Stephan Thober,     Jun 2014 - added timeStep_model_inputs and readPer
  !           Stephan Thober,     Jun 2014 - added perform_mpr, updated restart flags
  !           Cuntz M. & Mai J.,  Nov 2014 - LAI input from daily, monthly or yearly files
  !           Matthias Zink,      Dec 2014 - adopted inflow gauges to ignore headwater cells
  !           Matthias Zink,      Mar 2015 - added optional soil mositure readin: dirSoil_moisture, L1_sm
  !           Stephan Thober,     Aug 2015 - moved routing related variables to mRM
  !           Oldrich Rakovec,    Oct 2015 - added definition of Domain averaged TWS data
  !           Rohini Kumar,       Mar 2016 - new variables for handling different soil databases
  !           Johann Brenner,     Feb 2017 - added optional evapotranspiration readin: dirEvapotranspiration, L1_et
  !           Zink M. Demirel C., Mar 2017 - added Jarvis soil water stress variable for SM process(3)
  !           Demirel M.C.        May 2017 - added L1_petLAIcorFactor for PET correction
  !           O. Rakovec, R.Kumar Nov 2017 - added project description for the netcdf outputs
  !           Robert Schweppe,    Dec 2017 - expanded dimensions of effective parameters
  !           Robert Schweppe,    Dec 2017 - merged duplicated variables with mrm into common variables

  USE mo_kind, ONLY : i4, dp
  USE mo_constants, ONLY : YearMonths
  USE mo_mhm_constants, ONLY : nOutFlxState
  USE mo_optimization_types, ONLY : optidata
  use mo_common_variables, only : Grid

  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! DEFINE OUTPUTS
  ! -------------------------------------------------------------------
  integer(i4) :: output_deflate_level   !< deflate level in nc files
  logical :: output_double_precision    !< output precision in nc files
  integer(i4) :: timeStep_model_outputs !< timestep for writing model outputs
  logical, dimension(nOutFlxState) :: outputFlxState         !< Define model outputs see "mhm_outputs.nml"
                                                             !< dim1 = number of output variables to be written
  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mHM
  ! -------------------------------------------------------------------
  integer(i4), dimension(:), allocatable, public :: timeStep_model_inputs !< frequency for reading meteo input
  logical, public :: read_meteo_weights                                   !< read weights for tavg and pet
  character(256), public :: inputFormat_meteo_forcings                    !< format of meteo input data (nc)
  ! Optional data
  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nDomains
  character(256), dimension(:), allocatable, public :: dirPrecipitation   !< Directory where precipitation files are located
  character(256), dimension(:), allocatable, public :: dirTemperature     !< Directory where temperature files are located
  character(256), dimension(:), allocatable, public :: dirMinTemperature  !< Directory where minimum temp. files are located
  character(256), dimension(:), allocatable, public :: dirMaxTemperature  !< Directory where maximum temp. files are located
  character(256), dimension(:), allocatable, public :: dirNetRadiation    !< Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirabsVapPressure  !< Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirwindspeed       !< Directory where windspeed files are located
  character(256), dimension(:), allocatable, public :: dirReferenceET     !< Directory where reference-ET files are located
  ! riv-temp releated
  character(256), dimension(:), allocatable, public :: dirRadiation       !< Directory where short/long-wave rad. files are located

  ! ------------------------------------------------------------------
  ! CONSTANT
  ! ------------------------------------------------------------------
  integer(i4), public, parameter :: routingStates = 2  !< [-]   Routing states (2=current, 1=past)

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type(Grid), dimension(:), allocatable, public :: level2       !< Reference of the metereological variables

  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! Forcings
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  real(dp), public, dimension(:, :, :), allocatable :: L1_temp_weights  !< hourly temperature weights for daily values
  real(dp), public, dimension(:, :, :), allocatable :: L1_pet_weights   !< hourly pet weights for daily values
  real(dp), public, dimension(:, :, :), allocatable :: L1_pre_weights   !< hourly pre weights for daily values
  real(dp), public, dimension(:, :), allocatable :: L1_pre           !< [mm]    Precipitation
  real(dp), public, dimension(:, :), allocatable :: L1_temp          !< [degC]  Air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_pet           !< [mm TS-1] Potential evapotranspiration
  real(dp), public, dimension(:, :), allocatable :: L1_tmin          !< [degC]  minimum daily air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_tmax          !< [degC]  maximum daily air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_netrad        !< [W m2]  net radiation
  real(dp), public, dimension(:, :), allocatable :: L1_absvappress   !< [Pa]    absolute vapour pressure
  real(dp), public, dimension(:, :), allocatable :: L1_windspeed     !< [m s-1] windspeed
  ! riv-temp related
  real(dp), public, dimension(:, :), allocatable :: L1_ssrd          !< [W m2]  short wave radiation
  real(dp), public, dimension(:, :), allocatable :: L1_strd          !< [W m2]  long wave radiation
  real(dp), public, dimension(:, :), allocatable :: L1_tann          !< [degC]  annual mean air temperature


  ! soil moisture
  real(dp), public, dimension(:, :), allocatable :: L1_sm                  !< [-] soil moisture input for optimization
  logical, public, dimension(:, :), allocatable :: L1_sm_mask             !< [-] mask for valid data in L1_sm
  ! neutrons
  real(dp), public, dimension(:, :), allocatable :: L1_neutronsdata            !< [cph] ground albedo neutrons input
  logical, public, dimension(:, :), allocatable :: L1_neutronsdata_mask       !< [cph] mask for valid data in L1_neutrons

  ! soil moisture
  integer(i4) :: nSoilHorizons_sm_input ! No. of mhm soil horizons equivalent to sm input

  type(optidata), public, dimension(:), allocatable :: L1_smObs
  ! neutrons
  type(optidata), public, dimension(:), allocatable :: L1_neutronsObs
  ! evapotranspiration
  type(optidata), public, dimension(:), allocatable :: L1_etObs
  ! tws
  type(optidata), public, dimension(:), allocatable :: L1_twsaObs !< this stores L1_tws, the mask, the directory of the
                                                              !< observerd data, and the
                                                              !< timestepInput of the simulated data
                                                              ! ToDo: add unit
  logical, public                             :: BFI_calc     !< calculate observed BFI from gauges with Eckhardt filter
  real(dp), public, dimension(:), allocatable :: BFI_obs      !< given base-flow index per domain
  real(dp), public, dimension(:), allocatable :: BFI_qBF_sum  !< q2 weighted sum for each domain
  real(dp), public, dimension(:), allocatable :: BFI_qT_sum   !< q2 weighted sum for each domain

  ! State variables
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable :: L1_inter        !< [mm]  Interception
  real(dp), public, dimension(:), allocatable :: L1_snowPack     !< [mm]  Snowpack
  real(dp), public, dimension(:), allocatable :: L1_sealSTW      !< [mm]  Retention storage of impervious areas
  real(dp), public, dimension(:, :), allocatable :: L1_soilMoist !< [mm]  Soil moisture of each horizon
  real(dp), public, dimension(:), allocatable :: L1_unsatSTW     !< [mm]  upper soil storage
  real(dp), public, dimension(:), allocatable :: L1_satSTW       !< [mm]  groundwater storage
  real(dp), public, dimension(:), allocatable :: L1_neutrons     !< [mm]  Ground Albedo Neutrons


  ! Fluxes
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable :: L1_pet_calc     !< [mm TS-1] estimated/corrected potential evapotranspiration
  real(dp), public, dimension(:), allocatable :: L1_temp_calc    !< [degC] temperature for current time step
  real(dp), public, dimension(:), allocatable :: L1_prec_calc    !< [mm TS-1] precipitation for current time step
  real(dp), public, dimension(:, :), allocatable :: L1_aETSoil   !< [mm TS-1] Actual ET from soil layers
  real(dp), public, dimension(:), allocatable :: L1_aETCanopy    !< [mm TS-1] Real evaporation intensity from canopy
  real(dp), public, dimension(:), allocatable :: L1_aETSealed    !< [mm TS-1] Real evap. from free water surfaces
  real(dp), public, dimension(:), allocatable :: L1_baseflow     !< [mm TS-1] Baseflow
  real(dp), public, dimension(:, :), allocatable :: L1_infilSoil !< [mm TS-1] Infiltration intensity each soil horizon
  real(dp), public, dimension(:), allocatable :: L1_fastRunoff   !< [mm TS-1] Fast runoff component
  real(dp), public, dimension(:), allocatable :: L1_melt         !< [mm TS-1] Melting snow depth.
  real(dp), public, dimension(:), allocatable :: L1_percol       !< [mm TS-1] Percolation.
  real(dp), public, dimension(:), allocatable :: L1_preEffect    !< [mm TS-1] Effective precip. depth (snow melt + rain)
  real(dp), public, dimension(:), allocatable :: L1_rain         !< [mm TS-1] Rain precipitation depth
  real(dp), public, dimension(:), allocatable :: L1_runoffSeal   !< [mm TS-1] Direct runoff from impervious areas
  real(dp), public, dimension(:), allocatable :: L1_slowRunoff   !< [mm TS-1] Slow runoff component
  real(dp), public, dimension(:), allocatable :: L1_snow         !< [mm TS-1] Snow precipitation depth
  real(dp), public, dimension(:), allocatable :: L1_Throughfall  !< [mm TS-1] Throughfall.
  real(dp), public, dimension(:), allocatable :: L1_total_runoff !< [m3 TS-1] Generated runoff

  ! -------------------------------------------------------------------
  ! Monthly day/night variation of Meteorological variables
  ! for temporal disaggregation
  ! -------------------------------------------------------------------
  ! dim1 = number of months in a year
  real(dp), public, dimension(int(YearMonths, i4)) :: evap_coeff     !< [-] Evap. coef. for free-water surfaces
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_prec      !< [-] Day ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_prec    !< [-] Night ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_pet       !< [-] Day ratio PET  < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_pet     !< [-] Night ratio PET  < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_temp      !< [-] Day factor mean temp
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_temp    !< [-] Night factor mean temp
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_ssrd      !< [-] Day factor short-wave rad.
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_ssrd    !< [-] Night factor short-wave rad.
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_strd      !< [-] Day factor long-wave rad.
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_strd    !< [-] Night factor long-wave rad.

  ! -------------------------------------------------------------------
  ! AUXILIARY VARIABLES
  ! -------------------------------------------------------------------
  !

  real(dp), public, dimension(:), allocatable :: neutron_integral_AFast !< pre-calculated integrand for
  ! vertical projection of isotropic neutron flux

END MODULE mo_global_variables
