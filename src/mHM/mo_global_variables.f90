!> \file mo_global_variables.f90

!> \brief Global variables ONLY used in reading, writing and startup.

!> \details

!> \authors Luis Samaniego
!> \date Dec 2012

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
  !           Matthias Zink,      Nov 2013 - added "InflowGauge" and inflow gauge variabels in basin
  !           Rohini Kumar,       May 2014 - added options for the model run cordinate system
  !           Stephan Thober,     Jun 2014 - added timeStep_model_inputs and readPer
  !           Stephan Thober,     Jun 2014 - added perform_mpr, updated restart flags
  !           Cuntz M. & Mai J.,  Nov 2014 - LAI input from daily, monthly or yearly files
  !           Matthias Zink,      Dec 2014 - adopted inflow gauges to ignore headwater cells
  !           Matthias Zink,      Mar 2015 - added optional soil mositure readin: dirSoil_moisture, L1_sm
  !           Stephan Thober,     Aug 2015 - moved routing related variables to mRM
  !           Oldrich Rakovec,    Oct 2015 - added definition of basin averaged TWS data
  !           Rohini Kumar,       Mar 2016 - new variables for handling different soil databases
  !           Johann Brenner,     Feb 2017 - added optional evapotranspiration readin: dirEvapotranspiration, L1_et
  !           Zink M. Demirel C., Mar 2017 - added Jarvis soil water stress variable for SM process(3) 
  !           Demirel M.C.        May 2017 - added L1_petLAIcorFactor for PET correction
  !           O. Rakovec, R.Kumar Nov 2017 - added project description for the netcdf outputs
  !           Robert Schweppe,    Dec 2017 - expanded dimensions of effective parameters
  !           Robert Schweppe,    Dec 2017 - merged duplicated variables with mrm into common variables

  USE mo_kind, ONLY : i4, i8, dp
  USE mo_mhm_constants, ONLY : nOutFlxState
  USE mo_common_constants, ONLY : YearMonths, maxNoBasins, maxNLCovers
  use mo_common_variables, only : Grid

  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! DEFINE OUTPUTS
  ! -------------------------------------------------------------------
  integer(i4) :: timeStep_model_outputs ! timestep for writing model outputs
  logical, dimension(nOutFlxState) :: outputFlxState         ! Define model outputs see "mhm_outputs.nml"
  !                                                            dim1 = number of output variables to be written
  ! -------------------------------------------------------------------
  ! INPUT variables for configuration of mHM
  ! -------------------------------------------------------------------
  integer(i4), dimension(:), allocatable, public :: timeStep_model_inputs      ! frequency for reading meteo input
  logical, public :: read_meteo_weights         ! read weights for tavg and pet
  character(256), public :: inputFormat_meteo_forcings ! format of meteo input data (nc)
  ! Optional data
  integer(i4), public :: timeStep_sm_input          ! time step of optional data: soil moisture sm
  integer(i4), public :: timeStep_neutrons_input    ! time step of optional data: neutron counts
  integer(i4), public :: timeStep_et_input          ! time step of optional data: evapotransp. et
  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nBasins
  character(256), dimension(:), allocatable, public :: dirPrecipitation   ! Directory where precipitation files are located
  character(256), dimension(:), allocatable, public :: dirTemperature     ! Directory where temperature files are located
  character(256), dimension(:), allocatable, public :: dirMinTemperature  ! Directory where minimum temp. files are located
  character(256), dimension(:), allocatable, public :: dirMaxTemperature  ! Directory where maximum temp. files are located
  character(256), dimension(:), allocatable, public :: dirNetRadiation    ! Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirabsVapPressure  ! Directory where abs. vap. pressure files are located
  character(256), dimension(:), allocatable, public :: dirwindspeed       ! Directory where windspeed files are located
  character(256), dimension(:), allocatable, public :: dirReferenceET     ! Directory where reference-ET files are located
  character(256), dimension(:), allocatable, public :: dirSoil_moisture        ! File of monthly soil moisture
  character(256), dimension(:), allocatable, public :: fileTWS                 ! File of tws data
  character(256), dimension(:), allocatable, public :: dirNeutrons             ! File of spatio-temporal neutron data
  character(256), dimension(:), allocatable, public :: dirEvapotranspiration   ! File of monthly soil moisture

  ! ------------------------------------------------------------------
  ! CONSTANT
  ! ------------------------------------------------------------------
  integer(i4), public, parameter :: routingStates = 2  ! [-]   Routing states (2=current, 1=past)

  ! ------------------------------------------------------------------
  ! BASIN AVERAGED TOTAL WATER STORAGE DATA
  ! ------------------------------------------------------------------
  type TWSstructure
    integer(i4), dimension(:), allocatable :: basinId            ! Basin Id
    character(256), dimension(:), allocatable :: fname              ! file name
    real(dp), dimension(:, :), allocatable :: TWS                ! [mm]
  end type TWSstructure
  type(TWSstructure), public :: basin_avg_TWS_obs   ! [mm] basin average TWS observational data

  real(dp), public, dimension(:, :), allocatable :: basin_avg_TWS_sim  ! variable containing basin average TWS for each basin
  integer(i4), public :: nMeasPerDay_TWS    ! Number of WTS observations per day,
  !                                                                      ! e.g. 24 -> hourly, 1 -> daily

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type(Grid), dimension(:), allocatable, public :: level2       ! Reference of the metereological variables

  ! -------------------------------------------------------------------
  ! L1 DOMAIN description
  ! -------------------------------------------------------------------
  ! Forcings
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  real(dp), public, dimension(:, :, :), allocatable :: L1_temp_weights  ! hourly temperature weights for daily values
  real(dp), public, dimension(:, :, :), allocatable :: L1_pet_weights   ! hourly pet weights for daily values
  real(dp), public, dimension(:, :, :), allocatable :: L1_pre_weights   ! hourly pre weights for daily values
  real(dp), public, dimension(:, :), allocatable :: L1_pre           ! [mm]    Precipitation
  real(dp), public, dimension(:, :), allocatable :: L1_temp          ! [degC]  Air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_pet           ! [mm TST-1] Potential evapotranspiration
  real(dp), public, dimension(:, :), allocatable :: L1_tmin          ! [degC]  minimum daily air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_tmax          ! [degC]  maximum daily air temperature
  real(dp), public, dimension(:, :), allocatable :: L1_netrad        ! [W m2]  net radiation
  real(dp), public, dimension(:, :), allocatable :: L1_absvappress   ! [Pa]    absolute vapour pressure
  real(dp), public, dimension(:, :), allocatable :: L1_windspeed     ! [m s-1] windspeed

  ! optional data
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  ! soil moisture
  real(dp), public, dimension(:, :), allocatable :: L1_sm                  ! [-] soil moisture input for optimization
  logical, public, dimension(:, :), allocatable :: L1_sm_mask             ! [-] mask for valid data in L1_sm
  integer(i4) :: nTimeSteps_L1_sm       ! [-] number of time steps in L1_sm_mask
  integer(i4) :: nSoilHorizons_sm_input ! No. of mhm soil horizons equivalent to sm input
  ! neutrons
  real(dp), public, dimension(:, :), allocatable :: L1_neutronsdata            ! [cph] ground albedo neutrons input
  logical, public, dimension(:, :), allocatable :: L1_neutronsdata_mask       ! [cph] mask for valid data in L1_neutrons
  integer(i4) :: nTimeSteps_L1_neutrons     ! [-] number of time steps in L1_neutrons_mask
  ! evapotranspiration
  real(dp), public, dimension(:, :), allocatable :: L1_et                 ! [mm] Evapotranspiration input for optimization
  logical, public, dimension(:, :), allocatable :: L1_et_mask            ! [mm] mask for valid data in L1_neutrons
  integer(i4) :: nTimeSteps_L1_et      ! [-] number of time steps in L1_sm_mask

  ! State variables
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable :: L1_inter        ! [mm]  Interception
  real(dp), public, dimension(:), allocatable :: L1_snowPack     ! [mm]  Snowpack
  real(dp), public, dimension(:), allocatable :: L1_sealSTW      ! [mm]  Retention storage of impervious areas
  real(dp), public, dimension(:, :), allocatable :: L1_soilMoist    ! [mm]  Soil moisture of each horizon
  real(dp), public, dimension(:), allocatable :: L1_unsatSTW     ! [mm]  upper soil storage
  real(dp), public, dimension(:), allocatable :: L1_satSTW       ! [mm]  groundwater storage
  real(dp), public, dimension(:), allocatable :: L1_neutrons     ! [mm]  Ground Albedo Neutrons

  ! Fluxes
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons
  real(dp), public, dimension(:), allocatable :: L1_pet_calc     ! [mm TS-1] estimated/corrected potential evapotranspiration
  real(dp), public, dimension(:, :), allocatable :: L1_aETSoil      ! [mm TS-1] Actual ET from soil layers
  real(dp), public, dimension(:), allocatable :: L1_aETCanopy    ! [mm TS-1] Real evaporation intensity from canopy
  real(dp), public, dimension(:), allocatable :: L1_aETSealed    ! [mm TS-1] Real evap. from free water surfaces
  real(dp), public, dimension(:), allocatable :: L1_baseflow     ! [mm TS-1] Baseflow
  real(dp), public, dimension(:, :), allocatable :: L1_infilSoil    ! [mm TS-1] Infiltration intensity each soil horizon
  real(dp), public, dimension(:), allocatable :: L1_fastRunoff   ! [mm TS-1] Fast runoff component
  real(dp), public, dimension(:), allocatable :: L1_melt         ! [mm TS-1] Melting snow depth.
  real(dp), public, dimension(:), allocatable :: L1_percol       ! [mm TS-1] Percolation.
  real(dp), public, dimension(:), allocatable :: L1_preEffect    ! [mm TS-1] Effective precip. depth (snow melt + rain)
  real(dp), public, dimension(:), allocatable :: L1_rain         ! [mm TS-1] Rain precipitation depth
  real(dp), public, dimension(:), allocatable :: L1_runoffSeal   ! [mm TS-1] Direct runoff from impervious areas
  real(dp), public, dimension(:), allocatable :: L1_slowRunoff   ! [mm TS-1] Slow runoff component
  real(dp), public, dimension(:), allocatable :: L1_snow         ! [mm TS-1] Snow precipitation depth
  real(dp), public, dimension(:), allocatable :: L1_Throughfall  ! [mm TS-1] Throughfall.
  real(dp), public, dimension(:), allocatable :: L1_total_runoff ! [m3 TS-1] Generated runoff

  ! -------------------------------------------------------------------
  ! Monthly day/night variation of Meteorological variables
  ! for temporal disaggregation
  ! -------------------------------------------------------------------
  ! dim1 = number of months in a year
  real(dp), public, dimension(int(YearMonths, i4)) :: evap_coeff     ! [-] Evap. coef. for free-water surfaces
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_prec      ! [-] Day ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_prec    ! [-] Night ratio precipitation < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_pet       ! [-] Day ratio PET  < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_pet     ! [-] Night ratio PET  < 1
  real(dp), public, dimension(int(YearMonths, i4)) :: fday_temp      ! [-] Day factor mean temp
  real(dp), public, dimension(int(YearMonths, i4)) :: fnight_temp    ! [-] Night factor mean temp

  ! -------------------------------------------------------------------
  ! AUXILIARY VARIABLES
  ! -------------------------------------------------------------------
  !

  real(dp), public, dimension(:), allocatable :: neutron_integral_AFast ! pre-calculated integrand for
  ! vertical projection of isotropic neutron flux

END MODULE mo_global_variables
