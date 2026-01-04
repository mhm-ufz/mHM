!> \file    mo_namelists.f90
!> \copydoc mo_namelists

!> \brief   Module containing all namelists representations.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jul 2022
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_namelists

  use mo_kind, only : i4, i8, dp
  use mo_nml, only : position_nml
  use mo_constants, only : YearMonths
  use mo_mhm_constants, only : nOutFlxState
  use mo_common_constants, only : maxNLcovers, maxNoDomains, nColPars, nodata_dp, nodata_i4
  use mo_common_variables, only : nProcesses
  use mo_common_types, only : period
  use mo_common_mHM_mRM_variables, only : nerror_model
  use mo_mpr_constants, only : maxGeoUnit, maxNoSoilHorizons
  use mo_mrm_constants, only : maxNoGauges, mrm_nOutFlxState => nOutFlxState
  use mo_string_utils, only : num2str
  use mo_sentinel, only : set_sentinel

  implicit none

  !######## mo_common_read_config

  ! namelist /project_description/ &
  !   project_details, &
  !   setup_description, &
  !   simulation_type, &
  !   Conventions, &
  !   contact, &
  !   mHM_details, &
  !   history
  !
  !> \class   nml_project_description_t
  !> \brief   'project_description' namelist content
  type, public :: nml_project_description_t
    character(19) :: name = "project_description" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    character(1024) :: project_details !< project including funding instituion., PI, etc.
    character(1024) :: setup_description !< any specific description of simulation
    character(1024) :: simulation_type !< e.g. seasonal forecast, climate projection, ...
    character(256) :: Conventions !< convention used for dataset
    character(1024) :: contact !< contact details, incl. PI name
    character(1024) :: mHM_details !< developing institution, specific mHM revision
    character(1024) :: history !< details on version/creation date
  contains
    procedure, public :: read => read_project_description
  end type nml_project_description_t
  !> 'project_description' namelist content
  type(nml_project_description_t), public :: nml_project_description

  ! namelist /directories_general/ &
  !   dirConfigOut, &
  !   dirCommonFiles, &
  !   dir_Morpho, &
  !   dir_LCover, &
  !   dir_Out, &
  !   mhm_file_RestartOut, &
  !   mrm_file_RestartOut, &
  !   file_LatLon
  !
  !> \class   nml_directories_general_t
  !> \brief   'directories_general' namelist content
  type, public :: nml_directories_general_t
    character(19) :: name = "directories_general" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    character(256) :: dirConfigOut !< directory for config file output
    character(256) :: dirCommonFiles !< directory where common input files should be located
    character(256), dimension(maxNoDomains) :: mhm_file_RestartOut !< Directory where mhm output of restart is written
    character(256), dimension(maxNoDomains) :: mrm_file_RestartOut !< Directory where mrm output of restart is written
    character(256), dimension(maxNoDomains) :: dir_Morpho !< Directory where morphological files are located
    character(256), dimension(maxNoDomains) :: dir_LCover !< Directory where land cover files are located
    character(256), dimension(maxNoDomains) :: dir_Out !< Directory where output is written to
    character(256), dimension(maxNoDomains) :: file_LatLon !< Directory where the Lat Lon Files are located
  contains
    procedure, public :: read => read_directories_general
  end type nml_directories_general_t
  !> 'directories_general' namelist content
  type(nml_directories_general_t), public :: nml_directories_general

  ! namelist /mainconfig/ &
  !   iFlag_cordinate_sys, &
  !   resolution_Hydrology, &
  !   nDomains, &
  !   L0Domain, &
  !   write_restart, &
  !   read_opt_domain_data
  !
  !> \class   nml_mainconfig_t
  !> \brief   'mainconfig' namelist content
  type, public :: nml_mainconfig_t
    character(10) :: name = "mainconfig" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: iFlag_cordinate_sys !< options model for the run cordinate system
    real(dp), dimension(maxNoDomains) :: resolution_Hydrology !< [m or degree] resolution of hydrology - Level 1
    integer(i4) :: nDomains !< number of domains
    integer(i4), dimension(maxNoDomains) :: L0Domain !< specify same index for domains to share L0_data to save memory
    logical :: write_restart !< flag to write restart
    integer(i4), dimension(maxNoDomains) :: read_opt_domain_data !< read domain specific optional data
    contains
    procedure, public :: read => read_mainconfig
  end type nml_mainconfig_t
  !> 'mainconfig' namelist content
  type(nml_mainconfig_t), public :: nml_mainconfig

  ! namelist /processSelection/ &
  !   processCase
  !
  !> \class   nml_processselection_t
  !> \brief   'processSelection' namelist content
  type, public :: nml_processselection_t
    character(16) :: name = "processselection" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4), dimension(nProcesses) :: processCase !< ! Choosen process description number
  contains
    procedure, public :: read => read_processselection
  end type nml_processselection_t
  !> 'processSelection' namelist content
  type(nml_processselection_t), public :: nml_processselection

  ! namelist /LCover/ &
  !   nLcoverScene, &
  !   LCoverYearStart, &
  !   LCoverYearEnd, &
  !   LCoverfName
  !
  !> \class   nml_lcover_t
  !> \brief   'LCover' namelist content
  type, public :: nml_lcover_t
    character(6) :: name = "lcover" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: nLCoverScene !< Number of land cover scene (lcs)
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart !< starting year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd !< ending year LCover
    character(256), dimension(maxNLCovers) :: LCoverfName !< filename of Lcover file
  contains
    procedure, public :: read => read_lcover
  end type nml_lcover_t
  !> 'LCover' namelist content
  type(nml_lcover_t), public :: nml_lcover

  !######## mo_mHM_mRM_read_config

  ! namelist /mainconfig_mhm_mrm/ &
  !   timestep, &
  !   resolution_Routing, &
  !   optimize, &
  !   optimize_restart, &
  !   opti_method, &
  !   opti_function, &
  !   read_restart, &
  !   mrm_read_river_network, &
  !   read_old_style_restart_bounds, &
  !   mhm_file_RestartIn, &
  !   mrm_file_RestartIn
  !
  !> \class   nml_mainconfig_mhm_mrm_t
  !> \brief   'mainconfig_mhm_mrm' namelist content
  type, public :: nml_mainconfig_mhm_mrm_t
    character(18) :: name = "mainconfig_mhm_mrm" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: timeStep !< [h] simulation time step (= TS) in [h] either 1, 2, 3, 4, 6, 12 or 24
    real(dp), dimension(maxNoDomains) :: resolution_Routing !< resolution of Level-11 discharge routing [m or degree] per domain
    logical :: optimize !< Optimization (.true.) or Evaluation run (.false.)
    logical :: optimize_restart !< Optimization will be restarted from mo_<opti_method>.restart file (.true.)
    integer(i4) :: opti_method !< Optimization algorithm: 1 - DDS; 2 - Simulated Annealing; 3 - SCE
    integer(i4) :: opti_function !< Objective function
    logical :: read_restart !< flag for reading restart output
    logical :: mrm_read_river_network !< flag to read the river network for mRM (read_restart = .True. forces .True.)
    logical :: read_old_style_restart_bounds !< flag to use an old-style restart file created by mhm<=v5.11
    logical :: restart_reset_fluxes_states !< flag to reset fluxes and states read from restart to default values
    character(256), dimension(maxNoDomains) :: mhm_file_RestartIn !< mhm restart file paths
    character(256), dimension(maxNoDomains) :: mrm_file_RestartIn !< mrm restart file paths
  contains
    procedure, public :: read => read_mainconfig_mhm_mrm
  end type nml_mainconfig_mhm_mrm_t
  !> 'mainconfig_mhm_mrm' namelist content
  type(nml_mainconfig_mhm_mrm_t), public :: nml_mainconfig_mhm_mrm

  ! namelist /Optimization/ &
  !   nIterations, &
  !   seed, &
  !   dds_r, &
  !   sa_temp, &
  !   sce_ngs, &
  !   sce_npg, &
  !   sce_nps, &
  !   mcmc_opti, &
  !   mcmc_error_params
  !
  !> \class   nml_optimization_t
  !> \brief   'optimization' namelist content
  type, public :: nml_optimization_t
    character(12) :: name = "optimization" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: nIterations !< number of iterations for optimization
    integer(i8) :: seed !< seed used for optimization, default: -9 --> system time
    real(dp) :: dds_r !< DDS: perturbation rate, default: 0.2
    real(dp) :: sa_temp !< SA:  initial temperature, default: -9.0 --> estimated
    integer(i4) :: sce_ngs !< SCE: # of complexes, default: 2
    integer(i4) :: sce_npg !< SCE: # of points per complex,default: -9 --> 2n+1
    integer(i4) :: sce_nps !< SCE: # of points per subcomplex,default: -9 --> n+1
    logical :: mcmc_opti !< MCMC: optimization (.true.) or only parameter uncertainty (.false.)
    real(dp), dimension(nerror_model) :: mcmc_error_params !< error model para (mcmc_opti=.false.) e.g. for opti_function=8: .01, .3
  contains
    procedure, public :: read => read_optimization
  end type nml_optimization_t
  !> 'optimization' namelist content
  type(nml_optimization_t), public :: nml_optimization

  ! namelist /time_periods/ &
  !   warming_Days, &
  !   eval_Per
  !
  !> \class   nml_time_periods_t
  !> \brief   'time_periods' namelist content
  type, public :: nml_time_periods_t
    character(12) :: name = "time_periods" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4), dimension(maxNoDomains) :: warming_Days !< number of days for warm up period
    type(period), dimension(maxNoDomains) :: eval_Per !< time period for model evaluation
  contains
    procedure, public :: read => read_time_periods
  end type nml_time_periods_t
  !> 'time_periods' namelist content
  type(nml_time_periods_t), public :: nml_time_periods

  !######## mo_mhm_read_config

  ! namelist /directories_mhm/ &
  !   inputFormat_meteo_forcings, &
  !   dir_Precipitation, &
  !   dir_Temperature, &
  !   dir_ReferenceET, &
  !   dir_MinTemperature, &
  !   dir_MaxTemperature, &
  !   dir_absVapPressure, &
  !   dir_windspeed, &
  !   dir_NetRadiation, &
  !   dir_Radiation, &
  !   time_step_model_inputs
  !
  !> \class   nml_directories_mhm_t
  !> \brief   'directories_mhm' namelist content
  type, public :: nml_directories_mhm_t
    character(15) :: name = "directories_mhm" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    !> .FALSE. to only warn about bound (lower, upper) violations in meteo files, default = .TRUE. - raise an error
    logical :: bound_error = .true.
    character(256), public :: inputFormat_meteo_forcings !< format of meteo input data (nc)
    character(256), dimension(maxNoDomains) :: dir_meteo_header !< Directory where the meteo header file is located
    character(256), dimension(maxNoDomains) :: dir_Precipitation !< Directory where precipitation files are located
    character(256), dimension(maxNoDomains) :: dir_Temperature !< Directory where temperature files are located
    character(256), dimension(maxNoDomains) :: dir_ReferenceET !< Directory where reference-ET files are located
    character(256), dimension(maxNoDomains) :: dir_MinTemperature !< Directory where minimum temp. files are located
    character(256), dimension(maxNoDomains) :: dir_MaxTemperature !< Directory where maximum temp. files are located
    character(256), dimension(maxNoDomains) :: dir_absVapPressure !< Directory where abs. vap. pressure files are located
    character(256), dimension(maxNoDomains) :: dir_windspeed !< Directory where windspeed files are located
    character(256), dimension(maxNoDomains) :: dir_NetRadiation !< Directory where abs. vap. pressure files are located
    character(256), dimension(maxNoDomains) :: dir_Radiation !< riv-temp related: directory of (long/short-wave)radiation
    integer(i4), dimension(maxNoDomains) :: time_step_model_inputs !< frequency for reading meteo input
  contains
    procedure, public :: read => read_directories_mhm
  end type nml_directories_mhm_t
  !> 'directories_mhm' namelist content
  type(nml_directories_mhm_t), public :: nml_directories_mhm

  ! namelist /optional_data/ &
  !   nSoilHorizons_sm_input, &
  !   dir_soil_moisture, &
  !   dir_neutrons, &
  !   dir_evapotranspiration, &
  !   dir_TWS, &
  !   timeStep_sm_input, &
  !   timeStep_neutrons_input, &
  !   timeStep_et_input, &
  !   timeStep_tws_input
  !
  !> \class   nml_optional_data_t
  !> \brief   'optional_data' namelist content
  type, public :: nml_optional_data_t
    character(13) :: name = "optional_data" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: nSoilHorizons_sm_input !< No. of mhm soil horizons equivalent to sm input
    character(256), dimension(maxNoDomains) :: dir_soil_moisture !< soil moisture input
    character(256), dimension(maxNoDomains) :: dir_neutrons !< ground albedo neutron input
    character(256), dimension(maxNoDomains) :: dir_evapotranspiration !< evapotranspiration input
    character(256), dimension(maxNoDomains) :: dir_TWS !< tws input
    character(256), dimension(maxNoDomains) :: dir_spf !< spf input
    integer(i4) :: timeStep_sm_input !< time step of optional data: sm
    integer(i4) :: timeStep_neutrons_input !< time step of optional data: neutrons
    integer(i4) :: timeStep_et_input !< time step of optional data: et
    integer(i4) :: timeStep_tws_input !< time step of optional data: tws
    integer(i4) :: timeStep_spf_input !< time step of optional data: spf
    integer(i4) :: weight_for_optional_data   ! weight of optional data in OF 67 to 69
    integer(i4) :: snow_water_equivalent_threshold_for_spf ! threshold swe to convert to spf for OF 75, 76


  contains
    procedure, public :: read => read_optional_data
  end type nml_optional_data_t
  !> 'optional_data' namelist content
  type(nml_optional_data_t), public :: nml_optional_data

  ! namelist /panevapo/ &
  !   evap_coeff
  !
  !> \class   nml_panevapo_t
  !> \brief   'panevapo' namelist content
  type, public :: nml_panevapo_t
    character(8) :: name = "panevapo" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(int(YearMonths, i4)) :: evap_coeff !< [-] Evap. coef. for free-water surfaces
  contains
    procedure, public :: read => read_panevapo
  end type nml_panevapo_t
  !> 'panevapo' namelist content
  type(nml_panevapo_t), public :: nml_panevapo

  ! namelist /nightdayratio/ &
  !   read_meteo_weights, &
  !   fnight_prec, &
  !   fnight_pet, &
  !   fnight_temp, &
  !   fnight_ssrd, &
  !   fnight_strd
  !
  !> \class   nml_nightdayratio_t
  !> \brief   'nightdayratio' namelist content
  type, public :: nml_nightdayratio_t
    character(13) :: name = "nightdayratio" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    logical :: read_meteo_weights !< read weights for meteo data
    real(dp), dimension(int(YearMonths, i4)) :: fnight_prec !< [-] Night ratio precipitation < 1
    real(dp), dimension(int(YearMonths, i4)) :: fnight_pet !< [-] Night ratio PET  < 1
    real(dp), dimension(int(YearMonths, i4)) :: fnight_temp !< [-] Night factor mean temp
    real(dp), dimension(int(YearMonths, i4)) :: fnight_ssrd !< [-] Night factor short-wave rad.
    real(dp), dimension(int(YearMonths, i4)) :: fnight_strd !< [-] Night factor long-wave rad.
  contains
    procedure, public :: read => read_nightdayratio
  end type nml_nightdayratio_t
  !> 'nightdayratio' namelist content
  type(nml_nightdayratio_t), public :: nml_nightdayratio

  ! namelist /nloutputresults/ &
  !   output_deflate_level, &
  !   output_double_precision, &
  !   timeStep_model_outputs, &
  !   outputFlxState
  !
  !> \class   nml_nloutputresults_t
  !> \brief   'nloutputresults' namelist content
  type, public :: nml_nloutputresults_t
    character(15) :: name = "nloutputresults" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: output_deflate_level !< deflate level in nc files
    logical :: output_double_precision !< output precision in nc files
    integer(i4) :: timeStep_model_outputs !< timestep for writing model outputs
    integer(i4) :: output_time_reference !< time reference point location in output nc files
    logical, dimension(nOutFlxState) :: outputFlxState !< Define model outputs see "mhm_outputs.nml"
  contains
    procedure, public :: read => read_nloutputresults
  end type nml_nloutputresults_t
  !> 'nloutputresults' namelist content
  type(nml_nloutputresults_t), public :: nml_nloutputresults

  ! namelist /baseflow_config/ &
  !   BFI_calc, &
  !   BFI_obs
  !
  !> \class   nml_baseflow_config_t
  !> \brief   'baseflow_config' namelist content
  type, public :: nml_baseflow_config_t
    character(15) :: name = "baseflow_config" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    logical :: BFI_calc !< calculate observed BFI from gauges with Eckhardt filter
    real(dp), dimension(maxNoDomains) :: BFI_obs !< given base-flow index per domain
contains
    procedure, public :: read => read_baseflow_config
  end type nml_baseflow_config_t
  !> 'baseflow_config' namelist content
  type(nml_baseflow_config_t), public :: nml_baseflow_config

  !######## mo_mpr_read_config
  ! namelist /directories_MPR/ &
  !   dir_gridded_LAI
  !
  !> \class   nml_directories_mpr_t
  !> \brief   'directories_mpr' namelist content
  type, public :: nml_directories_mpr_t
    character(15) :: name = "directories_mpr" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    character(256), dimension(maxNoDomains) :: dir_gridded_LAI !< directory of gridded LAI data, used when timeStep_LAI_input<0
  contains
    procedure, public :: read => read_directories_mpr
  end type nml_directories_mpr_t
  !> 'directories_mpr' namelist content
  type(nml_directories_mpr_t), public :: nml_directories_mpr

  ! namelist /soildata/ &
  !   iFlag_soilDB, &
  !   tillageDepth, &
  !   nSoilHorizons_mHM, &
  !   soil_Depth
  !
  !> \class   nml_soildata_t
  !> \brief   'soildata' namelist content
  type, public :: nml_soildata_t
    character(8) :: name = "soildata" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: iFlag_soilDB !< options to handle different soil databases
    real(dp) :: tillageDepth !< [mm] Soil depth down to which organic
    integer(i4) :: nSoilHorizons_mHM !< Number of horizons to model
    real(dp), dimension(maxNoSoilHorizons) :: soil_Depth !< depth of the single horizons
  contains
    procedure, public :: read => read_soildata
  end type nml_soildata_t
  !> 'soildata' namelist content
  type(nml_soildata_t), public :: nml_soildata

  ! namelist /LAI_data_information/ &
  !   inputFormat_gridded_LAI, &
  !   timeStep_LAI_input
  !
  !> \class   nml_lai_data_information_t
  !> \brief   'lai_data_information' namelist content
  type, public :: nml_lai_data_information_t
    character(20) :: name = "lai_data_information" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    character(256) :: inputFormat_gridded_LAI !< format of gridded LAI data (nc only)
    integer(i4) :: timeStep_LAI_input !< time step of gridded LAI input
  contains
    procedure, public :: read => read_lai_data_information
  end type nml_lai_data_information_t
  !> 'lai_data_information' namelist content
  type(nml_lai_data_information_t), public :: nml_lai_data_information

  ! namelist /LCover_MPR/ &
  !   fracSealed_cityArea
  !
  !> \class   nml_lcover_mpr_t
  !> \brief   'lcover_mpr' namelist content
  type, public :: nml_lcover_mpr_t
    character(10) :: name = "lcover_mpr" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp) :: fracSealed_cityArea !< fraction of area within city assumed to be perfectly sealed [0-1]
  contains
    procedure, public :: read => read_lcover_mpr
  end type nml_lcover_mpr_t
  !> 'lcover_mpr' namelist content
  type(nml_lcover_mpr_t), public :: nml_lcover_mpr

  ! namelist /interception1/ &
  !   canopyInterceptionFactor
  !
  !> \class   nml_interception1_t
  !> \brief   'interception1' namelist content
  type, public :: nml_interception1_t
    character(13) :: name = "interception1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: canopyInterceptionFactor !< multiplier to relate LAI to interception storage [-]
  contains
    procedure, public :: read => read_interception1
  end type nml_interception1_t
  !> 'interception1' namelist content
  type(nml_interception1_t), public :: nml_interception1

  ! namelist /snow1/ &
  !   snowTreshholdTemperature, &
  !   degreeDayFactor_forest, &
  !   degreeDayFactor_impervious, &
  !   degreeDayFactor_pervious, &
  !   increaseDegreeDayFactorByPrecip, &
  !   maxDegreeDayFactor_forest, &
  !   maxDegreeDayFactor_impervious, &
  !   maxDegreeDayFactor_pervious
  !
  !> \class   nml_snow1_t
  !> \brief   'snow1' namelist content
  type, public :: nml_snow1_t
    character(5) :: name = "snow1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: snowTreshholdTemperature !< Threshold for rain/snow partitioning [degC]
    real(dp), dimension(nColPars) :: degreeDayFactor_forest !< forest: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: degreeDayFactor_impervious !< impervious: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: degreeDayFactor_pervious !< pervious: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: increaseDegreeDayFactorByPrecip !< increase of deg day factor in case of precipitation [degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_forest !< forest: maximum values for degree day factor [m degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_impervious !< impervious: maximum values for degree day factor [m degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_pervious !< pervious: maximum values for degree day factor [m degC-1]
  contains
    procedure, public :: read => read_snow1
  end type nml_snow1_t
  !> 'snow1' namelist content
  type(nml_snow1_t), public :: nml_snow1

  ! namelist /soilmoisture1/ &
  !   orgMatterContent_forest, &
  !   orgMatterContent_impervious, &
  !   orgMatterContent_pervious, &
  !   PTF_lower66_5_constant, &
  !   PTF_lower66_5_clay, &
  !   PTF_lower66_5_Db, &
  !   PTF_higher66_5_constant, &
  !   PTF_higher66_5_clay, &
  !   PTF_higher66_5_Db, &
  !   PTF_Ks_constant, &
  !   PTF_Ks_sand, &
  !   PTF_Ks_clay, &
  !   PTF_Ks_curveSlope, &
  !   rootFractionCoefficient_forest, &
  !   rootFractionCoefficient_impervious, &
  !   rootFractionCoefficient_pervious, &
  !   infiltrationShapeFactor
  !
  !> \class   nml_soilmoisture1_t
  !> \brief   'soilmoisture1' namelist content
  type, public :: nml_soilmoisture1_t
    character(13) :: name = "soilmoisture1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
  contains
    procedure, public :: read => read_soilmoisture1
  end type nml_soilmoisture1_t
  !> 'soilmoisture1' namelist content
  type(nml_soilmoisture1_t), public :: nml_soilmoisture1

  ! namelist /soilmoisture2/ &
  !   orgMatterContent_forest, &
  !   orgMatterContent_impervious, &
  !   orgMatterContent_pervious, &
  !   PTF_lower66_5_constant, &
  !   PTF_lower66_5_clay, &
  !   PTF_lower66_5_Db, &
  !   PTF_higher66_5_constant, &
  !   PTF_higher66_5_clay, &
  !   PTF_higher66_5_Db, &
  !   PTF_Ks_constant, &
  !   PTF_Ks_sand, &
  !   PTF_Ks_clay, &
  !   PTF_Ks_curveSlope, &
  !   rootFractionCoefficient_forest, &
  !   rootFractionCoefficient_impervious, &
  !   rootFractionCoefficient_pervious, &
  !   infiltrationShapeFactor, &
  !   jarvis_sm_threshold_c1
  !
  !> \class   nml_soilmoisture2_t
  !> \brief   'soilmoisture2' namelist content
  type, public :: nml_soilmoisture2_t
    character(13) :: name = "soilmoisture2" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1 !< soil moisture threshod for jarvis model
  contains
    procedure, public :: read => read_soilmoisture2
  end type nml_soilmoisture2_t
  !> 'soilmoisture2' namelist content
  type(nml_soilmoisture2_t), public :: nml_soilmoisture2

  ! namelist /soilmoisture3/ &
  !   orgMatterContent_forest, &
  !   orgMatterContent_impervious, &
  !   orgMatterContent_pervious, &
  !   PTF_lower66_5_constant, &
  !   PTF_lower66_5_clay, &
  !   PTF_lower66_5_Db, &
  !   PTF_higher66_5_constant, &
  !   PTF_higher66_5_clay, &
  !   PTF_higher66_5_Db, &
  !   PTF_Ks_constant, &
  !   PTF_Ks_sand, &
  !   PTF_Ks_clay, &
  !   PTF_Ks_curveSlope, &
  !   rootFractionCoefficient_forest, &
  !   rootFractionCoefficient_impervious, &
  !   rootFractionCoefficient_pervious, &
  !   infiltrationShapeFactor, &
  !   rootFractionCoefficient_sand, &
  !   rootFractionCoefficient_clay, &
  !   FCmin_glob, &
  !   FCdelta_glob, &
  !   jarvis_sm_threshold_c1
  !
  !> \class   nml_soilmoisture3_t
  !> \brief   'soilmoisture3' namelist content
  type, public :: nml_soilmoisture3_t
    character(13) :: name = "soilmoisture3" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: FCmin_glob !< global field capacity minimum
    real(dp), dimension(nColPars) :: FCdelta_glob !< difference between global field capacity minimum and maximum
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand !< threshold for actual ET reduction for sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay !< threshold for actual ET reduction for clay
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1 !< soil moisture threshod for jarvis model
  contains
    procedure, public :: read => read_soilmoisture3
  end type nml_soilmoisture3_t
  !> 'soilmoisture3' namelist content
  type(nml_soilmoisture3_t), public :: nml_soilmoisture3

  ! namelist /soilmoisture4/ &
  !   orgMatterContent_forest, &
  !   orgMatterContent_impervious, &
  !   orgMatterContent_pervious, &
  !   PTF_lower66_5_constant, &
  !   PTF_lower66_5_clay, &
  !   PTF_lower66_5_Db, &
  !   PTF_higher66_5_constant, &
  !   PTF_higher66_5_clay, &
  !   PTF_higher66_5_Db, &
  !   PTF_Ks_constant, &
  !   PTF_Ks_sand, &
  !   PTF_Ks_clay, &
  !   PTF_Ks_curveSlope, &
  !   rootFractionCoefficient_forest, &
  !   rootFractionCoefficient_impervious, &
  !   rootFractionCoefficient_pervious, &
  !   infiltrationShapeFactor, &
  !   rootFractionCoefficient_sand, &
  !   rootFractionCoefficient_clay, &
  !   FCmin_glob, &
  !   FCdelta_glob, &
  !
  !> \class   nml_soilmoisture4_t
  !> \brief   'soilmoisture4' namelist content
  type, public :: nml_soilmoisture4_t
    character(13) :: name = "soilmoisture4" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: FCmin_glob !< global field capacity minimum
    real(dp), dimension(nColPars) :: FCdelta_glob !< difference between global field capacity minimum and maximum
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand !< threshold for actual ET reduction for sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay !< threshold for actual ET reduction for clay
  contains
    procedure, public :: read => read_soilmoisture4
  end type nml_soilmoisture4_t
  !> 'soilmoisture4' namelist content
  type(nml_soilmoisture4_t), public :: nml_soilmoisture4

  ! namelist /directRunoff1/ &
  !   imperviousStorageCapacity
  !
  !> \class   nml_directrunoff1_t
  !> \brief   'directrunoff1' namelist content
  type, public :: nml_directrunoff1_t
    character(13) :: name = "directrunoff1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: imperviousStorageCapacity !< direct Runoff: Sealed Area storage capacity
  contains
    procedure, public :: read => read_directrunoff1
  end type nml_directrunoff1_t
  !> 'directrunoff1' namelist content
  type(nml_directrunoff1_t), public :: nml_directrunoff1

  ! namelist /PETminus1/  &
  !   PET_a_forest, &
  !   PET_a_impervious, &
  !   PET_a_pervious, &
  !   PET_b, &
  !   PET_c
  !
  !> \class   nml_petminus1_t
  !> \brief   'petminus1' namelist content
  !> \details PET is input, LAI driven correction
  type, public :: nml_petminus1_t
    character(9) :: name = "petminus1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: PET_a_forest !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_a_impervious !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_a_pervious !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_b !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_c !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
  contains
    procedure, public :: read => read_petminus1
  end type nml_petminus1_t
  !> 'petminus1' namelist content
  type(nml_petminus1_t), public :: nml_petminus1

  ! namelist /PET0/ &
  !   minCorrectionFactorPET, &
  !   maxCorrectionFactorPET, &
  !   aspectTresholdPET
  !
  !> \class   nml_pet0_t
  !> \brief   'pet0' namelist content
  !> \details PET is input, aspect driven correction
  type, public :: nml_pet0_t
    character(4) :: name = "pet0" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: minCorrectionFactorPET !< minimum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: maxCorrectionFactorPET !< maximum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: aspectTresholdPET !< aspect threshold for PET correction with aspect
  contains
    procedure, public :: read => read_pet0
  end type nml_pet0_t
  !> 'pet0' namelist content
  type(nml_pet0_t), public :: nml_pet0

  ! namelist /PET1/ &
  !   minCorrectionFactorPET, &
  !   maxCorrectionFactorPET, &
  !   aspectTresholdPET, &
  !   HargreavesSamaniCoeff
  !
  !> \class   nml_pet1_t
  !> \brief   'pet1' namelist content
  !> \details PET - Hargreaves Samani
  type, public :: nml_pet1_t
    character(4) :: name = "pet1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: minCorrectionFactorPET !< minimum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: maxCorrectionFactorPET !< maximum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: aspectTresholdPET !< aspect threshold for PET correction with aspect
    real(dp), dimension(nColPars) :: HargreavesSamaniCoeff !< coefficient for Hargreaves Samani
  contains
    procedure, public :: read => read_pet1
  end type nml_pet1_t
  !> 'pet1' namelist content
  type(nml_pet1_t), public :: nml_pet1

  ! namelist /PET2/ &
  !   PriestleyTaylorCoeff, &
  !   PriestleyTaylorLAIcorr
  !
  !> \class   nml_pet2_t
  !> \brief   'pet2' namelist content
  !> \details PET - Priestley Taylor
  type, public :: nml_pet2_t
    character(4) :: name = "pet2" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: PriestleyTaylorCoeff !< Priestley-Taylor coefficient
    real(dp), dimension(nColPars) :: PriestleyTaylorLAIcorr !< Priestley-Taylor LAI correction factor
  contains
    procedure, public :: read => read_pet2
  end type nml_pet2_t
  !> 'pet2' namelist content
  type(nml_pet2_t), public :: nml_pet2

  ! namelist /PET3/ &
  !   canopyheigth_forest, &
  !   canopyheigth_impervious, &
  !   canopyheigth_pervious, &
  !   displacementheight_coeff, &
  !   roughnesslength_momentum_coeff, &
  !   roughnesslength_heat_coeff, &
  !   stomatal_resistance
  !
  !> \class   nml_pet3_t
  !> \brief   'pet3' namelist content
  !> \details PET - Penman Monteith
  type, public :: nml_pet3_t
    character(4) :: name = "pet3" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: canopyheigth_forest !< canopy height for foreset
    real(dp), dimension(nColPars) :: canopyheigth_impervious !< canopy height for impervious
    real(dp), dimension(nColPars) :: canopyheigth_pervious !< canopy height for pervious
    real(dp), dimension(nColPars) :: displacementheight_coeff !< displacement height coefficient
    real(dp), dimension(nColPars) :: roughnesslength_momentum_coeff !< roughness length momentum coefficient
    real(dp), dimension(nColPars) :: roughnesslength_heat_coeff !< roughness length heat coefficient
    real(dp), dimension(nColPars) :: stomatal_resistance !< stomatal resistance
  contains
    procedure, public :: read => read_pet3
  end type nml_pet3_t
  !> 'pet3' namelist content
  type(nml_pet3_t), public :: nml_pet3

  ! namelist /interflow1/ &
  !   interflowStorageCapacityFactor, &
  !   interflowRecession_slope, &
  !   fastInterflowRecession_forest, &
  !   slowInterflowRecession_Ks, &
  !   exponentSlowInterflow
  !
  !> \class   nml_interflow1_t
  !> \brief   'interflow1' namelist content
  type, public :: nml_interflow1_t
    character(10) :: name = "interflow1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: interflowStorageCapacityFactor !< interflow storage capacity factor
    real(dp), dimension(nColPars) :: interflowRecession_slope !< multiplier for slope to derive interflow recession constant
    !> multiplier to derive fast interflow recession constant for forest
    real(dp), dimension(nColPars) :: fastInterflowRecession_forest
    !> multiplier for variability of saturated hydraulic conductivity to derive slow interflow recession constant
    real(dp), dimension(nColPars) :: slowInterflowRecession_Ks
    !> multiplier for variability of saturated hydraulic conductivity to derive slow interflow exponent
    real(dp), dimension(nColPars) :: exponentSlowInterflow
  contains
    procedure, public :: read => read_interflow1
  end type nml_interflow1_t
  !> 'interflow1' namelist content
  type(nml_interflow1_t), public :: nml_interflow1

  ! namelist /percolation1/ &
  !   rechargeCoefficient, &
  !   rechargeFactor_karstic, &
  !   gain_loss_GWreservoir_karstic
  !
  !> \class   nml_percolation1_t
  !> \brief   'percolation1' namelist content
  type, public :: nml_percolation1_t
    character(12) :: name = "percolation1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: rechargeCoefficient !< recharge coefficient
    real(dp), dimension(nColPars) :: rechargeFactor_karstic !< recharge factor for karstic percolation
    real(dp), dimension(nColPars) :: gain_loss_GWreservoir_karstic !< gain loss in ground water reservoir for karstic
  contains
    procedure, public :: read => read_percolation1
  end type nml_percolation1_t
  !> 'percolation1' namelist content
  type(nml_percolation1_t), public :: nml_percolation1

  ! namelist /neutrons1/ &
  !   Desilets_N0, &
  !   Desilets_LW0, &
  !   Desilets_LW1
  !
  !> \class   nml_neutrons1_t
  !> \brief   'neutrons1' namelist content
  type, public :: nml_neutrons1_t
    character(9) :: name = "neutrons1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: Desilets_N0 !< Desilets N0 parameter
    real(dp), dimension(nColPars) :: Desilets_LW0 !< Desilets LW0 parameter
    real(dp), dimension(nColPars) :: Desilets_LW1 !< Desilets LW1 parameter
  contains
    procedure, public :: read => read_neutrons1
  end type nml_neutrons1_t
  !> 'neutrons1' namelist content
  type(nml_neutrons1_t), public :: nml_neutrons1

  ! namelist /neutrons2/ &
  !   COSMIC_N0, &
  !   COSMIC_N1, &
  !   COSMIC_N2, &
  !   COSMIC_alpha0, &
  !   COSMIC_alpha1, &
  !   COSMIC_L30, &
  !   COSMIC_L31, &
  !   COSMIC_LW0, &
  !   COSMIC_LW1
  !
  !> \class   nml_neutrons2_t
  !> \brief   'neutrons2' namelist content
  type, public :: nml_neutrons2_t
    character(9) :: name = "neutrons2" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: COSMIC_N0 !< cosmic N0 parameter
    real(dp), dimension(nColPars) :: COSMIC_N1 !< cosmic N1 parameter
    real(dp), dimension(nColPars) :: COSMIC_N2 !< cosmic N2 parameter
    real(dp), dimension(nColPars) :: COSMIC_alpha0 !< cosmic alpha0 parameter
    real(dp), dimension(nColPars) :: COSMIC_alpha1 !< cosmic alpha1 parameter
    real(dp), dimension(nColPars) :: COSMIC_L30 !< cosmic L30 parameter
    real(dp), dimension(nColPars) :: COSMIC_L31 !< cosmic L31 parameter
    real(dp), dimension(nColPars) :: COSMIC_LW0 !< cosmic LW0 parameter
    real(dp), dimension(nColPars) :: COSMIC_LW1 !< cosmic LW1 parameter
  contains
    procedure, public :: read => read_neutrons2
  end type nml_neutrons2_t
  !> 'neutrons2' namelist content
  type(nml_neutrons2_t), public :: nml_neutrons2

  ! namelist /geoparameter/ &
  !   GeoParam
  !
  !> \class   nml_geoparameter_t
  !> \brief   'geoparameter' namelist content
  type, public :: nml_geoparameter_t
    character(12) :: name = "geoparameter" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    !> geological parameters (ordering according to file 'geology_classdefinition.txt')
    real(dp), dimension(maxGeoUnit, nColPars) :: GeoParam
  contains
    procedure, public :: read => read_geoparameter
  end type nml_geoparameter_t
  !> 'geoparameter' namelist content
  type(nml_geoparameter_t), public :: nml_geoparameter

  !######## mo_mrm_read_config
  ! namelist /mainconfig_mrm/ &
  !   ALMA_convention, &
  !   filenameTotalRunoff, &
  !   varnameTotalRunoff, &
  !   gw_coupling
  !
  !> \class   nml_mainconfig_mrm_t
  !> \brief   'mainconfig_mrm' namelist content
  type, public :: nml_mainconfig_mrm_t
    character(14) :: name = "mainconfig_mrm" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    logical :: ALMA_convention !< flag for ALMA convention (see http://www.lmd.jussieu.fr/~polcher/ALMA/convention_3.html)
    character(256) :: filenameTotalRunoff !< Filename of simulated total runoff file
    character(256) :: varnameTotalRunoff !< variable name of total runoff
    logical :: gw_coupling !< switch to enable ground water coupling
  contains
    procedure, public :: read => read_mainconfig_mrm
  end type nml_mainconfig_mrm_t
  !> 'mainconfig_mrm' namelist content
  type(nml_mainconfig_mrm_t), public :: nml_mainconfig_mrm

  ! namelist /directories_mRM/ &
  !   dir_Gauges, &
  !   dir_Total_Runoff, &
  !   dir_Bankfull_Runoff
  !
  !> \class   nml_directories_mrm_t
  !> \brief   'directories_mrm' namelist content
  type, public :: nml_directories_mrm_t
    character(15) :: name = "directories_mrm" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    character(256), dimension(maxNoDomains) :: dir_Gauges !< directory containing gauge time series
    character(256), dimension(maxNoDomains) :: dir_Total_Runoff !< directory where simulated runoff can be found
    character(256), dimension(maxNoDomains) :: dir_Bankfull_Runoff !< directory where runoff at bankfull conditions can be found
  contains
    procedure, public :: read => read_directories_mrm
  end type nml_directories_mrm_t
  !> 'directories_mrm' namelist content
  type(nml_directories_mrm_t), public :: nml_directories_mrm

  ! namelist /evaluation_gauges/ &
  !   nGaugesTotal, &
  !   NoGauges_domain, &
  !   Gauge_id, &
  !   gauge_filename
  !
  !> \class   nml_evaluation_gauges_t
  !> \brief   'evaluation_gauges' namelist content
  type, public :: nml_evaluation_gauges_t
    character(17) :: name = "evaluation_gauges" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: nGaugesTotal !< Number of evaluation gauges for all domains
    integer(i4), dimension(maxNoDomains) :: NoGauges_domain !< number of gauges per domain
    integer(i4), dimension(maxNoDomains, maxNoGauges) :: Gauge_id !< gauge ID for each gauge
    character(256), dimension(maxNoDomains, maxNoGauges) :: Gauge_filename !< filename for each gauge time series
  contains
    procedure, public :: read => read_evaluation_gauges
  end type nml_evaluation_gauges_t
  !> 'evaluation_gauges' namelist content
  type(nml_evaluation_gauges_t), public :: nml_evaluation_gauges

  ! namelist /inflow_gauges/ &
  !   nInflowGaugesTotal, &
  !   NoInflowGauges_domain, &
  !   InflowGauge_id, &
  !   InflowGauge_filename, &
  !   InflowGauge_Headwater
  !
  !> \class   nml_inflow_gauges_t
  !> \brief   'inflow_gauges' namelist content
  type, public :: nml_inflow_gauges_t
    character(13) :: name = "inflow_gauges" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: nInflowGaugesTotal !< Number of evaluation gauges for all domains
    integer(i4), dimension(maxNoDomains) :: NoInflowGauges_domain !< number of gauges for subdomain (1)
    integer(i4), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_id !< id of inflow gauge(1) for subdomain(1) --> (1,1)
    !> name of file with timeseries of inflow gauge(1) for subdomain(1) --> (1,1)
    character(256), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_filename
    !> consider flows from upstream/headwater cells of inflow gauge(1) for subdomain(1) --> (1,1)
    logical, dimension(maxNoDomains, maxNoGauges) :: InflowGauge_Headwater
  contains
    procedure, public :: read => read_inflow_gauges
  end type nml_inflow_gauges_t
  !> 'inflow_gauges' namelist content
  type(nml_inflow_gauges_t), public :: nml_inflow_gauges

  ! namelist /nloutputresults/ &
  !   output_deflate_level_mrm, &
  !   output_double_precision_mrm, &
  !   timeStep_model_outputs_mrm, &
  !   outputFlxState_mrm
  !
  !> \class   nml_mrm_outputs_t
  !> \brief   'mrm_outputs' namelist content
  type, public :: nml_mrm_outputs_t
    character(15) :: name = "nloutputresults" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: output_deflate_level_mrm !< netcdf deflate level
    logical :: output_double_precision_mrm !< switch to enable double precision in netcdf
    integer(i4) :: output_time_reference_mrm !< time reference point location in output nc files
    integer(i4) :: timeStep_model_outputs_mrm !< timestep for writing model outputs
    logical, dimension(mrm_nOutFlxState) :: outputFlxState_mrm !< Define model outputs see "mhm_outputs.nml"
  contains
    procedure, public :: read => read_mrm_outputs
  end type nml_mrm_outputs_t
  !> 'mrm_outputs' namelist content
  type(nml_mrm_outputs_t), public :: nml_mrm_outputs

  ! namelist /routing1/ &
  !   muskingumTravelTime_constant, &
  !   muskingumTravelTime_riverLength, &
  !   muskingumTravelTime_riverSlope, &
  !   muskingumTravelTime_impervious, &
  !   muskingumAttenuation_riverSlope
  !
  !> \class   nml_routing1_t
  !> \brief   'routing1' namelist content
  type, public :: nml_routing1_t
    character(8) :: name = "routing1" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: muskingumTravelTime_constant !< muskingum parameter constant
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength !< muskingum parameter river length
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope !< muskingum parameter river slope
    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious !< muskingum parameter impervious
    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope !< muskingum parameter attenuation river slope
  contains
    procedure, public :: read => read_routing1
  end type nml_routing1_t
  !> 'routing1' namelist content
  type(nml_routing1_t), public :: nml_routing1

  ! namelist /routing2/ &
  !   streamflow_celerity
  !
  !> \class   nml_routing2_t
  !> \brief   'routing2' namelist content
  type, public :: nml_routing2_t
    character(8) :: name = "routing2" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: streamflow_celerity !< streamflow celerity
  contains
    procedure, public :: read => read_routing2
  end type nml_routing2_t
  !> 'routing2' namelist content
  type(nml_routing2_t), public :: nml_routing2

  ! namelist /routing3/ &
  !   slope_factor
  !
  !> \class   nml_routing3_t
  !> \brief   'routing3' namelist content
  type, public :: nml_routing3_t
    character(8) :: name = "routing3" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp), dimension(nColPars) :: slope_factor !< slope factor
  contains
    procedure, public :: read => read_routing3
  end type nml_routing3_t
  !> 'routing3' namelist content
  type(nml_routing3_t), public :: nml_routing3

  !######## mo_mrm_riv_temp_class
  ! namelist /config_riv_temp/ &
  !   albedo_water, &
  !   pt_a_water, &
  !   emissivity_water, &
  !   turb_heat_ex_coeff, &
  !   max_iter, &
  !   delta_iter, &
  !   step_iter, &
  !   riv_widths_file, &
  !   riv_widths_name, &
  !   dir_riv_widths
  !
  !> \class   nml_config_riv_temp_t
  !> \brief   'config_riv_temp' namelist content
  type, public :: nml_config_riv_temp_t
    character(15) :: name = "config_riv_temp" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    real(dp) :: albedo_water !< albedo of open water
    real(dp) :: pt_a_water !< priestley taylor alpha parameter for PET on open water
    real(dp) :: emissivity_water !< emissivity of water
    real(dp) :: turb_heat_ex_coeff !< lateral heat exchange coefficient water <-> air
    integer(i4) :: max_iter !< maximum number of iterations
    real(dp) :: delta_iter !< convergence delta
    real(dp) :: step_iter !< step size for iterative solver
    character(256) :: riv_widths_file !< file name for river widths
    character(256) :: riv_widths_name !< variable name for river widths
    character(256), dimension(maxNoDomains) :: dir_riv_widths !< files for river widths
  contains
    procedure, public :: read => read_config_riv_temp
  end type nml_config_riv_temp_t
  !> 'config_riv_temp' namelist content
  type(nml_config_riv_temp_t), public :: nml_config_riv_temp

  !######## mo_coupling_type
  ! namelist /coupling/ &
  !   case, &
  !   meteo_timestep, &
  !   meteo_time_ref_endpoint, &
  !   meteo_expect_pre, &
  !   meteo_expect_temp, &
  !   meteo_expect_pet, &
  !   meteo_expect_tmin, &
  !   meteo_expect_tmax, &
  !   meteo_expect_netrad, &
  !   meteo_expect_absvappress, &
  !   meteo_expect_windspeed, &
  !   meteo_expect_ssrd, &
  !   meteo_expect_strd, &
  !   meteo_expect_tann
  !
  !> \class   nml_coupling_t
  !> \brief   'coupling' namelist content
  type, public :: nml_coupling_t
    character(8) :: name = "coupling" !< namelist name
    logical :: read_from_file = .true. !< whether the associated variables are already set by interfaces
    integer(i4) :: case !< coupling case
    integer(i4) :: meteo_timestep !< timestep for meteo-data from coupling
    logical :: meteo_time_ref_endpoint !< expect meteo has time reference point at end of associated time interval
    logical :: meteo_expect_pre !< expect meteo from coupling: [mm]      Precipitation
    logical :: meteo_expect_temp !< expect meteo from coupling: [degC]    Air temperature
    logical :: meteo_expect_pet !< expect meteo from coupling: [mm TS-1] Potential evapotranspiration
    logical :: meteo_expect_tmin !< expect meteo from coupling: [degC]    minimum daily air temperature
    logical :: meteo_expect_tmax !< expect meteo from coupling: [degC]    maximum daily air temperature
    logical :: meteo_expect_netrad !< expect meteo from coupling: [W m2]    net radiation
    logical :: meteo_expect_absvappress !< expect meteo from coupling: [Pa]      absolute vapour pressure
    logical :: meteo_expect_windspeed !< expect meteo from coupling: [m s-1]   windspeed
    logical :: meteo_expect_ssrd !< expect meteo from coupling: [W m2]    short wave radiation
    logical :: meteo_expect_strd !< expect meteo from coupling: [W m2]    long wave radiation
    logical :: meteo_expect_tann !< expect meteo from coupling: [degC]    annual mean air temperature
  contains
    procedure, public :: read => read_coupling
  end type nml_coupling_t
  !> 'coupling' namelist content
  type(nml_coupling_t), public :: nml_coupling

contains

  !> \brief Open namelist file and generate a new unit.
  subroutine open_new_nml(file, unit)
    use mo_message, only: error_message
    character(len = *), intent(in) :: file
    integer, intent(out) :: unit
    integer :: stat
    open(newunit=unit, file=file, iostat=stat, status='old', action='read', delim='apostrophe')
    if (stat .ne. 0) call error_message('open_new_nml: could not open namelist file ', trim(file))
  end subroutine open_new_nml

  !> \brief Close namelist file.
  subroutine close_nml(unit)
    use mo_message, only: error_message
    integer, intent(in) :: unit
    integer :: stat
    close(unit, iostat=stat)
    if (stat .ne. 0) call error_message('close_nml: could not close namelist file.')
  end subroutine close_nml

  !> \brief Read 'project_description' namelist content.
  subroutine read_project_description(self, file)
    implicit none
    class(nml_project_description_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(1024) :: project_details !< project including funding instituion., PI, etc.
    character(1024) :: setup_description !< any specific description of simulation
    character(1024) :: simulation_type !< e.g. seasonal forecast, climate projection, ...
    character(256) :: Conventions !< convention used for dataset
    character(1024) :: contact !< contact details, incl. PI name
    character(1024) :: mHM_details !< developing institution, specific mHM revision
    character(1024) :: history !< details on version/creation date

    namelist /project_description/ &
      project_details, &
      setup_description, &
      simulation_type, &
      Conventions, &
      contact, &
      mHM_details, &
      history

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=project_description)
      call close_nml(unit)
      self%project_details = project_details
      self%setup_description = setup_description
      self%simulation_type = simulation_type
      self%Conventions = Conventions
      self%contact = contact
      self%mHM_details = mHM_details
      self%history = history
      self%read_from_file = .false.
    end if
  end subroutine read_project_description

  !> \brief Read 'directories_general' namelist content.
  subroutine read_directories_general(self, file)
    implicit none
    class(nml_directories_general_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(256) :: dirConfigOut !< directory for config file output
    character(256) :: dirCommonFiles !< directory where common input files should be located
    character(256), dimension(maxNoDomains) :: mhm_file_RestartOut !< Directory where mhm output of restart is written
    character(256), dimension(maxNoDomains) :: mrm_file_RestartOut !< Directory where mrm output of restart is written
    character(256), dimension(maxNoDomains) :: dir_Morpho !< Directory where morphological files are located
    character(256), dimension(maxNoDomains) :: dir_LCover !< Directory where land cover files are located
    character(256), dimension(maxNoDomains) :: dir_Out !< Directory where output is written to
    character(256), dimension(maxNoDomains) :: file_LatLon !< Directory where the Lat Lon Files are located

    namelist /directories_general/ &
      dirConfigOut, &
      dirCommonFiles, &
      dir_Morpho, &
      dir_LCover, &
      dir_Out, &
      mhm_file_RestartOut, &
      mrm_file_RestartOut, &
      file_LatLon

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=directories_general)
      call close_nml(unit)
      self%dirConfigOut = dirConfigOut
      self%dirCommonFiles = dirCommonFiles
      self%dir_Morpho = dir_Morpho
      self%dir_LCover = dir_LCover
      self%dir_Out = dir_Out
      self%mhm_file_RestartOut = mhm_file_RestartOut
      self%mrm_file_RestartOut = mrm_file_RestartOut
      self%file_LatLon = file_LatLon
      self%read_from_file = .false.
    end if
  end subroutine read_directories_general

  !> \brief Read 'mainconfig' namelist content.
  subroutine read_mainconfig(self, file)
    implicit none
    class(nml_mainconfig_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: iFlag_cordinate_sys !< options model for the run cordinate system
    real(dp), dimension(maxNoDomains) :: resolution_Hydrology !< [m or degree] resolution of hydrology - Level 1
    integer(i4) :: nDomains !< number of domains
    integer(i4), dimension(maxNoDomains) :: L0Domain !< specify same index for domains to share L0_data to save memory
    logical :: write_restart !< flag to write restart
    integer(i4), dimension(maxNoDomains) :: read_opt_domain_data !< read domain specific optional data

    namelist /mainconfig/ &
      iFlag_cordinate_sys, &
      resolution_Hydrology, &
      nDomains, &
      L0Domain, &
      write_restart, &
      read_opt_domain_data

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=mainconfig)
      call close_nml(unit)
      self%iFlag_cordinate_sys = iFlag_cordinate_sys
      self%resolution_Hydrology = resolution_Hydrology
      self%nDomains = nDomains
      self%L0Domain = L0Domain
      self%write_restart = write_restart
      self%read_opt_domain_data = read_opt_domain_data
      self%read_from_file = .false.
    end if
  end subroutine read_mainconfig

  !> \brief Read 'processSelection' namelist content.
  subroutine read_processSelection(self, file)
    implicit none
    class(nml_processSelection_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4), dimension(nProcesses) :: processCase !< ! Choosen process description number

    namelist /processSelection/ &
      processCase

    if ( self%read_from_file ) then
      ! init the processCase matrix to 0 to be backward compatible
      ! if cases were added later (then there would be no values if not init here)
      processCase = 0_i4
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=processSelection)
      call close_nml(unit)
      self%processCase = processCase
      self%read_from_file = .false.
    end if
  end subroutine read_processSelection

  !> \brief Read 'LCover' namelist content.
  subroutine read_LCover(self, file)
    implicit none
    class(nml_LCover_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: nLCoverScene !< Number of land cover scene (lcs)
    integer(i4), dimension(maxNLCovers) :: LCoverYearStart !< starting year LCover
    integer(i4), dimension(maxNLCovers) :: LCoverYearEnd !< ending year LCover
    character(256), dimension(maxNLCovers) :: LCoverfName !< filename of Lcover file

    namelist /LCover/ &
      nLcoverScene, &
      LCoverYearStart, &
      LCoverYearEnd, &
      LCoverfName

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=LCover)
      call close_nml(unit)
      self%nLcoverScene = nLcoverScene
      self%LCoverYearStart = LCoverYearStart
      self%LCoverYearEnd = LCoverYearEnd
      self%LCoverfName = LCoverfName
      self%read_from_file = .false.
    end if
  end subroutine read_LCover

  !> \brief Read 'mainconfig_mhm_mrm' namelist content.
  subroutine read_mainconfig_mhm_mrm(self, file)
    implicit none
    class(nml_mainconfig_mhm_mrm_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: timeStep !< [h] simulation time step (= TS) in [h] either 1, 2, 3, 4, 6, 12 or 24
    real(dp), dimension(maxNoDomains) :: resolution_Routing !< resolution of Level-11 discharge routing [m or degree] per domain
    logical :: optimize !< Optimization (.true.) or Evaluation run (.false.)
    logical :: optimize_restart !< Optimization will be restarted from mo_<opti_method>.restart file (.true.)
    integer(i4) :: opti_method !< Optimization algorithm: 1 - DDS; 2 - Simulated Annealing; 3 - SCE
    integer(i4) :: opti_function !< Objective function
    logical :: read_restart !< flag for reading restart output
    logical :: mrm_read_river_network !< flag to read the river network for mRM (read_restart = .True. forces .True.)
    logical :: read_old_style_restart_bounds !< flag to use an old-style restart file created by mhm<=v5.11
    logical :: restart_reset_fluxes_states !< flag to reset fluxes and states read from restart to default values
    character(256), dimension(maxNoDomains) :: mhm_file_RestartIn !< mhm restart file paths
    character(256), dimension(maxNoDomains) :: mrm_file_RestartIn !< mrm restart file paths

    namelist /mainconfig_mhm_mrm/ &
      timestep, &
      resolution_Routing, &
      optimize, &
      optimize_restart, &
      opti_method, &
      opti_function, &
      read_restart, &
      mrm_read_river_network, &
      read_old_style_restart_bounds, &
      restart_reset_fluxes_states, &
      mhm_file_RestartIn, &
      mrm_file_RestartIn

    if ( self%read_from_file ) then
      ! set default values for optional arguments
      mrm_read_river_network = .false.
      read_old_style_restart_bounds = .false.
      restart_reset_fluxes_states = .false.
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=mainconfig_mhm_mrm)
      call close_nml(unit)
      self%timestep = timestep
      self%resolution_Routing = resolution_Routing
      self%optimize = optimize
      self%optimize_restart = optimize_restart
      self%opti_method = opti_method
      self%opti_function = opti_function
      self%read_restart = read_restart
      self%mrm_read_river_network = mrm_read_river_network
      self%read_old_style_restart_bounds = read_old_style_restart_bounds
      self%restart_reset_fluxes_states = restart_reset_fluxes_states
      self%mhm_file_RestartIn = mhm_file_RestartIn
      self%mrm_file_RestartIn = mrm_file_RestartIn
      self%read_from_file = .false.
    end if
  end subroutine read_mainconfig_mhm_mrm

  !> \brief Read 'optimization' namelist content.
  subroutine read_optimization(self, file)
    implicit none
    class(nml_optimization_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit, status
    integer(i4) :: nIterations !< number of iterations for optimization
    integer(i8) :: seed !< seed used for optimization, default: -9 --> system time
    real(dp) :: dds_r !< DDS: perturbation rate, default: 0.2
    real(dp) :: sa_temp !< SA:  initial temperature, default: -9.0 --> estimated
    integer(i4) :: sce_ngs !< SCE: # of complexes, default: 2
    integer(i4) :: sce_npg !< SCE: # of points per complex,default: -9 --> 2n+1
    integer(i4) :: sce_nps !< SCE: # of points per subcomplex,default: -9 --> n+1
    logical :: mcmc_opti !< MCMC: optimization (.true.) or only parameter uncertainty (.false.)
    real(dp), dimension(nerror_model) :: mcmc_error_params !< error model para (mcmc_opti=.false.) e.g. for opti_function=8: .01, .3

    namelist /optimization/ &
      nIterations, &
      seed, &
      dds_r, &
      sa_temp, &
      sce_ngs, &
      sce_npg, &
      sce_nps, &
      mcmc_opti, &
      mcmc_error_params

    if ( self%read_from_file ) then
      nIterations = 0_i4
      seed = -9_i8
      dds_r = 0.2_dp
      sa_temp = -9.0_dp
      sce_ngs = 2_i4
      sce_npg = -9_i4
      sce_nps = -9_i4
      mcmc_opti = .true.
      ! mcmc_error_params -> no defaults
      call open_new_nml(file, unit)
      call position_nml(self%name, unit, status=status)
      if (status == 0) read(unit, nml=optimization)
      call close_nml(unit)
      self%nIterations = nIterations
      self%seed = seed
      self%dds_r = dds_r
      self%sa_temp = sa_temp
      self%sce_ngs = sce_ngs
      self%sce_npg = sce_npg
      self%sce_nps = sce_nps
      self%mcmc_opti = mcmc_opti
      self%mcmc_error_params = mcmc_error_params
      self%read_from_file = .false.
    end if
  end subroutine read_optimization

  !> \brief Read 'time_periods' namelist content.
  subroutine read_time_periods(self, file)
    implicit none
    class(nml_time_periods_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4), dimension(maxNoDomains) :: warming_Days !< number of days for warm up period
    type(period), dimension(maxNoDomains) :: eval_Per !< time period for model evaluation

    namelist /time_periods/ &
      warming_Days, &
      eval_Per

      if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=time_periods)
      call close_nml(unit)
      self%warming_Days = warming_Days
      self%eval_Per = eval_Per
      self%read_from_file = .false.
    end if
  end subroutine read_time_periods

  !> \brief Read 'directories_mhm' namelist content.
  subroutine read_directories_mhm(self, file)
    implicit none
    class(nml_directories_mhm_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(256) :: inputFormat_meteo_forcings !< format of meteo input data (nc)
    !> .FALSE. to only warn about bound (lower, upper) violations in meteo files, default = .TRUE. - raise an error
    logical :: bound_error
    character(256), dimension(maxNoDomains) :: dir_meteo_header !< Directory where the meteo header file is located
    character(256), dimension(maxNoDomains) :: dir_Precipitation !< Directory where precipitation files are located
    character(256), dimension(maxNoDomains) :: dir_Temperature !< Directory where temperature files are located
    character(256), dimension(maxNoDomains) :: dir_ReferenceET !< Directory where reference-ET files are located
    character(256), dimension(maxNoDomains) :: dir_MinTemperature !< Directory where minimum temp. files are located
    character(256), dimension(maxNoDomains) :: dir_MaxTemperature !< Directory where maximum temp. files are located
    character(256), dimension(maxNoDomains) :: dir_absVapPressure !< Directory where abs. vap. pressure files are located
    character(256), dimension(maxNoDomains) :: dir_windspeed !< Directory where windspeed files are located
    character(256), dimension(maxNoDomains) :: dir_NetRadiation !< Directory where abs. vap. pressure files are located
    character(256), dimension(maxNoDomains) :: dir_Radiation !< riv-temp related: directory of (long/short-wave)radiation
    integer(i4), dimension(maxNoDomains) :: time_step_model_inputs !< frequency for reading meteo input

    namelist /directories_mhm/ &
      inputFormat_meteo_forcings, &
      bound_error, &
      dir_meteo_header, &
      dir_Precipitation, &
      dir_Temperature, &
      dir_ReferenceET, &
      dir_MinTemperature, &
      dir_MaxTemperature, &
      dir_absVapPressure, &
      dir_windspeed, &
      dir_NetRadiation, &
      dir_Radiation, &
      time_step_model_inputs

    if ( self%read_from_file ) then
      call set_sentinel(dir_meteo_header) ! set sentinal to check reading
      inputFormat_meteo_forcings = "nc"
      bound_error = .TRUE.
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=directories_mhm)
      call close_nml(unit)
      self%inputFormat_meteo_forcings = inputFormat_meteo_forcings
      self%bound_error = bound_error
      self%dir_meteo_header = dir_meteo_header
      self%dir_Precipitation = dir_Precipitation
      self%dir_Temperature = dir_Temperature
      self%dir_ReferenceET = dir_ReferenceET
      self%dir_MinTemperature = dir_MinTemperature
      self%dir_MaxTemperature = dir_MaxTemperature
      self%dir_absVapPressure = dir_absVapPressure
      self%dir_windspeed = dir_windspeed
      self%dir_NetRadiation = dir_NetRadiation
      self%dir_Radiation = dir_Radiation
      self%time_step_model_inputs = time_step_model_inputs
      self%read_from_file = .false.
    end if
  end subroutine read_directories_mhm

  !> \brief Read 'optional_data' namelist content.
  subroutine read_optional_data(self, file)
    implicit none
    class(nml_optional_data_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: nSoilHorizons_sm_input !< No. of mhm soil horizons equivalent to sm input
    character(256), dimension(maxNoDomains) :: dir_soil_moisture !< soil moisture input
    character(256), dimension(maxNoDomains) :: dir_neutrons !< ground albedo neutron input
    character(256), dimension(maxNoDomains) :: dir_evapotranspiration !< evapotranspiration input
    character(256), dimension(maxNoDomains) :: dir_TWS !< tws input
    character(256), dimension(maxNoDomains) :: dir_spf !< spf input
    integer(i4) :: timeStep_sm_input !< time step of optional data: sm
    integer(i4) :: timeStep_neutrons_input !< time step of optional data: neutrons
    integer(i4) :: timeStep_et_input !< time step of optional data: et
    integer(i4) :: timeStep_tws_input !< time step of optional data: tws
    integer(i4) :: timeStep_spf_input !< time step of optional data: tws
    integer(i4) :: weight_for_optional_data   ! weight of optional data in OF 67 to 69
    integer(i4) :: snow_water_equivalent_threshold_for_spf ! threshold swe to convert to spf for OF 75, 76


    namelist /optional_data/ &
      nSoilHorizons_sm_input, &
      dir_soil_moisture, &
      dir_neutrons, &
      dir_evapotranspiration, &
      dir_TWS, &
      dir_spf, &
      timeStep_sm_input, &
      timeStep_neutrons_input, &
      timeStep_et_input, &
      timeStep_tws_input, &
      timeStep_spf_input, &
      weight_for_optional_data, &
      snow_water_equivalent_threshold_for_spf

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=optional_data)
      call close_nml(unit)
      self%nSoilHorizons_sm_input = nSoilHorizons_sm_input
      self%dir_soil_moisture = dir_soil_moisture
      self%dir_neutrons = dir_neutrons
      self%dir_evapotranspiration = dir_evapotranspiration
      self%dir_TWS = dir_TWS
      self%dir_spf = dir_spf
      self%timeStep_sm_input = timeStep_sm_input
      self%timeStep_neutrons_input = timeStep_neutrons_input
      self%timeStep_et_input = timeStep_et_input
      self%timeStep_tws_input = timeStep_tws_input
      self%timeStep_spf_input = timeStep_spf_input
      self%weight_for_optional_data = weight_for_optional_data
      self%snow_water_equivalent_threshold_for_spf = snow_water_equivalent_threshold_for_spf
      self%read_from_file = .false.
    end if
  end subroutine read_optional_data

  !> \brief Read 'panevapo' namelist content.
  subroutine read_panevapo(self, file)
    implicit none
    class(nml_panevapo_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(int(YearMonths, i4)) :: evap_coeff !< [-] Evap. coef. for free-water surfaces

    namelist /panevapo/ &
      evap_coeff

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=panevapo)
      call close_nml(unit)
      self%evap_coeff = evap_coeff
      self%read_from_file = .false.
    end if
  end subroutine read_panevapo

  !> \brief Read 'nightdayratio' namelist content.
  subroutine read_nightdayratio(self, file)
    implicit none
    class(nml_nightdayratio_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    logical :: read_meteo_weights !< read weights for meteo data
    real(dp), dimension(int(YearMonths, i4)) :: fnight_prec !< [-] Night ratio precipitation < 1
    real(dp), dimension(int(YearMonths, i4)) :: fnight_pet !< [-] Night ratio PET  < 1
    real(dp), dimension(int(YearMonths, i4)) :: fnight_temp !< [-] Night factor mean temp
    real(dp), dimension(int(YearMonths, i4)) :: fnight_ssrd !< [-] Night factor short-wave rad.
    real(dp), dimension(int(YearMonths, i4)) :: fnight_strd !< [-] Night factor long-wave rad.

    namelist /nightdayratio/ &
      read_meteo_weights, &
      fnight_prec, &
      fnight_pet, &
      fnight_temp, &
      fnight_ssrd, &
      fnight_strd

    if ( self%read_from_file ) then
      ! default values for long/shortwave rad.
      fnight_ssrd = 0.0_dp
      fnight_strd = 0.45_dp
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=nightdayratio)
      call close_nml(unit)
      self%read_meteo_weights = read_meteo_weights
      self%fnight_prec = fnight_prec
      self%fnight_pet = fnight_pet
      self%fnight_temp = fnight_temp
      self%fnight_ssrd = fnight_ssrd
      self%fnight_strd = fnight_strd
      self%read_from_file = .false.
    end if
  end subroutine read_nightdayratio

  !> \brief Read 'nloutputresults' namelist content.
  subroutine read_nloutputresults(self, file)
    implicit none
    class(nml_nloutputresults_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: output_deflate_level !< deflate level in nc files
    logical :: output_double_precision !< output precision in nc files
    integer(i4) :: output_time_reference !< time reference point location in output nc files
    integer(i4) :: timeStep_model_outputs !< timestep for writing model outputs
    logical, dimension(nOutFlxState) :: outputFlxState !< Define model outputs see "mhm_outputs.nml"

    namelist /nloutputresults/ &
      output_deflate_level, &
      output_double_precision, &
      output_time_reference, &
      timeStep_model_outputs, &
      outputFlxState

    if ( self%read_from_file ) then
      ! default values
      output_deflate_level = 6
      output_double_precision = .true.
      output_time_reference = 0
      outputFlxState = .FALSE.
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=nloutputresults)
      call close_nml(unit)
      self%output_deflate_level = output_deflate_level
      self%output_double_precision = output_double_precision
      self%timeStep_model_outputs = timeStep_model_outputs
      self%output_time_reference = output_time_reference
      self%outputFlxState = outputFlxState
      self%read_from_file = .false.
    end if
  end subroutine read_nloutputresults

  !> \brief Read 'baseflow_config' namelist content.
  subroutine read_baseflow_config(self, file)
    implicit none
    class(nml_baseflow_config_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    logical :: BFI_calc !< calculate observed BFI from gauges with Eckhardt filter
    real(dp), dimension(maxNoDomains) :: BFI_obs !< given base-flow index per domain

    namelist /baseflow_config/ &
      BFI_calc, &
      BFI_obs

    if ( self%read_from_file ) then
      BFI_calc = .false. ! default value
      BFI_obs = -1.0_dp  ! negative value to flag missing values
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=baseflow_config)
      call close_nml(unit)
      self%BFI_calc = BFI_calc
      self%BFI_obs = BFI_obs
      self%read_from_file = .false.
    end if
  end subroutine read_baseflow_config

  !> \brief Read 'directories_mpr' namelist content.
  subroutine read_directories_mpr(self, file)
    implicit none
    class(nml_directories_mpr_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(256), dimension(maxNoDomains) :: dir_gridded_LAI !< directory of gridded LAI data, used when timeStep_LAI_input<0

    namelist /directories_mpr/ &
      dir_gridded_LAI

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=directories_mpr)
      call close_nml(unit)
      self%dir_gridded_LAI = dir_gridded_LAI
      self%read_from_file = .false.
    end if
  end subroutine read_directories_mpr

  !> \brief Read 'soildata' namelist content.
  subroutine read_soildata(self, file)
    implicit none
    class(nml_soildata_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: iFlag_soilDB !< options to handle different soil databases
    real(dp) :: tillageDepth !< [mm] Soil depth down to which organic
    integer(i4) :: nSoilHorizons_mHM !< Number of horizons to model
    real(dp), dimension(maxNoSoilHorizons) :: soil_Depth !< depth of the single horizons

    namelist /soildata/ &
      iFlag_soilDB, &
      tillageDepth, &
      nSoilHorizons_mHM, &
      soil_Depth

    if ( self%read_from_file ) then
      soil_Depth = 0.0_dp ! default soil depth
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=soildata)
      call close_nml(unit)
      self%iFlag_soilDB = iFlag_soilDB
      self%tillageDepth = tillageDepth
      self%nSoilHorizons_mHM = nSoilHorizons_mHM
      self%soil_Depth = soil_Depth
      self%read_from_file = .false.
    end if
  end subroutine read_soildata

  !> \brief Read 'lai_data_information' namelist content.
  subroutine read_lai_data_information(self, file)
    implicit none
    class(nml_lai_data_information_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(256) :: inputFormat_gridded_LAI !< format of gridded LAI data (nc only)
    integer(i4) :: timeStep_LAI_input !< time step of gridded LAI input

    namelist /lai_data_information/ &
      inputFormat_gridded_LAI, &
      timeStep_LAI_input

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=lai_data_information)
      call close_nml(unit)
      self%inputFormat_gridded_LAI = inputFormat_gridded_LAI
      self%timeStep_LAI_input = timeStep_LAI_input
      self%read_from_file = .false.
    end if
  end subroutine read_lai_data_information

  !> \brief Read 'lcover_mpr' namelist content.
  subroutine read_lcover_mpr(self, file)
    implicit none
    class(nml_lcover_mpr_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp) :: fracSealed_cityArea !< fraction of area within city assumed to be perfectly sealed [0-1]

    namelist /lcover_mpr/ &
      fracSealed_cityArea

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=lcover_mpr)
      call close_nml(unit)
      self%fracSealed_cityArea = fracSealed_cityArea
      self%read_from_file = .false.
    end if
  end subroutine read_lcover_mpr

  !> \brief Read 'interception1' namelist content.
  subroutine read_interception1(self, file)
    implicit none
    class(nml_interception1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: canopyInterceptionFactor !< multiplier to relate LAI to interception storage [-]

    namelist /interception1/ &
      canopyInterceptionFactor

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=interception1)
      call close_nml(unit)
      self%canopyInterceptionFactor = canopyInterceptionFactor
      self%read_from_file = .false.
    end if
  end subroutine read_interception1

  !> \brief Read 'snow1' namelist content.
  subroutine read_snow1(self, file)
    implicit none
    class(nml_snow1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: snowTreshholdTemperature !< Threshold for rain/snow partitioning [degC]
    real(dp), dimension(nColPars) :: degreeDayFactor_forest !< forest: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: degreeDayFactor_impervious !< impervious: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: degreeDayFactor_pervious !< pervious: deg day factors to determine melting flux [m degC-1]
    real(dp), dimension(nColPars) :: increaseDegreeDayFactorByPrecip !< increase of deg day factor in case of precipitation [degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_forest !< forest: maximum values for degree day factor [m degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_impervious !< impervious: maximum values for degree day factor [m degC-1]
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_pervious !< pervious: maximum values for degree day factor [m degC-1]

    namelist /snow1/ &
      snowTreshholdTemperature, &
      degreeDayFactor_forest, &
      degreeDayFactor_impervious, &
      degreeDayFactor_pervious, &
      increaseDegreeDayFactorByPrecip, &
      maxDegreeDayFactor_forest, &
      maxDegreeDayFactor_impervious, &
      maxDegreeDayFactor_pervious

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=snow1)
      call close_nml(unit)
      self%snowTreshholdTemperature = snowTreshholdTemperature
      self%degreeDayFactor_forest = degreeDayFactor_forest
      self%degreeDayFactor_impervious = degreeDayFactor_impervious
      self%degreeDayFactor_pervious = degreeDayFactor_pervious
      self%increaseDegreeDayFactorByPrecip = increaseDegreeDayFactorByPrecip
      self%maxDegreeDayFactor_forest = maxDegreeDayFactor_forest
      self%maxDegreeDayFactor_impervious = maxDegreeDayFactor_impervious
      self%maxDegreeDayFactor_pervious = maxDegreeDayFactor_pervious
      self%read_from_file = .false.
    end if
  end subroutine read_snow1

  !> \brief Read 'soilmoisture1' namelist content.
  subroutine read_soilmoisture1(self, file)
    implicit none
    class(nml_soilmoisture1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor

    namelist /soilmoisture1/ &
      orgMatterContent_forest, &
      orgMatterContent_impervious, &
      orgMatterContent_pervious, &
      PTF_lower66_5_constant, &
      PTF_lower66_5_clay, &
      PTF_lower66_5_Db, &
      PTF_higher66_5_constant, &
      PTF_higher66_5_clay, &
      PTF_higher66_5_Db, &
      PTF_Ks_constant, &
      PTF_Ks_sand, &
      PTF_Ks_clay, &
      PTF_Ks_curveSlope, &
      rootFractionCoefficient_forest, &
      rootFractionCoefficient_impervious, &
      rootFractionCoefficient_pervious, &
      infiltrationShapeFactor

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=soilmoisture1)
      call close_nml(unit)
      self%orgMatterContent_forest = orgMatterContent_forest
      self%orgMatterContent_impervious = orgMatterContent_impervious
      self%orgMatterContent_pervious = orgMatterContent_pervious
      self%PTF_lower66_5_constant = PTF_lower66_5_constant
      self%PTF_lower66_5_clay = PTF_lower66_5_clay
      self%PTF_lower66_5_Db = PTF_lower66_5_Db
      self%PTF_higher66_5_constant = PTF_higher66_5_constant
      self%PTF_higher66_5_clay = PTF_higher66_5_clay
      self%PTF_higher66_5_Db = PTF_higher66_5_Db
      self%PTF_Ks_constant = PTF_Ks_constant
      self%PTF_Ks_sand = PTF_Ks_sand
      self%PTF_Ks_clay = PTF_Ks_clay
      self%PTF_Ks_curveSlope = PTF_Ks_curveSlope
      self%rootFractionCoefficient_forest = rootFractionCoefficient_forest
      self%rootFractionCoefficient_impervious = rootFractionCoefficient_impervious
      self%rootFractionCoefficient_pervious = rootFractionCoefficient_pervious
      self%infiltrationShapeFactor = infiltrationShapeFactor
      self%read_from_file = .false.
    end if
  end subroutine read_soilmoisture1

  !> \brief Read 'soilmoisture2' namelist content.
  subroutine read_soilmoisture2(self, file)
    implicit none
    class(nml_soilmoisture2_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1 !< soil moisture threshod for jarvis model

    namelist /soilmoisture2/ &
      orgMatterContent_forest, &
      orgMatterContent_impervious, &
      orgMatterContent_pervious, &
      PTF_lower66_5_constant, &
      PTF_lower66_5_clay, &
      PTF_lower66_5_Db, &
      PTF_higher66_5_constant, &
      PTF_higher66_5_clay, &
      PTF_higher66_5_Db, &
      PTF_Ks_constant, &
      PTF_Ks_sand, &
      PTF_Ks_clay, &
      PTF_Ks_curveSlope, &
      rootFractionCoefficient_forest, &
      rootFractionCoefficient_impervious, &
      rootFractionCoefficient_pervious, &
      infiltrationShapeFactor, &
      jarvis_sm_threshold_c1

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=soilmoisture2)
      call close_nml(unit)
      self%orgMatterContent_forest = orgMatterContent_forest
      self%orgMatterContent_impervious = orgMatterContent_impervious
      self%orgMatterContent_pervious = orgMatterContent_pervious
      self%PTF_lower66_5_constant = PTF_lower66_5_constant
      self%PTF_lower66_5_clay = PTF_lower66_5_clay
      self%PTF_lower66_5_Db = PTF_lower66_5_Db
      self%PTF_higher66_5_constant = PTF_higher66_5_constant
      self%PTF_higher66_5_clay = PTF_higher66_5_clay
      self%PTF_higher66_5_Db = PTF_higher66_5_Db
      self%PTF_Ks_constant = PTF_Ks_constant
      self%PTF_Ks_sand = PTF_Ks_sand
      self%PTF_Ks_clay = PTF_Ks_clay
      self%PTF_Ks_curveSlope = PTF_Ks_curveSlope
      self%rootFractionCoefficient_forest = rootFractionCoefficient_forest
      self%rootFractionCoefficient_impervious = rootFractionCoefficient_impervious
      self%rootFractionCoefficient_pervious = rootFractionCoefficient_pervious
      self%infiltrationShapeFactor = infiltrationShapeFactor
      self%jarvis_sm_threshold_c1 = jarvis_sm_threshold_c1
      self%read_from_file = .false.
    end if
  end subroutine read_soilmoisture2

  !> \brief Read 'soilmoisture3' namelist content.
  subroutine read_soilmoisture3(self, file)
    implicit none
    class(nml_soilmoisture3_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: FCmin_glob !< global field capacity minimum
    real(dp), dimension(nColPars) :: FCdelta_glob !< difference between global field capacity minimum and maximum
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand !< threshold for actual ET reduction for sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay !< threshold for actual ET reduction for clay
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1 !< soil moisture threshod for jarvis model

    namelist /soilmoisture3/ &
      orgMatterContent_forest, &
      orgMatterContent_impervious, &
      orgMatterContent_pervious, &
      PTF_lower66_5_constant, &
      PTF_lower66_5_clay, &
      PTF_lower66_5_Db, &
      PTF_higher66_5_constant, &
      PTF_higher66_5_clay, &
      PTF_higher66_5_Db, &
      PTF_Ks_constant, &
      PTF_Ks_sand, &
      PTF_Ks_clay, &
      PTF_Ks_curveSlope, &
      rootFractionCoefficient_forest, &
      rootFractionCoefficient_impervious, &
      rootFractionCoefficient_pervious, &
      infiltrationShapeFactor, &
      rootFractionCoefficient_sand, &
      rootFractionCoefficient_clay, &
      FCmin_glob, &
      FCdelta_glob, &
      jarvis_sm_threshold_c1

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=soilmoisture3)
      call close_nml(unit)
      self%orgMatterContent_forest = orgMatterContent_forest
      self%orgMatterContent_impervious = orgMatterContent_impervious
      self%orgMatterContent_pervious = orgMatterContent_pervious
      self%PTF_lower66_5_constant = PTF_lower66_5_constant
      self%PTF_lower66_5_clay = PTF_lower66_5_clay
      self%PTF_lower66_5_Db = PTF_lower66_5_Db
      self%PTF_higher66_5_constant = PTF_higher66_5_constant
      self%PTF_higher66_5_clay = PTF_higher66_5_clay
      self%PTF_higher66_5_Db = PTF_higher66_5_Db
      self%PTF_Ks_constant = PTF_Ks_constant
      self%PTF_Ks_sand = PTF_Ks_sand
      self%PTF_Ks_clay = PTF_Ks_clay
      self%PTF_Ks_curveSlope = PTF_Ks_curveSlope
      self%rootFractionCoefficient_forest = rootFractionCoefficient_forest
      self%rootFractionCoefficient_impervious = rootFractionCoefficient_impervious
      self%rootFractionCoefficient_pervious = rootFractionCoefficient_pervious
      self%infiltrationShapeFactor = infiltrationShapeFactor
      self%rootFractionCoefficient_sand = rootFractionCoefficient_sand
      self%rootFractionCoefficient_clay = rootFractionCoefficient_clay
      self%FCmin_glob = FCmin_glob
      self%FCdelta_glob = FCdelta_glob
      self%jarvis_sm_threshold_c1 = jarvis_sm_threshold_c1
      self%read_from_file = .false.
    end if
  end subroutine read_soilmoisture3

  !> \brief Read 'soilmoisture4' namelist content.
  subroutine read_soilmoisture4(self, file)
    implicit none
    class(nml_soilmoisture4_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: orgMatterContent_forest !< organic matter content [%] for forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious !< organic matter content [%] for impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious !< organic matter content [%] for pervious
    !> Zacharias PTF parameters below 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay !< multiplier for clay constant (see PTF_lower66_5_constant)
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db !< multiplier for mineral bulk density (see PTF_lower66_5_constant)
    !> Zacharias PTF parameters above 66.5 % sand content (Zacharias et al., 2007, doi:10.2136/sssaj2006.0098)
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay !< multiplier for clay constant (see PTF_higher66_5_constant)
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db !< multiplier for mineral bulk density (see PTF_higher66_5_constant)
    !> PTF parameters for saturated hydraulic conductivity after Cosby et al. (1984)
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand !< multiplier for sand (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_clay !< multiplier for clay (see PTF_Ks_constant)
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope !< unit conversion factor from inch/h to cm/d
    !> shape factor for root distribution with depth, which follows an exponential function [-] for forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    !> shape factor for root distribution with depth, which follows an exponential function [-] for impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    !> shape factor for root distribution with depth, which follows an exponential function [-] for pervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    !> shape factor for partitioning effective precipitation into runoff and infiltration based on soil wetness [-]
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: FCmin_glob !< global field capacity minimum
    real(dp), dimension(nColPars) :: FCdelta_glob !< difference between global field capacity minimum and maximum
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand !< threshold for actual ET reduction for sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay !< threshold for actual ET reduction for clay

    namelist /soilmoisture4/ &
      orgMatterContent_forest, &
      orgMatterContent_impervious, &
      orgMatterContent_pervious, &
      PTF_lower66_5_constant, &
      PTF_lower66_5_clay, &
      PTF_lower66_5_Db, &
      PTF_higher66_5_constant, &
      PTF_higher66_5_clay, &
      PTF_higher66_5_Db, &
      PTF_Ks_constant, &
      PTF_Ks_sand, &
      PTF_Ks_clay, &
      PTF_Ks_curveSlope, &
      rootFractionCoefficient_forest, &
      rootFractionCoefficient_impervious, &
      rootFractionCoefficient_pervious, &
      infiltrationShapeFactor, &
      rootFractionCoefficient_sand, &
      rootFractionCoefficient_clay, &
      FCmin_glob, &
      FCdelta_glob

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=soilmoisture4)
      call close_nml(unit)
      self%orgMatterContent_forest = orgMatterContent_forest
      self%orgMatterContent_impervious = orgMatterContent_impervious
      self%orgMatterContent_pervious = orgMatterContent_pervious
      self%PTF_lower66_5_constant = PTF_lower66_5_constant
      self%PTF_lower66_5_clay = PTF_lower66_5_clay
      self%PTF_lower66_5_Db = PTF_lower66_5_Db
      self%PTF_higher66_5_constant = PTF_higher66_5_constant
      self%PTF_higher66_5_clay = PTF_higher66_5_clay
      self%PTF_higher66_5_Db = PTF_higher66_5_Db
      self%PTF_Ks_constant = PTF_Ks_constant
      self%PTF_Ks_sand = PTF_Ks_sand
      self%PTF_Ks_clay = PTF_Ks_clay
      self%PTF_Ks_curveSlope = PTF_Ks_curveSlope
      self%rootFractionCoefficient_forest = rootFractionCoefficient_forest
      self%rootFractionCoefficient_impervious = rootFractionCoefficient_impervious
      self%rootFractionCoefficient_pervious = rootFractionCoefficient_pervious
      self%infiltrationShapeFactor = infiltrationShapeFactor
      self%rootFractionCoefficient_sand = rootFractionCoefficient_sand
      self%rootFractionCoefficient_clay = rootFractionCoefficient_clay
      self%FCmin_glob = FCmin_glob
      self%FCdelta_glob = FCdelta_glob
      self%read_from_file = .false.
    end if
  end subroutine read_soilmoisture4

  !> \brief Read 'directrunoff1' namelist content.
  subroutine read_directrunoff1(self, file)
    implicit none
    class(nml_directrunoff1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: imperviousStorageCapacity !< direct Runoff: Sealed Area storage capacity

    namelist /directrunoff1/ &
      imperviousStorageCapacity

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=directrunoff1)
      call close_nml(unit)
      self%imperviousStorageCapacity = imperviousStorageCapacity
      self%read_from_file = .false.
    end if
  end subroutine read_directrunoff1

  !> \brief Read 'petminus1' namelist content.
  subroutine read_petminus1(self, file)
    implicit none
    class(nml_petminus1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: PET_a_forest !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_a_impervious !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_a_pervious !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_b !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    real(dp), dimension(nColPars) :: PET_c !< DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET

    namelist /petminus1/ &
      PET_a_forest, &
      PET_a_impervious, &
      PET_a_pervious, &
      PET_b, &
      PET_c

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=petminus1)
      call close_nml(unit)
      self%PET_a_forest = PET_a_forest
      self%PET_a_impervious = PET_a_impervious
      self%PET_a_pervious = PET_a_pervious
      self%PET_b = PET_b
      self%PET_c = PET_c
      self%read_from_file = .false.
    end if
  end subroutine read_petminus1

  !> \brief Read 'pet0' namelist content.
  subroutine read_pet0(self, file)
    implicit none
    class(nml_pet0_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: minCorrectionFactorPET !< minimum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: maxCorrectionFactorPET !< maximum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: aspectTresholdPET !< aspect threshold for PET correction with aspect

    namelist /pet0/ &
      minCorrectionFactorPET, &
      maxCorrectionFactorPET, &
      aspectTresholdPET

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=pet0)
      call close_nml(unit)
      self%minCorrectionFactorPET = minCorrectionFactorPET
      self%maxCorrectionFactorPET = maxCorrectionFactorPET
      self%aspectTresholdPET = aspectTresholdPET
      self%read_from_file = .false.
    end if
  end subroutine read_pet0

  !> \brief Read 'pet1' namelist content.
  subroutine read_pet1(self, file)
    implicit none
    class(nml_pet1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: minCorrectionFactorPET !< minimum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: maxCorrectionFactorPET !< maximum factor for PET correction with aspect
    real(dp), dimension(nColPars) :: aspectTresholdPET !< aspect threshold for PET correction with aspect
    real(dp), dimension(nColPars) :: HargreavesSamaniCoeff !< coefficient for Hargreaves Samani

    namelist /pet1/ &
    minCorrectionFactorPET, &
    maxCorrectionFactorPET, &
    aspectTresholdPET, &
    HargreavesSamaniCoeff

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=pet1)
      call close_nml(unit)
      self%minCorrectionFactorPET = minCorrectionFactorPET
      self%maxCorrectionFactorPET = maxCorrectionFactorPET
      self%aspectTresholdPET = aspectTresholdPET
      self%HargreavesSamaniCoeff = HargreavesSamaniCoeff
      self%read_from_file = .false.
    end if
  end subroutine read_pet1

  !> \brief Read 'pet2' namelist content.
  subroutine read_pet2(self, file)
    implicit none
    class(nml_pet2_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: PriestleyTaylorCoeff !< Priestley-Taylor coefficient
    real(dp), dimension(nColPars) :: PriestleyTaylorLAIcorr !< Priestley-Taylor LAI correction factor

    namelist /pet2/ &
      PriestleyTaylorCoeff, &
      PriestleyTaylorLAIcorr

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=pet2)
      call close_nml(unit)
      self%PriestleyTaylorCoeff = PriestleyTaylorCoeff
      self%PriestleyTaylorLAIcorr = PriestleyTaylorLAIcorr
      self%read_from_file = .false.
    end if
  end subroutine read_pet2

  !> \brief Read 'pet3' namelist content.
  subroutine read_pet3(self, file)
    implicit none
    class(nml_pet3_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: canopyheigth_forest !< canopy height for foreset
    real(dp), dimension(nColPars) :: canopyheigth_impervious !< canopy height for impervious
    real(dp), dimension(nColPars) :: canopyheigth_pervious !< canopy height for pervious
    real(dp), dimension(nColPars) :: displacementheight_coeff !< displacement height coefficient
    real(dp), dimension(nColPars) :: roughnesslength_momentum_coeff !< roughness length momentum coefficient
    real(dp), dimension(nColPars) :: roughnesslength_heat_coeff !< roughness length heat coefficient
    real(dp), dimension(nColPars) :: stomatal_resistance !< stomatal resistance

    namelist /pet3/ &
      canopyheigth_forest, &
      canopyheigth_impervious, &
      canopyheigth_pervious, &
      displacementheight_coeff, &
      roughnesslength_momentum_coeff, &
      roughnesslength_heat_coeff, &
      stomatal_resistance

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=pet3)
      call close_nml(unit)
      self%canopyheigth_forest = canopyheigth_forest
      self%canopyheigth_impervious = canopyheigth_impervious
      self%canopyheigth_pervious = canopyheigth_pervious
      self%displacementheight_coeff = displacementheight_coeff
      self%roughnesslength_momentum_coeff = roughnesslength_momentum_coeff
      self%roughnesslength_heat_coeff = roughnesslength_heat_coeff
      self%stomatal_resistance = stomatal_resistance
      self%read_from_file = .false.
    end if
  end subroutine read_pet3

  !> \brief Read 'interflow1' namelist content.
  subroutine read_interflow1(self, file)
    implicit none
    class(nml_interflow1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: interflowStorageCapacityFactor !< interflow storage capacity factor
    real(dp), dimension(nColPars) :: interflowRecession_slope !< multiplier for slope to derive interflow recession constant
    !> multiplier to derive fast interflow recession constant for forest
    real(dp), dimension(nColPars) :: fastInterflowRecession_forest
    !> multiplier for variability of saturated hydraulic conductivity to derive slow interflow recession constant
    real(dp), dimension(nColPars) :: slowInterflowRecession_Ks
    !> multiplier for variability of saturated hydraulic conductivity to derive slow interflow exponent
    real(dp), dimension(nColPars) :: exponentSlowInterflow

    namelist /interflow1/ &
      interflowStorageCapacityFactor, &
      interflowRecession_slope, &
      fastInterflowRecession_forest, &
      slowInterflowRecession_Ks, &
      exponentSlowInterflow

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=interflow1)
      call close_nml(unit)
      self%interflowStorageCapacityFactor = interflowStorageCapacityFactor
      self%interflowRecession_slope = interflowRecession_slope
      self%fastInterflowRecession_forest = fastInterflowRecession_forest
      self%slowInterflowRecession_Ks = slowInterflowRecession_Ks
      self%exponentSlowInterflow = exponentSlowInterflow
      self%read_from_file = .false.
    end if
  end subroutine read_interflow1

  !> \brief Read 'percolation1' namelist content.
  subroutine read_percolation1(self, file)
    implicit none
    class(nml_percolation1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: rechargeCoefficient !< recharge coefficient
    real(dp), dimension(nColPars) :: rechargeFactor_karstic !< recharge factor for karstic percolation
    real(dp), dimension(nColPars) :: gain_loss_GWreservoir_karstic !< gain loss in ground water reservoir for karstic

    namelist /percolation1/ &
      rechargeCoefficient, &
      rechargeFactor_karstic, &
      gain_loss_GWreservoir_karstic

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=percolation1)
      call close_nml(unit)
      self%rechargeCoefficient = rechargeCoefficient
      self%rechargeFactor_karstic = rechargeFactor_karstic
      self%gain_loss_GWreservoir_karstic = gain_loss_GWreservoir_karstic
      self%read_from_file = .false.
    end if
  end subroutine read_percolation1

  !> \brief Read 'neutrons1' namelist content.
  subroutine read_neutrons1(self, file)
    implicit none
    class(nml_neutrons1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: Desilets_N0 !< Desilets N0 parameter
    real(dp), dimension(nColPars) :: Desilets_LW0 !< Desilets LW0 parameter
    real(dp), dimension(nColPars) :: Desilets_LW1 !< Desilets LW1 parameter

    namelist /neutrons1/ &
      Desilets_N0, &
      Desilets_LW0, &
      Desilets_LW1

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=neutrons1)
      call close_nml(unit)
      self%Desilets_N0 = Desilets_N0
      self%Desilets_LW0 = Desilets_LW0
      self%Desilets_LW1 = Desilets_LW1
      self%read_from_file = .false.
    end if
  end subroutine read_neutrons1

  !> \brief Read 'neutrons2' namelist content.
  subroutine read_neutrons2(self, file)
    implicit none
    class(nml_neutrons2_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: COSMIC_N0 !< cosmic N0 parameter
    real(dp), dimension(nColPars) :: COSMIC_N1 !< cosmic N1 parameter
    real(dp), dimension(nColPars) :: COSMIC_N2 !< cosmic N2 parameter
    real(dp), dimension(nColPars) :: COSMIC_alpha0 !< cosmic alpha0 parameter
    real(dp), dimension(nColPars) :: COSMIC_alpha1 !< cosmic alpha1 parameter
    real(dp), dimension(nColPars) :: COSMIC_L30 !< cosmic L30 parameter
    real(dp), dimension(nColPars) :: COSMIC_L31 !< cosmic L31 parameter
    real(dp), dimension(nColPars) :: COSMIC_LW0 !< cosmic LW0 parameter
    real(dp), dimension(nColPars) :: COSMIC_LW1 !< cosmic LW1 parameter

    namelist /neutrons2/ &
      COSMIC_N0, &
      COSMIC_N1, &
      COSMIC_N2, &
      COSMIC_alpha0, &
      COSMIC_alpha1, &
      COSMIC_L30, &
      COSMIC_L31, &
      COSMIC_LW0, &
      COSMIC_LW1

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=neutrons2)
      call close_nml(unit)
      self%COSMIC_N0 = COSMIC_N0
      self%COSMIC_N1 = COSMIC_N1
      self%COSMIC_N2 = COSMIC_N2
      self%COSMIC_alpha0 = COSMIC_alpha0
      self%COSMIC_alpha1 = COSMIC_alpha1
      self%COSMIC_L30 = COSMIC_L30
      self%COSMIC_L31 = COSMIC_L31
      self%COSMIC_LW0 = COSMIC_LW0
      self%COSMIC_LW1 = COSMIC_LW1
      self%read_from_file = .false.
    end if
  end subroutine read_neutrons2

  !> \brief Read 'geoparameter' namelist content.
  subroutine read_geoparameter(self, file)
    implicit none
    class(nml_geoparameter_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    !> geological parameters (ordering according to file 'geology_classdefinition.txt')
    real(dp), dimension(maxGeoUnit, nColPars) :: GeoParam

    namelist /geoparameter/ &
      GeoParam

    if ( self%read_from_file ) then
      GeoParam = nodata_dp
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=geoparameter)
      call close_nml(unit)
      self%GeoParam = GeoParam
      self%read_from_file = .false.
    end if
  end subroutine read_geoparameter

  !> \brief Read 'mainconfig_mrm' namelist content.
  subroutine read_mainconfig_mrm(self, file)
    implicit none
    class(nml_mainconfig_mrm_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    logical :: ALMA_convention !< flag for ALMA convention (see http://www.lmd.jussieu.fr/~polcher/ALMA/convention_3.html)
    character(256) :: filenameTotalRunoff !< Filename of simulated total runoff file
    character(256) :: varnameTotalRunoff !< variable name of total runoff
    logical :: gw_coupling !< switch to enable ground water coupling

    namelist /mainconfig_mrm/ &
      ALMA_convention, &
      filenameTotalRunoff, &
      varnameTotalRunoff, &
      gw_coupling

    if ( self%read_from_file ) then
      ALMA_convention = .false.
      filenameTotalRunoff = 'total_runoff'
      varnameTotalRunoff = 'total_runoff'
      gw_coupling = .false.
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=mainconfig_mrm)
      call close_nml(unit)
      self%ALMA_convention = ALMA_convention
      self%filenameTotalRunoff = filenameTotalRunoff
      self%varnameTotalRunoff = varnameTotalRunoff
      self%gw_coupling = gw_coupling
      self%read_from_file = .false.
    end if
  end subroutine read_mainconfig_mrm

  !> \brief Read 'directories_mrm' namelist content.
  subroutine read_directories_mrm(self, file)
    implicit none
    class(nml_directories_mrm_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    character(256), dimension(maxNoDomains) :: dir_Gauges !< directory containing gauge time series
    character(256), dimension(maxNoDomains) :: dir_Total_Runoff !< directory where simulated runoff can be found
    character(256), dimension(maxNoDomains) :: dir_Bankfull_Runoff !< directory where runoff at bankfull conditions can be found

    namelist /directories_mrm/ &
      dir_Gauges, &
      dir_Total_Runoff, &
      dir_Bankfull_Runoff

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=directories_mrm)
      call close_nml(unit)
      self%dir_Gauges = dir_Gauges
      self%dir_Total_Runoff = dir_Total_Runoff
      self%dir_Bankfull_Runoff = dir_Bankfull_Runoff
      self%read_from_file = .false.
    end if
  end subroutine read_directories_mrm

  !> \brief Read 'evaluation_gauges' namelist content.
  subroutine read_evaluation_gauges(self, file)
    implicit none
    class(nml_evaluation_gauges_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: nGaugesTotal !< Number of evaluation gauges for all domains
    integer(i4), dimension(maxNoDomains) :: NoGauges_domain !< number of gauges per domain
    integer(i4), dimension(maxNoDomains, maxNoGauges) :: Gauge_id !< gauge ID for each gauge
    character(256), dimension(maxNoDomains, maxNoGauges) :: Gauge_filename !< filename for each gauge time series

    namelist /evaluation_gauges/ &
      nGaugesTotal, &
      NoGauges_domain, &
      Gauge_id, &
      gauge_filename

    if ( self%read_from_file ) then
      nGaugesTotal = nodata_i4
      NoGauges_domain = nodata_i4
      Gauge_id = nodata_i4
      gauge_filename = num2str(nodata_i4)
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=evaluation_gauges)
      call close_nml(unit)
      self%nGaugesTotal = nGaugesTotal
      self%NoGauges_domain = NoGauges_domain
      self%Gauge_id = Gauge_id
      self%gauge_filename = gauge_filename
      self%read_from_file = .false.
    end if
  end subroutine read_evaluation_gauges

  !> \brief Read 'inflow_gauges' namelist content.
  subroutine read_inflow_gauges(self, file)
    implicit none
    class(nml_inflow_gauges_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: nInflowGaugesTotal !< Number of evaluation gauges for all domains
    integer(i4), dimension(maxNoDomains) :: NoInflowGauges_domain !< number of gauges for subdomain (1)
    integer(i4), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_id !< id of inflow gauge(1) for subdomain(1) --> (1,1)
    !> name of file with timeseries of inflow gauge(1) for subdomain(1) --> (1,1)
    character(256), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_filename
    !> consider flows from upstream/headwater cells of inflow gauge(1) for subdomain(1) --> (1,1)
    logical, dimension(maxNoDomains, maxNoGauges) :: InflowGauge_Headwater

    namelist /inflow_gauges/ &
      nInflowGaugesTotal, &
      NoInflowGauges_domain, &
      InflowGauge_id, &
      InflowGauge_filename, &
      InflowGauge_Headwater

    if ( self%read_from_file ) then
      nInflowGaugesTotal = 0
      NoInflowGauges_domain = 0
      InflowGauge_id = nodata_i4
      InflowGauge_filename = num2str(nodata_i4)
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=inflow_gauges)
      call close_nml(unit)
      self%nInflowGaugesTotal = nInflowGaugesTotal
      self%NoInflowGauges_domain = NoInflowGauges_domain
      self%InflowGauge_id = InflowGauge_id
      self%InflowGauge_filename = InflowGauge_filename
      self%InflowGauge_Headwater = InflowGauge_Headwater
      self%read_from_file = .false.
    end if
  end subroutine read_inflow_gauges

  !> \brief Read 'mrm_outputs' namelist content.
  subroutine read_mrm_outputs(self, file)
    use mo_message, only : message
    implicit none
    class(nml_mrm_outputs_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    integer(i4) :: output_deflate_level_mrm !< netcdf deflate level
    logical :: output_double_precision_mrm !< switch to enable double precision in netcdf
    integer(i4) :: output_time_reference_mrm !< time reference point location in output nc files
    integer(i4) :: timeStep_model_outputs_mrm !< timestep for writing model outputs
    logical, dimension(mrm_nOutFlxState) :: outputFlxState_mrm !< Define model outputs see "mhm_outputs.nml"

    logical :: file_exists

    namelist /nloutputresults/ &
      output_deflate_level_mrm, &
      output_double_precision_mrm, &
      output_time_reference_mrm, &
      timeStep_model_outputs_mrm, &
      outputFlxState_mrm

    if ( self%read_from_file ) then
      output_deflate_level_mrm = 6
      output_double_precision_mrm = .true.
      output_time_reference_mrm = 0
      outputFlxState_mrm = .FALSE.
      timeStep_model_outputs_mrm = -2
      inquire(file = file, exist = file_exists)
      if (file_exists) then
        call open_new_nml(file, unit)
        call position_nml(self%name, unit)
        read(unit, nml=nloutputresults)
        call close_nml(unit)
      else
        call message('***Warning: No file specifying mRM output fluxes exists')
      end if
      self%output_deflate_level_mrm = output_deflate_level_mrm
      self%output_double_precision_mrm = output_double_precision_mrm
      self%output_time_reference_mrm = output_time_reference_mrm
      self%timeStep_model_outputs_mrm = timeStep_model_outputs_mrm
      self%outputFlxState_mrm = outputFlxState_mrm
      self%read_from_file = .false.
    end if
  end subroutine read_mrm_outputs

  !> \brief Read 'routing1' namelist content.
  subroutine read_routing1(self, file)
    implicit none
    class(nml_routing1_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: muskingumTravelTime_constant !< muskingum parameter constant
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength !< muskingum parameter river length
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope !< muskingum parameter river slope
    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious !< muskingum parameter impervious
    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope !< muskingum parameter attenuation river slope

    namelist /routing1/ &
      muskingumTravelTime_constant, &
      muskingumTravelTime_riverLength, &
      muskingumTravelTime_riverSlope, &
      muskingumTravelTime_impervious, &
      muskingumAttenuation_riverSlope

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=routing1)
      call close_nml(unit)
      self%muskingumTravelTime_constant = muskingumTravelTime_constant
      self%muskingumTravelTime_riverLength = muskingumTravelTime_riverLength
      self%muskingumTravelTime_riverSlope = muskingumTravelTime_riverSlope
      self%muskingumTravelTime_impervious = muskingumTravelTime_impervious
      self%muskingumAttenuation_riverSlope = muskingumAttenuation_riverSlope
      self%read_from_file = .false.
    end if
  end subroutine read_routing1

  !> \brief Read 'routing2' namelist content.
  subroutine read_routing2(self, file)
    implicit none
    class(nml_routing2_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: streamflow_celerity !< streamflow celerity

    namelist /routing2/ &
      streamflow_celerity
    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=routing2)
      call close_nml(unit)
      self%streamflow_celerity = streamflow_celerity
      self%read_from_file = .false.
    end if
  end subroutine read_routing2

  !> \brief Read 'routing3' namelist content.
  subroutine read_routing3(self, file)
    implicit none
    class(nml_routing3_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp), dimension(nColPars) :: slope_factor !< slope factor

    namelist /routing3/ &
      slope_factor

    if ( self%read_from_file ) then
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=routing3)
      call close_nml(unit)
      self%slope_factor = slope_factor
      self%read_from_file = .false.
    end if
  end subroutine read_routing3

  !> \brief Read 'config_riv_temp' namelist content.
  subroutine read_config_riv_temp(self, file)
    implicit none
    class(nml_config_riv_temp_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit !< file unit to open the given file
    real(dp) :: albedo_water !< albedo of open water
    real(dp) :: pt_a_water !< priestley taylor alpha parameter for PET on open water
    real(dp) :: emissivity_water !< emissivity of water
    real(dp) :: turb_heat_ex_coeff !< lateral heat exchange coefficient water <-> air
    integer(i4) :: max_iter !< maximum number of iterations
    real(dp) :: delta_iter !< convergence delta
    real(dp) :: step_iter !< step size for iterative solver
    character(256) :: riv_widths_file !< file name for river widths
    character(256) :: riv_widths_name !< variable name for river widths
    character(256), dimension(maxNoDomains) :: dir_riv_widths !< files for river widths

    namelist /config_riv_temp/ &
      albedo_water, &
      pt_a_water, &
      emissivity_water, &
      turb_heat_ex_coeff, &
      max_iter, &
      delta_iter, &
      step_iter, &
      riv_widths_file, &
      riv_widths_name, &
      dir_riv_widths

    if ( self%read_from_file ) then
      albedo_water = 0.15_dp
      pt_a_water = 1.26_dp
      emissivity_water = 0.96_dp
      turb_heat_ex_coeff = 20.0_dp
      max_iter = 20_i4
      delta_iter = 1.0e-02_dp
      step_iter = 5.0_dp
      call open_new_nml(file, unit)
      call position_nml(self%name, unit)
      read(unit, nml=config_riv_temp)
      call close_nml(unit)
      self%albedo_water = albedo_water
      self%pt_a_water = pt_a_water
      self%emissivity_water = emissivity_water
      self%turb_heat_ex_coeff = turb_heat_ex_coeff
      self%max_iter = max_iter
      self%delta_iter = delta_iter
      self%step_iter = step_iter
      self%riv_widths_file = riv_widths_file
      self%riv_widths_name = riv_widths_name
      self%dir_riv_widths = dir_riv_widths
      self%read_from_file = .false.
    end if
  end subroutine read_config_riv_temp

  !> \brief Read 'coupling' namelist content.
  subroutine read_coupling(self, file)
    implicit none
    class(nml_coupling_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelist

    integer :: unit, status
    integer(i4) :: case !< coupling case
    integer(i4) :: meteo_timestep !< timestep for meteo-data from coupling
    logical :: meteo_time_ref_endpoint !< expect meteo has time reference point at end of associated time interval
    logical :: meteo_expect_pre !< expect meteo from coupling: [mm]      Precipitation
    logical :: meteo_expect_temp !< expect meteo from coupling: [degC]    Air temperature
    logical :: meteo_expect_pet !< expect meteo from coupling: [mm TS-1] Potential evapotranspiration
    logical :: meteo_expect_tmin !< expect meteo from coupling: [degC]    minimum daily air temperature
    logical :: meteo_expect_tmax !< expect meteo from coupling: [degC]    maximum daily air temperature
    logical :: meteo_expect_netrad !< expect meteo from coupling: [W m2]    net radiation
    logical :: meteo_expect_absvappress !< expect meteo from coupling: [Pa]      absolute vapour pressure
    logical :: meteo_expect_windspeed !< expect meteo from coupling: [m s-1]   windspeed
    logical :: meteo_expect_ssrd !< expect meteo from coupling: [W m2]    short wave radiation
    logical :: meteo_expect_strd !< expect meteo from coupling: [W m2]    long wave radiation
    logical :: meteo_expect_tann !< expect meteo from coupling: [degC]    annual mean air temperature

    namelist /coupling/ &
      case, &
      meteo_timestep, &
      meteo_time_ref_endpoint, &
      meteo_expect_pre, &
      meteo_expect_temp, &
      meteo_expect_pet, &
      meteo_expect_tmin, &
      meteo_expect_tmax, &
      meteo_expect_netrad, &
      meteo_expect_absvappress, &
      meteo_expect_windspeed, &
      meteo_expect_ssrd, &
      meteo_expect_strd, &
      meteo_expect_tann

    if ( self%read_from_file ) then
      case = 0_i4 ! no coupling by default
      meteo_timestep = 0_i4 ! only valid if no meteo expected
      meteo_time_ref_endpoint = .false. ! meteo data usually given at begin of time interval (i.e. 00:00 for current day)
      meteo_expect_pre = .false.
      meteo_expect_temp = .false.
      meteo_expect_pet = .false.
      meteo_expect_tmin = .false.
      meteo_expect_tmax = .false.
      meteo_expect_netrad = .false.
      meteo_expect_absvappress = .false.
      meteo_expect_windspeed = .false.
      meteo_expect_ssrd = .false.
      meteo_expect_strd = .false.
      meteo_expect_tann = .false.
      call open_new_nml(file, unit)
      call position_nml(self%name, unit, status=status)
      if (status == 0) read(unit, nml=coupling)
      call close_nml(unit)
      self%case = case
      self%meteo_timestep = meteo_timestep
      self%meteo_time_ref_endpoint = meteo_time_ref_endpoint
      self%meteo_expect_pre = meteo_expect_pre
      self%meteo_expect_temp = meteo_expect_temp
      self%meteo_expect_pet = meteo_expect_pet
      self%meteo_expect_tmin = meteo_expect_tmin
      self%meteo_expect_tmax = meteo_expect_tmax
      self%meteo_expect_netrad = meteo_expect_netrad
      self%meteo_expect_absvappress = meteo_expect_absvappress
      self%meteo_expect_windspeed = meteo_expect_windspeed
      self%meteo_expect_ssrd = meteo_expect_ssrd
      self%meteo_expect_strd = meteo_expect_strd
      self%meteo_expect_tann = meteo_expect_tann
    end if
  end subroutine read_coupling

end module mo_namelists
