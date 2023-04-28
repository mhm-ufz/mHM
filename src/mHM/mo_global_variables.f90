!> \dir src
!> \brief Source code of mHM.
!> \details All Fortran source files for mHM, mRM and MPR.

!> \dir mHM
!> \brief \copybrief f_mhm
!> \details \copydetails f_mhm

!> \defgroup   f_mhm mHM - Fortran modules
!> \brief      Core modules of mHM.
!> \details    These modules provide the core components of mHM.

!> \file mo_global_variables.f90
!> \brief \copybrief mo_global_variables
!> \details \copydetails mo_global_variables

!> \brief Main global variables for mHM.
!> \details Global variables ONLY used in reading, writing and startup.
!> \changelog
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!! - Luis Samaniego,     Feb 2013
!!   - new variable names, new modules, units
!! - Rohini Kumar,       Jul 2013
!!   - fraction of perfectly sealed area within city added
!! - Rohini Kumar,       Aug 2013
!!   - name changed "inputFormat" to "inputFormat_meteo_forcings"
!! - Rohini Kumar,       Aug 2013
!!   - name changed from "L0_LAI" to "L0_LCover_LAI"
!! - Rohini Kumar,       Aug 2013
!!   - added dirSoil_LUT and dirGeology_LUT
!! - Luis Samaniego,     Nov 2013
!!   - documentation of dimensions
!! - Matthias Zink,      Nov 2013
!!   - added "InflowGauge" and inflow gauge variabels in Domain
!! - Rohini Kumar,       May 2014
!!   - added options for the model run cordinate system
!! - Stephan Thober,     Jun 2014
!!   - added timeStep_model_inputs and readPer
!! - Stephan Thober,     Jun 2014
!!   - added perform_mpr, updated restart flags
!! - Cuntz M. & Mai J.,  Nov 2014
!!   - LAI input from daily, monthly or yearly files
!! - Matthias Zink,      Dec 2014
!!   - adopted inflow gauges to ignore headwater cells
!! - Matthias Zink,      Mar 2015
!!   - added optional soil mositure readin: dirSoil_moisture, L1_sm
!! - Stephan Thober,     Aug 2015
!!   - moved routing related variables to mRM
!! - Oldrich Rakovec,    Oct 2015
!!   - added definition of Domain averaged TWS data
!! - Rohini Kumar,       Mar 2016
!!   - new variables for handling different soil databases
!! - Johann Brenner,     Feb 2017
!!   - added optional evapotranspiration readin: dirEvapotranspiration, L1_et
!! - Zink M. Demirel C., Mar 2017
!!   - added Jarvis soil water stress variable for SM process(3)
!! - Demirel M.C.        May 2017
!!   - added L1_petLAIcorFactor for PET correction
!! - O. Rakovec, R.Kumar Nov 2017
!!   - added project description for the netcdf outputs
!! - Robert Schweppe,    Dec 2017
!!   - expanded dimensions of effective parameters
!! - Robert Schweppe,    Dec 2017
!!   - merged duplicated variables with mrm into common variables
!> \authors Luis Samaniego
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_global_variables

  USE mo_kind, ONLY : i4, dp
  USE mo_constants, ONLY : YearMonths
  USE mo_mhm_constants, ONLY : nOutFlxState
  USE mo_optimization_types, ONLY : optidata
  use mo_meteo_handler, only : meteo_handler_type
  use mo_coupling_type, only : couple_cfg_type

  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! COUPLING CONFIG
  ! -------------------------------------------------------------------
  type(couple_cfg_type), public :: couple_cfg !< coupling configuration class

  ! -------------------------------------------------------------------
  ! METEO HANDLER
  ! -------------------------------------------------------------------
  type(meteo_handler_type), public :: meteo_handler !< the meteo handler class

  ! -------------------------------------------------------------------
  ! DEFINE OUTPUTS
  ! -------------------------------------------------------------------
  integer(i4) :: output_deflate_level   !< deflate level in nc files
  integer(i4) :: output_time_reference   !< time reference point location in output nc files
  logical :: output_double_precision    !< output precision in nc files
  integer(i4) :: timeStep_model_outputs !< timestep for writing model outputs
  logical, dimension(nOutFlxState) :: outputFlxState         !< Define model outputs see "mhm_outputs.nml"
                                                             !< dim1 = number of output variables to be written

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
  ! disaggregated meteo forcings
  real(dp), public, dimension(:), allocatable :: L1_pet_calc     !< [mm TS-1] estimated/corrected potential evapotranspiration
  real(dp), public, dimension(:), allocatable :: L1_temp_calc    !< [degC] temperature for current time step
  real(dp), public, dimension(:), allocatable :: L1_prec_calc    !< [mm TS-1] precipitation for current time step
  ! dim2 = number model soil horizons
  ! states and fluxes
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

  ! -------------------------------------------------------------------
  ! AUXILIARY VARIABLES
  ! -------------------------------------------------------------------
  !

  real(dp), public, dimension(:), allocatable :: neutron_integral_AFast !< pre-calculated integrand for
  ! vertical projection of isotropic neutron flux

END MODULE mo_global_variables
