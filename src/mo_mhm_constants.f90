!> \file mo_mhm_constants.f90

!> \brief Provides mHM specific constants

!> \details Provides mHM specific constants such as flood plain elevation.

!> \author Matthias Cuntz
!> \date Nov 2011

MODULE mo_mhm_constants

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  real(dp),    public, parameter :: deltaH             = 5.000_dp  ! [m]     flood plain elevation, transept, above riverbed
  integer(i4), public, parameter :: nLCover_class      = 3_i4      ! [-]     Number of land cover class


  integer(i4), public, parameter :: nodata_i4          = -9999_i4  ! [-]     global no data value
  real(dp),    public, parameter :: nodata_dp          = -9999._dp ! [-]     global no data value

  integer(i4), public, parameter :: nRoutingStates     = 2         ! Dimension of the auxiliary vectors
  !                                                                ! which store current and past states of
  !                                                                ! incoming and outgoing of discharge at
  !                                                                ! a given node 
  !                                                                ! (1 - past)
  !                                                                ! (2 - current)
  integer(i4), public, parameter :: maxNoGauges        = 100_i4    ! maximal number of gauges allowed
  integer(i4), public, parameter :: nColPars           = 5_i4      ! number of properties of the global variables
  integer(i4), public, parameter :: maxNoSoilHorizons  = 10_i4     ! maximum number of allowed soil layers
  integer(i4), public, parameter :: maxNoBasins        = 50_i4     ! maximum number of allowed basins
  integer(i4), public, parameter :: maxNLcovers        = 50_i4     ! maximum number of allowed LCover scenes                 
  integer(i4), public, parameter :: maxGeoUnit         = 25_i4     ! maximum number of allowed geological classes

  ! default inital values for states and fluxes as well as parameter fields
  real(dp),    public, parameter :: P1_InitStateFluxes =    0.00_dp
  real(dp),    public, parameter :: P2_InitStateFluxes =   15.00_dp
  real(dp),    public, parameter :: P3_InitStateFluxes =   10.00_dp
  real(dp),    public, parameter :: P4_InitStateFluxes =   75.00_dp
  real(dp),    public, parameter :: P5_InitStateFluxes = 1500.00_dp
  real(dp),    public, parameter :: C1_InitStateSM     =    0.25_dp

  ! maximum number of outputs (fluxes states) for mHM
  integer(i4), public, parameter :: nOutFlxState       = 17_i4     ! max. number of outputs to write into a netcdf file

   ! Time constants
  real(dp),    public, parameter :: DayHours           =     24.0_dp  ! hours per day
  real(dp),    public, parameter :: HourSecs           =   3600.0_dp  ! seconds per hour
  real(dp),    public, parameter :: YearMonths         =     12.0_dp  ! months per year
  integer(i4), public, parameter :: YearMonths_i4      =     12       ! months per year
  real(dp),    public, parameter :: YearDays           =    365.0_dp  ! days in a year
  real(dp),    public, parameter :: DaySecs            =  86400.0_dp  ! sec in a day
 
  ! constant connected to soil paramterization (mo_mpr_soilmoist)
  ! organic matter constant for calculation of mineral bulk density following RAWL  
  real(dp),   public, parameter :: BulkDens_OrgMatter  = 0.224_dp     ! [g/cm3] from W.R. RAWLS
  ! constants for determinination of the field capacity following Twarakavi
  real(dp),   public, parameter :: field_cap_c1        = -0.60_dp     ! field capacity constant 1
  real(dp),   public, parameter :: field_cap_c2        =  2.0_dp      ! field capacity constant 2
  ! constants for determinination of the van Genuchten parameter n and sand treshold
  real(dp),   public, parameter :: vGenuchten_sandtresh = 66.5_dp     ! van Genuchten snad treshold
  real(dp),   public, parameter :: vGenuchtenN_c1       =  1.392_dp   ! constants for van Genuchten n 
  real(dp),   public, parameter :: vGenuchtenN_c2       =  0.418_dp
  real(dp),   public, parameter :: vGenuchtenN_c3       = -0.024_dp
  real(dp),   public, parameter :: vGenuchtenN_c4       =  1.212_dp
  real(dp),   public, parameter :: vGenuchtenN_c5       = -0.704_dp
  real(dp),   public, parameter :: vGenuchtenN_c6       = -0.648_dp
  real(dp),   public, parameter :: vGenuchtenN_c7       =  0.023_dp
  real(dp),   public, parameter :: vGenuchtenN_c8       =  0.044_dp
  real(dp),   public, parameter :: vGenuchtenN_c9       =  3.168_dp
  real(dp),   public, parameter :: vGenuchtenN_c10      = -2.562_dp
  real(dp),   public, parameter :: vGenuchtenN_c11      =  7.0E-9_dp
  real(dp),   public, parameter :: vGenuchtenN_c12      =  4.004_dp
  real(dp),   public, parameter :: vGenuchtenN_c13      =  3.750_dp
  real(dp),   public, parameter :: vGenuchtenN_c14      = -0.016_dp
  real(dp),   public, parameter :: vGenuchtenN_c15      = -4.197_dp
  real(dp),   public, parameter :: vGenuchtenN_c16      =  0.013_dp
  real(dp),   public, parameter :: vGenuchtenN_c17      =  0.076_dp
  real(dp),   public, parameter :: vGenuchtenN_c18      =  0.276_dp
  ! constants for determinination Ks
  real(dp),   public, parameter :: Ks_c                 = 10.0_dp
  ! constant for permanent wiltung point (PWP)
  real(dp),   public, parameter :: PWP_c                = 1.0_dp
  real(dp),   public, parameter :: PWP_matPot_ThetaR    = 15000.0_dp ! [hPa] matrix potential of -1500 kPa, assumed as thetaR=0


  !> Stefan-Boltzmann constant [W m^-2 K^-4] 
  real(dp),   public, parameter  :: StBoltzmann          = 5.67e-08_dp                
  !> Constant for Hargreaves-Samani ref. ET formula [deg C]
  real(dp),   public, parameter  :: HarSamConst          = 17.800_dp             
  !> assumed meteorol. measurement hight for estimation of aeroResist and surfResist
  real(dp),    public, parameter :: WindMeasHeight       = 10.0_dp 
  !> von karman constant 
  real(dp),    public, parameter :: karman               = 0.41_dp

  !> LAI factor for bulk surface resistance formulation
  real(dp),    public, parameter :: LAI_factor_surfResi  = 0.3_dp
  !> LAI offset for bulk surface resistance formulation
  real(dp),    public, parameter :: LAI_offset_surfResi  = 1.2_dp
  !> maximum bulk surface resistance
  real(dp),    public, parameter :: max_surfResist       = 250.0_dp


  ! constants in the Duffie formula for computing extraterrestrial radiation
  real(dp),    public, parameter :: DuffieDr             =   0.0330_dp
  real(dp),    public, parameter :: DuffieDelta1         =   0.4090_dp
  real(dp),    public, parameter :: DuffieDelta2         =   1.3900_dp

  !> constants for Tetens's formula to calculate saturated vapour pressure
  real(dp),   public, parameter  :: tetens_c1            = 0.6108_dp
  real(dp),   public, parameter  :: tetens_c2            = 17.270_dp
  real(dp),   public, parameter  :: tetens_c3            = 237.30_dp
  !> constant for calculating the slope of the saturation vapour pressure curve following Tetens
  real(dp),   public, parameter  :: satpressureslope1    = 4098.0_dp

END MODULE mo_mhm_constants
