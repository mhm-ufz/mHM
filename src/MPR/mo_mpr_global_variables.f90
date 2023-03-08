!> \dir MPR
!> \brief \copybrief f_mpr
!> \details \copydetails f_mpr

!> \defgroup   f_mpr MPR - Fortran modules
!> \brief      Core modules of MPR.
!> \details    These modules provide the core components of the Multiscale Parameter Regionalization scheme of mHM.

!> \file mo_mpr_global_variables.f90
!> \brief \copybrief mo_mpr_global_variables
!> \details \copydetails mo_mpr_global_variables

!> \brief Global variables for mpr only
!> \details Global variables used to run MPR for mHM.
!> \authors Robert Schweppe
!> \date Dec 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_global_variables

  use mo_kind, only : i4, dp
  use mo_common_types, only: period

  implicit none

  private

  ! ------------------------------------------------------------------
  ! SOIL DATA
  ! ------------------------------------------------------------------
  real(dp), public :: tillageDepth       ! [mm]  Soil depth down to which organic
  !                                                                               matter is possible
  integer(i4), public :: nSoilTypes         !       Number of soil types
  integer(i4), public :: iFlag_soilDB       ! options to handle different soil databases
  integer(i4), public :: nSoilHorizons_mHM  !       Number of horizons to model
  real(dp), dimension(:), allocatable, public :: HorizonDepth_mHM   ! [mm]  Horizon depth from surface,
  !                                                                               positive downwards

  type soilType
    ! dim1 =  nSoilType (e.g. i=1..72 for BUEK)
    ! dim2 =  the maximum of nHorizons
    ! dim3 =  land cover classes
    ! input data
    integer(i4), dimension(:), allocatable :: id                 !            Soil Id
    integer(i4), dimension(:), allocatable :: nHorizons          !            Number of horizons
    integer(i4), dimension(:), allocatable :: is_present         !            Wether this soil type is present in
    !                                                                !            this domain or not
    real(dp), dimension(:, :), allocatable :: UD                 ! [mm]       Upper Bound of depth
    real(dp), dimension(:, :), allocatable :: LD                 ! [mm]       Lower Bound of depth
    real(dp), dimension(:, :), allocatable :: clay               ! [%]        Clay content
    real(dp), dimension(:, :), allocatable :: sand               ! [%]        Sand content
    real(dp), dimension(:, :), allocatable :: DbM                ! [g/cm2]    Mineral Bulk density
    real(dp), dimension(:, :), allocatable :: depth              ! [mm]       Depth of the soil Horizon
    real(dp), dimension(:), allocatable :: RZdepth            ! [mm]       Total soil depth
    real(dp), dimension(:, :, :), allocatable :: Wd                 ! [1]        Weights of mHM Horizons according to
    !                                                                !            horizons provided in soil database
    integer(i4), dimension(:), allocatable :: nTillHorizons      ! [1]        Number of tillage horizons

    ! derived soil hydraulic properties
    real(dp), dimension(:, :, :), allocatable :: thetaS_Till        ! [1]        Saturated water content of soil horizons
    !                                                                !            tillage depth - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaS             ! [1]        Saturated water content of soil horizons
    !                                                                !            after tillage depth
    real(dp), dimension(:, :, :), allocatable :: Db                 ! [g/cm2]    Bulk density, LUC dependent
    !                                                                !            = f( OM, management)
    real(dp), dimension(:, :, :), allocatable :: thetaFC_Till       ! [1]        Field capacity of tillage layers;
    !                                                                !            LUC dependent - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaFC            ! [1]        Field capacity of deeper layers
    real(dp), dimension(:, :, :), allocatable :: thetaPW_Till       ! [1]        Permament wilting point of tillage layers;
    !                                                                !            LUC dependent - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaPW            ! [1]        Permanent wilting point of deeper layers
    real(dp), dimension(:, :, :), allocatable :: Ks                 ! [cm/d]     Saturated hydaulic conductivity
  end type soilType
  type(soilType), public :: soilDB             !            The soil database

  ! -----------------------------------------------------------------
  ! GEOLOGICAL FORMATION data
  ! -----------------------------------------------------------------
  integer(i4), public :: nGeoUnits   ! Number of geological formations
  integer(i4), dimension(:), allocatable, public :: GeoUnitList ! List of ids of each geological formations
  integer(i4), dimension(:), allocatable, public :: GeoUnitKar  ! Id of Karstic formation (0 == does not exist)

  ! -----------------------------------------------------------------
  ! Land cover, LAI LUT data
  ! -----------------------------------------------------------------
  character(256), public :: inputFormat_gridded_LAI    ! format of gridded LAI data (nc only)
  integer(i4), public :: timeStep_LAI_input         ! time step of gridded LAI input
  ! LAI data
  ! variables used when timeStep_LAI_input == 0
  integer(i4), public :: nLAIclass         ! Number of LAI classes
  integer(i4), public :: nLAI              ! Number of LAI slices (a.k.a timestep)
  real(dp), dimension(:), allocatable, public :: LAIBoundaries        !
  integer(i4), public, dimension(:), allocatable :: LAIUnitList       ! List of ids of each LAI class in LAILUT
  real(dp), public, dimension(:, :), allocatable :: LAILUT            ! [m2/m2] Leaf area index for LAIUnit
  !                                                                        ! dim1=land cover class, dim2=month of year
  type(period), dimension(:), allocatable, public :: LAIPer            ! time period for LAI_readin
  real(dp), public :: fracSealed_cityArea ! fraction of area within city assumed to be
  !                                                                          ! perfectly sealed [0-1]

  ! -------------------------------------------------------------------
  ! L0 DOMAIN description -> <only domain>
  ! -------------------------------------------------------------------
  ! mHM derived variables
  ! dim1 = number grid cells L0
  real(dp), public, dimension(:), allocatable :: L0_slope_emp         ! Empirical quantiles of slope
  !
  real(dp), public, dimension(:, :), allocatable :: L0_gridded_LAI       ! gridded LAI data used when timeStep_LAI_input<0 or==1
  !                                                                          ! dim1=number of gridcells, dim2=number LAI timesteps

  real(dp), public, dimension(:), allocatable :: L0_slope    ! [%]      Slope
  real(dp), public, dimension(:), allocatable :: L0_asp     ! [degree]  Aspect degree
  !  [dim1=number grid cells, dim2=Number of soil horizons] note: for iFlag_soilDB=0, dim2=1
  integer(i4), public, dimension(:, :), allocatable :: L0_soilId  !           soil id (iFlag_soilDB = 0)
  integer(i4), public, dimension(:), allocatable :: L0_geoUnit  !      Geologic formation (unit)

  ! ------------------------------------------------------------------
  ! DIRECTORIES
  ! ------------------------------------------------------------------
  ! has the dimension of nDomains
  character(256), dimension(:), allocatable, public :: dirgridded_LAI     ! Directory where gridded LAI is located
  ! used when timeStep_LAI_input < 0

  ! Effective parameters
  ! dim1 = number grid cells L1
  ! dim2 = number model soil horizons or YearMonths or other auxiliary dimension
  ! dim3 = number of LCscenes
  real(dp), public, dimension(:, :, :), allocatable :: L1_fSealed       ! [1]  Fraction of sealed area (nCells, nLCscenes)

  real(dp), public, dimension(:, :, :), allocatable :: L1_alpha               ! [1]            Exponent for the upper reservoir
  real(dp), public, dimension(:, :, :), allocatable :: L1_degDayInc           ! [d-1 degC-1]   Increase of the Degree-day factor
  !                                                                           !                per mm of increase in precipitation
  real(dp), public, dimension(:, :, :), allocatable :: L1_degDayMax           ! [mm-1 degC-1]  Maximum Degree-day factor
  real(dp), public, dimension(:, :, :), allocatable :: L1_degDayNoPre         ! [mm-1 degC-1]  Degree-day factor with no
                                                                              ! precipitation.
  real(dp), public, dimension(:, :, :), allocatable :: L1_degDay              ! [mm d-1degC-1] Degree-day factor.
  real(dp), public, dimension(:, :, :), allocatable :: L1_karstLoss           ! [1]    Karstic percolation loss
  real(dp), public, dimension(:, :, :), allocatable :: L1_fAsp                ! [1]    PET correction for aspect
  real(dp), public, dimension(:, :, :), allocatable :: L1_petLAIcorFactor     ! [-]   PET correction based on LAI (KC by GEUS.dk)

  real(dp), public, dimension(:, :, :), allocatable :: L1_HarSamCoeff         ! [1]    Hargreaves Samani coeffiecient
  real(dp), public, dimension(:, :, :), allocatable :: L1_PrieTayAlpha        ! [1]    Priestley Taylor coeffiecient
  real(dp), public, dimension(:, :, :), allocatable :: L1_aeroResist          ! [s m-1] aerodynamical resitance
  real(dp), public, dimension(:, :, :), allocatable :: L1_surfResist          ! [s m-1] bulk surface resitance
  real(dp), public, dimension(:, :, :), allocatable :: L1_fRoots              ! [1]    Fraction of roots in soil horizons
  real(dp), public, dimension(:, :, :), allocatable :: L1_maxInter            ! [mm]   Maximum interception

  real(dp), public, dimension(:, :, :), allocatable :: L1_kfastFlow           ! [d-1]  Fast interflow recession coefficient
  real(dp), public, dimension(:, :, :), allocatable :: L1_kSlowFlow           ! [d-1]  Slow interflow recession coefficient
  real(dp), public, dimension(:, :, :), allocatable :: L1_kBaseFlow           ! [d-1]  Baseflow recession coefficient
  real(dp), public, dimension(:, :, :), allocatable :: L1_kPerco              ! [d-1]  percolation coefficient
  real(dp), public, dimension(:, :, :), allocatable :: L1_soilMoistFC         ! [mm]   Soil moisture below which actual ET
  !                                                                           !        is reduced linearly till PWP
  real(dp), public, dimension(:, :, :), allocatable :: L1_soilMoistSat        ! [mm]   Saturation soil moisture for each horizon [mm]
  real(dp), public, dimension(:, :, :), allocatable :: L1_soilMoistExp        ! [1]    Exponential parameter to how non-linear
  !                                                                           !        is the soil water retention
  real(dp), public, dimension(:, :, :), allocatable :: L1_jarvis_thresh_c1    ![1] jarvis critical value for normalized soil
  !                                                                           !        water content
  real(dp), public, dimension(:, :, :), allocatable :: L1_tempThresh          ! [degC]   Threshold temperature for snow/rain
  real(dp), public, dimension(:, :, :), allocatable :: L1_unsatThresh         ! [mm]  Threshold waterdepth controlling fast interflow
  real(dp), public, dimension(:, :, :), allocatable :: L1_sealedThresh        ! [mm]  Threshold waterdepth for surface runoff
  !                                                                           !       in sealed surfaces
  real(dp), public, dimension(:, :, :), allocatable :: L1_wiltingPoint        ! [mm]  Permanent wilting point: below which neither
  !                                                                           !       plant can take water nor water can drain in
  ! >> COSMIC neutron count realated parameters -- only those which are regionlized
  !!   defined here others are treated as global parameters...
  real(dp), public, dimension(:,:,:), allocatable :: L1_No_Count     !   N0 COUNT      >> in Desilets and COSMIC routines
  real(dp), public, dimension(:,:,:), allocatable :: L1_bulkDens     !   Bulk density  >> in COSMIC routines
  real(dp), public, dimension(:,:,:), allocatable :: L1_latticeWater !   lattice water >> in COSMIC routines
  real(dp), public, dimension(:,:,:), allocatable :: L1_COSMICL3     !   !COSMIC L3    >> in COSMIC routines

end module mo_mpr_global_variables
