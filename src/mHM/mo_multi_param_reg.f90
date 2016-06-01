!> \file mo_multi_param_reg.f90

!> \brief   Multiscale parameter regionalization (MPR).

!> \details This module provides the routines for multiscale parameter regionalization (MPR).

!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012

!          created  Sa   16.02.2006
!           update  Sa   17.09.2007 betas new number
!           update  Sa   03.10.2007 new name, land cover state  
!           update  Ku   25.03.2008 all parameters are regionalised  
!           update  Ku   04.10.2010 vector version
!           update  Th   20.12.2012 modular version
!           update  MZ   27.11.2014 added parameterization of PET
!                   Ku   Mar 2016 - changes for handling multiple soil database options

MODULE mo_multi_param_reg

  use mo_kind, only: i4, dp

  implicit none

  private

  PUBLIC :: mpr                     ! calculates effective regionalised parameters
  PUBLIC :: canopy_intercept_param  ! estimate effective max. canopy interception

contains
  ! ---------------------------------------------------------------------------

  !      NAME
  !         mpr

  !>        \brief Regionalizing and Upscaling process parameters

  !>        \details calculating process parameters at L0 scale (Regionalization), like:\n
  !>                 - Baseflow recession parameter\n
  !>                 - Soil moisture parameters\n
  !>                 - PET correction for aspect\n
  !>
  !>                 and upscale these parameters to retrieve effective parameters at scale
  !>                 L1. \n 
  !>                 Further parameter regionalizations are done for:\n
  !>                 - snow accumulation and melting parameters\n
  !>                 - threshold parameter for runoff generation on impervious layer\n
  !>                 - karstic percolation loss\n
  !>                 - setting up the Regionalized Routing Parameters\n
  !>                 \n 

  !      INTENT(IN)
  !>       \param[in] "integer(i4) :: tt"       - simulation time step
  !>       \param[in] "integer(i4) :: newYear"  -  new year
  !>       \param[in] "integer(i4) :: LCyearId" -  mapping of landcover scenes
  !>       \param[in] "integer(i4) :: proc_flag(:,:)" - indicates which process shall be run
  !>                                            the shape is number of processes times 3, the
  !>                                            first column indicates which kind of process
  !>                                            is run, 0 indicating do not run, the second
  !>                                            column indicates the number of parameters
  !>                                            required for this process and the third column
  !>                                            indicates the position of the first parameter
  !>                                            for this process in the array param
  !>       \param[in] "real(dp)    :: param(:)"         - given global parameter array
  !>       \param[in] "integer(i4) :: iFlag_soil"       - flags for handling multiple soil databases
  !>       \param[in] "real(dp)    :: nodata"           - given nodata value
  !>       \param[in] "integer(i4) :: geoUnit0(:,:)"    - geological units at Level 0
  !>       \param[in] "integer(i4) :: geo_unit_list(:)" - index list of geological units
  !>       \param[in] "real(dp)    :: LAILUT(:,:)"      - [1] Leaf area index for LAIUnit
  !>       \param[in] "integer(i4) :: LAIUnitList(:)"   - [1] List of ids of each LAI class in LAILUT
  !>       \param[in] "integer(i4) :: is_present(:)"    - indicates whether soiltype exists
  !>       \param[in] "integer(i4) :: nHorizons(:)"     - Number of Horizons per soiltype
  !>       \param[in] "integer(i4) :: nTillHorizons(:)" - Number of Tillage Horizons
  !>       \param[in] "real(dp)    :: sand(:,:)"        - sand content
  !>       \param[in] "real(dp)    :: clay(:,:)"        - clay content
  !>       \param[in] "real(dp)    :: DbM(:,:)"         - mineral Bulk density
  !>       \param[in] "real(dp)    :: Wd(:,:,:)"        - weights of mHM horizons
  !>       \param[in] "real(dp)    :: RZdepth(:)"       - [mm] Total soil depth
  !>       \param[in] "integer(i4) :: nHorizons_mHM"    - Number of horizons in mHM
  !>       \param[in] "integer(i4) :: horizon_depth(:)" - depth of each horizon
  !>       \param[in] "real(dp)    :: c2TSTu"           - unit transformation coefficient
  !>       \param[in] "real(dp)    :: fForest1(:)"      - fraction of forest cover at scale L1
  !>       \param[in] "real(dp)    :: fIperm1(:)"       - fraction of sealed area at scale L1
  !>       \param[in] "real(dp)    :: fPerm1(:)"        - fraction of permeable area at scale L1
  !>       \param[in] "integer(i4) :: soilID0(:,:)"     - [1] soil IDs at level 0
  !>       \param[in] "real(dp)    :: Asp0(:,:)"        - [degree] Aspect at Level 0
  !>       \param[in] "real(dp)    :: LCover_LAI0(:)    - [1] land cover ID for LAI estimation
  !>       \param[in] "integer(i4) :: LCover0(:)"       - [1] land use cover at level 0 
  !>       \param[in] "real(dp)    :: length(:)"        - [m] total length
  !>       \param[in] "real(dp)    :: slope(:)"         - average slope
  !>       \param[in] "real(dp)    :: fFPimp(:)"        - fraction of the flood plain with
  !>                                                      impervious layer
  !>       \param[in] "real(dp)    :: TS"               - [h] time step in
  !>       \param[in] "integer(i4) :: cell_id0(:,:)"    - cell ids of high resolution field, 
  !>                                                      Number of rows times Number of columns of
  !>                                                      high resolution field
  !>       \param[in] "integer(i4) :: upp_row_L1(:)"    - Upper row id in high resolution field 
  !>                                                      (L0) of low resolution cell (L1 cell)
  !>       \param[in] "integer(i4) :: low_row_L1(:)"    - Lower row id in high resolution field 
  !>                                                      (L0) of low resolution cell (L1 cell)
  !>       \param[in] "integer(i4) :: lef_col_L1(:)"    - Left column id in high resolution field 
  !>                                                      (L0) of low resolution cell (L1 cell)
  !>       \param[in] "integer(i4) :: rig_col_L1(:)"    - Right column id in high resolution field 
  !>                                                      (L0) of low resolution cell (L1 cell)
  !>       \param[in] "integer(i4) :: nL0_in_L1(:)"     - Number of high resolution cells (L0) in 
  !>                                                      low resolution cell (L1 cell)

  !     INTENT(INOUT)
  !         None
  
  !      INTENT(OUT)
  !>       \param[in,out] "real(dp) :: k2_1(:,:)"          - baseflow recession parameter at L1
  !>       \param[in,out] "real(dp) :: KsVar_H0(:,:)"      - relative variability of saturated
  !>                                                         hydraulic cound. for Horizantal flow
  !>       \param[in,out] "real(dp) :: KsVar_V0(:,:)"      - relative variability of saturated
  !>                                                         hydraulic cound. for Vertical flow
  !>       \param[in,out] "real(dp) :: SMs_tot0(:,:)"      - total saturated soil moisture content
  !>       \param[in,out] "real(dp) :: SMs_FC0(:,:)"       - soil mositure deficit from
  !>                                                         field cap. w.r.t to saturation
  !>       \param[in,out] "real(dp) :: beta1(:,:)"         - Parameter that determines the
  !>                                                         relative contribution to SM, upscaled
  !>                                                         Bulk density. Number of cells at L1
  !>                                                         times number of horizons in mHM
  !>       \param[in,out] "real(dp) :: SMs1(:,:)"          - [10^-3 m] depth of saturated SM cont
  !>                                                         Number of cells at L1 times number
  !>                                                         of horizons in mHM
  !>       \param[in,out] "real(dp) :: FC1(:,:)"           - [10^-3 m] field capacity. Number
  !>                                                         of cells at L1 times number of horizons
  !>                                                         in mHM
  !>       \param[in,out] "real(dp) :: PW1(:,:)"           - [10^-3 m] permanent wilting point.
  !>                                                         Number of cells at L1 times number 
  !>                                                         of horizons in mHM
  !>       \param[in,out] "real(dp) :: fRoots1(:,:)"       - fraction of roots in soil horizons.
  !>                                                         Number of cells at L1 times number
  !>                                                         of horizons in mHM
  !>       \param[in,out] "real(dp) :: TT1(:) "            - threshold temperature for snow rain
  !>       \param[in,out] "real(dp) :: DD1(:) "            - Degree-day factor
  !>       \param[in,out] "real(dp) :: DDmax1(:)"          - Maximum Degree-day factor
  !>       \param[in,out] "real(dp) :: IDDP1(:)"           - increase of the degree-day factor per mm
  !>                                                         of increase in precipitation
  !>       \param[in,out] "real(dp) :: fAsp1(:)"           - [1]     PET correction for Aspect at level 1
  !>       \param[in,out] "real(dp) :: HarSamCoeff1(:)"    - [1]     PET Hargreaves Samani coefficient at level 1
  !>       \param[in,out] "real(dp) :: PrieTayAlpha1(:,:)" - [1]     PET Priestley Taylor coefficient at level 1
  !>       \param[in,out] "real(dp) :: aeroResist1(:,:)"   - [s m-1] PET aerodynamical resitance at level 1
  !>       \param[in,out] "real(dp) :: surfResist1(:,:)"   - [s m-1] PET bulk surface resitance at level 1
  !>       \param[in,out] "real(dp) :: HL3(:)"             - threshold parameter for runoff generation
  !>                                                         on impervious layer
  !>       \param[in,out] "real(dp) :: K(:)"               - [d] Muskingum travel time parameters
  !>       \param[in,out] "real(dp) :: xi(:)"              - [1] Muskingum diffusion parameter (attenuation)
  !>       \param[in,out] "real(dp) :: C1(:)"              - routing parameter C1 (Chow, 25-41)
  !>       \param[in,out] "real(dp) :: C2(:)"              - routing parameter C2 (")

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
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
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written Stephan Thober, Dec 2012
  !         Modified, Stephan Thober, Jan 2013 - updated calling sequence for upscaling operators
  !                   Luis Samaniego, Feb 2013 - calling sequence, initial CHECK, call mpr_runoff
  !                   Stephan Thober, Feb 2013 - added subroutine for karstic percolation loss
  !                                              removed L1_, L0_ in variable names
  !                   Stephan Thober, Aug 2015 - moved regionalization of routing to mRM
  !                  Rohini Kumar,    Mar 2016 - changes for handling multiple soil database options

  !TO DOS: all variable names have to be updated as in the mHM call and the sorted. Documentation has to be updated


  subroutine mpr( proc_Mat, & ! IN:    determines which regionalization shall be done
       iFlag_soil,          & ! IN:    flag to handle different soil database
       param,               & ! IN:    global parameter array
       nodata,              & ! IN:    no data value
       mask0,               & ! IN:    mask at Level 0
       geoUnit0,            & ! IN:    geological units at level 0
       geoUnitList,         & ! IN:    List of Ids for geological units
       GeoUnitKar,          & ! IN:    location karstic cells on level 0
       LAILUT,              & ! IN:    Leaf area index for LAIUnit
       LAIUnitList,         & ! IN:    List of ids of each LAI class in LAILUT
       SDB_is_present,      & ! IN:    indicates whether soiltype exists
       SDB_nHorizons,       & ! IN:    Number of Horizons per soiltype
       SDB_nTillHorizons,   & ! IN:    Number of Tillage Horizons
       SDB_sand,            & ! IN:    sand content
       SDB_clay,            & ! IN:    clay content
       SDB_DbM,             & ! IN:    mineral Bulk density
       SDB_Wd,              & ! IN:    weights of mHM
       SDB_RZdepth,         & ! IN:    soil depth
       nHorizons_mHM,       & ! IN:    Number of Horizons in mHM
       horizon_depth,       & ! IN:    Depth of each horizon
       c2TSTu,              & ! IN:    unit transformation coefficient
       fForest1,            & ! IN:    fraction of forest cover at scale L1
       fIperm1,             & ! IN:    fraction of sealed area at scale L1
       fPerm1,              & ! IN:    fraction of permeable area at scale L1
       soilId0,             & ! IN:    soil Ids at level 0
       Asp0,                & ! IN:    [degree] Aspect at Level 0
       LCover_LAI0,         & ! IN:    [1] land cover ID for LAI estimation
       LCover0,             & ! IN:    land use cover at level 0
       slope_emp0,          & ! IN:    Empirical quantiles of slope at Level 0
       cell_id0,            & ! IN:    cell Ids at level 0
       upp_row_L1,          & ! IN:    upper row of L0 block within L1 cell
       low_row_L1,          & ! IN:    lower row of L0 block within L1 cell
       lef_col_L1,          & ! IN:    left column of L0 block within L1 cell
       rig_col_L1,          & ! IN:    right column of L0 block within L1 cell
       nL0_in_L1,           & ! IN:    Number of L0 cells in L0 block within L1 cell
       latitude,            & ! IN:    latitude at level 0
       alpha1,              & ! INOUT: [1]       Exponent for the upper reservoir
       IDDP1,               & ! INOUT:           increase of the degree-day factor per mm of increase in precipitation
       DDmax1,              & ! INOUT:           Maximum Degree-day factor
       DD1,                 & ! INOUT:           Degree-day factor with no precipitation
       fAsp1,               & ! INOUT: [1]       PET correction for Aspect at level 1
       HarSamCoeff1,        & ! INOUT: [1]       PET Hargreaves Samani coefficient at level 1
       PrieTayAlpha1,       & ! INOUT: [1]       PET Priestley Taylor coefficient at level 1
       aeroResist1,         & ! INOUT: [s m-1]   PET aerodynamical resitance at level 1
       surfResist1,         & ! INOUT: [s m-1]   PET bulk surface resitance at level 1
       fRoots1,             & ! INOUT:           fraction of roots in soil horizons
       K0_1,                & ! INOUT: [10^-3 m] Recession coefficient of the upper reservoir, upper outlet
       K1_1,                & ! INOUT: [10^-3 m] Recession coefficient of the upper reservoir, lower outlet
       K2_1,                & ! INOUT:           baseflow recession parameter at level 1
       Kp1,                 & ! INOUT: [d-1]     percolation coefficient
       karstic_loss,        & ! INOUT:           karstic percolation loss parameter
       FC1,                 & ! INOUT: [10^-3 m] field capacity
       SMs1,                & ! INOUT: [10^-3 m] depth of saturated SM cont
       beta1,               & ! INOUT:           Parameter that determines the relative contribution to SM
       TT1,                 & ! INOUT:           threshold temperature for snow rain
       HL1_1,               & ! INOUT: [10^-3 m] Threshhold water depth in upper reservoir 
       !                      !                  (for Runoff  contribution)
       HL3,                 & ! INOUT:           threshold parameter for runoff generation on impervious Layer
       PW1                  & ! INOUT: [10^-3 m] permanent wilting point 
       )

    use mo_message,             only: message
    use mo_upscaling_operators, only: upscale_arithmetic_mean
    use mo_mpr_soilmoist,       only: mpr_sm
    use mo_mpr_SMhorizons,      only: mpr_SMhorizons
    use mo_mpr_runoff,          only: mpr_runoff
    
    implicit none

    ! Input ----------------------------------------------------------
    integer(i4), dimension(:,:),             intent(in)    :: proc_Mat          ! indicate processes
    integer(i4),                             intent(in)    :: iFlag_soil        ! flag to handle different soil database
    real(dp), dimension(:),                  intent(in)    :: param             ! array of global parameters
    real(dp),                                intent(in)    :: nodata            ! nodata value

    logical, dimension(:,:),                 intent(in)    :: mask0             ! mask at level 0 field

    ! baseflow recession
    integer(i4), dimension(:),               intent(in)    :: geoUnit0          ! L0 geological units
    integer(i4), dimension(:),               intent(in)    :: geoUnitList       ! List of geological units
    integer(i4), dimension(:),               intent(in)    :: GeoUnitKar        ! Id of Karstic formations

    real(dp),    dimension(:,:),             intent(in)    :: LAILUT            ! Leaf area index for LAIUnit
    integer(i4), dimension(:),               intent(in)    :: LAIUnitList       ! List of ids of each LAI class in LAILUT

    ! moisture parametrization
    integer(i4), dimension(:),               intent(in)    :: SDB_is_present    ! indicates whether soiltype exists
    integer(i4), dimension(:),               intent(in)    :: SDB_nHorizons     ! Number of Horizons per soiltype
    integer(i4), dimension(:),               intent(in)    :: SDB_nTillHorizons ! Number of Tillage Horizons
    real(dp),    dimension(:,:),             intent(in)    :: SDB_sand          ! sand content
    real(dp),    dimension(:,:),             intent(in)    :: SDB_clay          ! clay content
    real(dp),    dimension(:,:),             intent(in)    :: SDB_DbM           ! mineral Bulk density
    real(dp),    dimension(:,:,:),           intent(in)    :: SDB_Wd            ! weights of mHM horizons
    real(dp),    dimension(:),               intent(in)    :: SDB_RZdepth       ! [mm] Total soil depth
    integer(i4),                             intent(in)    :: nHorizons_mHM     ! Number of Horizons in mHM
    real(dp),    dimension(nHorizons_mHM),   intent(in)    :: horizon_depth     ! [10^-3 m] Depth of each horizon
    real(dp),                                intent(in)    :: c2TSTu            ! unit transformations
    real(dp),    dimension(:),               intent(in)    :: fForest1          ! [1] fraction of forest cover
    real(dp),    dimension(:),               intent(in)    :: fIperm1           ! [1] fraction of sealed area
    real(dp),    dimension(:),               intent(in)    :: fPerm1            ! [1] fraction of permeable area
    integer(i4), dimension(:,:),             intent(in)    :: soilId0           ! soil Ids at level 0
    real(dp),    dimension(:),               intent(in)    :: Asp0              ! [degree] Aspect at Level 0
    integer(i4), dimension(:),               intent(in)    :: LCover_LAI0       ! land cover ID for LAI estimation at level 0
    integer(i4), dimension(:),               intent(in)    :: LCOVER0           ! land cover at level 0

    ! Ids of L0 cells beneath L1 cell
    real(dp),    dimension(:),               intent(in)    :: slope_emp0        ! Empirical quantiles of slope
    integer(i4), dimension(:),               intent(in)    :: cell_id0          ! Cell ids at level 0
    integer(i4), dimension(:),               intent(in)    :: upp_row_L1        ! Upper row of hi res block
    integer(i4), dimension(:),               intent(in)    :: low_row_L1        ! Lower row of hi res block
    integer(i4), dimension(:),               intent(in)    :: lef_col_L1        ! Left column of hi res block
    integer(i4), dimension(:),               intent(in)    :: rig_col_L1        ! Right column of hi res block
    integer(i4), dimension(:),               intent(in)    :: nL0_in_L1         ! Number of L0 cells within a L1 cell

    ! Output of baseflow recession coefficient
    real(dp), dimension(size(upp_row_L1,1)), intent(inout) :: k2_1              ! Level 1 baseflow recession

    ! Output of soilmoisture parametrization
    real(dp), dimension(:,:),                intent(inout) :: beta1             ! Parameter that determines the rel.
    !                                                                           ! contribution to SM, upscal. Bulk den.
    real(dp), dimension(:,:),                intent(inout) :: SMs1              ! [10^-3 m] depth of saturated SM
    real(dp), dimension(:,:),                intent(inout) :: FC1               ! [10^-3 m] field capacity
    real(dp), dimension(:,:),                intent(inout) :: PW1               ! [10^-3 m] permanent wilting point
    real(dp), dimension(:,:),                intent(inout) :: fRoots1           ! fraction of roots in soil horizon
    real(dp), dimension(:),                  intent(inout) :: TT1               ! [degreeC] threshold temperature 
    !                                                                           ! for snow rain 
    real(dp), dimension(:),                  intent(inout) :: DD1               ! [mm-1 degreeC-1] Degree-day factor with
    !                                                                           ! no precipitation
    real(dp), dimension(:),                  intent(inout) :: DDmax1            ! [mm-1 degreeC-1] Maximum Degree-day factor
    real(dp), dimension(:),                  intent(inout) :: IDDP1             ! [d-1 degreeC-1]  Increase of the 
    !                                                                           ! Degree-day factor per mm of
    !                                                                           ! increase in precipitation

    ! Output for PET parameterization
    real(dp), dimension(:),                  intent(inout) :: fAsp1             ! [1]     PET correction for Aspect at level 1
    real(dp), dimension(:),                  intent(inout) :: HarSamCoeff1      ! [1]     PET Hargreaves Samani coeff. at level 1
    real(dp), dimension(:,:),                intent(inout) :: PrieTayAlpha1     ! [1]     PET Priestley Taylor coeff. at level 1
    real(dp), dimension(:,:),                intent(inout) :: aeroResist1       ! [s m-1] PET aerodynamical resitance at level 1
    real(dp), dimension(:,:),                intent(inout) :: surfResist1       ! [s m-1] PET bulk surface resitance at level 1

    ! Output for impervious layer threshold generation
    real(dp), dimension(:),                  intent(inout) :: HL3               ! threshold parameter

    ! Output of mpr runoff
    real(dp),    dimension(:),               intent(inout) :: HL1_1             ! [10^-3 m] Threshhold water depth
    !                                                                           ! in upper reservoir (for Runoff
    !                                                                           ! contribution)
    real(dp),    dimension(:),               intent(inout) :: K0_1              ! [10^-3 m] Recession coefficient
    !                                                                           ! of the upper reservoir, upper outlet
    real(dp),    dimension(:),               intent(inout) :: K1_1              ! [10^-3 m] Recession coefficient
    !                                                                           ! of the upper reservoir, lower outlet
    real(dp),    dimension(:),               intent(in)    :: latitude          ! latitude at level 0
    real(dp),    dimension(:),               intent(inout) :: alpha1            ! [1] Exponent for the upper reservoir
    real(dp),    dimension(:),               intent(inout) :: Kp1               ! [d-1] percolation coefficient

    ! Output of karstic percolation loss
    real(dp),    dimension(:),               intent(inout) :: karstic_loss

    ! Local Variables
    real(dp), dimension(:,:,:), allocatable :: thetaS_till
    real(dp), dimension(:,:,:), allocatable :: thetaFC_till
    real(dp), dimension(:,:,:), allocatable :: thetaPW_till
    real(dp), dimension(:,:,:), allocatable :: Ks       ! saturated hydraulic conductivity
    real(dp), dimension(:,:,:), allocatable :: Db       ! Bulk density
    real(dp), dimension(:,:), allocatable   :: thetaS
    real(dp), dimension(:,:), allocatable   :: thetaFC
    real(dp), dimension(:,:), allocatable   :: thetaPW
    real(dp), dimension(:),   allocatable   :: KsVar_H0 ! relative variability of saturated
    !                                                   ! hydraulic cound. for Horizantal flow
    real(dp), dimension(:),   allocatable   :: KsVar_V0 ! relative variability of saturated
    !                                                   ! hydraulic cound. for vertical flow
    real(dp), dimension(:),   allocatable   :: SMs_FC0  ! soil mositure deficit from
    !                                                   ! field cap. w.r.t to saturation
    real(dp), dimension(size(cell_id0,1))   :: k2_0     ! L0 baseflow parameter
    real(dp), dimension(size(cell_id0,1))   :: fAsp0    ! L0 Aspect
    integer(i4)                             :: mSoil    ! number of soil classes
    integer(i4)                             :: mTill    ! maximum of number of Tillage horizons
    integer(i4)                             :: mHor     ! maximum number of horizons
    integer(i4)                             :: mLC      ! number of Landcover classes
    integer(i4)                             :: iStart 
    integer(i4)                             :: iEnd

    ! ------------------------------------------------------------------
    ! snow parameters 
    ! ------------------------------------------------------------------
    select case( proc_Mat(2,1) )
    case(1)

       iStart = proc_Mat(2,3) - proc_Mat(2,2) + 1
       iEnd   = proc_Mat(2,3)

       call snow_acc_melt_param( param( iStart:iEnd ), c2TSTu, & ! intent(in)
            fForest1, fIperm1, fPerm1,                         & ! intent(in)
            TT1, DD1, IDDP1, DDmax1                            & ! intent(out)
            )
    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "snow pack" does not exist! mo_multi_param_reg')
       stop
    end select

    ! ------------------------------------------------------------------
    ! Soil moisture parametrization 
    ! ------------------------------------------------------------------
    select case( proc_Mat(3,1) )
    case(1)
       
       msoil =   size( SDB_is_present, 1 )
       mLC   = maxval( LCover0, ( LCover0 .ne. int(nodata,i4) )  )
       
       ! depending on which kind of soil database processing is to be performed
       if( iFlag_soil .eq. 0 )then
          mtill = maxval(SDB_nTillHorizons, ( SDB_nTillHorizons .ne. int(nodata,i4) )    )
          mHor  = maxval(SDB_nHorizons,     ( SDB_nHorizons     .ne. int(nodata,i4) )    )
       else if(iFlag_soil .eq. 1) then
          ! here for each soil type both till and non-till soil hydraulic properties are to be estimated
          ! since a given soil type can lie in any horizon (till or non-till ones)
          ! adopt it in a way that it do not break the consistency of iFlag_soil = 0
          ! ** NOTE: SDB_nTillHorizons and SDB_nHorizons are also assigned in
          !          this flag option (see mo_soildatabase.f90 file - read_soil_LUT).
          !          But we are not using those variables here since in this case we have not
          !          varying number of soil horizons or either tillage horizons. 
          !          So assigning them with a value = 1 is more than enough.   
          mtill = 1
          mHor  = 1
       end if
       
       allocate(  thetaS_till(msoil, mtill, mLC) ) 
       allocate( thetaFC_till(msoil, mtill, mLC) ) 
       allocate( thetaPW_till(msoil, mtill, mLC) ) 
       allocate(       thetaS(msoil, mHor      ) ) 
       allocate(      thetaFC(msoil, mHor      ) ) 
       allocate(      thetaPW(msoil, mHor      ) )
       allocate(           Ks(msoil, mHor, mLC ) )
       allocate(           Db(msoil, mHor, mLC ) )       

       ! earlier these variables were allocated with  size(soilId0,1)
       ! in which the variable "soilId0" changes according to the iFlag_soil
       ! so better to use other variable which is common to both soilDB (0 AND 1) flags
       allocate( KsVar_H0( size(cell_id0,1) ) )
       allocate( KsVar_V0( size(cell_id0,1) ) )
       allocate(  SMs_FC0( size(cell_id0,1) ) )
 
       ! first thirteen parameters go to this routine
       iStart = proc_Mat(3,3) - proc_Mat(3,2) + 1
       iEnd   = proc_Mat(3,3) - 4    

       call mpr_sm( param(iStart:iEnd), nodata, iFlag_soil,    &
            SDB_is_present, SDB_nHorizons, SDB_nTillHorizons,  &
            SDB_sand, SDB_clay, SDB_DbM,                       &
            cell_id0, soilId0, LCOVER0,                        &
            thetaS_till, thetaFC_till, thetaPW_till, thetaS,   &
            thetaFC, thetaPW, Ks, Db, KsVar_H0, KsVar_V0, SMs_FC0)

       ! next four parameters go here
       ! (the first three for the fRoots and the fourth one for the beta)
       iStart = proc_Mat(3,3) - 4 + 1
       iEnd   = proc_Mat(3,3)
       call mpr_SMhorizons( param(iStart:iEnd), nodata, iFlag_soil,    &
            nHorizons_mHM, horizon_depth, LCOVER0, soilId0,            &
            SDB_nHorizons, SDB_nTillHorizons,                          &
            thetaS_till,thetaFC_till, thetaPW_till,                    &
            thetaS, thetaFC, thetaPW, SDB_Wd, Db, SDB_DbM, SDB_RZdepth,&
            mask0, cell_id0,                                           &
            Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1, nL0_in_L1, &
            beta1, SMs1, FC1, PW1, fRoots1 )

       deallocate( thetaS_till ) 
       deallocate( thetaFC_till ) 
       deallocate( thetaPW_till ) 
       deallocate( thetaS  ) 
       deallocate( thetaFC ) 
       deallocate( thetaPW )
       deallocate( Ks )
       deallocate( Db )

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "soil moisture parametrization" does not exist! mo_multi_param_reg')
       stop
    END select

    ! ------------------------------------------------------------------
    ! sealed area threshold for runoff generation 
    ! ------------------------------------------------------------------
    select case( proc_Mat( 4, 1) )
    case (1)
       iStart = proc_Mat(4,3) - proc_Mat(4,2) + 1
       iEnd   = proc_Mat(4,3)
       call iper_thres_runoff( param( iStart : iEnd ), HL3 )
    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "runoff_generation" does not exist! mo_multi_param_reg')
       stop
    end select

    ! ------------------------------------------------------------------
    ! potential evapotranspiration (PET)
    ! ------------------------------------------------------------------
    select case( proc_Mat( 5,1 ) )
    case(0) ! aspect correction of input PET
       iStart = proc_Mat(5,3) - proc_Mat(5,2) + 1
       iEnd   = proc_Mat(5,3)    
       call pet_correct( cell_id0, latitude, Asp0, param( iStart : iEnd), nodata, fAsp0 )
       fAsp1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, fAsp0 )
    case(1) ! Hargreaves-Samani method
       iStart = proc_Mat(5,3) - proc_Mat(5,2) + 1
       iEnd   = proc_Mat(5,3)    
       call pet_correct( cell_id0, latitude, Asp0, param( iStart : iEnd - 1), nodata, fAsp0 )
       fAsp1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, fAsp0 )
       HarSamCoeff1 = param(iEnd)
    case(2) ! Priestley-Taylor Method
       iStart = proc_Mat(5,3) - proc_Mat(5,2) + 1
       iEnd   = proc_Mat(5,3)    
       call priestley_taylor_alpha(LCover_LAI0, LAILUT, LAIUnitList, param(iStart : iEnd),       & 
            mask0, nodata, cell_id0, nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1,  &
            PrieTayAlpha1)
    case(3) ! Penman-Monteith method
       iStart = proc_Mat(5,3) - proc_Mat(5,2) + 1
       iEnd   = proc_Mat(5,3) 
       call aerodynamical_resistance(LCover0, LAILUT, param(iStart : iEnd - 1), mask0,   & 
            nodata, cell_id0, nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1, &
            aeroResist1)
       call bulksurface_resistance(LCover_LAI0, LAILUT, LAIUnitList, param(iEnd), mask0,                             & 
            nodata, cell_id0, nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1, &
            surfResist1)
    case default
       call message()
       call message('***ERROR: Process description for process "pet correction" does not exist! mo_multi_param_reg')
       stop
    end select

    ! ------------------------------------------------------------------
    ! interflow
    ! ------------------------------------------------------------------
    select case( proc_Mat( 6, 1) )
    case (1)
       !
       iStart = proc_Mat(6,3) - proc_Mat(6,2) + 1
       iEnd   = proc_Mat(6,3)
       call mpr_runoff( LCOVER0, mask0, nodata, SMs_FC0, slope_emp0,        &
            KsVar_H0, param(iStart:iEnd), cell_id0, upp_row_L1, low_row_L1, &
            lef_col_L1, rig_col_L1, nL0_in_L1, c2TSTu, HL1_1, K0_1,         &
            K1_1 , alpha1 )
    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "interflow" does not exist! mo_multi_param_reg')
       stop
    END select

    ! ------------------------------------------------------------------
    ! percolation cofficient, karstic percolation loss
    ! ------------------------------------------------------------------
    select case( proc_Mat( 7, 1) )
    case(1)

       iStart = proc_Mat(7,3) - proc_Mat(7,2) + 1
       iEnd   = proc_Mat(7,3)
       call karstic_layer(                                                     & ! In
            param( iStart : iEnd ),                                            & ! In
            geoUnitKar,  geoUnit0, geoUnitList, mask0, nodata,                 & ! In
            SMs_FC0, KsVar_V0, cell_id0,                                       & ! In
            nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1, c2TSTu, & ! In
            karstic_loss, Kp1                                                  & ! Out
            )

       deallocate( KsVar_H0 )
       deallocate( KsVar_V0 )
       deallocate( SMs_FC0 )

    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "percolation" does not exist! mo_multi_param_reg')
       stop
    end select

    ! ------------------------------------------------------------------
    ! Regionalized routing parameters
    ! ------------------------------------------------------------------
    select case( proc_Mat( 8, 1) )
    case(0)  ! routing is off
    case(1)  !ST this part is moved to mRM
    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "routing" does not exist! mo_multi_param_reg')
       stop
    end select

    ! ------------------------------------------------------------------
    ! baseflow recession parameter 
    ! ------------------------------------------------------------------
    select case( proc_Mat(9,1) )
    case(1)

       ! the number of process parameters, so the number in proc_Mat(9,2) has
       ! to be equal to the size of geo_unit_list
       iStart = proc_Mat(9,3) - proc_Mat(9,2) + 1
       iEnd   = proc_Mat(9,3)      

       call baseflow_param( param( iStart : iEnd ),&
            geoUnit0, geoUnitList, nodata, k2_0)
       !
       ! Upscale by arithmetic mean
       k2_1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, k2_0 )
       !
       ! correction and unit conversion
       ! if percolation is ON: correct K2 such that it is at least k1
       if ( proc_Mat(7,1) .gt. 0 ) k2_1 = merge( K1_1, k2_1, k2_1 .lt. K1_1 )
       k2_1 = c2TSTu / k2_1
       !
    case DEFAULT
       call message()
       call message('***ERROR: Process description for process "baseflow Recession" does not exist! mo_multi_param_reg')
       stop
    end select

  end subroutine mpr

  ! ----------------------------------------------------------------------------

  !      NAME
  !        baseflow_param

  !>       \brief baseflow recession parameter

  !>       \details This subroutine calculates the baseflow recession parameter
  !>                based on the geological units at the Level 0 scale. For each level 0
  !>                cell, it assigns the value specified in the parameter array param for the
  !>                geological unit in this cell.\n
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>                    - param(1) = GeoParam(1,:) \n
  !>                    - param(2) = GeoParam(2,:) \n
  !>                    - ...\n

  !      INTENT(IN)
  !>       \param[in] "real(dp) :: param(:)" - array of global baseflow recession
  !>                                           parameters
  !>       \param[in] "integer(i4) :: geoUnit0(:,:)" - array of geological units
  !>                                           at Level 0
  !>       \param[in] "integer(i4) :: geoUnitList(:)" - array of indices for 
  !>                                           geological units.
  !>       \param[in] "real(dp) :: nodata"    - no data value
  
  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>       \param[out] "real(dp) :: k2_0" - baseflow recession parameter at Level 0

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None
  
  !      HISTORY
  !>       \author Stephan Thober, Rohini Kumar
  !>       \date Dec 2012
  !        Written  Stephan Thober, Dec 2012
  !        Modified Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)

  subroutine baseflow_param( &
       param,                & ! Intent(in)
       geoUnit0,             & ! Intent(in)
       geoUnitList,          & ! Intent(in)
       nodata,               & ! Intent(in)
       k2_0                  & ! Intent(out)
       )

    !$  use omp_lib

    implicit none

    ! Input
    real(dp),    dimension(:), intent(in) :: param         ! list of required parameters
    integer(i4), dimension(:), intent(in) :: geoUnit0      ! ids of geological units at L0
    integer(i4), dimension(:), intent(in) :: geoUnitList   ! list of geological units
    real(dp),                  intent(in) :: nodata        ! nodata value

    ! output
    real(dp),    dimension(:), intent(out):: k2_0          ! baseflow recession coefficient

    ! local variables
    integer(i4)                           :: ii            ! loop variable
    integer(i4), dimension(1)             :: gg            ! geo unit

    if ( size(param) .ne. size( geoUnitList ) ) &
         stop ' mo_multi_param_reg: baseflow_param: size mismatch, subroutine baseflow parameters '

    k2_0 = nodata

    !$OMP PARALLEL
    !$OMP DO PRIVATE(gg) SCHEDULE(STATIC)
    do ii = 1, size(k2_0)
       ! get parameter index in geoUnitList
       gg = minloc( abs( geoUnitList - geoUnit0(ii) ) )
       k2_0(ii) = param( gg(1) )
    end do
    !$OMP END DO
    !$OMP END PARALLEL


  end subroutine baseflow_param

  ! ----------------------------------------------------------------------------

  !      NAME
  !         snow_acc_melt_param

  !>        \brief Calculates the snow parameters.

  !>        \details This subroutine calculates the snow parameters
  !>        threshold temperature (TT), degree-day factor without precipitation (DD)
  !>        and maximum degree-day factor (DDmax) as well as increase of degree-day 
  !>        factor per mm of increase in precipitation (IDDP).\n
  !>
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param(1) = snowTreshholdTemperature        \n
  !>           - param(2) = degreeDayFactor_forest          \n
  !>           - param(3) = degreeDayFactor_impervious      \n
  !>           - param(4) = degreeDayFactor_pervious        \n
  !>           - param(5) = increaseDegreeDayFactorByPrecip \n
  !>           - param(6) = maxDegreeDayFactor_forest       \n
  !>           - param(7) = maxDegreeDayFactor_impervious   \n
  !>           - param(8) = maxDegreeDayFactor_pervious     \n

  !>     INTENT(IN)
  !>        \param[in] "real(dp) :: param(8)" - There are eight snow parameters required
  !>        \param[in] "real(dp) :: c2TSTu"   - unit transformation coefficient
  !>        \param[in] "real(dp) :: fForest1(:)" - fraction of forest cover at scale L1
  !>        \param[in] "real(dp) :: fIperm1(:)" - fraction of sealed area at scale L1
  !>        \param[in] "real(dp) :: fPerm1(:)" - fraction of permeable area at scale L1

  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: TT1(:) "   - threshold temperature for snow rain
  !>        \param[out] "real(dp) :: DD1(:) "   - Degree-day factor
  !>        \param[out] "real(dp) :: DDmax1(:)" - Maximum Degree-day factor
  !>        \param[out] "real(dp) :: IDDP1(:)"  - increase of the degree-day factor per mm of
  !>                                                increase in precipitation

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written   Stephan Thober, Dec 2012
  !         Modified, Juliane Mai,    Oct 2013 - OLD parametrization
  !                                                --> param(1) = snowTreshholdTemperature 
  !                                                --> param(2) = degreeDayFactor_forest
  !                                                --> param(3) = degreeDayFactor_impervious   
  !                                                --> param(4) = degreeDayFactor_pervious     
  !                                                --> param(5) = increaseDegreeDayFactorByPrecip 
  !                                                --> param(6) = maxDegreeDayFactor_forest     
  !                                                --> param(7) = maxDegreeDayFactor_impervious  
  !                                                --> param(8) = maxDegreeDayFactor_pervious   
  !                                             -------------------------------
  !                                             degreeDayFactor_impervious    = degreeDayFactor_forest + delta_1 + delta_2 
  !                                             degreeDayFactor_pervious      = degreeDayFactor_forest + delta_1
  !                                             maxDegreeDayFactor_forest     = degreeDayFactor_forest + delta_3
  !                                             maxDegreeDayFactor_impervious = degreeDayFactor_impervious + delta_5
  !                                                                           = degreeDayFactor_forest + delta_1 + delta_2 + delta_5
  !                                             maxDegreeDayFactor_pervious   = degreeDayFactor_pervious + delta_4
  !                                                                           = degreeDayFactor_forest + delta_1 + delta_4
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = snowTreshholdTemperature 
  !                                                --> param(2) = degreeDayFactor_forest
  !                                                --> param(3) = delta_2   
  !                                                --> param(4) = delta_1    
  !                                                --> param(5) = increaseDegreeDayFactorByPrecip 
  !                                                --> param(6) = delta_3     
  !                                                --> param(7) = delta_5
  !                                                --> param(8) = delta_4 
  !                    Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)

  subroutine snow_acc_melt_param(                &
       param, c2TSTu, fForest1, fIperm1, fPerm1, &  ! Intent(in)
       TT1, DD1, IDDP1, DDmax1                   &  ! Intent(out)
       )

    implicit none

    ! Input
    real(dp), dimension(8), intent(in)  :: param      ! eight global parameters
    real(dp),               intent(in)  :: c2TSTu     ! unit transformations
    real(dp), dimension(:), intent(in)  :: fForest1 ! [1] fraction of forest cover
    real(dp), dimension(:), intent(in)  :: fIperm1  ! [1] fraction of sealed area
    real(dp), dimension(:), intent(in)  :: fPerm1   ! [1] fraction of permeable area

    ! Output
    real(dp), dimension(:), intent(out) :: TT1   ! [degreeC] threshold temperature for snow rain 
    real(dp), dimension(:), intent(out) :: DD1   ! [mm-1 degreeC-1] Degree-day factor with
    ! no precipitation
    real(dp), dimension(:), intent(out) :: DDmax1! [mm-1 degreeC-1] Maximum Degree-day factor
    real(dp), dimension(:), intent(out) :: IDDP1 ! [d-1 degreeC-1]  Increase of the Degree-day 
    ! factor per mm of increase in precipitation

    ! local
    real(dp) :: tmp_degreeDayFactor_forest,    tmp_degreeDayFactor_impervious,    tmp_degreeDayFactor_pervious   
    real(dp) :: tmp_maxDegreeDayFactor_forest, tmp_maxDegreeDayFactor_impervious, tmp_maxDegreeDayFactor_pervious

    tmp_degreeDayFactor_forest         = param(2)                                    ! OLD: param(2)
    tmp_degreeDayFactor_impervious     = param(2) + param(4) + param(3)              ! OLD: param(3)
    tmp_degreeDayFactor_pervious       = param(2) + param(4)                         ! OLD: param(4)
    tmp_maxDegreeDayFactor_forest      = param(2)                       + param(6)   ! OLD: param(6)
    tmp_maxDegreeDayFactor_impervious  = param(2) + param(4) + param(3) + param(7)   ! OLD: param(7)
    tmp_maxDegreeDayFactor_pervious    = param(2) + param(4)            + param(8)   ! OLD: param(8)

    TT1    =  param(1)
    IDDP1  =  param(5)

    DD1    = ( &
         tmp_degreeDayFactor_forest     * fForest1 + &
         tmp_degreeDayFactor_impervious * fIperm1  + &
         tmp_degreeDayFactor_pervious   * fPerm1     ) * c2TSTu
    DDmax1 = ( &
         tmp_maxDegreeDayFactor_forest     * fForest1 + &
         tmp_maxDegreeDayFactor_impervious * fIperm1  + &
         tmp_maxDegreeDayFactor_pervious   * fPerm1     ) * c2TSTu

  end subroutine snow_acc_melt_param

  ! ----------------------------------------------------------------------------

  !      NAME
  !         PET correction due to aspect

  !>        \brief correction of PET

  !>        \details to be done by Kumar\n
  !>        ....
  !>        
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param(1) = minCorrectionFactorPET \n
  !>           - param(2) = maxCorrectionFactorPET \n
  !>           - param(3) = aspectTresholdPET      \n

  !      INTENT(IN)
  !>        \param[in] "integer(i4) :: id0(:,:)" - cell id at Level 0
  !>        \param[in] "real(dp) :: nodata"        - no data value
  !>        \param[in] "real(dp) :: param(3)"      - the three process parameters

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: fAsp0(:,:)" - [1] PET correction factor for aspect

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written  Stephan Thober, Dec 2012
  !         Modified Juliane Mai,    Oct 2013 - OLD parametrization
  !                                                --> param(1) = minCorrectionFactorPET
  !                                                --> param(2) = maxCorrectionFactorPET 
  !                                                --> param(3) = aspectTresholdPET
  !                                             -------------------------------
  !                                             maxCorrectionFactorPET = minCorrectionFactorPET + delta
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = minCorrectionFactorPET 
  !                                                --> param(2) = delta 
  !                                                --> param(3) = aspectTresholdPET
  !                  Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)
  !                  Stephan Thober, Sep 2015 - Mapping L1 to Lo, latitude on L0
  !                  Luis Samaniego, Sep 2015 - PET correction on the southern hemisphere

  subroutine pet_correct( Id0, latitude_l0, Asp0, param, nodata, fAsp0 )

    implicit none

    ! Input
    integer(i4), dimension(:), intent(in) :: id0      ! Level 0 cell id
    real(dp),    dimension(:), intent(in) :: latitude_l0 ! latitude on l0
    real(dp),                  intent(in) :: nodata   ! no data value
    real(dp),    dimension(3), intent(in) :: param    ! process parameters
    real(dp),    dimension(:), intent(in) :: Asp0     ! [degree] Aspect at Level 0

    ! Output
    real(dp),    dimension(:), intent(out):: fAsp0    ! PET correction for Aspect

    ! local
    real(dp), dimension(size(id0, 1)) :: fAsp0S       ! PET correction for Aspect, south

    logical,  dimension(size(id0, 1)) :: mask_north_hemisphere_l0
    
    real(dp)                          :: tmp_maxCorrectionFactorPET

    mask_north_hemisphere_l0 = merge(.TRUE.,.FALSE., latitude_l0 .gt. 0.0_dp)
       
    tmp_maxCorrectionFactorPET = param(1) + param(2)

    ! for cells on the northern hemisphere
    !$OMP PARALLEL
    fAsp0 = merge(  &
         param(1) + ( tmp_maxCorrectionFactorPET - param(1)) / param(3) * asp0, &
         param(1) + ( tmp_maxCorrectionFactorPET - param(1) ) / (360._dp - param(3)) * (360._dp - Asp0), &
         !         ( asp0 < param(3) ) .and. mask_north_hemisphere_l0  )
          asp0 < param(3)   )
    fAsp0 = merge( fAsp0, nodata, Id0 /= int(nodata, i4)  )
    !$OMP END PARALLEL

    ! for cells on the southern hemisphere
    !$OMP PARALLEL
    fAsp0S = merge( &
         param(1) + ( tmp_maxCorrectionFactorPET - param(1) ) / (360._dp - param(3)) * (360._dp - Asp0), &
         param(1) + ( tmp_maxCorrectionFactorPET - param(1)) / param(3) * asp0, &
         asp0 < param(3) )
    fAsp0S = merge( fAsp0S, nodata, Id0 /= int(nodata, i4) )
    !$OMP END PARALLEL

    !$OMP PARALLEL
    fAsp0 = merge( fAsp0, fAsp0S, mask_north_hemisphere_l0 )
    !$OMP END PARALLEL 

  end subroutine pet_correct

  ! ----------------------------------------------------------------------------

  !      NAME
  !         iper_thres_runoff

  !>        \brief sets the impervious layer threshold parameter for runoff generation

  !>        \details to be done by Kumar\n
  !>        ....
  !>        
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param(1) = imperviousStorageCapacity \n

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: param" - given threshold parameter

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: HL3(:)" - distributed parameter field

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written Stephan Thober, Dec 2012
  !         Modified Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)

  subroutine iper_thres_runoff( param, HL3 )

    implicit none

    ! Input
    real(dp), dimension(1), intent(in)  :: param ! threshold parameters

    ! Output
    real(dp), dimension(:), intent(out) :: HL3

    HL3 = param(1)

  end subroutine iper_thres_runoff

  ! ----------------------------------------------------------------------------

  !     NAME
  !         karstic_layer

  !>        \brief calculates the Karstic percolation loss

  !>        \details This subroutine calls first the karstic_fraction upscaling
  !>                 routine for determine the karstic fraction area for every Level 1
  !>                 cell. Then, the karstic percolation loss is estimated given two
  !>                 shape parameters by
  !>                 \f[ karstic_loss = 1 + ( fKarArea * param(1)) *( (-1)**INT(param(2),i4) ) \f]
  !>                 where \f$ karstic_loss \f$ is the karstic percolation loss and \f$ fKarArea \f$
  !>                 ist the fraction of karstic area at level 1\n
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>                    - param(1) = rechargeCoefficient           \n
  !>                    - param(2) = rechargeFactor_karstic        \n
  !>                    - param(3) = gain_loss_GWreservoir_karstic \n

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: nGeoUnits"      - number of geological formations
  !>        \param[in] "integer(i4) :: geoUnitKar(:)"  - number of Karstic formation
  !>        \param[in] "integer(i4) :: geoUnit0(:)"    - id of the Karstic formation
  !>        \param[in] "logical     :: mask0(:,:)"     - mask at level 0
  !>        \param[in] "real(dp)    :: nodata"         - given nodata value

  !>        \param[in] "integer(i4) :: nL0_in_L1(:)"   - number of l0 cells within a l1 cell
  !>        \param[in] "integer(i4) :: Upp_row_L1(:)"  - upper row of a l1 cell in l0 grid
  !>        \param[in] "integer(i4) :: Low_row_L1(:)"  - lower row of a l1 cell in l0 grid
  !>        \param[in] "integer(i4) :: Lef_col_L1(:)"  - left col of a l1 cell in l0 grid
  !>        \param[in] "integer(i4) :: Rig_col_L1(:)"  - right col of a l1 cell in l0 grid
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp)    :: karstic_loss(:)"  [-]    Karstic percolation loss
  !>        \param[out] "real(dp)    :: L1_Kp(:)"         [d-1] percolation coefficient

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
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
  !>        \author Rohini Kumar, Stephan Thober
  !>        \date Feb 2013
  !         Modified Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)
  !         Modified Stephan Thober, Dec 2013 - changed intent(inout) to intent(out)
  !
  subroutine karstic_layer( &
       param,        & ! IN:         parameters
       geoUnitKar,   & ! IN:         number of Karstic formation
       geoUnit0,     & ! IN:         id of the Karstic formation
       geoUnitlist,  & ! IN:         id of geo units
       mask0,        & ! IN:         mask at level 0
       nodata,       & ! IN:         given nodata value
       SMs_FC0,      & ! IN:  [-]    soil mositure deficit from field
       !               !             capacity w.r.t to saturation
       KsVar_V0,     & ! IN:  [-]    relative variability of saturated
       cell_id0,     & ! IN:         cell id at Level 0
       nL0_in_L1,    & ! IN:         number of l0 cells within a l1 cell
       Upp_row_L1,   & ! IN:         upper row of a l1 cell in l0 grid
       Low_row_L1,   & ! IN:         lower row of a l1 cell in l0 grid
       Lef_col_L1,   & ! IN:         left col of a l1 cell in l0 grid
       Rig_col_L1,   & ! IN:         right col of a l1 cell in l0 grid
       c2TSTu,       & ! IN:         unit transformations
       karstic_loss, & ! OUT: [-]    Karstic percolation loss
       L1_Kp         & ! OUT: [d-1]  percolation coefficient
       )

    use mo_kind,                only: i4, dp
    use mo_upscaling_operators, only: L0_fractionalCover_in_Lx, upscale_arithmetic_mean
    !$  use omp_lib

    implicit none

    ! Input
    real(dp),    dimension(3),   intent(in) :: param        ! parameters
    integer(i4), dimension(:),   intent(in) :: geoUnitKar   ! number of Karstic formation
    integer(i4), dimension(:),   intent(in) :: geoUnit0     ! id of the Karstic formation
    integer(i4), dimension(:),   intent(in) :: geoUnitlist  ! id of geo units

    logical,     dimension(:,:), intent(in) :: mask0        ! mask at level 0
    real(dp),                    intent(in) :: nodata       ! given nodata value
    real(dp),    dimension(:),   intent(in) :: SMs_FC0      ! [-] soil mositure deficit from field
    !                                                       ! capacity w.r.t to saturation
    real(dp),    dimension(:),   intent(in) :: KsVar_V0     ! [-] relative variability of saturated

    integer(i4), dimension(:),   intent(in) :: cell_id0     ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in) :: nL0_in_L1    ! number of l0 cells within a l1 cell
    integer(i4), dimension(:),   intent(in) :: Upp_row_L1   ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in) :: Low_row_L1   ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in) :: Lef_col_L1   ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in) :: Rig_col_L1   ! right col of a l1 cell in l0 grid
    real(dp),                    intent(in) :: c2TSTu       ! unit transformations
    
    ! Output
    real(dp), dimension(:),   intent(out)   :: karstic_loss ! [-]    Karstic percolation loss
    real(dp), dimension(:),   intent(out)   :: L1_Kp        ! [d-1] percolation coefficient

    ! Local variables
    real(dp), dimension(:), allocatable     :: fKarArea     ! fraction of karstic area
    real(dp), dimension(size(SMs_FC0,1))    :: tmp          ! temporal variable
    integer(i4)                             :: nGeoUnits
    integer(i4)                             :: i


    ! ------------------------------------------------------------------
    ! PERCOLATION; 1/Kp = f(Ks)
    ! ------------------------------------------------------------------
    ! Regionalise Kp with variability of last soil layer property
    !$OMP PARALLEL
    tmp = merge( param(1) * ( 1.0_dp + SMs_FC0 ) / ( 1.0_dp + KsVar_V0 ), &
         nodata, cell_id0 .ne. int(nodata, i4) )
    !$OMP END PARALLEL

    L1_Kp = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
         Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, tmp )

    ! minimum constrains
    L1_Kp = merge( 2.0_dp, L1_Kp, L1_Kp .lt. 2.0_dp )
    !
    ! Unit conversion
    L1_Kp = c2TSTu / L1_Kp

    nGeoUnits = size(geoUnitlist,1)

    ! 1st calculate fraction of Karstic area
    allocate( fKarArea( size(karstic_loss,1) ) )
    fKarArea = 0.0_dp

    do i = 1, nGeoUnits
       if(GeoUnitKar(i) .eq. 0) cycle
       fKarArea(:) = L0_fractionalCover_in_Lx( geoUnit0, geoUnitlist(i), mask0, &
            Upp_row_L1, Low_row_L1, Lef_col_L1, Rig_col_L1, nL0_in_L1 )
    end do

    ! 2nd calculate karstic_loss
    karstic_loss = 1.0_dp - ( fKarArea * param(2) )

    deallocate( fKarArea )

  end subroutine karstic_layer

  ! ----------------------------------------------------------------------------

  !      NAME
  !        canopy_intercept_param

  !>       \brief estimate effective maximum interception capacity at L1

  !>       \details estimate effective maximum interception capacity at L1 for a given 
  !>                Leaf Area Index field. \n
  !>                Global parameters needed (see mhm_parameter.nml):\n
  !>                Process Case 1:\n
  !>                   - param(1) = canopyInterceptionFactor \n

  !      INTENT(IN)
  !>       \param[in] "integer(i4)  :: proc_Mat(:,:)"  - process matrix
  !>       \param[in] "real(dp)     :: param(:)"       - array of global parameters
  !>       \param[in] "real(dp)     :: LAI0(:)"        - LAI at level-0
  !>       \param[in] "integer(i4)  :: nL0_in_L1 (:)"  - Number of L0 cells within a L1 cell
  !>       \param[in] "integer(i4)  :: upp_row_L1(:)"  - Upper row of high resolution block
  !>       \param[in] "integer(i4)  :: low_row_L1(:)"  - Lower row of high resolution block
  !>       \param[in] "integer(i4)  :: lef_col_L1(:)"  - Left column of high resolution block
  !>       \param[in] "integer(i4)  :: rig_col_L1(:)"  - Right column of high resolution block
  !>       \param[in] "integer(i4)  :: cell_id0  (:)"  - Cell ids at level 0
  !>       \param[in] "logical      :: mask0(:,:)"     - mask at level 0 field
  !>       \param[in] "real(dp)     :: nodata"         - nodata value

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>       \param[out] "real(dp) :: max_intercept1(:)" - maximum canopy interception at Level-1

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      EXAMPLE
  !         calling sequence
  !         call canopy_intercept_param(proc_Mat, param,             &    
  !				        LAI0, nL0_in_L1, upp_row_L1, &    
  !				        low_row_L1, lef_col_L1,      &   
  !				        rig_col_L1, cell_id0, mask0, &  
  !				        nodata, max_intercept1 ) 

  !      HISTORY
  !>         \author Rohini Kumar
  !>         \date Aug. 2013

  ! ------------------------------------------------------------------
  subroutine canopy_intercept_param(proc_Mat, param,             &
       LAI0, nL0_in_L1, upp_row_L1, &
       low_row_L1, lef_col_L1,      &
       rig_col_L1, cell_id0, mask0, &
       nodata, max_intercept1 ) 

    use mo_upscaling_operators, only: upscale_arithmetic_mean
    use mo_string_utils,        only: num2str
    use mo_message,             only: message

    implicit none

    ! input
    integer(i4), dimension(:,:), intent(in)  :: proc_Mat       ! indicate processes
    real(dp),    dimension(:),   intent(in)  :: param          ! array of global parameters
    real(dp),    dimension(:),   intent(in)  :: LAI0           ! LAI at level-0
    integer(i4), dimension(:),   intent(in)  :: nL0_in_L1      ! Number of L0 cells within a L1 cell
    integer(i4), dimension(:),   intent(in)  :: upp_row_L1     ! Upper row of high resolution block
    integer(i4), dimension(:),   intent(in)  :: low_row_L1     ! Lower row of high resolution block
    integer(i4), dimension(:),   intent(in)  :: lef_col_L1     ! Left column of high resolution block
    integer(i4), dimension(:),   intent(in)  :: rig_col_L1     ! Right column of high resolution block
    integer(i4), dimension(:),   intent(in)  :: cell_id0       ! Cell ids at level 0
    logical,     dimension(:,:), intent(in)  :: mask0          ! mask at level 0 field
    real(dp),                    intent(in)  :: nodata         ! nodata value
    ! output
    real(dp),    dimension(:),   intent(out) :: max_intercept1 ! max interception at level-1

    ! local variables
    integer(i4)                             :: iStart, iEnd
    real(dp), dimension(:), allocatable     :: max_intercept0
    real(dp), dimension(:), allocatable     :: gamma_intercept
    
    ! ------------------------------------------------------------------
    ! Maximum interception parameter 
    ! ------------------------------------------------------------------
    select case( proc_Mat(1,1))
    case(1)
       iStart = proc_Mat(1,3) - proc_Mat(1,2) + 1
       iEnd   = proc_Mat(1,3)      

       ! allocate space
       allocate( gamma_intercept(iEnd-iStart+1     ) )
       allocate( max_intercept0 (size(cell_id0, 1) ) )

       ! estimate max. intercept at Level-0
       gamma_intercept(:) = param(iStart:iEnd)
       !$OMP PARALLEL
       max_intercept0(:)  = LAI0(:) * gamma_intercept(1)
       !$OMP END PARALLEL

       ! Upscale by arithmetic mean
       max_intercept1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1,      &
            Rig_col_L1, cell_id0, mask0, nodata, max_intercept0 )

       deallocate( gamma_intercept)
       deallocate( max_intercept0 )        
    CASE DEFAULT
       call message('mo_multi_param_reg: This proc_Mat=',num2str(proc_Mat(1,1)),' is not implemented!')
       stop
    end select

  end subroutine canopy_intercept_param


  ! ----------------------------------------------------------------------------

  !      NAME
  !        aerodynamical_resistance

  !>       \brief Regionalization of aerodynamic resistance

  !>       \details estimation of aerodynamical resistance
  !>                Global parameters needed (see mhm_parameter.nml):\n
  !>                   - param(1) = canopyheigth_forest             \n
  !>                   - param(2) = canopyheigth_impervious         \n
  !>                   - param(3) = canopyheigth_pervious           \n
  !>                   - param(4) = displacementheight_coeff        \n
  !>                   - param(5) = roughnesslength_momentum_coeff  \n
  !>                   - param(6) = roughnesslength_heat_coeff      \n

  !      INTENT(IN)
  !>       \param[in] "integer(i4)  :: LCover0(:)"     - land cover at level 0
  !>       \param[in] "real(dp)     :: LAILUT(:)"      - LUT of LAi values
  !>       \param[in] "integer(i4)  :: LAIUnitList(:)" - List of ids of each LAI class in LAILUT
  !>       \param[in] "real(dp)     :: param(:)"       - vector with global parameters
  !>       \param[in] "logical      :: mask0(:,:)"     - mask at level 0 field
  !>       \param[in] "real(dp)     :: nodata"         - nodata value 
  !>       \param[in] "integer(i4)  :: cell_id0  (:)"  - Cell ids at level 0
  !>       \param[in] "integer(i4)  :: nL0_in_L1 (:)"  - Number of L0 cells within a L1 cell
  !>       \param[in] "integer(i4)  :: Upp_row_L1(:)"  - Upper row of high resolution block
  !>       \param[in] "integer(i4)  :: Low_row_L1(:)"  - Lower row of high resolution block
  !>       \param[in] "integer(i4)  :: Lef_col_L1(:)"  - Left column of high resolution block
  !>       \param[in] "integer(i4)  :: Rig_col_L1(:)"  - Right column of high resolution block

  !     INTENT(INOUT)
  !        None

  !     INTENT(OUT)
  !>       \param[out] "real(dp)    :: aerodyn_resistance1(:)" - [s m-1] aerodynamical resistance

  !     INTENT(IN), OPTIONAL
  !        None

  !     INTENT(INOUT), OPTIONAL
  !        None

  !     INTENT(OUT), OPTIONAL
  !        None

  !     RETURN
  !        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !        None

  !     HISTORY
  !>       \author Matthias Zink
  !>       \date   Apr 2013

  subroutine aerodynamical_resistance( &
       LCover0,                        & ! land cover at level 0
       LAILUT,                         & ! look up table for LAI
       param,                          & ! parameter values (size=6)
       mask0,                          & ! mask at level 0
       nodata,                         & ! given nodata value
       cell_id0,                       & ! cell id at Level 0
       nL0_in_L1,                      & ! number of l0 cells within a l1 cell
       Upp_row_L1,                     & ! upper row of a l1 cell in l0 grid
       Low_row_L1,                     & ! lower row of a l1 cell in l0 grid
       Lef_col_L1,                     & ! left col of a l1 cell in l0 grid
       Rig_col_L1,                     & ! right col of a l1 cell in l0 grid
       aerodyn_resistance1             & ! aerodynmaical resistance
       )

    use mo_upscaling_operators, only: upscale_arithmetic_mean
    use mo_mhm_constants,       only: YearMonths_i4, WindMeasHeight, karman
    use mo_constants,           only: eps_dp

    implicit none

    integer(i4), dimension(:),   intent(in)  :: LCover0    ! land cover field
    real(dp),    dimension(:,:), intent(in)  :: LAILUT     ! look up table for LAI
    real(dp),    dimension(6),   intent(in)  :: param      ! input parameter
    !                                                      ! dim1=land cover class, dim2=month of year
    logical,     dimension(:,:), intent(in)  :: mask0      ! mask at level 0
    real(dp),                    intent(in)  :: nodata     ! given nodata value
    integer(i4), dimension(:),   intent(in)  :: cell_id0   ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in)  :: nL0_in_L1  ! number of l0 cells within a l1 cell
    integer(i4), dimension(:),   intent(in)  :: Upp_row_L1 ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Low_row_L1 ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Lef_col_L1 ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Rig_col_L1 ! right col of a l1 cell in l0 grid
    ! Output
    real(dp),    dimension(:,:), intent(out) :: aerodyn_resistance1

    ! local
    integer(i4)                            :: iMon
    real(dp)                               :: maxLAI
    real(dp), dimension(:),   allocatable  :: zm
    real(dp), dimension(:),   allocatable  :: canopy_height0
    real(dp), dimension(:),   allocatable  :: zm_zero, zh_zero, displace
    real(dp), dimension(:,:), allocatable  :: aerodyn_resistance0        ! dim 1 = number of cells on level 0,
    !                                                                    ! dim 2 = number of months in year (12)
    !
    ! ID   LAI classes                 
    ! 1    Coniferous-forest        
    ! 2    Deciduous-forest         
    ! 3    Mixed-forest             
    ! 4    Sparsely-populated-forest
    ! 5    Sealed-Water-bodies      
    ! 6    Viniculture              
    ! 7    Intensive-orchards       
    ! 8    Pasture                  
    ! 9    Fields                   
    ! 10   Wetlands                 

    ! initialize some things
    allocate(zm                  (size(LCover0, dim=1)               )) ; zm                  = nodata
    allocate(zm_zero             (size(LCover0, dim=1)               )) ; zm_zero             = nodata
    allocate(zh_zero             (size(LCover0, dim=1)               )) ; zh_zero             = nodata
    allocate(displace            (size(LCover0, dim=1)               )) ; displace            = nodata
    allocate(canopy_height0      (size(LCover0, dim=1)               )) ; canopy_height0      = nodata
    allocate(aerodyn_resistance0 (size(LCover0, dim=1), YearMonths_i4)) ; aerodyn_resistance0 = nodata
    aerodyn_resistance1 = nodata
    !
    ! regionalization of canopy height
    ! substitute with canopy height
    canopy_height0 = merge(param(1), canopy_height0, LCover0 == 1)  ! forest
    canopy_height0 = merge(param(2), canopy_height0, LCover0 == 2)  ! impervious
    !
    maxLAI = MAXVAL(LAILUT(7,:))
    !
    do iMon = 1, YearMonths_i4
       !
       ! pervious canopy height is scaled with LAI
       canopy_height0 = merge( (param(3) * LAILUT(7,iMon) / maxLAI), canopy_height0, LCover0 == 3)  ! pervious

       ! estimation of the aerodynamic resistance on the lower level
       ! see FAO Irrigation and Draingae Paper No. 56 (p. 19 ff) for more information
       zm     = WindMeasHeight
       ! correction: if wind measurement height is below canopy height loagarithm becomes negative
       zm = merge(canopy_height0 + zm, zm, ((abs(zm - nodata) .GT. eps_dp) .AND. (zm .LT. canopy_height0)))
       !
       ! zh       = zm
       displace = param(4) * canopy_height0 
       zm_zero  = param(5) * canopy_height0 
       zh_zero  = param(6) * zm_zero
       !
       ! calculate aerodynamic resistance (changes monthly)
       aerodyn_resistance0(:,iMon) = log((zm - displace)/zm_zero) * log((zm - displace)/zh_zero)  / (karman**2.0_dp)
       aerodyn_resistance1(:,iMon) = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, aerodyn_resistance0(:,iMon))
       !
    end do
    !
  end subroutine aerodynamical_resistance


  ! ----------------------------------------------------------------------------

  !      NAME
  !        bulksurface_resistance

  !>       \brief Regionalization of bulk surface resistance

  !>       \details estimation of bulk surface resistance
  !>                Global parameters needed (see mhm_parameter.nml):\n
  !>                - param(1) = stomatal_resistance  \n

  !      INTENT(IN)
  !>       \param[in] "integer(i4)  :: LCover_LAI0(:)" - land cover id for LAI at level 0
  !>       \param[in] "real(dp)     :: LAILUT(:)"      - LUT of LAi values
  !>       \param[in] "integer(i4)  :: LAIUnitList(:)" - List of ids of each LAI class in LAILUT
  !>       \param[in] "real(dp)     :: param"          - global parameter
  !>       \param[in] "logical      :: mask0(:,:)"     - mask at level 0 field
  !>       \param[in] "real(dp)     :: nodata"         - nodata value 
  !>       \param[in] "integer(i4)  :: cell_id0 (:)"   - Cell ids at level 0
  !>       \param[in] "integer(i4)  :: nL0_in_L1 (:)"  - Number of L0 cells within a L1 cell
  !>       \param[in] "integer(i4)  :: Upp_row_L1(:)"  - Upper row of high resolution block
  !>       \param[in] "integer(i4)  :: Low_row_L1(:)"  - Lower row of high resolution block
  !>       \param[in] "integer(i4)  :: Lef_col_L1(:)"  - Left column of high resolution block
  !>       \param[in] "integer(i4)  :: Rig_col_L1(:)"  - Right column of high resolution block

  !     INTENT(INOUT)
  !        None

  !     INTENT(OUT)
  !>       \param[out] "real(dp)    :: bulksurface_resistance1(:)" - [s m-1] bulk surface resistance

  !     INTENT(IN), OPTIONAL
  !        None

  !     INTENT(INOUT), OPTIONAL
  !        None

  !     INTENT(OUT), OPTIONAL
  !        None

  !     RETURN
  !        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !>        McMahon et al, 2013: Estimating actual, potential, reference crop and pan evaporation using standard 
  !>        meteorological data: a pragmatic synthesis , HESS

  !     HISTORY
  !>       \author Matthias Zink
  !>       \date   Apr 2013

  subroutine bulksurface_resistance( &
       LCover_LAI0,                  & ! land cover id for LAI at level 0
       LAILUT,                       & ! look up table for LAI
       LAIUnitList,                  & ! List of ids of each LAI class in LAILUT
       param,                        & ! parameter values (size=1)
       mask0,                        & ! mask at level 0
       nodata,                       & ! given nodata value
       cell_id0,                     & ! cell id at Level 0
       nL0_in_L1,                    & ! number of l0 cells within a l1 cell
       Upp_row_L1,                   & ! upper row of a l1 cell in l0 grid
       Low_row_L1,                   & ! lower row of a l1 cell in l0 grid
       Lef_col_L1,                   & ! left col of a l1 cell in l0 grid
       Rig_col_L1,                   & ! right col of a l1 cell in l0 grid
       bulksurface_resistance1       & ! bulk surface resistance
       )

    use mo_upscaling_operators, only: upscale_arithmetic_mean
    use mo_mhm_constants,       only: YearMonths_i4, LAI_factor_surfResi, LAI_offset_surfResi, max_surfResist
    use mo_constants,           only: eps_dp

    implicit none

    integer(i4), dimension(:),   intent(in)  :: LCover_LAI0 ! land cover id for LAI
    real(dp),    dimension(:,:), intent(in)  :: LAILUT      ! look up table for LAI
    !                                                       ! dim1=land cover class, dim2=month of year
    integer(i4), dimension(:),   intent(in)  :: LAIUnitList ! List of ids of each LAI class in LAILUT
    real(dp),                    intent(in)  :: param       ! input parameter
    logical,     dimension(:,:), intent(in)  :: mask0       ! mask at level 0
    real(dp),                    intent(in)  :: nodata      ! given nodata value
    integer(i4), dimension(:),   intent(in)  :: cell_id0    ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in)  :: nL0_in_L1   ! number of l0 cells within a l1 cell
    integer(i4), dimension(:),   intent(in)  :: Upp_row_L1  ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Low_row_L1  ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Lef_col_L1  ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Rig_col_L1  ! right col of a l1 cell in l0 grid
    ! Output
    real(dp),    dimension(:,:), intent(out) :: bulksurface_resistance1

    ! local
    integer(i4)                            :: iMon, ll
    real(dp), dimension(:,:), allocatable  :: leafarea0
    real(dp), dimension(:,:), allocatable  :: bulksurface_resistance0    ! dim 1 = number of cells on level 0,
    !                                                                    ! dim 2 = number of months in year (12)

    ! initialize some things
    allocate(bulksurface_resistance0 (size(LCover_LAI0, dim=1), YearMonths_i4)) ; bulksurface_resistance0 = nodata
    allocate(leafarea0               (size(LCover_LAI0, dim=1), YearMonths_i4)) ; leafarea0               = nodata
    bulksurface_resistance1 = nodata
    !
    do iMon = 1, YearMonths_i4

       ! determine LAIs 
       do ll = 1, size(LAILUT, dim=1)
          leafarea0(:,iMon) = merge( LAILUT(ll, iMon),  leafarea0(:,iMon), LCover_LAI0(:) .EQ. LAIUnitList(ll))
       end do
       ! correction for 0 LAI values
       leafarea0(:,iMon) = merge( 1.00E-10_dp,  leafarea0(:,iMon), leafarea0(:,iMon) .LT. eps_dp)

       bulksurface_resistance0(:,iMon) = param / (  leafarea0(:,iMon) / &
            (LAI_factor_surfResi * leafarea0(:,iMon) + LAI_offset_surfResi))
       ! efeective LAI from McMahon et al ,2013 , HESS supplements

       ! since LAI may be very low, rs becomes very high 
       ! thus the values are restricted to maximum literaure values (i.e. McMahon et al ,2013 , HESS)
       bulksurface_resistance0(:,iMon) = merge(max_surfResist, bulksurface_resistance0(:,iMon), &
            bulksurface_resistance0(:,iMon) .GT. max_surfResist)

       bulksurface_resistance1(:,iMon) = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, bulksurface_resistance0(:,iMon))
    end do

  end subroutine bulksurface_resistance


  ! ----------------------------------------------------------------------------

  !      NAME
  !        priestley_taylor_alpha

  !>       \brief Regionalization of priestley taylor alpha

  !>       \details estimation of priestley taylor alpha
  !>                Global parameters needed (see mhm_parameter.nml):\n
  !>                   - param(1) = PriestleyTaylorCoeff    \n
  !>                   - param(2) = PriestleyTaylorLAIcorr  \n

  !      INTENT(IN)
  !>       \param[in] "integer(i4)  :: LCover_LAI0(:)" - land cover id for LAI at level 0
  !>       \param[in] "real(dp)     :: LAILUT(:)"      - LUT of LAi values
  !>       \param[in] "integer(i4)  :: LAIUnitList(:)" - List of ids of each LAI class in LAILUT
  !>       \param[in] "real(dp)     :: param(:)"       - global parameter
  !>       \param[in] "logical      :: mask0(:,:)"     - mask at level 0 field
  !>       \param[in] "real(dp)     :: nodata"         - nodata value 
  !>       \param[in] "integer(i4)  :: cell_id0 (:)"   - Cell ids at level 0
  !>       \param[in] "integer(i4)  :: nL0_in_L1 (:)"  - Number of L0 cells within a L1 cell
  !>       \param[in] "integer(i4)  :: Upp_row_L1(:)"  - Upper row of high resolution block
  !>       \param[in] "integer(i4)  :: Low_row_L1(:)"  - Lower row of high resolution block
  !>       \param[in] "integer(i4)  :: Lef_col_L1(:)"  - Left column of high resolution block
  !>       \param[in] "integer(i4)  :: Rig_col_L1(:)"  - Right column of high resolution block

  !     INTENT(INOUT)
  !        None

  !     INTENT(OUT)
  !>       \param[out] "real(dp)    :: priestley_taylor_alpha1(:)" - [s m-1] bulk surface resistance

  !     INTENT(IN), OPTIONAL
  !        None

  !     INTENT(INOUT), OPTIONAL
  !        None

  !     INTENT(OUT), OPTIONAL
  !        None

  !     RETURN
  !        None

  !     LITERATURE
  !        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     HISTORY
  !>       \author Matthias Zink
  !>       \date   Apr 2013

  subroutine priestley_taylor_alpha( &
       LCover_LAI0,            & ! land cover id for LAI at level 0
       LAILUT,                 & ! look up table for LAI
       LAIUnitList,            & ! List of ids of each LAI class in LAILUT
       param,                  & ! parameter values (size=2)
       mask0,                  & ! mask at level 0
       nodata,                 & ! given nodata value
       cell_id0,               & ! cell id at Level 0
       nL0_in_L1,              & ! number of l0 cells within a l1 cell
       Upp_row_L1,             & ! upper row of a l1 cell in l0 grid
       Low_row_L1,             & ! lower row of a l1 cell in l0 grid
       Lef_col_L1,             & ! left col of a l1 cell in l0 grid
       Rig_col_L1,             & ! right col of a l1 cell in l0 grid
       priestley_taylor_alpha1 & ! bulk surface resistance
       ) 

    use mo_upscaling_operators, only: upscale_arithmetic_mean
    use mo_mhm_constants,       only: YearMonths_i4
    use mo_constants,           only: eps_dp

    implicit none

    integer(i4), dimension(:),   intent(in)  :: LCover_LAI0 ! land cover id for LAI
    real(dp),    dimension(:,:), intent(in)  :: LAILUT      ! look up table for LAI
    !                                                       ! dim1=land cover class, dim2=month of year
    integer(i4), dimension(:),   intent(in)  :: LAIUnitList ! List of ids of each LAI class in LAILUT
    real(dp),    dimension(:),   intent(in)  :: param       ! input parameter
    logical,     dimension(:,:), intent(in)  :: mask0       ! mask at level 0
    real(dp),                    intent(in)  :: nodata      ! given nodata value
    integer(i4), dimension(:),   intent(in)  :: cell_id0    ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in)  :: nL0_in_L1   ! number of l0 cells within a l1 cell
    integer(i4), dimension(:),   intent(in)  :: Upp_row_L1  ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Low_row_L1  ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Lef_col_L1  ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:),   intent(in)  :: Rig_col_L1  ! right col of a l1 cell in l0 grid
    ! Output
    real(dp),    dimension(:,:), intent(out) :: priestley_taylor_alpha1

    ! local
    integer(i4)                            :: iMon, ll
    real(dp), dimension(:,:), allocatable  :: leafarea0
    real(dp), dimension(:,:), allocatable  :: priestley_taylor_alpha0    ! dim 1 = number of cells on level 0,
    !                                                                    ! dim 2 = number of months in year (12)

    ! initialize some things
    allocate(priestley_taylor_alpha0 (size(LCover_LAI0, dim=1), YearMonths_i4)) ; priestley_taylor_alpha0 = nodata
    allocate(leafarea0               (size(LCover_LAI0, dim=1), YearMonths_i4)) ; leafarea0               = nodata
    priestley_taylor_alpha1 = nodata
    !
    do iMon = 1, YearMonths_i4
       ! determine LAIs per month
       do ll = 1, size(LAILUT, dim=1)
          leafarea0(:,iMon) = merge( LAILUT(ll, iMon),  leafarea0(:,iMon), LCover_LAI0(:) .EQ. LAIUnitList(ll))
       end do
       ! correction for 0 LAI values to avoid numerical instabilities
       leafarea0(:,iMon) = merge( 1.00E-10_dp,  leafarea0(:,iMon), leafarea0(:,iMon) .LT. eps_dp)

       priestley_taylor_alpha0(:,iMon) = param(1) + param(2) * leafarea0(:,iMon)

       priestley_taylor_alpha1(:,iMon) = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, priestley_taylor_alpha0(:,iMon))
    end do

  end subroutine priestley_taylor_alpha


END MODULE mo_multi_param_reg
