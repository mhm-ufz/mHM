!> \file mo_mpr_pet.f90

!> \Dynamic scaling function for PET correction using LAI at level-0

!> \details This module sets up pet correction factor at level-1 based on LAI 

!> \author Mehmet Cuneyd Demirel, Simon Stisen
!> \date May 2017

module mo_mpr_pet

  use mo_kind, only: i4, dp

  implicit none

  PUBLIC :: pet_correctbyLAI          ! estimate PET correction factor with distributed LAI
  PUBLIC :: pet_correctbyASP          ! estimate PET correction factor with distributed Aspect
  PUBLIC :: priestley_taylor_alpha    ! factor (alpha) for Presley-Taylor ET estimation
  !PUBLIC :: aerodynamical_resistance  ! aerodynamical resistance (ra) for Penman-Monteith ET estimation
  PUBLIC :: bulksurface_resistance    ! bulk surface (stomatal) resistance (rs) for Penman-Monteith ET estimation

  private

contains

         
! ----------------------------------------------------------------------------

  !      NAME
  !        pet_correctbyLAI

  !>       \brief estimate PET correction factor based on LAI at L1

  !>       \details estimate PET correction factor based on LAI at L1 for a given 
  !>                Leaf Area Index field. \n
  !>                Global parameters needed (see mhm_parameter.nml):\n
  !>                Process Case 5:\n
  !>                   - param(1) = PET_a_forest \n
  !>                   - param(2) = PET_a_impervious \n
  !>                   - param(3) = PET_a_pervious \n
  !>                   - param(4) = PET_b \n
  !>                   - param(5) = PET_c \n
  
  !>     Example DSF=PET_a+PET_b*(1-exp(PET_c*LAI))

  !>     Similar to the crop coefficient concept Kc=a+b*(1-exp(c*LAI)) by Allen, R. G., L. S. Pereira, 
  !>     D. Raes, and M. Smith (1998), Crop evapotranspiration - Guidelines for computing crop water requirements., FAO
  !>     Irrigation and drainage paper 56. See Chapter 9, Equation 97  <http://www.fao.org/docrep/X0490E/x0490e0f.htm>

  !>     Date: 17/5/2017

  !      INTENT(IN)
  !>       \param[in] "integer(i4)  :: cell_id0(:)"                  - Cell ids at level 0
  !>       \param[in] "real(dp)     :: LCOVER0(:)"                   - Landcover at level-0
  !>       \param[in] "real(dp)     :: LAI0(:)"                      - LAI at level-0
  !>       \param[in] "real(dp)     :: param(:)"                     - array of global parameters
  !>       \param[in] "real(dp)     :: nodata"                       - nodata value

  !>        \param[in] "integer(i4) :: L0_cell_id(:,:)"  - cell ids of high resolution field, 
  !>                                                       Number of rows times Number of columns
  !>                                                       of high resolution field
  !>        \param[in] "integer(i4) :: upp_row_L1(:)"    - Upper row id in high resolution field
  !>                                                       (L0) of low resolution cell (L1 cell)
  !>        \param[in] "integer(i4) :: low_row_L1(:)"    - Lower row id in high resolution field
  !>                                                       (L0) of low resolution cell (L1 cell)
  !>        \param[in] "integer(i4) :: lef_col_L1(:)"    - Left column id in high resolution 
  !>                                                       field (L0) of low resolution cell
  !>        \param[in] "integer(i4) :: rig_col_L1(:)"    - Right column id in high resolution
  !>                                                       field (L0) of low resolution cell
  !>        \param[in] "integer(i4) :: nL0_in_L1(:)"     - Number of high resolution cells (L0)
  !>                                                       in low resolution cell (L1 cell)

  !      INTENT(OUT)
  !>       \param[in,out] "real(dp) :: L1_petLAIcorFactor(:)"   - PET correction using LAI.


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
  !         call pet_correctbyLAI(param,nodata,LCOVER0,LAI0,mask0,cell_id0,upp_row_L1,low_row_L1, &
  !                               lef_col_L1,rig_col_L1,nL0_in_L1,L1_petLAIcorFactor)
  
  !      HISTORY
  !>         \author M. Cuneyd Demirel and Simon Stisen from GEUS.dk
  !>         \date May. 2017

  subroutine pet_correctbyLAI( &
       ! Input -----------------------------------------------------------------
       param         , & ! global parameters, three are required
       nodata        , & ! no data value
       LCOVER0       , & ! land use cover at L0
       LAI0          , & ! LAI at L0
       mask0         , & ! mask at L0
       cell_id0      , & ! cell ids at L0
       upp_row_L1    , & ! upper row of L0 block within L1 cell
       low_row_L1    , & ! lower row of L0 block within L1 cell
       lef_col_L1    , & ! left column of L0 block within L1 cell
       rig_col_L1    , & ! right column of L0 block within L1 cell
       nL0_in_L1     , & ! Number of L0 cells in L0 block within L1 cell
       ! Output ----------------------------------------------------------------
       L1_petLAIcorFactor )       ! fraction of roots in soil horizons

    use mo_upscaling_operators, only: upscale_harmonic_mean
    !$  use omp_lib

    implicit none

    ! Input
    real(dp),    dimension(5),     intent(in) :: param         ! parameters
    real(dp),                      intent(in) :: nodata        ! no data value
    integer(i4), dimension(:),     intent(in) :: LCOVER0       ! Land cover at level 0
    real(dp),    dimension(:),     intent(in) :: LAI0          ! LAI at level-0
                                                         
    ! Ids of L0 cells beneath L1 cell
    logical,     dimension(:,:), intent(in)   :: mask0
    integer(i4), dimension(:),   intent(in)   :: cell_id0   ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in)   :: upp_row_L1 ! Upper row of hi res block
    integer(i4), dimension(:),   intent(in)   :: low_row_L1 ! Lower row of hi res block
    integer(i4), dimension(:),   intent(in)   :: lef_col_L1 ! Left column of hi res block
    integer(i4), dimension(:),   intent(in)   :: rig_col_L1 ! Right column of hi res block
    integer(i4), dimension(:),   intent(in)   :: nL0_in_L1  ! Number of L0 cells within a L1 cel

    ! Output
    ! The following five variables have the dimension: Number of cells at L1 times nHorizons_mHM
    real(dp),   dimension(:), intent(inout) :: L1_petLAIcorFactor ! pet cor factor at level-1

    ! Local Variables
    real(dp), dimension(size(LCOVER0,1))     :: petLAIcorFactor_0     ! pet cor factor at level-0

    ! local variables
    integer(i4)                             :: kk         ! loop index
    integer(i4)                             :: LL         ! loop index   
    
    ! ------------------------------------------------------------------
    ! Estimate DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    ! ------------------------------------------------------------------

    !$OMP PARALLEL
    !$OMP DO PRIVATE( LL ) SCHEDULE( STATIC )
    
    ! need to be done for every landcover to get DSF 
    do kk = 1, size(LCOVER0, 1)

       LL = LCOVER0(kk)
       
       select case(LL)
       case(1) ! forest
          petLAIcorFactor_0(kk) = param(1)+(param(4)*(1.0_dp-exp(param(5)*LAI0(kk))))
       case(2) ! impervious
          petLAIcorFactor_0(kk) = param(2)+(param(4)*(1.0_dp-exp(param(5)*LAI0(kk))))
       case(3) ! permeable 
          petLAIcorFactor_0(kk) = param(3)+(param(4)*(1.0_dp-exp(param(5)*LAI0(kk))))
       end select

    end do
    
    !$OMP END DO
    
    L1_petLAIcorFactor(:) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
         Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, petLAIcorFactor_0 )
    
    !$OMP END PARALLEL

  end subroutine pet_correctbyLAI

  ! ----------------------------------------------------------------------------

  !      NAME
  !         PET correction due to aspect

  !>        \brief correction of PET

  !>        \details Correction of PET based on L0 aspect data. \n
  !> 
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
  !                  Matthias Zink,  Jun 2017 - renamed and moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90

  subroutine pet_correctbyASP( Id0, latitude_l0, Asp0, param, nodata, fAsp0 )

    !$ use omp_lib
    
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

  end subroutine pet_correctbyASP

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
  !        Modified    Matthias Zink,   Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90

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
  !        Modified    Matthias Zink,   Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90

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
  
end module mo_mpr_pet
