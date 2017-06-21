!> \file mo_mpr_pet.f90

!> \Dynamic scaling function for PET correction using LAI at level-0

!> \details This module sets up pet correction factor at level-1 based on LAI 

!> \author Mehmet Cuneyd Demirel, Simon Stisen
!> \date May 2017

module mo_mpr_pet

  use mo_kind, only: i4, dp

  implicit none

  PUBLIC :: pet_correctbyLAI          ! estimate PET correction factor with distributed LAI

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
  !>     D. Raes, and M. Smith (1998), Crop evapotranspiration - Guidelines for computing crop water requirements., 
  !>  FAO Irrigation and drainage paper 56. See Chapter 9, Equation 97  <http://www.fao.org/docrep/X0490E/x0490e0f.htm>

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
  !         call pet_correctbyLAI(param,nodata,LCOVER0,LAI0,mask0,cell_id0,upp_row_L1,low_row_L1,lef_col_L1,rig_col_L1,nL0_in_L1,L1_petLAIcorFactor)
  
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

    use mo_upscaling_operators, only: upscale_harmonic_mean,upscale_arithmetic_mean,upscale_geometric_mean
    use mo_message,             only: message
    use mo_string_utils,        only: num2str

    
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

            !L1_petLAIcorFactor(:) = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            !                    Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, petLAIcorFactor_0 )            
            
            L1_petLAIcorFactor(:) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
                              Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, petLAIcorFactor_0 )


            !L1_petLAIcorFactor(:) = upscale_geometric_mean( Upp_row_L1, Low_row_L1, &
            !                  Lef_col_L1, Rig_col_L1, mask0, nodata, petLAIcorFactor_0 )        
         

   
            !$OMP END PARALLEL

  end subroutine pet_correctbyLAI

end module mo_mpr_pet
