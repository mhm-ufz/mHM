!> \file mo_mpr_runoff.f90

!> \brief multiscale parameter regionalization for runoff generation

!> \details This contains the routine for multiscale parameter regionalization of the runoff parametrization.

!> \author Stephan Thober, Rohini Kumar
!> \date Dec 2012

module mo_mpr_runoff

  use mo_kind, only: i4, dp

  implicit none

  private

  public :: mpr_runoff

  ! ----------------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------

  !      NAME
  !         mpr_runoff

  !>        \brief multiscale parameter regionalization for runoff parameters

  !>        \details Perform the multiscale parameter regionalization for runoff
  !>        global parameters (see mhm_parameter.nml). These are the following five
  !>        parameters:\n
  !>           - param(1) = interflowStorageCapacityFactor \n
  !>           - param(2) = interflowRecession_slope       \n
  !>           - param(3) = fastInterflowRecession_forest  \n
  !>           - param(4) = slowInterflowRecession_Ks      \n
  !>           - param(5) = exponentSlowInterflow          \n

  !     INTENT(IN)
  !>      \param[in] "real(dp), dimension(5) :: param"         global parameters
  !>      \param[in] "real(dp)               :: nodata"        nodata value
  !>      \param[in] "real(dp)               :: SMs_FC0(:)"    [-] soil mositure deficit from field
  !>                                                           capacity w.r.t to saturation
  !>      \param[in] "real(dp)               :: slope_emp0(:)" empirical quantile values F(slope)
  !>      \param[in] "real(dp)               :: KsVar_H0(:)"   [-] relative variability of saturated
  !>                                                         hydraulic counductivity for
  !>                                                         Horizantal flow at level 0
  !>      \param[in] "integer(i4)            :: LCOVER0(:)"    land cover at level 0
  !>      \param[in] "logical                :: mask0(:,:)"    mask at Level 0
  !>      \param[in] "integer(i4)            :: cell_id0(:)"   Cell ids of hi res field
  !>      \param[in] "integer(i4)            :: upp_row_L1(:)" Upper row of hi res block
  !>      \param[in] "integer(i4)            :: low_row_L1(:)" Lower row of hi res block
  !>      \param[in] "integer(i4)            :: lef_col_L1(:)" Left column of hi res block
  !>      \param[in] "integer(i4)            :: rig_col_L1(:)" Right column of hi res block
  !>      \param[in] "integer(i4)            :: nL0_in_L1(:)"  Number of L0 cells within a L1 cell
  !>      \param[in] "real(dp)               :: c2TSTu"        unit transformations      

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>      \param[out] "real(dp)              :: L1_HL1(:)"    [10^-3 m] Threshhold water depth
  !>                                                         in upper reservoir (for Runoff contribution)
  !>      \param[out] "real(dp)              :: L1_K0(:)"     [10^-3 m] Recession coefficient
  !>                                                         of the upper reservoir, upper outlet
  !>      \param[out] "real(dp)              :: L1_K1(:)"     [10^-3 m] Recession coefficient
  !>                                                         of the upper reservoir, lower outlet
  !>      \param[out] "real(dp)              :: L1_alpha(:)"  [1] Exponent for the upper reservoir

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !        None

  !     RESTRICTIONS
  !        None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !      HISTORY
  !         Modified, Stephan Thober, Jan 2013 - updated calling sequence for upscaling operators
  !                   Stephan Thober, Dec 2013 - made header conform with mo_template

  subroutine mpr_runoff( &
       ! Input -------------------------------------------------------
       LCOVER0,    &
       mask0,      &
       nodata,     &
       SMs_FC0,    &
       slope_emp0, &
       KsVar_H0,   &
       param,      &
       cell_id0,   &
       upp_row_L1, &
       low_row_L1, &
       lef_col_L1, &
       rig_col_L1, &
       nL0_in_L1,  &
       c2TSTu,     &
       ! Output ------------------------------------------------------
       L1_HL1,     &
       L1_K0,      &
       L1_K1,      &
       L1_alpha )

    use mo_upscaling_operators, only: upscale_arithmetic_mean

    implicit none

    ! Input
    real(dp),    dimension(5), intent(in)   :: param      ! global parameter
    real(dp),                  intent(in)   :: nodata     ! nodata value
    real(dp),    dimension(:), intent(in)   :: SMs_FC0    ! [-] soil mositure deficit from field
                                                          ! capacity w.r.t to saturation
    real(dp),    dimension(:), intent(in)   :: slope_emp0 ! empirical quantile values F(slope)
    real(dp),    dimension(:), intent(in)   :: KsVar_H0   ! [-] relative variability of saturated
                                                          ! hydraulic counductivity for
                                                          ! Horizantal flow at level 0
    integer(i4), dimension(:), intent(in)   :: LCOVER0    ! land cover at level 0

    ! Ids of L0 cells beneath L1 cell
    logical,     dimension(:,:), intent(in) :: mask0      ! mask at Level 0
    integer(i4), dimension(:),   intent(in) :: cell_id0   ! Cell ids of hi res field
    integer(i4), dimension(:),   intent(in) :: upp_row_L1 ! Upper row of hi res block
    integer(i4), dimension(:),   intent(in) :: low_row_L1 ! Lower row of hi res block
    integer(i4), dimension(:),   intent(in) :: lef_col_L1 ! Left column of hi res block
    integer(i4), dimension(:),   intent(in) :: rig_col_L1 ! Right column of hi res block
    integer(i4), dimension(:),   intent(in) :: nL0_in_L1  ! Number of L0 cells within a L1 cell
    real(dp),                    intent(in) :: c2TSTu     ! unit transformations

    ! Output
    real(dp),    dimension(:),   intent(out):: L1_HL1     ! [10^-3 m] Threshhold water depth
                                                          ! in upper reservoir (for Runoff
                                                          ! contribution)
    real(dp),    dimension(:),   intent(out):: L1_K0      ! [10^-3 m] Recession coefficient
                                                          ! of the upper reservoir, upper outlet
    real(dp),    dimension(:),   intent(out):: L1_K1      ! [10^-3 m] Recession coefficient
                                                          ! of the upper reservoir, lower outlet
    real(dp),    dimension(:),   intent(out):: L1_alpha   ! [1] Exponent for the upper reservoir
    
    ! Local Variables
    real(dp), dimension(size(SMs_FC0,1)) :: tmp      ! temporal variable

    !-----------------------------
    ! FAST INTERFLOW
    !-----------------------------
    ! HL1 = f(soil properties; No refrence found)

    ! Based on the saturation deficit from the field capacity status
    ! seems more reasonable and intutative.
    ! NOTE: This value for the sandy soils will have higher value of HL1, as compared to 
    !       to clayey soil and so these soils can hold larger amount of amount.

    tmp    = merge( param(1) * SMs_FC0, nodata, cell_id0 .ne. int(nodata, i4))
    L1_HL1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
             Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, tmp )

    ! 1/K0 = f(terrian slope) [Booij, et. al.(2005), JoH]
    ! Steeper slopes resists (1/K0) fast water flows lesser as 
    ! compared to that on the flater slope areas. 
    ! Assuming that above relationship holds for all kind of land cover classes

    ! In the forested area surface resistance to fast interflow is higher as compared
    ! to the permeable land surface

    tmp   = merge( param(2) * (2.0_dp - slope_emp0), nodata, cell_id0 .ne. int(nodata,i4) )
    tmp   = merge( tmp * param(3), tmp, LCOVER0 .eq. 1 )
    L1_K0 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, tmp )

    ! To avoid numerical error in fully impervious areas (K0 == 0)
    ! minimum value of K0 is 1-day 
    L1_K0 = merge( 1.0_dp, L1_K0, L1_K0 .lt. 1.0_dp )

    ! ------------------------------------------------------------------
    ! SLOW INTERFLOW
    ! ------------------------------------------------------------------
    !   K1 = f(terrian slope, Booij, et. al.(2005), JoH)
    !      = f(soil properties, LC & New modification)  
    !  K1  = K0 + K1(soil-Ks)

    tmp   = merge( param(2) * (2.0_dp - slope_emp0) + param(4) * (1.0_dp + KsVar_H0), &
            nodata, cell_id0 .ne. int(nodata,i4) )
    L1_K1 = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, tmp )

    ! minimum value of K1 is 1-day 
    L1_K1 = merge( 2.0_dp, L1_K1, L1_K1 .lt. 2.0_dp )

 
    ! alpha = f(soil type; variabitity of Ks)
    ! Lower the alpha (exponent of slow interflow) means lower amount of 
    ! water released from the storage to contribute for slow interflow. 
    ! For instance sandy soils will have lower value of alpha as comapred to
    ! the clayey soils.
    ! This assumption is quite realistic in physical sense... 
    tmp = merge( param(5) * (1.0_dp / KsVar_H0) * ( 1.0_dp / ( 1.0_dp + SMs_FC0) ),    &
          nodata, cell_id0 .ne. int(nodata,i4) )
    L1_alpha = upscale_arithmetic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, &
               Rig_col_L1, cell_id0, mask0, nodata, tmp )

    ! constraints and unit transformation
    L1_K0 = merge( L1_K1, L1_K0, L1_K0 .gt. L1_K1 )

    L1_K0 = c2TSTu / L1_K0
    L1_K1 = c2TSTu / L1_K1
    
  end subroutine mpr_runoff

end module mo_mpr_runoff
