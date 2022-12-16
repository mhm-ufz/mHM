!> \file mo_mpr_runoff.f90
!> \brief \copybrief mo_mpr_runoff
!> \details \copydetails mo_mpr_runoff

!> \brief multiscale parameter regionalization for runoff generation
!> \details This contains the routine for multiscale parameter regionalization of the runoff parametrization.
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_runoff

  use mo_kind, only : i4, dp

  implicit none

  private

  public :: mpr_runoff

  ! ----------------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------

  !    NAME
  !        mpr_runoff

  !    PURPOSE
  !>       \brief multiscale parameter regionalization for runoff parameters

  !>       \details Perform the multiscale parameter regionalization for runoff
  !>       global parameters (see mhm_parameter.nml). These are the following five
  !>       parameters:
  !>       - param(1) = interflowStorageCapacityFactor
  !>       - param(2) = interflowRecession_slope
  !>       - param(3) = fastInterflowRecession_forest
  !>       - param(4) = slowInterflowRecession_Ks
  !>       - param(5) = exponentSlowInterflow

  !    INTENT(IN)
  !>       \param[in] "integer(i4), dimension(:) :: LCOVER0"    land cover at level 0
  !>       \param[in] "logical, dimension(:, :) :: mask0"       mask at Level 0
  !>       \param[in] "real(dp), dimension(:) :: SMs_FC0"       [-] soil mositure deficit from field
  !>       \param[in] "real(dp), dimension(:) :: slope_emp0"    empirical quantile values F(slope)
  !>       \param[in] "real(dp), dimension(:) :: KsVar_H0"      [-] relative variability of saturated
  !>       \param[in] "real(dp), dimension(5) :: param"         global parameters
  !>       \param[in] "integer(i4), dimension(:) :: cell_id0"   Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: upp_row_L1" Upper row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: low_row_L1" Lower row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: lef_col_L1" Left column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: rig_col_L1" Right column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: nL0_in_L1"  Number of L0 cells within a L1 cell
  !>       \param[in] "real(dp) :: c2TSTu"                      unit transformations

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: L1_HL1"   [10^-3 m] Threshhold water depth
  !>       \param[out] "real(dp), dimension(:) :: L1_K0"    [10^-3 m] Recession coefficient
  !>       \param[out] "real(dp), dimension(:) :: L1_K1"    [10^-3 m] Recession coefficient
  !>       \param[out] "real(dp), dimension(:) :: L1_alpha" [1] Exponent for the upper reservoir

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober Jan 2013 - updated calling sequence for upscaling operators
  ! Stephan Thober Dec 2013 - made header conform with mo_template
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mpr_runoff(LCOVER0, mask0, SMs_FC0, slope_emp0, KsVar_H0, param, cell_id0, upp_row_L1, low_row_L1, &
                       lef_col_L1, rig_col_L1, nL0_in_L1, L1_HL1, L1_K0, L1_K1, L1_alpha)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    ! global parameters
    real(dp), dimension(5), intent(in) :: param

    ! [-] soil mositure deficit from field
    real(dp), dimension(:), intent(in) :: SMs_FC0

    ! empirical quantile values F(slope)
    real(dp), dimension(:), intent(in) :: slope_emp0

    ! [-] relative variability of saturated
    real(dp), dimension(:), intent(in) :: KsVar_H0

    ! land cover at level 0
    integer(i4), dimension(:), intent(in) :: LCOVER0

    ! mask at Level 0
    logical, dimension(:, :), intent(in) :: mask0

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: cell_id0

    ! Upper row of hi res block
    integer(i4), dimension(:), intent(in) :: upp_row_L1

    ! Lower row of hi res block
    integer(i4), dimension(:), intent(in) :: low_row_L1

    ! Left column of hi res block
    integer(i4), dimension(:), intent(in) :: lef_col_L1

    ! Right column of hi res block
    integer(i4), dimension(:), intent(in) :: rig_col_L1

    ! Number of L0 cells within a L1 cell
    integer(i4), dimension(:), intent(in) :: nL0_in_L1

    ! [10^-3 m] Threshhold water depth
    real(dp), dimension(:), intent(out) :: L1_HL1

    ! [10^-3 m] Recession coefficient
    real(dp), dimension(:), intent(out) :: L1_K0

    ! [10^-3 m] Recession coefficient
    real(dp), dimension(:), intent(out) :: L1_K1

    ! [1] Exponent for the upper reservoir
    real(dp), dimension(:), intent(out) :: L1_alpha

    ! temporal variable
    real(dp), dimension(size(SMs_FC0, 1)) :: tmp


    !-----------------------------
    ! FAST INTERFLOW
    !-----------------------------
    ! HL1 = f(soil properties; No reference found)
    ! Based on the saturation deficit from the field capacity status
    ! seems more reasonable and intutative.
    ! NOTE: This value for the sandy soils will have higher value of HL1, as compared to
    !       to clayey soil and so these soils can hold larger amount of amount.
    tmp = merge(param(1) * SMs_FC0, nodata_dp, cell_id0 .ne. nodata_i4)
    L1_HL1 = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, tmp)

    ! 1/K0 = f(terrain slope) [Booij, et. al.(2005), JoH]
    ! Steeper slopes resists (1/K0) fast water flows lesser as
    ! compared to that on the flater slope areas.
    ! Assuming that above relationship holds for all kind of land cover classes

    ! In the forested area surface resistance to fast interflow is higher as compared
    ! to the permeable land surface

    tmp = merge(param(2) * (2.0_dp - slope_emp0), nodata_dp, cell_id0 .ne. nodata_i4)
    tmp = merge(tmp * param(3), tmp, LCOVER0 .eq. 1)
    L1_K0 = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, tmp)

    ! To avoid numerical error in fully impervious areas (K0 == 0)
    ! minimum value of K0 is 1-day
    L1_K0 = merge(1.0_dp, L1_K0, L1_K0 .lt. 1.0_dp)

    ! ------------------------------------------------------------------
    ! SLOW INTERFLOW
    ! ------------------------------------------------------------------
    !   K1 = f(terrian slope, Booij, et. al.(2005), JoH)
    !      = f(soil properties, LC & New modification)
    !  K1  = K0 + K1(soil-Ks)

    tmp = merge(param(2) * (2.0_dp - slope_emp0) + param(4) * (1.0_dp + KsVar_H0), &
            nodata_dp, cell_id0 .ne. nodata_i4)
    L1_K1 = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
            Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, tmp)

    ! minimum value of K1 is 1-day
    L1_K1 = merge(2.0_dp, L1_K1, L1_K1 .lt. 2.0_dp)


    ! alpha = f(soil type; variabitity of Ks)
    ! Lower the alpha (exponent of slow interflow) means lower amount of
    ! water released from the storage to contribute for slow interflow.
    ! For instance sandy soils will have lower value of alpha as comapred to
    ! the clayey soils.
    ! This assumption is quite realistic in physical sense...
    tmp = merge(param(5) * (1.0_dp / KsVar_H0) * (1.0_dp / (1.0_dp + SMs_FC0)), &
            nodata_dp, cell_id0 .ne. nodata_i4)
    L1_alpha = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, &
            Rig_col_L1, cell_id0, mask0, nodata_dp, tmp)

    ! constraints and unit transformation
    L1_K0 = merge(L1_K1, L1_K0, L1_K0 .gt. L1_K1)


  end subroutine mpr_runoff

end module mo_mpr_runoff
