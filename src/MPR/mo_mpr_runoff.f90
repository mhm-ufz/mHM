!> \file mo_mpr_runoff.f90
!> \copydoc mo_mpr_runoff

!> \brief multiscale parameter regionalization for runoff generation
!> \details This contains the routine for multiscale parameter regionalization of the runoff parametrization.
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
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

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober Jan 2013 - updated calling sequence for upscaling operators
  ! Stephan Thober Dec 2013 - made header conform with mo_template
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mpr_runoff(LCOVER0, SMs_FC0, slope_emp0, KsVar_H0, param, upscaler, L1_HL1, L1_K0, L1_K1, L1_alpha)
    use mo_grid_scaler, only: scaler_t, up_a_mean
    use mo_message, only: error_message
    real(dp), intent(in) :: param(5) !< Process-specific interflow parameters in registry order.
    real(dp), intent(in) :: SMs_FC0(:) !< Packed fine-cell saturation deficit relative to field capacity.
    real(dp), intent(in) :: slope_emp0(:) !< Empirical fine-cell slope quantiles.
    real(dp), intent(in) :: KsVar_H0(:) !< Packed horizontal conductivity variability.
    integer(i4), intent(in) :: LCOVER0(:) !< Packed fine-cell land-cover classes.
    class(scaler_t), intent(inout) :: upscaler !< Connected fine-to-coarse grid mapping.
    real(dp), intent(out) :: L1_HL1(:) !< Unsaturated storage threshold [mm].
    real(dp), intent(out) :: L1_K0(:) !< Fast interflow recession coefficient [day].
    real(dp), intent(out) :: L1_K1(:) !< Slow interflow recession coefficient [day].
    real(dp), intent(out) :: L1_alpha(:) !< Slow-interflow exponent.
    real(dp), allocatable :: tmp(:)
    integer(i4) :: i

    if (size(SMs_FC0) /= size(LCOVER0) .or. size(slope_emp0) /= size(LCOVER0) .or. size(KsVar_H0) /= size(LCOVER0)) &
      call error_message('MPR: runoff L0 inputs must have matching packed sizes.')
    if (size(L1_HL1) /= size(L1_alpha) .or. size(L1_K0) /= size(L1_alpha) .or. size(L1_K1) /= size(L1_alpha)) &
      call error_message('MPR: runoff L1 outputs must have matching sizes.')
    allocate(tmp(size(LCOVER0)))
    ! Saturation-deficit-based fast interflow threshold.
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i4, size(tmp)
      tmp(i) = param(1) * SMs_FC0(i)
    end do
    !$omp end parallel do
    call upscaler%execute(tmp, L1_HL1, upscaling_operator=up_a_mean)

    ! Slope-dependent fast interflow, with the forest surface-resistance factor (Booij et al., 2005).
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i4, size(tmp)
      tmp(i) = param(2) * (2.0_dp - slope_emp0(i))
      if (LCOVER0(i) == 1_i4) tmp(i) = tmp(i) * param(3)
    end do
    !$omp end parallel do
    call upscaler%execute(tmp, L1_K0, upscaling_operator=up_a_mean)
    L1_K0 = merge(1.0_dp, L1_K0, L1_K0 < 1.0_dp)

    !$omp parallel do default(shared) schedule(static)
    do i = 1_i4, size(tmp)
      tmp(i) = param(2) * (2.0_dp - slope_emp0(i)) + param(4) * (1.0_dp + KsVar_H0(i))
    end do
    !$omp end parallel do
    call upscaler%execute(tmp, L1_K1, upscaling_operator=up_a_mean)
    L1_K1 = merge(2.0_dp, L1_K1, L1_K1 < 2.0_dp)

    !$omp parallel do default(shared) schedule(static)
    do i = 1_i4, size(tmp)
      tmp(i) = param(5) * (1.0_dp / KsVar_H0(i)) * (1.0_dp / (1.0_dp + SMs_FC0(i)))
    end do
    !$omp end parallel do
    call upscaler%execute(tmp, L1_alpha, upscaling_operator=up_a_mean)
    L1_K0 = merge(L1_K1, L1_K0, L1_K0 > L1_K1)
  end subroutine mpr_runoff

end module mo_mpr_runoff
