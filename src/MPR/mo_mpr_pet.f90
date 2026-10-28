!> \file mo_mpr_pet.f90
!> \copydoc mo_mpr_pet

!> \brief MPR routine for PET.
!> \details This module sets up pet correction factor at level-1 based on LAI
!> \authors Mehmet Cuneyd Demirel, Simon Stisen
!> \date May 2017
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_pet

  use mo_kind, only : i4, dp

  implicit none

  PUBLIC :: pet_correctbyASP          ! estimate PET correction factor with distributed Aspect
  PUBLIC :: priestley_taylor_alpha    ! factor (alpha) for Presley-Taylor ET estimation
  PUBLIC :: bulksurface_resistance    ! bulk surface (stomatal) resistance (rs) for Penman-Monteith ET estimation

  private

contains

  ! ----------------------------------------------------------------------------

  !    NAME
  !        pet_correctbyASP

  !    PURPOSE
  !>       \brief correction of PET

  !>       \details Correction of PET based on L0 aspect data.

  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = minCorrectionFactorPET
  !>       - param(2) = delta added to minCorrectionFactorPET
  !>       - param(3) = aspectTresholdPET

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Juliane Mai    Oct 2013 - OLD parametrization  --> param(1) = minCorrectionFactorPET
  !                                                --> param(2) = maxCorrectionFactorPET
  !                                                --> param(3) = aspectTresholdPET
  !                                             -------------------------------
  !                                             maxCorrectionFactorPET = minCorrectionFactorPET + delta
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = minCorrectionFactorPET
  !                                                --> param(2) = delta
  !                                                --> param(3) = aspectTresholdPET
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Stephan Thober Sep 2015 - Mapping L1 to Lo, latitude on L0
  ! Luis Samaniego Sep 2015 - PET correction on the southern hemisphere
  ! Matthias Zink  Jun 2017 - renamed and moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine pet_correctbyASP(latitude_l0, Asp0, param, fAsp0)
    real(dp), intent(in) :: latitude_l0(:) !< Physical latitude of packed fine cells [degree north].
    real(dp), intent(in) :: Asp0(:) !< Fine-cell aspect [degree].
    real(dp), intent(in) :: param(3) !< Minimum correction, correction delta, and aspect threshold.
    real(dp), intent(out) :: fAsp0(:) !< Packed fine-cell PET correction.
    integer(i4) :: i
    real(dp) :: maximum, north, south

    maximum = param(1) + param(2)
    !$omp parallel do default(shared) private(north,south) schedule(static)
    do i = 1_i4, size(fAsp0)
      north = merge(param(1) + (maximum - param(1)) / param(3) * Asp0(i), &
        param(1) + (maximum - param(1)) / (360.0_dp - param(3)) * (360.0_dp - Asp0(i)), Asp0(i) < param(3))
      south = merge(param(1) + (maximum - param(1)) / (360.0_dp - param(3)) * (360.0_dp - Asp0(i)), &
        param(1) + (maximum - param(1)) / param(3) * Asp0(i), Asp0(i) < param(3))
      fAsp0(i) = merge(north, south, latitude_l0(i) > 0.0_dp)
    end do
    !$omp end parallel do
  end subroutine pet_correctbyASP

  ! ----------------------------------------------------------------------------

  !    NAME
  !        priestley_taylor_alpha

  !    PURPOSE
  !>       \brief Regionalization of priestley taylor alpha

  !>       \details estimation of priestley taylor alpha
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = PriestleyTaylorCoeff
  !>       - param(2) = PriestleyTaylorLAIcorr

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine priestley_taylor_alpha(LAI0, param, upscaler, priestley_taylor_alpha1)
    use mo_grid_scaler, only: scaler_t, up_a_mean
    use mo_message, only: error_message
    real(dp), intent(in) :: LAI0(:, :) !< Packed fine-cell LAI by temporal slice.
    real(dp), intent(in) :: param(2) !< Priestley-Taylor intercept and LAI coefficient.
    class(scaler_t), intent(inout) :: upscaler !< Connected fine-to-coarse grid mapping.
    real(dp), intent(out) :: priestley_taylor_alpha1(:, :) !< Coarse-cell coefficients by temporal slice.
    real(dp), allocatable :: alpha0(:)
    integer(i4) :: tt, i

    if (size(priestley_taylor_alpha1, 2) /= size(LAI0, 2)) &
      call error_message('MPR: priestley_taylor_alpha output LAI dimension does not match the input.')
    allocate(alpha0(size(LAI0, 1)))
    do tt = 1_i4, size(LAI0, 2)
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i4, size(alpha0)
        alpha0(i) = param(1) + param(2) * LAI0(i, tt)
      end do
      !$omp end parallel do
      call upscaler%execute(alpha0, priestley_taylor_alpha1(:, tt), upscaling_operator=up_a_mean)
    end do
  end subroutine priestley_taylor_alpha

  ! ----------------------------------------------------------------------------

  !    NAME
  !        bulksurface_resistance

  !    PURPOSE
  !>       \brief Regionalization of bulk surface resistance

  !>       \details estimation of bulk surface resistance
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = stomatal_resistance

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine bulksurface_resistance(LAI0, param, upscaler, bulksurface_resistance1)
    use mo_mpr_constants, only: LAI_factor_surfResi, LAI_offset_surfResi, max_surfResist
    use mo_grid_scaler, only: scaler_t, up_a_mean
    use mo_message, only: error_message
    real(dp), intent(in) :: LAI0(:, :) !< Packed fine-cell LAI by temporal slice.
    real(dp), intent(in) :: param !< Bulk surface resistance parameter.
    class(scaler_t), intent(inout) :: upscaler !< Connected fine-to-coarse grid mapping.
    real(dp), intent(out) :: bulksurface_resistance1(:, :) !< Coarse-cell resistance by temporal slice.
    real(dp), allocatable :: resistance0(:)
    integer(i4) :: tt, i

    if (size(bulksurface_resistance1, 2) /= size(LAI0, 2)) &
      call error_message('MPR: bulksurface_resistance output LAI dimension does not match the input.')
    allocate(resistance0(size(LAI0, 1)))
    do tt = 1_i4, size(LAI0, 2)
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i4, size(resistance0)
        resistance0(i) = param / (LAI0(i, tt) / (LAI_factor_surfResi * LAI0(i, tt) + LAI_offset_surfResi))
        ! Effective LAI and resistance ceiling follow McMahon et al. (2013), HESS supplements.
        if (resistance0(i) > max_surfResist) resistance0(i) = max_surfResist
      end do
      !$omp end parallel do
      call upscaler%execute(resistance0, bulksurface_resistance1(:, tt), upscaling_operator=up_a_mean)
    end do
  end subroutine bulksurface_resistance

end module mo_mpr_pet
