!> \file mo_mpr_aerodynamics.f90
!> \copydoc mo_mpr_aerodynamics

!> \brief Regionalize aerodynamic resistance from LAI and land cover.
!> \details Extracted from mo_multi_param_reg without its global orchestration dependencies.
!> \changelog
!! - Matthias Zink, Apr 2013: implemented aerodynamic resistance regionalization.
!! - Matthias Zink, Jun 2017: moved the PET regionalization helper.
!! - Robert Schweppe, Jun 2018: refactored and reformatted the routine.
!> \authors Matthias Zink
!> \authors Robert Schweppe
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_aerodynamics
  use mo_kind, only: i4, dp
  implicit none
  private
  public :: aerodynamical_resistance
contains
  !> \brief Regionalize aerodynamic resistance from LAI and land cover.
  subroutine aerodynamical_resistance(LAI0, LCover0, param, upscaler, aerodyn_resistance1)
    use mo_mpr_constants, only: WindMeasHeight, karman
    use mo_grid_scaler, only: scaler_t, up_a_mean
    use mo_message, only: error_message
    real(dp), intent(in) :: LAI0(:, :) !< Packed fine-cell LAI by temporal slice.
    integer(i4), intent(in) :: LCover0(:) !< Packed fine-cell land-cover class (forest, impervious, pervious).
    real(dp), intent(in) :: param(6) !< Canopy heights and displacement/roughness coefficients.
    class(scaler_t), intent(inout) :: upscaler !< Connected fine-to-coarse grid mapping.
    real(dp), intent(out) :: aerodyn_resistance1(:, :) !< Coarse-cell aerodynamic resistance by temporal slice.
    integer(i4) :: tt, i
    real(dp) :: canopy_height, zm, zm_zero, zh_zero, displace
    real(dp), allocatable :: maxLAI(:), resistance0(:)

    if (size(LAI0, 1) /= size(LCover0)) call error_message('MPR: aerodynamic LAI and land-cover sizes differ.')
    if (size(aerodyn_resistance1, 2) /= size(LAI0, 2)) &
      call error_message('MPR: aerodynamic resistance output LAI dimension does not match the input.')
    if (any(LCover0 < 1_i4) .or. any(LCover0 > 3_i4)) &
      call error_message('MPR: aerodynamic resistance requires land-cover classes 1 through 3.')
    allocate(maxLAI(size(LCover0)), resistance0(size(LCover0)))
    maxLAI = maxval(LAI0, dim=2)
    do tt = 1_i4, size(LAI0, 2)
      !$omp parallel do default(shared) private(canopy_height,zm,zm_zero,zh_zero,displace) schedule(static)
      do i = 1_i4, size(LCover0)
        select case (LCover0(i))
        case (1_i4)
          canopy_height = param(1)
        case (2_i4)
          canopy_height = param(2)
        case (3_i4)
          canopy_height = param(3) * LAI0(i, tt) / maxLAI(i)
        end select
        ! FAO Irrigation and Drainage Paper 56, pp. 19 ff.: measurement height must exceed the canopy.
        zm = WindMeasHeight
        if (zm < canopy_height) zm = canopy_height + zm
        displace = param(4) * canopy_height
        zm_zero = param(5) * canopy_height
        zh_zero = param(6) * zm_zero
        resistance0(i) = log((zm - displace) / zm_zero) * log((zm - displace) / zh_zero) / (karman**2.0_dp)
      end do
      !$omp end parallel do
      call upscaler%execute(resistance0, aerodyn_resistance1(:, tt), upscaling_operator=up_a_mean)
    end do
  end subroutine aerodynamical_resistance
end module mo_mpr_aerodynamics
