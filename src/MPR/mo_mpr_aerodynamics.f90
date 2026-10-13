!> \file mo_mpr_aerodynamics.f90
!> \brief \copybrief mo_mpr_aerodynamics
!> \details \copydetails mo_mpr_aerodynamics

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
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_aerodynamics
  use mo_kind, only: i4, dp
  use mo_constants, only: nodata_dp
  implicit none
  private
  public :: aerodynamical_resistance
contains
  !> \brief Regionalize aerodynamic resistance from LAI and land cover.
  subroutine aerodynamical_resistance(LAI0, LCover0, param, mask0, Id0, n_subcells1, upper_bound1, lower_bound1, &
       left_bound1, right_bound1, aerodyn_resistance1)

    use mo_constants, only : eps_dp
    use mo_mpr_constants, only : WindMeasHeight, karman
    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    real(dp), dimension(:, :), intent(in) :: LAI0 !< LAI at level-0

    integer(i4), dimension(:), intent(in) :: LCover0 !< land cover field

    real(dp), dimension(6), intent(in) :: param !< input parameter

    logical, dimension(:, :), intent(in) :: mask0 !< mask at level 0

    integer(i4), dimension(:), intent(in) :: Id0 !< Cell ids of hi res field

    integer(i4), dimension(:), intent(in) :: n_subcells1 !< number of l0 cells within a l1 cell

    integer(i4), dimension(:), intent(in) :: upper_bound1 !< upper row of a l1 cell in l0 grid

    integer(i4), dimension(:), intent(in) :: lower_bound1 !< lower row of a l1 cell in l0 grid

    integer(i4), dimension(:), intent(in) :: left_bound1 !< left col of a l1 cell in l0 grid

    integer(i4), dimension(:), intent(in) :: right_bound1 !< right col of a l1 cell in l0 grid

    real(dp), dimension(:, :), intent(out) :: aerodyn_resistance1 !< aerodynmaical resistance

    integer(i4) :: tt

    real(dp), dimension(:), allocatable :: maxLAI

    real(dp), dimension(:), allocatable :: zm

    real(dp), dimension(:), allocatable :: canopy_height0

    real(dp), dimension(:), allocatable :: zm_zero, zh_zero, displace

    ! dim 1 = number of cells on level 0,
    ! dim2=month of year
    real(dp), dimension(:, :), allocatable :: aerodyn_resistance0


    ! initialize some things
    allocate(zm                  (size(LCover0, dim = 1))) ; zm = nodata_dp
    allocate(zm_zero             (size(LCover0, dim = 1))) ; zm_zero = nodata_dp
    allocate(zh_zero             (size(LCover0, dim = 1))) ; zh_zero = nodata_dp
    allocate(displace            (size(LCover0, dim = 1))) ; displace = nodata_dp
    allocate(canopy_height0      (size(LCover0, dim = 1))) ; canopy_height0 = nodata_dp
    allocate(aerodyn_resistance0 (size(LCover0, dim = 1), size(LAI0, 2))) ; aerodyn_resistance0 = nodata_dp
    allocate(maxLAI              (size(LCover0, dim = 1))) ; maxLAI = nodata_dp

    ! regionalization of canopy height
    ! substitute with canopy height
    canopy_height0 = merge(param(1), canopy_height0, LCover0 == 1)  ! forest
    canopy_height0 = merge(param(2), canopy_height0, LCover0 == 2)  ! impervious

    ! old implementation used values from LUT statically for all cells (Jan-Dec, 12 values):
    ! 7 Intensive-orchards: 2.0, 2.0,  2.0,  2.0,  3.0, 3.5,  4.0,  4.0,  4.0,  2.5,  2.0,  2.0
    maxLAI = MAXVAL(LAI0, dim=2)

    do tt = 1, size(LAI0, 2)

       ! pervious canopy height is scaled with LAI
       canopy_height0 = merge((param(3) * LAI0(:, tt) / maxLAI), canopy_height0, LCover0 == 3)  ! pervious

       ! estimation of the aerodynamic resistance on the lower level
       ! see FAO Irrigation and Drainage Paper No. 56 (p. 19 ff) for more information
       zm = WindMeasHeight
       ! correction: if wind measurement height is below canopy height loagarithm becomes negative
       zm = merge(canopy_height0 + zm, zm, ((abs(zm - nodata_dp) .GT. eps_dp) .AND. (zm .LT. canopy_height0)))

       ! zh       = zm
       displace = param(4) * canopy_height0
       zm_zero = param(5) * canopy_height0
       zh_zero = param(6) * zm_zero
       !
       ! calculate aerodynamic resistance (changes monthly)
       aerodyn_resistance0(:, tt) = log((zm - displace) / zm_zero) * log((zm - displace) / zh_zero) / (karman**2.0_dp)
       aerodyn_resistance1(:, tt) = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
            left_bound1, right_bound1, Id0, mask0, nodata_dp, aerodyn_resistance0(:, tt))

    end do

  end subroutine aerodynamical_resistance
end module mo_mpr_aerodynamics
