!> \file mo_meteo_constants.f90
!> \copydoc mo_meteo_constants

!> \brief Provides meteorological constants
!> \details Provides coefficients for PET, extraterrestrial radiation, and vapour pressure.
! Modifications:
! Robert Schweppe Jun 2018 - refactoring and reformatting
!> \authors Matthias Cuntz
!> \date Nov 2011
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_meteo
module mo_meteo_constants
  use mo_kind, only: dp
  implicit none
  private

  !> Hargreaves-Samani ref. ET formula [deg C]
  real(dp), public, parameter :: HarSamConst = 17.800_dp

  ! Duffie formula for computing extraterrestrial radiation
  real(dp), public, parameter :: DuffieDr = 0.0330_dp
  real(dp), public, parameter :: DuffieDelta1 = 0.4090_dp
  real(dp), public, parameter :: DuffieDelta2 = 1.3900_dp

  !> Tetens's formula to calculate saturated vapour pressure
  real(dp), public, parameter :: tetens_c1 = 0.6108_dp
  real(dp), public, parameter :: tetens_c2 = 17.270_dp
  real(dp), public, parameter :: tetens_c3 = 237.30_dp
  !> calculation of the slope of the saturation vapour pressure curve following Tetens
  real(dp), public, parameter :: satpressureslope1 = 4098.0_dp

end module mo_meteo_constants
