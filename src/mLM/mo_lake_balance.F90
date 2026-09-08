!> \dir mLM
!> \copydoc f_mLM

!> \defgroup   f_mLM mLM - Fortran modules
!> \brief      Modules for lake-process topology, forcing, and balance calculations.
!> \details    Reusable lake logic shared by the mLM and meteorology exchange containers.

!> \file    mo_lake_balance.f90
!> \copydoc mo_lake_balance

!> \brief   Stateless lake water-balance helpers.
!> \details Shared lake-process calculations independent of exchange lifecycle.
!> \version 0.1
!> \changelog
!! - Pallav Shrestha (2018-2023): lake coupling concepts.
!! - Sebastian Mueller (2026): minimal lake water-balance extraction.
!> \authors Pallav Shrestha, Sebastian Mueller
!> \date    September 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
module mo_lake_balance
  use mo_kind, only: dp
  implicit none
  private
  public :: lake_balance_outflow
contains
  !> \brief Apply one model-step atmospheric balance to routed lake inflow.
  pure function lake_balance_outflow(inflow, area, pre, pet, step_seconds) result(outflow)
    real(dp), intent(in) :: inflow(:), area(:), pre(:), pet(:)
    real(dp), intent(in) :: step_seconds
    real(dp) :: outflow(size(inflow))
    outflow = max(0.0_dp, inflow + area * (pre - pet) * 1.0e-3_dp / step_seconds)
  end function lake_balance_outflow
end module mo_lake_balance
