!> \file mo_common_functions.f90
!> \brief \copybrief mo_common_functions
!> \details \copydetails mo_common_functions

!> \brief Provides small utility functions used by multiple parts of the code (mHM, mRM, MPR)
!> \details Provides the functions in_bound used to check global_parameter ranges
!> \changelog
!!  - Robert Schweppe Dec 2017
!!    - refactoring
!> \authors Robert Schweppe
!> \date Dec 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_functions
  use mo_kind, only : dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: in_bound ! check parameter bounds

contains

  !    NAME
  !        in_bound

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: params" parameter:
  !>       col_1=Lower bound,
  !>       col_2=Upper bound
  !>       col_3=initial

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  function in_bound(params)
    implicit none

    ! parameter:
    ! col_1=Lower bound,
    ! col_2=Upper bound
    ! col_3=initial
    real(dp), dimension(:, :), intent(in) :: params

    logical :: in_bound


    if (any(params(:, 3) .lt. params(:, 1)) .or. any(params(:, 3) .gt. params(:, 2))) then
      in_bound = .false.
    else
      in_bound = .true.
    end if

  end function in_bound
end module mo_common_functions