!> \file mo_cli_support.f90
!> \copydoc mo_cli_support

!> \brief Verbosity support for retained scientific message routines.
!> \details Extracted from mo_mhm_cli for the standalone v6 driver. Sets the
!! message, warning, and error visibility used by scientific message routines.
!> \authors Sebastian Mueller
!> \date Oct 2021
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_cli_support
  use mo_kind, only: i4
  implicit none
  private
  public :: set_verbosity_level
contains
  !> \brief Set the message verbosity for the application.
  subroutine set_verbosity_level(level)
    use mo_message, only: SHOW_MSG, SHOW_ERR, SHOW_WARN
    implicit none
    integer(i4), intent(in), optional :: level !< verbosity level (0: no output, 1: errors, 2: warnings, 3: all (default))
    integer(i4) :: level_
    SHOW_MSG = .false.
    SHOW_WARN = .false.
    SHOW_ERR = .false.
    level_ = 3_i4
    if ( present(level) ) level_ = level
    if ( level_ > 0_i4 ) SHOW_ERR = .true.
    if ( level_ > 1_i4 ) SHOW_WARN = .true.
    if ( level_ > 2_i4 ) SHOW_MSG = .true.
  end subroutine set_verbosity_level
end module mo_cli_support
