!> \file    mo_check.f90
!> \brief   \copybrief mo_check
!> \details \copydetails mo_check

!> \brief   Input checking routines
!> \details This module provides sanity checks for the input data.
!> \authors Sebastian Mueller
!> \date    Nov 2020
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_check

  USE mo_kind, ONLY : i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: check_dir

CONTAINS

  !>       \brief   Check if a given directory exists.
  !>       \details Check if a given directory exists and write out a message about it.
  !!                Will also give potential information about prefixes given with the path
  !>       \authors Sebastian Mueller
  !>       \date    Nov 2020

  subroutine check_dir(path, text, raise, tab, text_length)

    use mo_constants, ONLY : nout, nerr
    use mo_message, only : message, show_msg, show_err
    use mo_os, only : path_split, path_isdir

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)           :: path        !< input path to check
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text        !< text to write out
    LOGICAL, INTENT(IN), OPTIONAL          :: raise       !< whether to throw an error if folder does not exist
    integer(i4), INTENT(in), OPTIONAL      :: tab         !< tab-depth
    integer(i4), INTENT(in), OPTIONAL      :: text_length !< maximal text length (for aligning)

    LOGICAL            :: raise_
    integer(i4)        :: tab_
    integer(i4)        :: text_length_, uni
    LOGICAL            :: is_dir, error, show
    CHARACTER(len=255) :: head, tail, info, prefix_info, ws, text_

    ! set standard values
    prefix_info = ""
    ws = " " ! this should hold 255 whitespaces
    text_ = "Directory:"
    raise_ = .false.
    tab_ = 0
    if (present(text)) text_ = text
    if (present(raise)) raise_ = raise
    if (present(tab)) tab_ = tab
    text_length_ = len_trim(text_)
    if (present(text_length)) text_length_ = text_length

    ! split path to retrieve potential prefix to output files
    call path_split(path, head, tail)
    ! check if base directory exists
    is_dir = path_isdir(head) ! allow file prefix as path tail

    if ( is_dir ) then
      info = trim(head) // " (found)"
    else
      info = trim(head) // " (not found)"
    end if
    if ( len_trim(tail) > 0 ) prefix_info = "added file prefix: " // trim(tail)

    error = .not. is_dir .and. raise_
    show = show_msg
    uni = nout
    if ( error ) then
      show = show_err
      uni = nerr
    end if

    call message( &
      ws(1:tab_), &
      trim(text_), &
      ws(1:max(0, text_length_-len_trim(text_))), &
      trim(info), &
      " ", &
      trim(prefix_info), &
      show=show, &
      uni=uni &
    )
    ! throw error if wanted
    if ( error ) stop 1

  end subroutine check_dir

END MODULE mo_check
