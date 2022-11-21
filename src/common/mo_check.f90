!>       \file    mo_check.f90
!>       \copydoc mo_check

!>       \brief   Input checking routines
!>       \details This module provides sanity checks for the input data.
!>       \authors Sebastian Mueller
!>       \date    Nov 2020
MODULE mo_check

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: check_dir

CONTAINS

  !>       \brief   Check if a given directory exists.
  !>       \details Check if a given directory exists and write out a message about it.
  !!                Will also give potential information about prefixes given with the path
  !>       \authors Sebastian Mueller
  !>       \date    Nov 2020

  subroutine check_dir(path, text_, throwError_, tab_, text_length_)

    use mo_message, only : message
    use mo_os, only : path_split, path_isdir

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)           :: path         !< input path to check
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text_        !< text to write out
    LOGICAL, INTENT(IN), OPTIONAL          :: throwError_  !< wheather to throw an error if folder not existing
    integer(i4), INTENT(in), OPTIONAL      :: tab_         !< tab-depth
    integer(i4), INTENT(in), OPTIONAL      :: text_length_ !< maximal text length (for aligning)

    LOGICAL            :: throwError = .false.
    integer(i4)        :: tab = 0
    integer(i4)        :: text_length
    LOGICAL            :: is_dir
    CHARACTER(len=255) :: head, tail, info, prefix_info, ws, text

    ! set standard values
    prefix_info = ""
    ws = " " ! this should hold 255 whitespaces
    text = "Directory:"
    if (present(text_)) text = text_
    if (present(throwError_)) throwError = throwError_
    if (present(tab_)) tab = tab_
    text_length = len_trim(text)
    if (present(text_length_)) text_length = text_length_

    ! split path to retrieve potential prefix to output files
    call path_split(path, head, tail)
    ! check if base directory exists
    call path_isdir(head, quiet_=.true., result_=is_dir) ! allow file prefix as path tail

    if ( is_dir ) then
      info = trim(head) // " (found)"
    else
      info = trim(head) // " (not found)"
    end if
    if ( len_trim(tail) > 0 ) prefix_info = "added file prefix: " // trim(tail)

    call message( &
      ws(1:tab), &
      trim(text), &
      ws(1:max(0, text_length-len_trim(text))), &
      trim(info), &
      " ", &
      trim(prefix_info) &
    )
    ! throw error if wanted
    if (.not. is_dir .and. throwError) stop 1

  end subroutine check_dir

END MODULE mo_check
