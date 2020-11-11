MODULE mo_os

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  USE mo_message, ONLY: message

  IMPLICIT NONE

  PUBLIC :: path_exists
  PUBLIC :: path_isfile
  PUBLIC :: path_isdir
  PUBLIC :: path_splitext
  PUBLIC :: path_split

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         path_exists
  !
  !     PURPOSE
  !>        \brief Existence of a path
  !
  !>        \details Checks whether a given path exists.
  !
  !     INTENT(IN)
  !>        \paran[in] "character(len=*) :: path"        To be checked path
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: quiet_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  no message\n
  !>                                                            .false.: write out message if path is not found
  !>       \param[in] "logical, optional :: throwError_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  quit program and throw error message if path is not found\n
  !>                                                            .false.: continue program and no message
  !>       \param[in] "logical, optional :: dirOrFile_"   Be verbose or not (default: potential error message with imprecise information)\n
  !>                                                            .true.:  if the given path does not exist, the potential error message relate to a directory\n
  !>                                                            .false.: if the given path does not exist, the potential error message relate to a file
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>       \param[out] "logical, optional :: result_"   Be verbose or not. Contains the result of checking the path.
  !
  !     RETURN
  !        None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         call path_exists(path, quiet)
  !         -> see also example in pf_test directory
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Nicola Doering
  !>        \date Aug 2020

  SUBROUTINE path_exists(path, quiet_, throwError_, dirOrFile_, result_)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: path
    LOGICAL, INTENT(IN), OPTIONAL :: quiet_
    LOGICAL, INTENT(IN), OPTIONAL :: throwError_
    LOGICAL, INTENT(IN), OPTIONAL :: dirOrFile_
    LOGICAL, INTENT(OUT), OPTIONAL :: result_

    LOGICAL :: quiet = .false.
    LOGICAL :: throwError = .false.
    LOGICAL :: exists
    CHARACTER(LEN=40) :: messagetext

    inquire (file=path, exist=exists)

#ifdef INTEL
    if (.not. exists) inquire (directory=path, exist=exists)
#endif

    if (present(quiet_)) quiet = quiet_
    if (present(throwError_)) throwError = throwError_
    if (.not. present(dirOrFile_)) then
      messagetext = "Cannot find the file or directory: "
    else if (dirOrFile_) then
      messagetext = "Cannot find the directory: "
    else
      messagetext = "Cannot find the file: "
    endif

    if (.not. exists) then
      if (.not. quiet .or. throwError) call message(trim(messagetext), path)
      if (throwError) stop 1
    endif

    if (present(result_)) result_ = exists

  END SUBROUTINE path_exists

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         path_isfile
  !
  !     PURPOSE
  !>        \brief Whether the path describes a file.
  !
  !>        \details Checks whether a given path exists and describes a file.
  !
  !     INTENT(IN)
  !>        \paran[in] "character(len=*) :: path"        To be checked path
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: quiet_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  no message\n
  !>                                                            .false.: write out message if path is not found or does not desrcibes a file
  !>       \param[in] "logical, optional :: throwError_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  quit program and throw error message if path is not found or does not desrcibes a file\n
  !>                                                            .false.: continue program and no message
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "logical, optional :: result_"   Be verbose or not. Contains the result of checking the path.
  !
  !     RETURN
  !               None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         call path_isfile(path, quiet, throwError, result)
  !         -> see also example in pf_test directory
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Nicola Doering
  !>        \date Aug 2020

  SUBROUTINE path_isfile(path, quiet_, throwError_, result_)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: path
    LOGICAL, INTENT(IN), OPTIONAL :: quiet_
    LOGICAL, INTENT(IN), OPTIONAL :: throwError_
    LOGICAL, INTENT(OUT), OPTIONAL :: result_

    LOGICAL :: quiet
    LOGICAL :: throwError
    LOGICAL :: exists
    LOGICAL :: isfile

    quiet = .false.
    throwError = .false.
    isfile = .true.

    if (present(quiet_)) quiet = quiet_
    if (present(throwError_)) throwError = throwError_

    call path_exists(path, quiet, throwError, .false., exists)
    if (exists) then
      !checking whether the path is ending with '/' or '/.' which would indicates a directory
      if (path(len(path):len(path)) .eq. '/' .or. path(len(path) - 1:len(path)) .eq. '/.') then
        isfile = .false.
      else
        !checking whether would still exist if '/.' is added to the end, in this case it is a directory
#ifdef INTEL
        inquire (directory=path//'/.', exist=exists)
#else
        inquire (file=path//'/.', exist=exists)
#endif
        if (exists) then
          isfile = .false.
          if (.not. quiet .or. throwError) call message('The following path describes a directory and not a file: ', path)
          if (throwError) stop 1
        endif
      endif
    else
      isfile = .false.
    endif

    if (present(result_)) result_ = isfile

  END SUBROUTINE path_isfile

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         path_isdir
  !
  !     PURPOSE
  !>        \brief Whether the path describes a directory.
  !
  !>        \details Checks whether a given path exists and describes a directory.
  !
  !     INTENT(IN)
  !>        \paran[in] "character(len=*) :: path"                To be checked path
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: quiet_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  no message\n
  !>                                                            .false.: write out message if path is not found or does not desrcibes a directory
  !>       \param[in] "logical, optional :: throwError_"   Be verbose or not (default: .false.)\n
  !>                                                            .true.:  quit program and throw error message if path is not found or does not desrcibes a directory\n
  !>                                                            .false.: continue program and no message
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>       \param[out] "logical, optional :: result_"   Be verbose or not. Contains the result of checking the path.
  !
  !     RETURN
  !                                        None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         call path_isdir(path, quiet, throwError, result)
  !         -> see also example in pf_test directory
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Nicola Doering
  !>        \date Aug 2020

  SUBROUTINE path_isdir(path, quiet_, throwError_, result_)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)        :: path
    LOGICAL, INTENT(IN), OPTIONAL :: quiet_
    LOGICAL, INTENT(IN), OPTIONAL :: throwError_
    LOGICAL, INTENT(OUT), OPTIONAL :: result_

    LOGICAL :: quiet
    LOGICAL :: throwError
    LOGICAL :: isdir
    LOGICAL :: exists

    quiet = .false.
    throwError = .false.
    isdir = .true.

    if (present(quiet_)) quiet = quiet_
    if (present(throwError_)) throwError = throwError_

    call path_exists(path, quiet, throwError, .true., exists)
    if (exists) then
      call path_isfile(path, .true., .false., exists)
      if (exists) then
        isdir = .false.
        if (.not. quiet .or. throwError) call message('The following path describes a file and not a directory: ', path)
        if (throwError) stop 1
      endif
    else
      isdir = .false.
    endif

    if (present(result_)) result_ = isdir

  END SUBROUTINE path_isdir

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         path_splitext
  !
  !     PURPOSE
  !>        \brief Splitting the path into root and ext
  !
  !>        \details Splitting the path name into a pair root and ext.
  !>                 Here, ext stands for extension and has the extension portion
  !>                 of the specified path while root is everything except ext part.
  !>                 If the path describes a directory or there is no extension
  !>                 portion ext is returned empty.
  !
  !     INTENT(IN)
  !>        \paran[in] "character(len=*) :: path"       To be checked path
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         \paran[out] "character(len=*) :: root"      Path without the extension
  !         \paran[out] "character(len=*) :: ext"       Extension of path
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !               None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         call path_splitext(path, root, ext)
  !         -> see also example in pf_test directory
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Nicola Doering
  !>        \date Aug 2020

  SUBROUTINE path_splitext(path, root, ext)

    CHARACTER(LEN=*), INTENT(IN)  :: path
    CHARACTER(LEN=*), INTENT(OUT) :: root
    CHARACTER(LEN=*), INTENT(OUT) :: ext

    INTEGER   :: i
    CHARACTER :: c
    LOGICAL :: isdir

    i = len(path) - 1
    c = path(len(path):len(path))

    !Checking, whether the path describes a directory so it cannot ends with an extension.
    call path_isdir(path, .true., .false., isdir)
    if (isdir) then
      i = len(path)
    else
      !running through the path, beginning at the end until a point is found that probably indicates
      !the seperation of a file name and its extension or a '/' occurs what means that the rest of the
      !path is consisting of directories
      do while (.not. (c .eq. '.' .or. c .eq. '/' .or. i .eq. 0))
        c = path(i:i)
        i = i - 1
      end do
      !checking whether the last symbol of the path is a point or the while-loop run through the whole path
      !without finding a point or ended at a '/'. In any case it is not possible to seperate an extension.
      if (i .eq. len(path) - 1 .or. i .eq. 0 .or. c .eq. '/') then
        i = len(path)
      endif
    endif

    root = path(1:i)
    ext = path(i + 1:len(path))
    return

  END SUBROUTINE path_splitext

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         path_split
  !
  !     PURPOSE
  !>        \brief Splitting the path into head and tail
  !
  !>        \details Splitting the path name into a pair head and tail.
  !>                 Here, tail is the last path name component and head is
  !>                 everything leading up to that.
  !>                 If the path ends with an '/' tail is returned empty and
  !>                 if there is no '/' in path head is returned empty.
  !
  !     INTENT(IN)
  !>        \paran[in] "character(len=*) :: path"       To be checked path
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         \paran[out] "character(len=*) :: head"      Everything leading up to the last path component
  !         \paran[out] "character(len=*) :: tail"      Last path name component
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !               None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         call path_split(path, head, tail)
  !         -> see also example in pf_test directory
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Nicola Doering
  !>        \date Aug 2020

  SUBROUTINE path_split(path, head, tail)

    CHARACTER(LEN=*), INTENT(IN)  :: path
    CHARACTER(LEN=*), INTENT(OUT) :: head
    CHARACTER(LEN=*), INTENT(OUT) :: tail

    INTEGER   :: i
    CHARACTER :: c
    LOGICAL :: isdir

    i = len(path) - 1
    c = path(len(path):len(path))

    !running through the path, beginning at the end until a point is found that probably indicates
    !the seperation of a file name and its extension or a '/' occurs what means that the rest of the
    !path is consisting of directories
    do while (.not. (c .eq. '/' .or. i .eq. 0))
      c = path(i:i)
      i = i - 1
    end do
    !checking whether the while-loop run through the whole path without finding a '/'
    if (i .eq. 0) then
      head = ''
      tail = path
    else
      head = path(1:i + 1)
      tail = path(i + 2:len(path))
    endif

    return

  END SUBROUTINE path_split

  ! ------------------------------------------------------------------

END MODULE mo_os
