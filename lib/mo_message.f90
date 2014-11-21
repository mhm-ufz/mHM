!> \file mo_message.f90

!> \brief Write out concatenated strings

!> \details Write out several strings concatenated on standard out or a given unit, either advancing or not.

!> \author Matthias Cuntz
!> \date Jul 2011

MODULE mo_message

  ! This module supplies routines to write out text

  ! Written Jul 2011, Matthias Cuntz - Inspired from Echam5 mo_exception.f90

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

  ! Copyright 2011 Matthias Cuntz

  USE mo_constants, ONLY: nout

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text    ! dummy string to use in subroutines
  PUBLIC :: message         ! versatile routine to write out strings in file or on screen

  CHARACTER(len=1024) :: message_text = ''

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         message

  !     PURPOSE
  !>        \brief Write out several string concatenated either on screen or in a file.

  !     CALLING SEQUENCE
  !         call message(t01=t01, t02=t02, t03=t03, t04=t04, t05=t05, t06=t06, t07=t07, &
  !                      t08=t08, t09=t09, t10=t10, unit=unit, advance=advance)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optional :: t01"        1st string
  !>        \param[in] "character(len=*), optional :: t02"        2nd string
  !>        \param[in] "character(len=*), optional :: t03"        3rd string
  !>        \param[in] "character(len=*), optional :: t04"        4th string
  !>        \param[in] "character(len=*), optional :: t05"        5th string
  !>        \param[in] "character(len=*), optional :: t06"        6th string
  !>        \param[in] "character(len=*), optional :: t07"        7th string
  !>        \param[in] "character(len=*), optional :: t08"        8th string
  !>        \param[in] "character(len=*), optional :: t09"        9th string
  !>        \param[in] "character(len=*), optional :: t10"        10th string
  !>        \param[in] "integer         , optional :: unit"       Unit to write to (default: nout)
  !>        \param[in] "character(len=*), optional :: advance"    WRITE advance keyword (default: 'yes')\n
  !>                                       yes: newline will be written after message\n
  !>                                       no:  no newline at end of message

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         "Only" 10 input strings

  !     EXAMPLE
  !         call message('A=',advance='no')
  !         call message(num2str(a))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011

  SUBROUTINE message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, uni, advance)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t01
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t02
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t03
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t04
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t05
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t06
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t07
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t08
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t09
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t10
    INTEGER,          INTENT(IN), OPTIONAL :: uni
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: advance

    INTEGER              :: iout
    CHARACTER(len=32000) :: out
    CHARACTER(len=3)     :: iadv
#ifdef GFORTRAN
    CHARACTER(len=32000) :: nold
#endif

    if (present(uni)) then
       iout = uni
    else
       iout = nout
    end if
    if (present(advance)) then
       iadv = ''
       iadv(1:min(len(advance),3)) = advance(1:min(len(advance),3))
    else
       iadv = 'yes'
    end if

    out = ''
    ! start from back so that trim does not remove user desired blanks
#ifdef GFORTRAN
    ! GFORTRAN has problems with concatenation operator //
    ! It is also weird in write:
    !    write(out,'(A,A)') t10, trim(out)
    ! writes t10 twice into out.
    nold = out
    if (present(t10)) write(out,'(A,A)') t10, trim(nold)
    nold = out
    if (present(t09)) write(out,'(A,A)') t09, trim(nold)
    nold = out
    if (present(t08)) write(out,'(A,A)') t08, trim(nold)
    nold = out
    if (present(t07)) write(out,'(A,A)') t07, trim(nold)
    nold = out
    if (present(t06)) write(out,'(A,A)') t06, trim(nold)
    nold = out
    if (present(t05)) write(out,'(A,A)') t05, trim(nold)
    nold = out
    if (present(t04)) write(out,'(A,A)') t04, trim(nold)
    nold = out
    if (present(t03)) write(out,'(A,A)') t03, trim(nold)
    nold = out
    if (present(t02)) write(out,'(A,A)') t02, trim(nold)
    nold = out
    if (present(t01)) write(out,'(A,A)') t01, trim(nold)
    ! output at least one space otherwise some compilers get confused on Mac (empty assembler statement)
    if ((lle(trim(out),'') .and. lge(trim(out),''))) then
       nold = out
       write(out,'(A,A)') trim(nold), ' '
    endif
    write(iout,'(a)',advance=iadv) trim(out)
#else
    if (present(t10)) out = t10//trim(out)
    if (present(t09)) out = t09//trim(out)
    if (present(t08)) out = t08//trim(out)
    if (present(t07)) out = t07//trim(out)
    if (present(t06)) out = t06//trim(out)
    if (present(t05)) out = t05//trim(out)
    if (present(t04)) out = t04//trim(out)
    if (present(t03)) out = t03//trim(out)
    if (present(t02)) out = t02//trim(out)
    if (present(t01)) out = t01//trim(out)
    ! output at least one space otherwise some compilers get confused on Mac (empty assembler statement)
    if ((lle(trim(out),'') .and. lge(trim(out),''))) then
       write(iout,'(a)',advance=iadv) trim(out)//' '
    else
       write(iout,'(a)',advance=iadv) trim(out)
    endif
#endif

  END SUBROUTINE message

END MODULE mo_message
