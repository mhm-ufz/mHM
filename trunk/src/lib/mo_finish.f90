!> \file mo_finish.f90

!> \brief Finish a program gracefully

!> \details This module supplies a routine that writes out final comments and then stops.

!> \authors Original of Echam5, (C) MPI-MET, Hamburg, Germany.\n
!> Modified Matthias Cuntz
!> \date Jan 2011

MODULE mo_finish

  ! This module supplies a routine to finish a program gracefully.

  ! Written Jan 2011, Matthias Cuntz - Adapted from Echam5, (C) MPI-MET, Hamburg, Germany

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

  USE mo_string_utils, ONLY: separator

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finish     ! Write out error message and stop program

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         finish

  !     PURPOSE
  !>        \brief Finish a program gracefully

  !>        \details Stop a program but writing out a message first that is separated
  !>        from earlier output by -------------- (i.e. the separator of mo_string_utils)

  !     CALLING SEQUENCE
  !         call finish(name, text=text, unit=unit)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: name"         First string separated from otional second by :

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optional :: text"         Second string separated by :
  !>        \param[in] "integer, optional          :: unit"         File unit for write (default: *)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call finish('main','End gracefully')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011

  SUBROUTINE finish(name, text, unit)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN)           :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text
    INTEGER,          INTENT(IN), OPTIONAL :: unit

    IF (PRESENT(unit)) THEN
       WRITE (unit,'(a)') separator

       IF (PRESENT(text)) THEN
          WRITE (unit,'(a,a,a)') name, ': ', text
       ELSE
          WRITE (unit,'(a)') name
       END IF

       WRITE (unit,'(a)') separator
    ELSE
       WRITE (*,'(a)') separator

       IF (PRESENT(text)) THEN
          WRITE (*,'(a,a,a)') name, ': ', text
       ELSE
          WRITE (*,'(a)') name
       END IF

       WRITE (*,'(a)') separator
    ENDIF

    STOP

  END SUBROUTINE finish

END MODULE mo_finish
