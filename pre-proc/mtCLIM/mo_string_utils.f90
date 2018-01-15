!> \file mo_string_utils.f90

!> \brief String utilities

!> \details This module provides string conversion and checking utilities.

!> \authors Matthias Cuntz, Matthias Zink, Giovanni Dalmasso, David Schaefer
!> \date Dec 2011

MODULE mo_string_utils

  ! This module holds string conversion utilities

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

  ! Copyright 2011-2015 Matthias Cuntz

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: compress      ! Conversion   : 'A b C x Y z' -> 'AbCxYz'
#ifndef ABSOFT
  PUBLIC :: divide_string ! split string in substring with the help of delimiter
#endif
  PUBLIC :: equalStrings  ! compares two strings
  PUBLIC :: nonull        ! Check if string is still NULL
  PUBLIC :: num2str       ! Convert a number to a string
  PUBLIC :: separator     ! Format string: '-----...-----'
  PUBLIC :: splitString   ! splits string at given delimiter
  PUBLIC :: startsWith    ! checks if string starts with a certain string
  PUBLIC :: str2num       ! Converts string into an array of its numerical representation
  PUBLIC :: tolower       ! Conversion   : 'ABCXYZ' -> 'abcxyz'
  PUBLIC :: toupper       ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  
  ! public :: numarray2str

  ! ------------------------------------------------------------------

  !     NAME
  !         num2str

  !     PURPOSE
  !>        \brief Convert to string.

  !>        \details Convert a number or logical to a string with an optional format.

  !     CALLING SEQUENCE
  !         str = num2str(num,form=form)

  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8)/real(sp/dp)/logical :: num"    Number or logical

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optinal :: form"   Format string\n
  !>                                                         Defaults are:\n
  !>                                                         i4    - '(I10)'\n
  !>                                                         i8    - '(I20)'\n
  !>                                                         sp/dp - '(G32.5)'\n
  !>                                                         log   - '(L10)'

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return character(len=X) :: str &mdash; String of formatted input number or logical\n
  !>                                                Ouput length X is:\n
  !>                                                i4    - 10\n
  !>                                                i8    - 20\n
  !>                                                sp/dp - 32\n
  !>                                                log   - 10

  !     RESTRICTIONS
  !         Uses WRITE to write into string. Recursive write is not permitted before Fortran 2003
  !         so that one cannot use
  !             write(*,*) 'A='//num2str(a)
  !         Use 'call message' from mo_messages.f90
  !             use mo_messages, only message
  !             call message('A=', trim(num2str(a)))
  !         or write into another string first:
  !             str = 'A='//num2str(a)
  !             write(*,*) trim(str)

  !     EXAMPLE
  !         str = num2str(3.1415217_i4,'(F3.1)')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011
  INTERFACE num2str
     MODULE PROCEDURE i42str, i82str, sp2str, dp2str, log2str
  END INTERFACE num2str


  ! ------------------------------------------------------------------

  !     NAME
  !         numarray2str

  !     PURPOSE
  !>        \brief Convert to string.

  !>        \details Convert a array of numbers or logicals to a string.

  !     CALLING SEQUENCE
  !         str = numarray2str(num)

  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8)/real(sp/dp)/logical :: num(:)"    Array of numbers or logicals

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None
  
  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return character(len=X) :: str &mdash; String of formatted input number or logical\n
  !

  !     RESTRICTIONS
  !        Unknown

  !     EXAMPLE
  !        None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011

  INTERFACE numarray2str
     MODULE PROCEDURE i4array2str
  END INTERFACE numarray2str

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

  CHARACTER(len=*), PARAMETER :: separator = repeat('-',70)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         compress

  !     PURPOSE
  !         \brief Remove white spaces

  !         \details Return a copy of an input string with all whitespace (spaces and tabs) removed

  !     CALLING SEQUENCE
  !         noSpaces = compress(whiteSpaces)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: whiteSpaces"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         \param[out] "integer(i4) :: n"  Integer

  !     RETURN
  !         \return character(len = len(whiteSpaces)) :: compress;  String where all all whitespace (spaces and tabs) are removed

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'Hallo'
  !         whiteSpaces = compress('H a l l o')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author Giovanni Dalmasso - modified from Paul van Delst, CIMSS/SSEC 18-Oct-1999
  !         \date Jan 2013

  function compress( whiteSpaces, n )

        use mo_kind,    only : i4

        implicit none

        character(len=*),               intent(in)  :: whiteSpaces
        integer(i4),        optional,   intent(out) :: n

        character(len(whiteSpaces))                 ::  compress

        ! Local parameters
        integer(i4),    parameter                   :: iachar_space = 32_i4
        integer(i4),    parameter                   :: iachar_tab   = 9_i4

        ! Local variables
        integer(i4)                                 :: i, j
        integer(i4)                                 :: iachar_character

        ! Setup

        ! Initialise compress
        compress = ' '
        ! Initialise counter
        j = 0_i4

        ! Loop over string
        do i = 1, len(whiteSpaces)
            ! Convert the current character to its position
            iachar_character = iachar(whiteSpaces(i:i))

            ! If the character is NOT a space ' ' or a tab '->|' copy it to the output string.
            if ( iachar_character .ne. iachar_space .and. iachar_character .ne. iachar_tab )    then
                j = j + 1
                compress(j:j) = whiteSpaces(i:i)
            end if
        end do

        ! Save the non-whitespace count
        if ( present(n) ) n = j

    end function compress

#ifndef ABSOFT
  ! ------------------------------------------------------------------

  !     NAME
  !         divide_string

  !     PURPOSE
  !>        \brief Divide string in substrings.

  !>        \details Divides a string in several substrings (array of strings) with the help of a user
  !>        specified delimiter.

  !     CALLING SEQUENCE
  !         call divide_string(string, delim, strArr(:))

  !     INTENT(IN)
  !>        \param[in] "CHARACTER(len=*), INTENT(IN) :: string"     - string to be divided
  !>        \param[in] "CHARACTER(len=*), INTENT(IN) :: delim"      - delimiter specifying places for division

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "CHARACTER(len=*), DIMENSION(:), ALLOCATABLE,  INTENT(OUT) :: strArr"
  !>                     Array of substrings, has to be allocateable and is handed to the routine unallocated

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         only character types allowed
  !         output array should be allocateable array, which is unallocated handed to the subroutine
  !             allocation is done in in devide_string

  !     EXAMPLE
  !        divide_string('I want to test this routine!', ' ', strArr(:))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Oct 2012

  SUBROUTINE divide_string(string, delim, strArr)

    IMPLICIT NONE

    CHARACTER(len=*)             , INTENT(IN)        :: string
    CHARACTER(len=*)             , INTENT(IN)        :: delim
    CHARACTER(len=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: strArr

    CHARACTER(256)                                   :: stringDummy   ! string in fisrt place but cutted in pieces
    CHARACTER(256), DIMENSION(:) , ALLOCATABLE       :: strDummyArr   ! Dummy arr until number of substrings is known
    INTEGER(i4)                                      :: pos           ! position of dilimiter
    INTEGER(i4)                                      :: nosubstr      ! number of substrings in string

    stringDummy = string

    allocate(strDummyArr(len_trim(stringDummy)))
    pos=999_i4
    nosubstr=0_i4
    ! search for substrings and theirs count
    do
       pos = index(trim(adjustl(stringDummy)), delim)
       ! exit if no more delimiter is find and save the last part of the string
       if (pos .EQ. 0_i4) then
          nosubstr = nosubstr + 1_i4
          StrDummyArr(nosubstr) = trim(stringDummy)
          exit
       end if

       nosubstr = nosubstr + 1_i4
       strDummyArr(nosubstr) = stringDummy(1:pos-1)
       stringDummy = stringDummy(pos+1:len_trim(stringDummy))
    end do
    ! hand over results to strArr
    if (nosubstr .EQ. 0_i4) then
       print*, '***WARNING: string does not contain delimiter. There are no substrings. (subroutine DIVIDE_STRING)'
       return
    else
       allocate(strArr(nosubstr))
       strArr = StrDummyArr(1:nosubstr)
    end if

    deallocate(strDummyArr)

  END SUBROUTINE divide_string
#endif

    ! ------------------------------------------------------------------

  !     NAME
  !         equalStrings

  !     PURPOSE
  !         \brief Checks if two string are equal

  !         \details Returns true if the given string arguments are equal

  !     CALLING SEQUENCE
  !         isequal = equalString(string1,string2)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: string1"    String
  !         \param[in] "character(len=*) :: string2"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return logical

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author David Schaefer
  !         \date Mar 2015

  function equalStrings(string1,string2)
    implicit none

    character(len=*), intent(in)     :: string1, string2
    integer(i4),      allocatable    :: array1(:), array2(:)
    integer(i4)                      :: i
    logical                          :: equalStrings

    array1 = str2num(trim(string1))
    array2 = str2num(trim(string2))
    equalStrings = .false.

    if (size(array1) .eq. size(array2)) then
       equalStrings = .true.
       do i=1, size(array1)
          if (array1(i) .ne. array2(i)) then
             equalStrings = .false.
             exit 
          end if
       end do
    end if

  end function equalStrings

  ! ------------------------------------------------------------------

  !     NAME
  !         nonull

  !     PURPOSE
  !>        \brief Checks if string was already used

  !>        \details Checks if string was already used, i.e. does not contain NULL character anymore.

  !     CALLING SEQUENCE
  !         used = nonull(str)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: str"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return logical :: used &mdash;  .true.: string was already set; .false.: string still in initialised state

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         if (nonull(str)) write(*,*) trim(str)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Jan 2012

  FUNCTION nonull(str)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in) :: str
    LOGICAL                      :: nonull

    if (scan(str, achar(0)) == 0) then
       nonull = .true.
    else
       nonull = .false.
    endif

  END FUNCTION nonull

  ! ------------------------------------------------------------------

  !     NAME
  !         splitString

  !     PURPOSE
  !         \brief split string at delimiter

  !         \details Split string at delimiter an return an array of strings

  !     CALLING SEQUENCE
  !         string_parts = splitString(string,delim)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: string"    String
  !         \param[in] "character(len=*) :: delim"     String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return character(len=245) :: out(:)

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author David Schaefer
  !         \date Mar 2015

  function splitString(string,delim) result(out)
    
    use mo_append, only : append    
    implicit none
    
    character(len=*),   intent(in)        :: string
    character(len=*),   intent(in)        :: delim
    character(len=256), allocatable       :: out(:)
    integer(i4),        allocatable       :: string_array(:), delim_array(:)
    integer(i4)                           :: i, start
    !
    if (allocated(out)) deallocate(out)
    string_array = str2num(string//delim)
    delim_array = str2num(delim)
    start = 1

    do i=1, size(string_array) - size(delim_array) + 1
       if (all(string_array(i:i+size(delim_array)-1) .eq. delim_array)) then
          call append(out, numarray2str(string_array(start:i-1)))
          start = i + size(delim_array)
       end if
    end do
    !
  end function splitString

  ! ------------------------------------------------------------------

  !     NAME
  !         startsWith

  !     PURPOSE
  !         \brief Checks if string starts with character(s)

  !         \details Returns true if string starts with characters, flase otherwise

  !     CALLING SEQUENCE
  !         starts_with = startsWith(string,start)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: string"    String
  !         \param[in] "character(len=*) :: start"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return logical

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author David Schaefer
  !         \date Mar 2015

  function startsWith(string, start)
    
    implicit none

    character(len=*), intent(in)     :: string, start
    integer(i4), allocatable         :: string_array(:), start_array(:)
    logical                          :: startsWith

    string_array = str2num(string)
    start_array = str2num(start)

    startsWith = .false.
    if (all(string_array(1:1+size(start_array)-1) .eq. start_array)) then 
       startsWith = .true.
    end if
  end function startsWith

  ! ------------------------------------------------------------------

  !     NAME
  !         tolower

  !     PURPOSE
  !>        \brief Convert to lower case

  !>        \details Convert all upper case letters in string to lower case letters.

  !     CALLING SEQUENCE
  !         low = tolower(upper)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: upper"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return character(len=len_trim(upper)) :: low  &mdash;  String where all uppercase in input is converted to lowercase

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'hallo'
  !         low = tolower('Hallo')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011

  FUNCTION tolower(upper)

    IMPLICIT NONE

    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
       IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
            ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
          tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
       ELSE
          tolower(i:i) = upper(i:i)
       END IF
    END DO

  END FUNCTION tolower

  ! ------------------------------------------------------------------

  !     NAME
  !         toupper

  !     PURPOSE
  !         \brief Convert to upper case

  !         \details Convert all lower case letters in string to upper case letters.

  !     CALLING SEQUENCE
  !         up = toupper(lower)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: lower"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return character(len=len_trim(lower)) :: up  &mdash;  String where all lowercase in input is converted to uppercase

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'HALLO'
  !         up = toupper('Hallo')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !         \date Dec 2011

  pure function toupper (lower)

    implicit none

    character(len=*)              ,intent(in) :: lower
    character(len=len_trim(lower))            :: toupper

    integer            :: i
    integer, parameter :: idel = ichar('A')-ichar('a')

    do i=1, len_trim(lower)
       if (ichar(lower(i:i)) >= ichar('a') .and. &
            ichar(lower(i:i)) <= ichar('z')) then
          toupper(i:i) = char( ichar(lower(i:i)) + idel )
       else
          toupper(i:i) = lower(i:i)
       end if
    end do

  end function toupper


  ! -----------------------------------------------------------
  ! PRIVATE ROUTINES
  ! (no "template" documentation required)
  ! -----------------------------------------------------------

  PURE FUNCTION i42str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i4),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: i42str

    if (present(form)) then
       write(i42str,form) nn
    else
       write(i42str,'(I10)') nn
    end if
    !i42str = adjustl(i42str)

  END FUNCTION i42str


  PURE FUNCTION i82str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i8),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=20) :: i82str

    if (present(form)) then
       write(i82str,form) nn
    else
       write(i82str,'(I20)') nn
    end if
    !i82str = adjustl(i82str)

  END FUNCTION i82str


  PURE FUNCTION sp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(sp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: sp2str

    if (present(form)) then
       write(sp2str,form) rr
    else
       write(sp2str,'(G32.5)') rr
    end if
    !sp2str = adjustl(sp2str)

  END FUNCTION sp2str


  PURE FUNCTION dp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(dp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: dp2str

    if (present(form)) then
       write(dp2str,form) rr
    else
       write(dp2str,'(G32.5)') rr
    end if
    !dp2str = adjustl(dp2str)

  END FUNCTION dp2str


  PURE FUNCTION log2str(ll,form)
    ! returns logical ll as a string (often needed in printing messages)
    IMPLICIT NONE
    LOGICAL,          INTENT(in)           :: ll
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: log2str

    if (present(form)) then
       write(log2str,form) ll
    else
       write(log2str,'(L10)') ll
    end if
    !log2str = adjustl(log2str)

  END FUNCTION log2str

  function i4array2str(arr) result(out)

    integer(i4), intent(in)     :: arr(:)
    integer(i4)                 :: ii
    character(len=size(arr))    :: out

    out = " "
    do ii=1,size(arr)
       out(ii:ii) = char(arr(ii))
    end do

  end function i4array2str

    ! ------------------------------------------------------------------

  !     NAME
  !         str2num

  !     PURPOSE
  !         \brief Converts string into an array of its numerical representation

  !         \details Converts string into an integer array of the numerical values of the letters

  !     CALLING SEQUENCE
  !         str2num = startsWith(string)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: string"    String

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return integer  :: out(:)

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author David Schaefer
  !         \date Mar 2015

  function str2num(string) result(out)
  
    implicit none

    character(len=*), intent(in)       :: string
    integer(i4), allocatable           :: out(:)
    integer(i4)                        :: i  

    if (allocated(out)) deallocate(out)
    allocate(out(len(string)))

    do i=1,len(string)
       out(i) = ichar(string(i:i))
    end do

  end function str2num


END MODULE mo_string_utils
