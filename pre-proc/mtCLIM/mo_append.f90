!> \file mo_append.f90

!> \brief Append values on existing arrays.

!> \details Provides routines to append (rows) and paste (columns) scalars, vectors, 
!>          and matrixes onto existing arrays.

!> \author Juliane Mai
!> \date Aug 2012 

MODULE mo_append

  ! This module is appending and pasting scalars, vectors, and matrixes into one.
  ! and is part of the UFZ CHS Fortran library.


  ! Written  Juliane Mai,    Aug 2012
  ! Modified Juliane Mai,    Aug 2012 - character append & paste
  ! Modified Matthias Cuntz, Jan 2013 - removed 256 character restriction
  ! Modified Matthias Cuntz, Feb 2013 - logical append and paste

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

  ! Copyright 2012-2013 Juliane Mai, Matthias Cuntz

  USE mo_kind, only: i4, i8, sp, dp

  IMPLICIT NONE

#ifndef ABSOFT
  PUBLIC :: append    ! Returns input1 appended (on rows) with input2.  (like Unix cat)
  PUBLIC :: paste     ! Returns input1 pasted (on columns) with input2. (like Unix paste)

  ! ------------------------------------------------------------------

  !     NAME
  !         append

  !     PURPOSE
  !>        \brief Append (rows) scalars, vectors, and matrixes onto existing array.
  
  !>        \details Appends one input to the rows of another, i.e. append
  !>        on the first dimension.\n
  !>        The input might be a scalar, a vector or a matrix.\n
  !>        Possibilities are:\n
  !>        (1)     append scalar to vector\n
  !>        (2)     append vector to vector\n
  !>        (3)     append matrix to matrix

  !     CALLING SEQUENCE
  !         input1 = (/ 1.0_dp , 2.0_dp /)
  !         input2 = 3.0_dp

  !         call append(input1, input2)
  !         --> input1 = (/ 1.0_dp , 2.0_dp, 3.0_dp /)

  !         See also test folder for a detailed example.


  !     INTENT(IN)
  !>        \param[in] "input2" values to append. Can be INTEGER(I4/I8), REAL(SP/DP), or CHARACTER(len=*)
  !>                            and also scalar, DIMENSION(:), or DIMENSION(:,:)\n
  !>                            If not scalar then the columns have to agree with input1

  !     INTENT(INOUT)
  !>        \param[in,out] "allocatable :: input1" array to be appended to. Can be INTEGER(I4/I8), REAL(SP/DP),
  !>                            or CHARACTER(len=*). Must be DIMENSION(:) or DIMENSION(:,:), and allocatable.\n
  !>                            If input2 is not scalar then it must be size(input1,2) = size(input2,2).

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Size of input1 and input2 have to fit together,
  !         i.e. number of columns input1 = number of columns input2

  !     EXAMPLE
  !         see test/test_mo_append/

  !     LITERATURE

  !     HISTORY
  !>       \author Juliane Mai
  !>       \date Aug 2012
  !        Modified Matthias Cuntz, Jan 2013 - removed 256 character restriction
  !        Modified Matthias Cuntz, Feb 2013 - logical append and paste
  !        Modified Matthias Zink,  Feb 2015 - added optional 'fill_value' for logical append

  
  INTERFACE append
     MODULE PROCEDURE append_i4_v_s, append_i4_v_v, append_i4_m_m, &
          append_i8_v_s, append_i8_v_v, append_i8_m_m, &
          append_sp_v_s, append_sp_v_v, append_sp_m_m, &
          append_dp_v_s, append_dp_v_v, append_dp_m_m, &
          append_char_v_s, append_char_v_v, append_char_m_m, &
          append_lgt_v_s, append_lgt_v_v, append_lgt_m_m

  END INTERFACE append

  ! ------------------------------------------------------------------

  !     NAME
  !         paste

  !     PURPOSE
  !>        \brief Paste (columns) scalars, vectors, and matrixes onto existing array.
  
  !>        \details Pastes one input to the columns of another, i.e. append
  !>        on the second dimension.\n
  !>        The input might be a scalar, a vector or a matrix.\n
  !>        Possibilities are:\n
  !>        (1)     paste scalar to one-line matrix\n
  !>        (3)     paste vector to a matrix\n
  !>        (5)     paste matrix to matrix

  !     CALLING SEQUENCE
  !         input1 = (/ 1.0_dp , 2.0_dp /)
  !         input2 = (/ 3.0_dp , 4.0_dp /)

  !         call paste(input1, input2)
  !         --> input1(1,:) = (/ 1.0_dp , 3.0_dp /)
  !             input1(2,:) = (/ 2.0_dp , 4.0_dp /)

  !         See also test folder for a detailed example.

  !     INTENT(IN)
  !>        \param[in] "input2" values to paste. Can be INTEGER(I4/I8), REAL(SP/DP), or CHARACTER(len=*)
  !>                            and also scalar, DIMENSION(:), or DIMENSION(:,:)\n
  !>                            If not scalar then the rows have to agree with input1

  !     INTENT(INOUT)
  !>        \param[in,out] "allocatable :: input1" array to be pasted to. Can be INTEGER(I4/I8), REAL(SP/DP),
  !>                            or CHARACTER(len=*). Must be DIMENSION(:) or DIMENSION(:,:), and allocatable.\n
  !>                            If input2 is not scalar then it must be size(input1,1) = size(input2,1).

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Size of input1 and input2 have to fit together,
  !         i.e. number of rows input1 = number of rows input2

  !     EXAMPLE
  !         see test/test_mo_append/

  !     LITERATURE

  !     HISTORY
  !>       \author Juliane Mai
  !>       \date Aug 2012
  !        Modified Matthias Cuntz, Jan 2013 - removed 256 character restriction
  !        Modified Matthias Cuntz, Feb 2013 - logical append and paste

  INTERFACE paste
     MODULE PROCEDURE paste_i4_m_s, paste_i4_m_v, paste_i4_m_m, &
          paste_i8_m_s, paste_i8_m_v, paste_i8_m_m, &
          paste_sp_m_s, paste_sp_m_v, paste_sp_m_m, &
          paste_dp_m_s, paste_dp_m_v, paste_dp_m_m, &
          paste_char_m_s, paste_char_m_v, paste_char_m_m, &
          paste_lgt_m_s, paste_lgt_m_v, paste_lgt_m_m

  END INTERFACE paste

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  SUBROUTINE append_i4_v_s(vec1, sca2)

    implicit none

    integer(i4), dimension(:), allocatable, intent(inout)   :: vec1
    integer(i4),                            intent(in)      :: sca2

    ! local variables
    integer(i4)                             :: n1, n2
    integer(i4), dimension(:), allocatable  :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_i4_v_s

  SUBROUTINE append_i4_v_v(vec1, vec2)

    implicit none

    integer(i4), dimension(:), allocatable, intent(inout)   :: vec1
    integer(i4), dimension(:), intent(in)                   :: vec2

    ! local variables
    integer(i4)                             :: n1, n2    ! length of vectors
    integer(i4), dimension(:), allocatable  :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_i4_v_v

  SUBROUTINE append_i4_m_m(mat1, mat2, fill_value)

    implicit none

    integer(i4), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i4), dimension(:,:), intent(in)                   :: mat2
    integer(i4), optional,       intent(in)                   :: fill_value

    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                               :: n1, n2    ! dim2 of matrixes: columns
    integer(i4), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)    ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ((n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if

       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if
          
   else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_i4_m_m

  SUBROUTINE append_i8_v_s(vec1, sca2)

    implicit none

    integer(i8), dimension(:), allocatable, intent(inout)   :: vec1
    integer(i8),                            intent(in)      :: sca2

    ! local variables
    integer(i4)                             :: n1, n2
    integer(i8), dimension(:), allocatable  :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_i8_v_s

  SUBROUTINE append_i8_v_v(vec1, vec2)

    implicit none

    integer(i8), dimension(:), allocatable, intent(inout)   :: vec1
    integer(i8), dimension(:), intent(in)                   :: vec2

    ! local variables
    integer(i4)                             :: n1, n2    ! length of vectors
    integer(i8), dimension(:), allocatable  :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_i8_v_v

  SUBROUTINE append_i8_m_m(mat1, mat2, fill_value)

    implicit none

    integer(i8), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i8), dimension(:,:), intent(in)                   :: mat2
    integer(i8), optional,       intent(in)                   :: fill_value
    
    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                               :: n1, n2    ! dim2 of matrixes: columns
    integer(i8), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)    ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ((n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if

       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if

    else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_i8_m_m

  SUBROUTINE append_sp_v_s(vec1, sca2)

    implicit none

    real(sp), dimension(:), allocatable, intent(inout)   :: vec1
    real(sp),                            intent(in)      :: sca2

    ! local variables
    integer(i4)                             :: n1, n2
    real(sp), dimension(:), allocatable     :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_sp_v_s

  SUBROUTINE append_sp_v_v(vec1, vec2)

    implicit none

    real(sp), dimension(:), allocatable, intent(inout)   :: vec1
    real(sp), dimension(:),              intent(in)      :: vec2

    ! local variables
    integer(i4)                             :: n1, n2    ! length of vectors
    real(sp), dimension(:), allocatable     :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_sp_v_v

  SUBROUTINE append_sp_m_m(mat1, mat2, fill_value)

    implicit none

    real(sp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(sp), dimension(:,:), intent(in)                   :: mat2
    real(sp), optional,       intent(in)                   :: fill_value

    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                               :: n1, n2    ! dim2 of matrixes: columns
    real(sp), dimension(:,:), allocatable     :: tmp

    m2 = size(mat2,1)    ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ((n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

      if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if

       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if

    else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_sp_m_m

  SUBROUTINE append_dp_v_s(vec1, sca2)

    implicit none

    real(dp), dimension(:), allocatable, intent(inout)   :: vec1
    real(dp),                            intent(in)      :: sca2

    ! local variables
    integer(i4)                             :: n1, n2
    real(dp), dimension(:), allocatable     :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_dp_v_s

  SUBROUTINE append_dp_v_v(vec1, vec2)

    implicit none

    real(dp), dimension(:), allocatable, intent(inout)   :: vec1
    real(dp), dimension(:), intent(in)                   :: vec2

    ! local variables
    integer(i4)                             :: n1, n2    ! length of vectors
    real(dp), dimension(:), allocatable     :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_dp_v_v

  SUBROUTINE append_dp_m_m(mat1, mat2, fill_value)

    implicit none

    real(dp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(dp), dimension(:,:), intent(in)                   :: mat2
    real(dp), optional,       intent(in)                   :: fill_value

    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                               :: n1, n2    ! dim2 of matrixes: columns
    real(dp), dimension(:,:), allocatable     :: tmp

    m2 = size(mat2,1)    ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ((n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if

       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if

    else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_dp_m_m

  SUBROUTINE append_char_v_s(vec1, sca2)

    implicit none

    character(len=*), dimension(:), allocatable, intent(inout)   :: vec1
    character(len=*),                            intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: n1, n2
    character(len(vec1)), dimension(:), allocatable :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_char_v_s

  SUBROUTINE append_char_v_v(vec1, vec2)

    character(len=*), dimension(:), allocatable, intent(inout)   :: vec1
    character(len=*), dimension(:),              intent(in)      :: vec2

    ! local variables
    integer(i4)                               :: n1, n2
    character(len(vec1)), dimension(:), allocatable :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_char_v_v

  SUBROUTINE append_char_m_m(mat1, mat2, fill_value)

    implicit none

    character(len=*), dimension(:,:), allocatable, intent(inout)   :: mat1
    character(len=*), dimension(:,:),              intent(in)      :: mat2
    character(len=*), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                                 :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                                 :: n1, n2    ! dim2 of matrixes: columns
    character(len(mat1)), dimension(:,:), allocatable :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ((n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if

       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if

   else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_char_m_m

  SUBROUTINE append_lgt_v_s(vec1, sca2)

    implicit none

    logical, dimension(:), allocatable, intent(inout)   :: vec1
    logical,                            intent(in)      :: sca2

    ! local variables
    integer(i4)                             :: n1, n2
    logical, dimension(:), allocatable  :: tmp

    n2 = 1_i4

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4)       = sca2
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(1_i4) = sca2
    end if

  END SUBROUTINE append_lgt_v_s

  SUBROUTINE append_lgt_v_v(vec1, vec2)

    implicit none

    logical, dimension(:), allocatable, intent(inout)   :: vec1
    logical, dimension(:), intent(in)                   :: vec2

    ! local variables
    integer(i4)                             :: n1, n2    ! length of vectors
    logical, dimension(:), allocatable  :: tmp

    n2 = size(vec2)

    if (allocated(vec1)) then
       n1 = size(vec1)
       ! save vec1
       allocate(tmp(n1))
       tmp=vec1
       deallocate(vec1)

       allocate(vec1(n1+n2))
       vec1(1:n1)          = tmp(1:n1)
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    else
       n1 = 0_i4

       allocate(vec1(n2))
       vec1(n1+1_i4:n1+n2) = vec2(1:n2)
    end if

  END SUBROUTINE append_lgt_v_v

  SUBROUTINE append_lgt_m_m(mat1, mat2, fill_value)

    implicit none

    logical, dimension(:,:), allocatable, intent(inout)   :: mat1
    logical, dimension(:,:), intent(in)                   :: mat2
    logical, optional,       intent(in)                   :: fill_value
    
    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes: rows
    integer(i4)                               :: n1, n2    ! dim2 of matrixes: columns
    logical, dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)    ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns

       if ( (n1 .ne. n2) .and. .not. present(fill_value) ) then
          print*, 'append: columns of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if

       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( n1 .eq. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)          = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,:) = mat2(1:m2,:)
       end if

       if ( n1 .gt. n2 ) then
          allocate(mat1(m1+m2,n1))
          mat1(1:m1,:)                = tmp(1:m1,:)
          mat1(m1+1_i4:m1+m2,   1:n2) = mat2(1:m2,:)
          mat1(m1+1_i4:m1+m2,n2+1:n1) = fill_value
       end if
       
       if ( n1 .lt. n2 ) then
          allocate(mat1(m1+m2,n2))
          mat1(      1:m1,      1:n1) = tmp(1:m1,:)
          mat1(      1:m1,   n1+1:n2) = fill_value
          mat1(m1+1_i4:m1+m2,    :  ) = mat2(1:m2,:)
       end if
       
    else
       n1 = 0_i4

       allocate(mat1(m2,n2))
       mat1 = mat2
    end if

  END SUBROUTINE append_lgt_m_m

  ! ------------------------------------------------------------------

  SUBROUTINE paste_i4_m_s(mat1, sca2)

    implicit none

    integer(i4), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i4),                              intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: m1    ! dim1 of matrix
    integer(i4)                               :: n1    ! dim2 of matrix
    integer(i4), dimension(:,:), allocatable  :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_i4_m_s

  SUBROUTINE paste_i4_m_v(mat1, vec2, fill_value)

    implicit none

    integer(i4), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i4), dimension(:),                intent(in)      :: vec2
    integer(i4), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    integer(i4), dimension(:,:), allocatable  :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(1:m1,1:n1)       = tmp(:,1:n1)
          mat1(1:m2,n1+n2)      = vec2(1:m2)
       end if
       
       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
          mat1(m2+1:m1,n1+n2)   = fill_value
       end if

       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(m1+1:m2,1:n1)    = fill_value
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_i4_m_v

  SUBROUTINE paste_i4_m_m(mat1, mat2, fill_value)

    implicit none

    integer(i4), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i4), dimension(:,:),              intent(in)      :: mat2
    integer(i4), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    integer(i4), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(:,1:n1)          = tmp(:,1:n1)
          mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(    :  ,1:n1)          = tmp(:,1:n1)
          mat1(   1:m2,n1+1_i4:n1+n2) = mat2(:,1:n2)
          mat1(m2+1:m1,n1+1_i4:n1+n2) = fill_value
       end if
       
       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,      1:n1   ) = tmp(:,1:n1)
          mat1(m1+1:m2,      1:n1   ) = fill_value
          mat1(    :  ,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if
          
    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_i4_m_m

  SUBROUTINE paste_i8_m_s(mat1, sca2)

    implicit none

    integer(i8), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i8),                              intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: m1    ! dim1 of matrix
    integer(i4)                               :: n1    ! dim2 of matrix
    integer(i8), dimension(:,:), allocatable  :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_i8_m_s

  SUBROUTINE paste_i8_m_v(mat1, vec2, fill_value)

    implicit none

    integer(i8), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i8), dimension(:),                intent(in)      :: vec2
    integer(i8), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                               :: m1, m2    ! dim1 of matrixes
    integer(i4)                               :: n1, n2    ! dim2 of matrixes
    integer(i8), dimension(:,:), allocatable  :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(1:m1,1:n1)       = tmp(:,1:n1)
          mat1(1:m2,n1+n2)      = vec2(1:m2)
       end if
       
       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
          mat1(m2+1:m1,n1+n2)   = fill_value
       end if

       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(m1+1:m2,1:n1)    = fill_value
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_i8_m_v

  SUBROUTINE paste_i8_m_m(mat1, mat2, fill_value)

    implicit none

    integer(i8), dimension(:,:), allocatable, intent(inout)   :: mat1
    integer(i8), dimension(:,:),              intent(in)      :: mat2
    integer(i8), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    integer(i8), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(:,1:n1)          = tmp(:,1:n1)
          mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(    :  ,1:n1)          = tmp(:,1:n1)
          mat1(   1:m2,n1+1_i4:n1+n2) = mat2(:,1:n2)
          mat1(m2+1:m1,n1+1_i4:n1+n2) = fill_value
       end if
       
       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,      1:n1   ) = tmp(:,1:n1)
          mat1(m1+1:m2,      1:n1   ) = fill_value
          mat1(    :  ,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if
          
   else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_i8_m_m

  SUBROUTINE paste_sp_m_s(mat1, sca2)

    implicit none

    real(sp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(sp),                              intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: m1    ! dim1 of matrix
    integer(i4)                               :: n1    ! dim2 of matrix
    real(sp), dimension(:,:), allocatable  :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_sp_m_s

  SUBROUTINE paste_sp_m_v(mat1, vec2, fill_value)

    implicit none

    real(sp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(sp), dimension(:),                intent(in)      :: vec2
    real(sp), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    real(sp), dimension(:,:), allocatable  :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(1:m1,1:n1)       = tmp(:,1:n1)
          mat1(1:m2,n1+n2)      = vec2(1:m2)
       end if
       
       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
          mat1(m2+1:m1,n1+n2)   = fill_value
       end if

       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(m1+1:m2,1:n1)    = fill_value
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_sp_m_v

  SUBROUTINE paste_sp_m_m(mat1, mat2, fill_value)

    implicit none

    real(sp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(sp), dimension(:,:),              intent(in)      :: mat2
    real(sp), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    real(sp), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(:,1:n1)          = tmp(:,1:n1)
          mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(    :  ,1:n1)          = tmp(:,1:n1)
          mat1(   1:m2,n1+1_i4:n1+n2) = mat2(:,1:n2)
          mat1(m2+1:m1,n1+1_i4:n1+n2) = fill_value
       end if
       
       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,      1:n1   ) = tmp(:,1:n1)
          mat1(m1+1:m2,      1:n1   ) = fill_value
          mat1(    :  ,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if


    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_sp_m_m

  SUBROUTINE paste_dp_m_s(mat1, sca2)

    implicit none

    real(dp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(dp),                              intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: m1    ! dim1 of matrix
    integer(i4)                               :: n1    ! dim2 of matrix
    real(dp), dimension(:,:), allocatable  :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_dp_m_s

  SUBROUTINE paste_dp_m_v(mat1, vec2, fill_value)

    implicit none

    real(dp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(dp), dimension(:),                intent(in)      :: vec2
    real(dp), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    real(dp), dimension(:,:), allocatable  :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(1:m1,1:n1)       = tmp(:,1:n1)
          mat1(1:m2,n1+n2)      = vec2(1:m2)
       end if
       
       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
          mat1(m2+1:m1,n1+n2)   = fill_value
       end if

       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(m1+1:m2,1:n1)    = fill_value
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_dp_m_v

  SUBROUTINE paste_dp_m_m(mat1, mat2, fill_value)

    implicit none

    real(dp), dimension(:,:), allocatable, intent(inout)   :: mat1
    real(dp), dimension(:,:),              intent(in)      :: mat2
    real(dp), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    real(dp), dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(:,1:n1)          = tmp(:,1:n1)
          mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(    :  ,1:n1)          = tmp(:,1:n1)
          mat1(   1:m2,n1+1_i4:n1+n2) = mat2(:,1:n2)
          mat1(m2+1:m1,n1+1_i4:n1+n2) = fill_value
       end if
       
       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,      1:n1   ) = tmp(:,1:n1)
          mat1(m1+1:m2,      1:n1   ) = fill_value
          mat1(    :  ,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_dp_m_m

  SUBROUTINE paste_char_m_s(mat1, sca2)

    implicit none

    character(len=*), dimension(:,:), allocatable, intent(inout)   :: mat1
    character(len=*),                              intent(in)      :: sca2

    ! local variables
    integer(i4)                                  :: m1    ! dim1 of matrix
    integer(i4)                                  :: n1    ! dim2 of matrix
    character(len(mat1)), dimension(:,:), allocatable :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_char_m_s

  SUBROUTINE paste_char_m_v(mat1, vec2, fill_value)

    implicit none

    character(len=*), dimension(:,:), allocatable, intent(inout)   :: mat1
    character(len=*), dimension(:),                intent(in)      :: vec2
    character(len=*), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                                  :: m1, m2    ! dim1 of matrixes
    integer(i4)                                  :: n1, n2    ! dim2 of matrixes
    character(len(mat1)), dimension(:,:), allocatable :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(1:m1,1:n1)       = tmp(:,1:n1)
          mat1(1:m2,n1+n2)      = vec2(1:m2)
       end if
       
       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
          mat1(m2+1:m1,n1+n2)   = fill_value
       end if

       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,1:n1)    = tmp(:,1:n1)
          mat1(m1+1:m2,1:n1)    = fill_value
          mat1(   1:m2,n1+n2)   = vec2(1:m2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_char_m_v

  SUBROUTINE paste_char_m_m(mat1, mat2, fill_value)

    implicit none

    character(len=*), dimension(:,:), allocatable, intent(inout)   :: mat1
    character(len=*), dimension(:,:),              intent(in)      :: mat2
    character(len=*), optional,                    intent(in)      :: fill_value

    ! local variables
    integer(i4)                                  :: m1, m2    ! dim1 of matrixes
    integer(i4)                                  :: n1, n2    ! dim2 of matrixes
    character(len(mat1)), dimension(:,:), allocatable :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if ( (m1 .ne. m2) .and. .not. present( fill_value ) ) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       if ( m1 .eq. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(:,1:n1)          = tmp(:,1:n1)
          mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

       if ( m1 .gt. m2 ) then
          allocate(mat1(m1,n1+n2))
          mat1(    :  ,1:n1)          = tmp(:,1:n1)
          mat1(   1:m2,n1+1_i4:n1+n2) = mat2(:,1:n2)
          mat1(m2+1:m1,n1+1_i4:n1+n2) = fill_value
       end if
       
       if ( m1 .lt. m2 ) then
          allocate(mat1(m2,n1+n2))
          mat1(   1:m1,      1:n1   ) = tmp(:,1:n1)
          mat1(m1+1:m2,      1:n1   ) = fill_value
          mat1(    :  ,n1+1_i4:n1+n2) = mat2(:,1:n2)
       end if

    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_char_m_m

  SUBROUTINE paste_lgt_m_s(mat1, sca2)

    implicit none

    logical, dimension(:,:), allocatable, intent(inout)   :: mat1
    logical,                              intent(in)      :: sca2

    ! local variables
    integer(i4)                               :: m1    ! dim1 of matrix
    integer(i4)                               :: n1    ! dim2 of matrix
    logical, dimension(:,:), allocatable  :: tmp

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. 1_i4) then
          print*, 'paste: scalar paste to matrix only works with one-line matrix'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(1_i4,n1+1_i4))
       mat1(1,1:n1)          = tmp(1,1:n1)
       mat1(1,n1+1_i4)       = sca2
    else
       allocate(mat1(1_i4,1_i4))
       mat1(1,1) = sca2
    end if

  END SUBROUTINE paste_lgt_m_s

  SUBROUTINE paste_lgt_m_v(mat1, vec2)

    implicit none

    logical, dimension(:,:), allocatable, intent(inout)   :: mat1
    logical, dimension(:),                intent(in)      :: vec2

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    logical, dimension(:,:), allocatable  :: tmp

    m2 = size(vec2,1)   ! rows
    n2 = 1_i4           ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. m2) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(m1,n1+n2))
       mat1(:,1:n1)          = tmp(:,1:n1)
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(1:m2,n1+n2)      = vec2(1:m2)
    end if

  END SUBROUTINE paste_lgt_m_v

  SUBROUTINE paste_lgt_m_m(mat1, mat2)

    implicit none

    logical, dimension(:,:), allocatable, intent(inout)   :: mat1
    logical, dimension(:,:),              intent(in)      :: mat2

    ! local variables
    integer(i4)                             :: m1, m2    ! dim1 of matrixes
    integer(i4)                             :: n1, n2    ! dim2 of matrixes
    logical, dimension(:,:), allocatable  :: tmp

    m2 = size(mat2,1)   ! rows
    n2 = size(mat2,2)   ! columns

    if (allocated(mat1)) then
       m1 = size(mat1,1)   ! rows
       n1 = size(mat1,2)   ! columns
       if (m1 .ne. m2) then
          print*, 'paste: rows of matrix1 and matrix2 are unequal : (',m1,',',n1,')  and  (',m2,',',n2,')'
          STOP
       end if
       ! save mat1
       allocate(tmp(m1,n1))
       tmp=mat1
       deallocate(mat1)

       allocate(mat1(m1,n1+n2))
       mat1(:,1:n1)          = tmp(:,1:n1)
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    else
       n1 = 0_i4
       m1 = m2

       allocate(mat1(m2,n2))
       mat1(:,n1+1_i4:n1+n2) = mat2(:,1:n2)
    end if

  END SUBROUTINE paste_lgt_m_m
#endif


END MODULE mo_append
