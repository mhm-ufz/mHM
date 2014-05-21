!> \file mo_sort.f90

!> \brief Quicksort array or indices

!> \details This module contains a sorting routines from numerical recipes
!> 1. quicksort one array
!> 2. return indeces that would sort the array

!> \authors Matthias Cuntz - modified Numerical Recipes routines sort and indexx
!> \date March 2011

MODULE mo_sort

  ! This module contains a sorting routines from numerical recipes:
  ! 1. quicksort one array
  ! 2. return indeces that would sort the array

  ! Usage:
  !   USE mo_sort, ONLY: sort
  !   call sort(arr)
  ! Or if you want to sort more than one array
  !   USE mo_sort, ONLY: sort_index
  !   integer(i4), dimension(size(arr)) :: ii
  !   ii = sort_index(arr)
  !   arr = arr(ii)

  ! Written March 2011, Matthias Cuntz
  !   - modified numerical recipes: sort and indexx

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
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011 Matthias Cuntz


  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred.

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  USE mo_kind,   ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: sort       ! Quicksorts one array
  PUBLIC :: sort_index ! Returns indeces that sort the input array

  ! ------------------------------------------------------------------

  !     NAME
  !         sort

  !     PURPOSE
  !>        \brief Sorts the input array in ascending order.

  !     CALLING SEQUENCE
  !         call sort(arr)
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !>        \param[in,out] "real(sp/dp) :: arr(:)"     On input:  unsorted 1D-array\n
  !>                                                   On output: ascendingly sorted input array

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         No mask or undefined value provided.

  !     EXAMPLE
  !         vec = (/ 3., 2, 1., 4., 6., 5. /)
  !         call sort(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz - adapted routine sort from numerical recipes
  !>        \date Nov 2011
  INTERFACE sort
     MODULE PROCEDURE sort_i4, sort_sp, sort_dp
  END INTERFACE sort

  ! ------------------------------------------------------------------

  !     NAME
  !         sort_index

  !     PURPOSE
  !>        \brief Gives the indeces that would sort an array in ascending order.

  !     CALLING SEQUENCE
  !         ii = sort_index(arr)
  
  !     INDENT(IN)
  !>        \param[in] "real(sp/dp) :: arr(:)"     Unsorted 1D-array

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return integer(i4) :: indices(:) &mdash; Indices that would sort arr in ascending order

  !     RESTRICTIONS
  !         No mask or undefined value provided.

  !     EXAMPLE
  !         vec = (/ 3., 2, 1., 4., 6., 5. /)
  !         ii = sort_index(arr)
  !         arr = arr(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz - adapted routine indexx from numerical recipes
  !>        \date Nov 2011
  INTERFACE sort_index
     MODULE PROCEDURE sort_index_i4, sort_index_sp, sort_index_dp
  END INTERFACE sort_index

  PRIVATE

  ! Local routine for indeces exchange, from numerical recipes
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE arth
  INTERFACE icomp_xchg
     MODULE PROCEDURE icomp_xchg_i4, icomp_xchg_sp, icomp_xchg_dp
  END INTERFACE icomp_xchg
  INTERFACE swap
     MODULE PROCEDURE &
          swap_i,  swap_r,  swap_d, &  ! swap_c,  swap_z, &
          ! swap_iv, swap_rv, swap_dv, swap_cv, swap_zv, &
          ! swap_im, swap_rm, swap_dm, swap_cm, swap_zm, &
          masked_swap_is, & !masked_swap_iv, masked_swap_im, &
          masked_swap_rs, & !masked_swap_rv, masked_swap_rm, &
          masked_swap_ds !, & !masked_swap_dv, masked_swap_dm, &
          !  masked_swap_cs, masked_swap_cv, masked_swap_cm, &
          ! masked_swap_zs, masked_swap_zv, masked_swap_zm
  END INTERFACE swap

  INTEGER(I4), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! Array function returning an arithmetic progression.
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r

  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d

  FUNCTION arth_i(first,increment,n)
    INTEGER(I4), INTENT(IN) :: first,increment,n
    INTEGER(I4), DIMENSION(n) :: arth_i
    INTEGER(I4) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Swap or don''t swap integer arguments, depending on the ordering of their corresponding
  ! elements in an array arr.
  SUBROUTINE icomp_xchg_i4(arr, i, j)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_i4


  SUBROUTINE icomp_xchg_sp(arr, i, j)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_sp


  SUBROUTINE icomp_xchg_dp(arr, i, j)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_dp

! ------------------------------------------------------------------

  SUBROUTINE sort_i4(arr)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    INTEGER(i4) :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE sort_i4

  ! ------------------------------------------------------------------

  ! quicksort
  SUBROUTINE sort_sp(arr)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(sp)    :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE sort_sp


  SUBROUTINE sort_dp(arr)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(dp)    :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_dp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE sort_dp


  ! ------------------------------------------------------------------

  ! quicksort but return indices
  FUNCTION sort_index_sp(arr)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sort_index_sp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(sp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    sort_index_sp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_sp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(sort_index_sp(i)) <= a) exit
                sort_index_sp(i+1) = sort_index_sp(i)
             end do
             sort_index_sp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_sp(k), sort_index_sp(l+1))
          call icomp_xchg(arr, sort_index_sp(l),   sort_index_sp(r))
          call icomp_xchg(arr, sort_index_sp(l+1), sort_index_sp(r))
          call icomp_xchg(arr, sort_index_sp(l),   sort_index_sp(l+1))
          i = l+1
          j = r
          indext = sort_index_sp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(sort_index_sp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(sort_index_sp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_sp(i), sort_index_sp(j))
          end do
          sort_index_sp(l+1) = sort_index_sp(j)
          sort_index_sp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l=i
          end if
       end if
    end do

  END FUNCTION sort_index_sp


  FUNCTION sort_index_dp(arr)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sort_index_dp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(dp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    sort_index_dp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_dp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(sort_index_dp(i)) <= a) exit
                sort_index_dp(i+1) = sort_index_dp(i)
             end do
             sort_index_dp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_dp(k), sort_index_dp(l+1))
          call icomp_xchg(arr, sort_index_dp(l),   sort_index_dp(r))
          call icomp_xchg(arr, sort_index_dp(l+1), sort_index_dp(r))
          call icomp_xchg(arr, sort_index_dp(l),   sort_index_dp(l+1))
          i = l+1
          j = r
          indext = sort_index_dp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(sort_index_dp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(sort_index_dp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_dp(i), sort_index_dp(j))
          end do
          sort_index_dp(l+1) = sort_index_dp(j)
          sort_index_dp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_dp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l=i
          end if
       end if
    end do

  END FUNCTION sort_index_dp


  FUNCTION sort_index_i4(iarr)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(i4), DIMENSION(size(iarr))    :: sort_index_i4

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    INTEGER(i4) :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(iarr)
    sort_index_i4 = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_i4(j)
             a = iarr(indext)
             ! Corrected from Numerical Recipes Forum
             ! do i=j-1, 1, -1
             do i=j-1, l, -1
                if (iarr(sort_index_i4(i)) <= a) exit
                sort_index_i4(i+1) = sort_index_i4(i)
             end do
             sort_index_i4(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_i4(k), sort_index_i4(l+1))
          call icomp_xchg(iarr, sort_index_i4(l),   sort_index_i4(r))
          call icomp_xchg(iarr, sort_index_i4(l+1), sort_index_i4(r))
          call icomp_xchg(iarr, sort_index_i4(l),   sort_index_i4(l+1))
          i = l+1
          j = r
          indext = sort_index_i4(l+1)
          a = iarr(indext)
          do
             do
                i = i+1
                if (iarr(sort_index_i4(i)) >= a) exit
             end do
             do
                j = j-1
                if (iarr(sort_index_i4(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_i4(i), sort_index_i4(j))
          end do
          sort_index_i4(l+1) = sort_index_i4(j)
          sort_index_i4(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_i4: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END FUNCTION sort_index_i4

  ! ------------------------------------------------------------------

  ! ! Swap the contents of a and b.
  SUBROUTINE swap_i(a,b)
    INTEGER(I4), INTENT(INOUT) :: a,b
    INTEGER(I4) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_d(a,b)
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_d

  ! SUBROUTINE swap_c(a,b)
  !   COMPLEX(SPC), INTENT(INOUT) :: a,b
  !   COMPLEX(SPC) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_c

  ! SUBROUTINE swap_z(a,b)
  !   COMPLEX(DPC), INTENT(INOUT) :: a,b
  !   COMPLEX(DPC) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_z

  ! SUBROUTINE swap_iv(a,b)
  !   INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: a,b
  !   INTEGER(I4), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_iv

  ! SUBROUTINE swap_rv(a,b)
  !   REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
  !   REAL(SP), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_rv

  ! SUBROUTINE swap_dv(a,b)
  !   REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
  !   REAL(DP), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_dv

  ! SUBROUTINE swap_cv(a,b)
  !   COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
  !   COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_cv

  ! SUBROUTINE swap_zv(a,b)
  !   COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
  !   COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_zv

  ! SUBROUTINE swap_im(a,b)
  !   INTEGER(I4), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   INTEGER(I4), DIMENSION(size(a,1),size(a,2)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_im

  ! SUBROUTINE swap_rm(a,b)
  !   REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   REAL(SP), DIMENSION(size(a,1),size(a,2)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_rm

  ! SUBROUTINE swap_dm(a,b)
  !   REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   REAL(DP), DIMENSION(size(a,1),size(a,2)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_dm

  ! SUBROUTINE swap_cm(a,b)
  !   COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_cm

  ! SUBROUTINE swap_zm(a,b)
  !   COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
  !   dum=a
  !   a=b
  !   b=dum
  ! END SUBROUTINE swap_zm

  SUBROUTINE masked_swap_is(a,b,mask)
    INTEGER(I4), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    INTEGER(I4) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_is

  ! SUBROUTINE masked_swap_iv(a,b,mask)
  !   INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   INTEGER(I4), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_iv

  ! SUBROUTINE masked_swap_im(a,b,mask)
  !   INTEGER(I4), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   INTEGER(I4), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_im

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  ! SUBROUTINE masked_swap_rv(a,b,mask)
  !   REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   REAL(SP), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_rv

  ! SUBROUTINE masked_swap_rm(a,b,mask)
  !   REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_rm

  SUBROUTINE masked_swap_ds(a,b,mask)
    REAL(DP), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(DP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_ds

  ! SUBROUTINE masked_swap_dv(a,b,mask)
  !   REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   REAL(DP), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_dv

  ! SUBROUTINE masked_swap_dm(a,b,mask)
  !   REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   REAL(DP), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_dm

  ! SUBROUTINE masked_swap_cs(a,b,mask)
  !   COMPLEX(SPC), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   COMPLEX(SPC) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_cs

  ! SUBROUTINE masked_swap_cv(a,b,mask)
  !   COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   COMPLEX(SPC), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_cv

  ! SUBROUTINE masked_swap_cm(a,b,mask)
  !   COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_cm

  ! SUBROUTINE masked_swap_zs(a,b,mask)
  !   COMPLEX(DPC), INTENT(INOUT) :: a,b
  !   LOGICAL, INTENT(IN) :: mask
  !   COMPLEX(DPC) :: swp
  !   if (mask) then
  !      swp=a
  !      a=b
  !      b=swp
  !   end if
  ! END SUBROUTINE masked_swap_zs

  ! SUBROUTINE masked_swap_zv(a,b,mask)
  !   COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  !   COMPLEX(DPC), DIMENSION(size(a)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_zv

  ! SUBROUTINE masked_swap_zm(a,b,mask)
  !   COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
  !   LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
  !   COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: swp
  !   where (mask)
  !      swp=a
  !      a=b
  !      b=swp
  !   end where
  ! END SUBROUTINE masked_swap_zm

  ! ------------------------------------------------------------------

END MODULE mo_sort
