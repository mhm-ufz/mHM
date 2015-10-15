!> \file mo_global_structures.f90

!> \brief Provides structures needed by mHM and mRM.

!> \details Provides the global structure period that is used
!>     by both mHM and mRM.

!> \author Stephan Thober
!> \date Sep 2015
module mo_global_structures
  use mo_kind, only: i4
  implicit none
  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  type period
      integer(i4) :: dStart      ! first day
      integer(i4) :: mStart      ! first month
      integer(i4) :: yStart      ! first year
      integer(i4) :: dEnd        ! last  day
      integer(i4) :: mEnd        ! last  month
      integer(i4) :: yEnd        ! last  year
      integer(i4) :: julStart    ! first julian day 
      integer(i4) :: julEnd      ! last  julian day 
      integer(i4) :: nObs        ! total number of observations
  end type period

end module mo_global_structures
