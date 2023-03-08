!> \file mo_common_types.F90
!> \brief \copybrief mo_common_types
!> \details \copydetails mo_common_types

!> \brief Provides common types needed by mHM, mRM and/or mpr.
!> \changelog
!! - Stephan Thober  Nov 2016
!!   - moved processdescription from mo_global_variables to here
!! - Robert Schweppe Dec 2017
!!   - merged more duplicated variables from mhm and mrm global variables
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!! - Sebastian Müller Mar 2023
!!   - moving types to a separate module
!> \authors Stephan Thober
!> \date Sep 2015
!> \authors Sebastian Müller
!> \date Mar 2023
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_types

  use mo_kind, only : i4, dp
#ifdef MPI
  USE mpi_f08
#endif
  implicit none

  ! -------------------------------------------------------------------
  ! PERIOD description
  ! -------------------------------------------------------------------
  !> \class   period
  !> \brief   period description
  type, public :: period
    integer(i4) :: dStart      !< first day
    integer(i4) :: mStart      !< first month
    integer(i4) :: yStart      !< first year
    integer(i4) :: dEnd        !< last  day
    integer(i4) :: mEnd        !< last  month
    integer(i4) :: yEnd        !< last  year
    integer(i4) :: julStart    !< first julian day
    integer(i4) :: julEnd      !< last  julian day
    integer(i4) :: nObs        !< total number of observations
  end type period

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  !> \class   grid
  !> \brief   grid description
  type, public :: grid
    ! general domain information
    integer(i4) :: ncols     !< Number of columns
    integer(i4) :: nrows     !< Number of rows
    integer(i4) :: nCells    !< Number of cells in mask
    real(dp) :: xllcorner    !< x coordinate of the lowerleft corner
    real(dp) :: yllcorner    !< y coordinate of the lowerleft corner
    real(dp) :: cellsize     !< Cellsize x = cellsize y
    real(dp) :: nodata_value !< Code to define the mask
    real(dp), dimension(:, :), allocatable :: x  !< 2d longitude array (unmasked version is needed for output anyway)
    real(dp), dimension(:, :), allocatable :: y  !< 2d latitude  array (unmasked version is needed for output anyway)
    logical, dimension(:, :), allocatable :: mask  !< the mask for valid cells in the original grid (nrows*ncols)
    ! for referencing values in the nValidCells vector
    integer(i4) :: iStart          !< Starting cell index of a given domain
    integer(i4) :: iEnd            !< Ending cell index of a given domain
    ! dimension(nCells, (x,y) )
    integer(i4), dimension(:, :), allocatable :: CellCoor  !< this is only used for mRM
    real(dp), dimension(:), allocatable :: CellArea  !< area of the cell in sq m
    integer(i4), dimension(:), allocatable :: Id !< id

  end type grid

  !> \class   gridremapper
  !> \brief   grid remapper
  type, public :: gridremapper
    type(Grid), pointer :: high_res_grid !< high resolution grid
    type(Grid), pointer :: low_res_grid !< low resolution grid

    ! dimension nCells
    integer(i4), dimension(:), allocatable :: lower_bound  !< 1d index of lower side subgrid
    integer(i4), dimension(:), allocatable :: upper_bound  !< 1d index of upper side subgrid
    integer(i4), dimension(:), allocatable :: left_bound   !< 1d index of left side subgrid
    integer(i4), dimension(:), allocatable :: right_bound  !< 1d index of right side subgrid
    integer(i4), dimension(:), allocatable :: n_subcells   !< 1d numberof valid subgrid cells
    integer(i4), dimension(:, :), allocatable :: lowres_id_on_highres   !< 2d index array of lowres id

  end type gridremapper

  ! -------------------------------------------------------------------
  ! DOMAIN general description
  ! -------------------------------------------------------------------
  !> \class   domain_meta
  !> \brief   DOMAIN general description
  type, public :: domain_meta
    integer(i4)                            :: nDomains !< number of domains
    integer(i4)                            :: overallNumberOfDomains  !< Number of domains for multi-domain optimization
    integer(i4), dimension(:), allocatable :: indices !< indices
    integer(i4), dimension(:), allocatable :: L0DataFrom !< index of associated level-0 domain
    ! optidata saves for each domain which optional data is assigned to it
    ! (0) default: the program decides. If you are confused, choose 0
    ! (1) runoff
    ! (2) sm
    ! (3) tws
    ! (4) neutons
    ! (5) et
    ! (6) et & tws
    integer(i4), dimension(:), allocatable :: optidata !< optidata flag (0-6)
    logical,     dimension(:), allocatable :: doRouting !< flag to indicate if routing is done
#ifdef MPI
    logical                                :: isMasterInComLocal  !< true if the process is master proc in comLocal
    !> the communicater the domains are using to send messages to each other here are all processes wich have rank 0 in comLocal
    type(MPI_Comm)                         :: comMaster
    type(MPI_Comm)                         :: comLocal  !< the communicater the domain internal communication takes place
#endif
  end type domain_meta

end module mo_common_types
