!> \file mo_mrm_tools.f90

!> \brief Provide utility routines used within mRM.

!> \details This module contains subroutines that are used frequently
!> to obtain basin and grid properties.

!> \author Luis Samaniego
!> \date Dec 2012
!  Modified
!       Stephan Thober, Aug 2015 - adapted to mRM
module mo_mrm_tools
  implicit none
  public :: get_basin_info_mrm
  public :: calculate_grid_properties
  private
contains
  ! ------------------------------------------------------------------

  !      NAME
  !          get_basin_info_mrm

  !>        \brief Get basic basin information (e.g., nrows, ncols, indices, mask)
  
  !>        \details Get basic basin information (e.g., nrows, ncols, indices, mask) for
  !>                 different levels (L0, L1, L11, and L110).
  !
  !     CALLING SEQUENCE
  !         call get_basin_info_mrm(iBasin, iLevel,nrows,ncols, ncells, iStart, iEnd, &
  !                             iStartMask, iEndMask, mask, xllcorner, yllcorner, cellsize) 

  !     INTENT(IN)
  !>        \param[in] "integer(i4)             :: iBasin"    basin id
  !>        \param[in] "integer(i4)             :: iLevel"    level id (e.g., 0, 1, 11, 2)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)             :: nRows"      no. of rows
  !>        \param[out] "integer(i4)             :: nCols"      no. of coloums
  !>        \param[out] "integer(i4)             :: ncells"     no. of cells
  !>        \param[out] "integer(i4)             :: iStart"     start cell index of a given basin at a given level 
  !>        \param[out] "integer(i4)             :: iEnd"       end cell index of a given basin at a given level
  !>        \param[out] "integer(i4)             :: iStartMask" start cell index of mask a given basin at a given level
  !>        \param[out] "integer(i4)             :: iEndMask"   end cell index of mask a given basin at a given level

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL                   
  !>        \param[out] "logical, optional       :: mask"       Mask at a given level
  !>        \param[out] "real(dp), optional      :: xllcorner"  x coordinate of the lowerleft corner at a given level
  !>        \param[out] "real(dp), optional      :: yllcorner"  y coordinate of the lowerleft corner at a given level
  !>        \param[out] "real(dp), optional      :: cellsize"   cell size at a given level

  !     RETURN

  !     RESTRICTIONS
  !>    \note This subroutine cannot be called for level 2, that does not exist within mRM.

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         \authors  Rohini Kumar, Luis Samaniego
  !         \date     Jan 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template
  !                   Stephan Thober, Aug 2015 - adapted for mRM

  subroutine get_basin_info_mrm(iBasin, iLevel, nrows, ncols, ncells, iStart, iEnd, &
                            iStartMask, iEndMask, mask, xllcorner, yllcorner, cellsize) 
    use mo_kind, only: i4, dp
    use mo_message, only: message
    use mo_mrm_global_variables, only: level0, level1, level11, basin_mrm
    implicit none

    integer(i4), intent(in)                                      :: iBasin
    integer(i4), intent(in)                                      :: iLevel
    integer(i4), intent(out)                                     :: nrows, ncols
    integer(i4), optional, intent(out)                           :: ncells
    integer(i4), optional, intent(out)                           :: iStart, iEnd
    integer(i4), optional, intent(out)                           :: iStartMask, iEndMask
    logical, optional, dimension(:,:), allocatable,  intent(out) :: mask
    real(dp), optional, intent(out)                              :: xllcorner, yllcorner, cellsize

    ! level information
    select case (iLevel)
    case (0)
       nrows = level0%nrows(iBasin)
       ncols = level0%ncols(iBasin)
       if (present(ncells)) ncells = basin_mrm%L0_iEnd(iBasin) - basin_mrm%L0_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin_mrm%L0_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin_mrm%L0_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin_mrm%L0_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin_mrm%L0_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin_mrm%L0_Mask( basin_mrm%L0_iStartMask(iBasin): basin_mrm%L0_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level0%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level0%yllcorner(iBasin)
       if (present(cellsize))  cellsize  = level0%cellsize(iBasin)

    case (1)
       nrows = level1%nrows(iBasin)
       ncols = level1%ncols(iBasin)
       if (present(ncells)) ncells = basin_mrm%L1_iEnd(iBasin) - basin_mrm%L1_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin_mrm%L1_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin_mrm%L1_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin_mrm%L1_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin_mrm%L1_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin_mrm%L1_Mask( basin_mrm%L1_iStartMask(iBasin): basin_mrm%L1_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level1%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level1%yllcorner(iBasin)
       if (present(cellsize))  cellsize = level1%cellsize(iBasin)

    case (11)
       nrows = level11%nrows(iBasin)
       ncols = level11%ncols(iBasin)
       if (present(ncells)) ncells = basin_mrm%L11_iEnd(iBasin) - basin_mrm%L11_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin_mrm%L11_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin_mrm%L11_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin_mrm%L11_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin_mrm%L11_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin_mrm%L11_Mask( basin_mrm%L11_iStartMask(iBasin): basin_mrm%L11_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level11%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level11%yllcorner(iBasin)
       if (present(cellsize)) cellsize   = level11%cellsize(iBasin)

    case (110)
       if (present(iStart)) iStart = basin_mrm%L110_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin_mrm%L110_iEnd(iBasin)

    case (2)
       call message('***ERROR: get_basin_info_mrm has been called for level2 that does not exist within mRM')
       call message('***ERROR: use get_basin_info from mhm instead')
       stop
       
    end select

  end subroutine get_basin_info_mrm

  ! ------------------------------------------------------------------

  !      NAME
  !         calculate_grid_properties

  !     PURPOSE
  !>        \brief Calculates basic grid properties at a required coarser level using 
  !>              information of a given finer level.

  !>        \brief Calculates basic grid properties at a required coarser level (e.g., L11) using 
  !>              information of a given finer level (e.g., L0). Basic grid properties such as 
  !>              nrows, ncols, xllcorner, yllcorner cellsize are estimated in this
  !>              routine. 

  !     CALLING SEQUENCE
  !         call calculate_grid_properties( nrowsIn, ncolsIn,  xllcornerIn,                     &
  !                                         yllcornerIn,  cellsizeIn, nodata_valueIn,           &
  !                                         aimingResolution, nrowsOut, ncolsOut, xllcornerOut, &
  !                                         yllcornerOut, cellsizeOut, nodata_valueOut ) 
  !     INTENT(IN)
  !>        \param[in] "integer(i4)             :: nrowsIn"           no. of rows at an input level
  !>        \param[in] "integer(i4)             :: ncolsIn"           no. of cols at an input level
  !>        \param[in] "real(dp)                :: xllcornerIn"       xllcorner at an input level
  !>        \param[in] "real(dp)                :: yllcornerIn"       yllcorner at an input level
  !>        \param[in] "real(dp)                :: cellsizeIn"        cell size at an input level
  !>        \param[in] "real(dp)                :: nodata_valueIn"    nodata value at an input level
  !>        \param[in] "real(dp)                :: aimingResolution"  resolution of an output level

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)             :: nrowsOut"         no. of rows at an output level
  !>        \param[out] "integer(i4)             :: ncolsOut"         no. of cols at an output level
  !>        \param[out] "real(dp)                :: xllcornerOut"      xllcorner at an output level
  !>        \param[out] "real(dp)                :: yllcornerOut"      yllcorner at an output level
  !>        \param[out] "real(dp)                :: cellsizeOut"       cell size at an output level
  !>        \param[out] "real(dp)                :: nodata_valueOut"   nodata value at an output level

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note resolutions of input and output levels should confirm each other.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink & Rohini Kumar
  !>        \date Feb 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template

  subroutine calculate_grid_properties( &
       nrowsIn, ncolsIn,  xllcornerIn,  yllcornerIn,  cellsizeIn, nodata_valueIn,  &
       aimingResolution,                                                           &
       nrowsOut, ncolsOut, xllcornerOut, yllcornerOut, cellsizeOut, nodata_valueOut ) 
    use mo_kind, only: i4, dp
    use mo_message,      only: message       ! for print out
    use mo_string_utils, only: num2str

    implicit none

    integer(i4), intent(in) :: nrowsIn
    integer(i4), intent(in) :: ncolsIn
    real(dp), intent(in)    :: xllcornerIn
    real(dp), intent(in)    :: yllcornerIn
    real(dp), intent(in)    :: cellsizeIn
    real(dp), intent(in)    :: nodata_valueIn
    real(dp), intent(in)    :: aimingResolution

    integer(i4), intent(out) :: nrowsOut
    integer(i4), intent(out) :: ncolsOut
    real(dp), intent(out)    :: xllcornerOut
    real(dp), intent(out)    :: yllcornerOut
    real(dp), intent(out)    :: cellsizeOut    
    real(dp), intent(out)    :: nodata_valueOut

    ! local variables
    real(dp)                 :: cellfactor

    cellFactor = aimingResolution / cellsizeIn

    if ( nint(mod(aimingResolution, cellsizeIn)) .ne. 0) then
       call message()
       call message('***ERROR: Two resolutions size do not confirm: ',   &
            trim(adjustl(num2str(nint(AimingResolution)))),              &
            trim(adjustl(num2str(nint(cellsizeIn))))         )
       stop
    end if

    cellsizeOut     = cellsizeIn * cellFactor
    ncolsOut        = ceiling( real(ncolsIn, dp) / cellFactor)
    nrowsOut        = ceiling( real(nrowsIn, dp) / cellFactor)
    xllcornerOut    = xllcornerIn + real(ncolsIn,dp) * cellsizeIn - real(ncolsOut,dp) * cellsizeOut  
    yllcornerOut    = yllcornerIn + real(nrowsIn,dp) * cellsizeIn - real(nrowsOut,dp) * cellsizeOut  
    nodata_valueOut  = nodata_valueIn

  end subroutine calculate_grid_properties
  
end module mo_mrm_tools
