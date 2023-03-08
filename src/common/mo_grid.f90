!> \file mo_grid.f90
!> \brief \copybrief mo_grid
!> \details \copydetails mo_grid

!> \brief gridding tools
!> \details Common tools to deal with grids in mHM.
!> \authors Robert Schweppe
!> \date Jun 2018
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_grid
  use mo_kind, only : dp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_lowres_level, set_domain_indices, L0_grid_setup, &
          mapCoordinates, geoCoordinates
contains
  ! ------------------------------------------------------------------

  !    NAME
  !        init_lowres_level

  !    PURPOSE
  !>       \brief Level-1 variable initialization

  !>       \details following tasks are performed for L1 datasets
  !>       -  cell id & numbering
  !>       -  mask creation
  !>       -  storage of cell cordinates (row and coloum id)
  !>       -  sorage of four corner L0 cordinates
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine config_variables_set in
  !>       module mo_restart and in the subroutine set_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "type(Grid) :: highres"
  !>       \param[in] "real(dp) :: target_resolution"

  !    INTENT(INOUT)
  !>       \param[inout] "type(Grid) :: lowres"

  !    INTENT(INOUT), OPTIONAL
  !>       \param[inout] "type(GridRemapper), optional :: highres_lowres_remap"

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Jan 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine init_lowres_level(highres, target_resolution, lowres, highres_lowres_remap)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_types, only : Grid, GridRemapper

    implicit none

    type(Grid), target, intent(in) :: highres

    real(dp), intent(in) :: target_resolution

    type(Grid), target, intent(inout) :: lowres

    type(GridRemapper), intent(inout), optional :: highres_lowres_remap

    real(dp), dimension(:, :), allocatable :: areaCell0_2D

    real(dp) :: cellFactor

    integer(i4) :: iup, idown

    integer(i4) :: jl, jr

    integer(i4) :: i, j, k, ic, jc

    ! STEPS ::

    !--------------------------------------------------------
    ! 1) Estimate each variable locally for a given domain
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! grid properties
    if (.not. allocated(lowres%mask)) then
      call calculate_grid_properties(highres%nrows, highres%ncols, &
              highres%xllcorner, highres%yllcorner, highres%cellsize, &
              target_resolution, &
              lowres%nrows, lowres%ncols, &
              lowres%xllcorner, lowres%yllcorner, lowres%cellsize)
      ! cellfactor = leve1-1 / level-0
      cellFactor = anint(lowres%cellsize / highres%cellsize)

      ! allocation and initalization of mask at level-1
      allocate(lowres%mask(lowres%nrows, lowres%ncols))
      lowres%mask(:, :) = .FALSE.

      ! create mask at level-1
      do j = 1, highres%ncols
        jc = ceiling(real(j, dp) / cellFactor)
        do i = 1, highres%nrows
          if (.NOT. highres%mask(i, j)) cycle
          ic = ceiling(real(i, dp) / cellFactor)
          lowres%mask(ic, jc) = .TRUE.
        end do
      end do

      ! estimate ncells and initalize related variables
      lowres%nCells = count(lowres%mask)
      ! allocate and initalize cell1 related variables
      allocate(lowres%Id        (lowres%nCells))
      lowres%Id = (/ (k, k = 1, lowres%nCells) /)
    end if

    if (present(highres_lowres_remap)) then
      ! cellfactor = leve1-1 / level-0, set again in case not yet initialized
      cellFactor = lowres%cellsize / highres%cellsize

      ! lowres additional properties
      allocate(areaCell0_2D(highres%nrows, highres%ncols))
      areaCell0_2D(:, :) = UNPACK(highres%CellArea, highres%mask, nodata_dp)

      if (.not. allocated(lowres%CellCoor)) then
        allocate(lowres%CellCoor  (lowres%nCells, 2))
        allocate(lowres%CellArea  (lowres%nCells))
      end if

      allocate(highres_lowres_remap%lower_bound(lowres%nCells))
      allocate(highres_lowres_remap%upper_bound(lowres%nCells))
      allocate(highres_lowres_remap%left_bound (lowres%nCells))
      allocate(highres_lowres_remap%right_bound(lowres%nCells))
      allocate(highres_lowres_remap%n_subcells (lowres%nCells))
      allocate(highres_lowres_remap%lowres_id_on_highres (highres%nrows, highres%ncols))
      highres_lowres_remap%lowres_id_on_highres = nodata_i4

      highres_lowres_remap%high_res_grid => highres
      highres_lowres_remap%low_res_grid => lowres

      k = 0
      do jc = 1, lowres%ncols
        do ic = 1, lowres%nrows
          if (.NOT. lowres%mask(ic, jc)) cycle
          k = k + 1

          lowres%CellCoor(k, 1) = ic
          lowres%CellCoor(k, 2) = jc

          ! coord. of all corners -> of finer scale level-0
          iup = (ic - 1) * nint(cellFactor, i4) + 1
          idown = ic * nint(cellFactor, i4)
          jl = (jc - 1) * nint(cellFactor, i4) + 1
          jr = jc * nint(cellFactor, i4)

          ! constrain the range of up, down, left, and right boundaries
          if(iup   < 1) iup = 1
          if(idown > highres%nrows) idown = highres%nrows
          if(jl    < 1) jl = 1
          if(jr    > highres%ncols) jr = highres%ncols

          highres_lowres_remap%upper_bound   (k) = iup
          highres_lowres_remap%lower_bound (k) = idown
          highres_lowres_remap%left_bound (k) = jl
          highres_lowres_remap%right_bound(k) = jr

          ! effective area [km2] & total no. of L0 cells within a given L1 cell
          lowres%CellArea(k) = sum(areacell0_2D(iup : idown, jl : jr), highres%mask(iup : idown, jl : jr))
          highres_lowres_remap%n_subcells(k) = count(highres%mask(iup : idown, jl : jr))
          ! Delimitation of level-11 cells on level-0
          highres_lowres_remap%lowres_id_on_highres(iup : idown, jl : jr) = k
        end do
      end do

      ! free space
      deallocate(areaCell0_2D)

    end if

  end subroutine init_lowres_level

  !    NAME
  !        set_domain_indices

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(INOUT)
  !>       \param[inout] "type(Grid), dimension(:) :: grids"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  !        Stephan Thober, Aug 2019 - added optional indices for L0 data because L0 data can be shared among domains

  subroutine set_domain_indices(grids, indices)

    use mo_common_types, only: Grid

    implicit none

    type(Grid), intent(inout), dimension(:) :: grids
    integer(i4),   intent(in), dimension(:), optional :: indices

    integer(i4) :: iDomain


    do iDomain = 1, size(grids)
      ! Saving indices of mask and packed data
      if(iDomain .eq. 1_i4) then
        grids(iDomain)%iStart = 1_i4
      else
        if (present(indices)) then
          grids(iDomain)%iStart = grids(indices(iDomain - 1))%iEnd + 1_i4
        else
          grids(iDomain)%iStart = grids(iDomain - 1_i4)%iEnd + 1_i4
        end if
      end if
      grids(iDomain)%iEnd = grids(iDomain)%iStart + grids(iDomain)%nCells - 1_i4
    end do

  end subroutine set_domain_indices

  ! ------------------------------------------------------------------

  !    NAME
  !        L0_grid_setup

  !    PURPOSE
  !>       \brief level 0 variable initialization

  !>       \details following tasks are performed for L0 data sets
  !>       -  cell id & numbering
  !>       -  storage of cell cordinates (row and coloum id)
  !>       -  empirical dist. of terrain slope
  !>       -  flag to determine the presence of a particular soil id
  !>       in this configuration of the model run
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine config_variables_set in
  !>       module mo_restart and in the subroutine set_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(INOUT)
  !>       \param[inout] "type(Grid) :: new_grid"

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Jan 2013

  ! Modifications:
  ! Rohini Kumar & Matthias Cuntz  May 2014 - cell area calulation based on a regular lat-lon grid or
  !                                           on a regular X-Y coordinate system
  ! Matthias Cuntz                 May 2014 - changed empirical distribution function so that doubles get the same value
  ! Matthias Zink & Matthias Cuntz Feb 2016 - code speed up due to reformulation of CDF calculation
  ! Rohini Kumar                   Mar 2016 - changes for handling multiple soil database options
  ! Robert Schweppe                Jun 2018 - refactoring and reformatting

  subroutine L0_grid_setup(new_grid)

    use mo_common_types, only: Grid
    use mo_common_variables, only : iFlag_cordinate_sys
    use mo_constants, only : RadiusEarth_dp, TWOPI_dp

    implicit none

    type(Grid), intent(inout) :: new_grid

    real(dp), dimension(:, :), allocatable :: areaCell_2D

    integer(i4) :: i, j, k

    real(dp) :: rdum, degree_to_radian, degree_to_metre

    ! STEPS ::


    !--------------------------------------------------------
    ! 1) Estimate each variable locally for a given domain
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! level-0 information
    new_grid%nCells = count(new_grid%mask)

    allocate(new_grid%CellCoor(new_grid%nCells, 2))
    allocate(new_grid%Id(new_grid%nCells))
    allocate(new_grid%CellArea(new_grid%nCells))
    allocate(areaCell_2D(new_grid%nrows, new_grid%ncols))

    new_grid%Id = (/ (k, k = 1, new_grid%nCells) /)

    !------------------------------------------------
    ! start looping for cell cordinates and ids
    !------------------------------------------------
    k = 0
    do j = 1, new_grid%ncols
      do i = 1, new_grid%nrows
        if (.NOT. new_grid%mask(i, j)) cycle
        k = k + 1
        new_grid%cellCoor(k, 1) = i
        new_grid%cellCoor(k, 2) = j
      end do
    end do

    ! ESTIMATE AREA [m2]

    ! regular X-Y coordinate system
    if(iFlag_cordinate_sys .eq. 0) then
      new_grid%CellArea(:) = new_grid%cellsize * new_grid%cellsize

      ! regular lat-lon coordinate system
    else if(iFlag_cordinate_sys .eq. 1) then

      degree_to_radian = TWOPI_dp / 360.0_dp
      degree_to_metre = RadiusEarth_dp * TWOPI_dp / 360.0_dp
      do i = new_grid%ncols, 1, -1
        j = new_grid%ncols - i + 1
        ! get latitude in degrees
        rdum = new_grid%yllcorner + (real(j, dp) - 0.5_dp) * new_grid%cellsize
        ! convert to radians
        rdum = rdum * degree_to_radian
        !    AREA [m2]
        areaCell_2D(:, i) = (new_grid%cellsize * cos(rdum) * degree_to_metre) * (new_grid%cellsize * degree_to_metre)
      end do
      new_grid%CellArea(:) = pack(areaCell_2D(:, :), new_grid%mask)

    end if

    ! free space
    deallocate(areaCell_2D)

  end subroutine L0_grid_setup


  !------------------------------------------------------------------
  !    NAME
  !        mapCoordinates

  !    PURPOSE
  !>       \brief Generate map coordinates

  !>       \details Generate map coordinate arrays for given domain and level

  !    INTENT(IN)
  !>       \param[in] "type(Grid) :: level" -> grid reference

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: x, y"
  !>       \param[out] "real(dp), dimension(:) :: x, y"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Stephan Thober Nov 2013 - removed fproj dependency
  ! David Schaefer Jun 2015 - refactored the former subroutine CoordSystem
  ! Robert Schweppe Jun 2018 - refactoring and reformatting


  subroutine mapCoordinates(level, y, x)

    use mo_common_types, only: Grid

    implicit none

    ! -> grid reference
    type(Grid), intent(in) :: level

    real(dp), intent(out), allocatable, dimension(:) :: x, y

    integer(i4) :: ii, ncols, nrows

    real(dp) :: cellsize


    cellsize = level%cellsize
    nrows = level%nrows
    ncols = level%ncols

    allocate(x(nrows), y(ncols))

    x(1) = level%xllcorner + 0.5_dp * cellsize
    do ii = 2, nrows
      x(ii) = x(ii - 1) + cellsize
    end do

    ! inverse for Panoply, ncview display
    y(ncols) = level%yllcorner + 0.5_dp * cellsize
    do ii = ncols - 1, 1, -1
      y(ii) = y(ii + 1) + cellsize
    end do

  end subroutine mapCoordinates

  !------------------------------------------------------------------
  !    NAME
  !        geoCoordinates

  !    PURPOSE
  !>       \brief Generate geographic coordinates

  !>       \details Generate geographic coordinate arrays for given domain and level

  !    INTENT(IN)
  !>       \param[in] "type(Grid) :: level" -> grid reference

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: lat, lon"
  !>       \param[out] "real(dp), dimension(:, :) :: lat, lon"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Stephan Thober  Nov 2013 - removed fproj dependency
  ! David Schaefer  Jun 2015 - refactored the former subroutine CoordSystem
  ! Stephan Thober  Sep 2015 - using mask to unpack coordinates
  ! Stephan Thober  Oct 2015 - writing full lat/lon again
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine geoCoordinates(level, lat, lon)

    use mo_common_types, only: Grid

    implicit none

    ! -> grid reference
    type(Grid), intent(in) :: level

    real(dp), intent(out), allocatable, dimension(:, :) :: lat, lon


    lat = level%y
    lon = level%x

  end subroutine geoCoordinates

  ! ------------------------------------------------------------------

  !    NAME
  !        calculate_grid_properties

  !    PURPOSE
  !>       \brief Calculates basic grid properties at a required coarser level using
  !>       information of a given finer level.
  !>       Calculates basic grid properties at a required coarser level (e.g., L11) using
  !>       information of a given finer level (e.g., L0). Basic grid properties such as
  !>       nrows, ncols, xllcorner, yllcorner cellsize are estimated in this
  !>       routine.

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: nrowsIn"       no. of rows at an input level
  !>       \param[in] "integer(i4) :: ncolsIn"       no. of cols at an input level
  !>       \param[in] "real(dp) :: xllcornerIn"      xllcorner at an input level
  !>       \param[in] "real(dp) :: yllcornerIn"      yllcorner at an input level
  !>       \param[in] "real(dp) :: cellsizeIn"       cell size at an input level
  !>       \param[in] "real(dp) :: aimingResolution" resolution of an output level

  !    INTENT(OUT)
  !>       \param[out] "integer(i4) :: nrowsOut"  no. of rows at an output level
  !>       \param[out] "integer(i4) :: ncolsOut"  no. of cols at an output level
  !>       \param[out] "real(dp) :: xllcornerOut" xllcorner at an output level
  !>       \param[out] "real(dp) :: yllcornerOut" yllcorner at an output level
  !>       \param[out] "real(dp) :: cellsizeOut"  cell size at an output level

  !    HISTORY
  !>       \authors Matthias Zink & Rohini Kumar

  !>       \date Feb 2013

  ! Modifications:
  ! R. Kumar        Sep 2013 - documentation added according to the template
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine calculate_grid_properties(nrowsIn, ncolsIn, xllcornerIn, yllcornerIn, cellsizeIn, aimingResolution, &
                                      nrowsOut, ncolsOut, xllcornerOut, yllcornerOut, cellsizeOut)

    use mo_message, only : error_message
    use mo_string_utils, only : num2str

    implicit none

    ! no. of rows at an input level
    integer(i4), intent(in) :: nrowsIn

    ! no. of cols at an input level
    integer(i4), intent(in) :: ncolsIn

    ! xllcorner at an input level
    real(dp), intent(in) :: xllcornerIn

    ! yllcorner at an input level
    real(dp), intent(in) :: yllcornerIn

    ! cell size at an input level
    real(dp), intent(in) :: cellsizeIn

    ! resolution of an output level
    real(dp), intent(in) :: aimingResolution

    ! no. of rows at an output level
    integer(i4), intent(out) :: nrowsOut

    ! no. of cols at an output level
    integer(i4), intent(out) :: ncolsOut

    ! xllcorner at an output level
    real(dp), intent(out) :: xllcornerOut

    ! yllcorner at an output level
    real(dp), intent(out) :: yllcornerOut

    ! cell size at an output level
    real(dp), intent(out) :: cellsizeOut

    real(dp) :: cellFactor, rounded
    integer(i4) :: rounded_int


    cellFactor = aimingResolution / cellsizeIn
    rounded = anint(cellFactor)
    rounded_int = nint(cellFactor)

    if (abs(rounded - cellFactor) > 1.e-7_dp) then
      call error_message( &
        '***ERROR: Two resolutions size do not confirm: ', &
        trim(adjustl(num2str(nint(AimingResolution)))), &
        trim(adjustl(num2str(nint(cellsizeIn)))))
    end if

    cellsizeOut = aimingResolution
    ncolsOut = nint(real(ncolsIn, dp) / cellFactor)
    nrowsOut = nint(real(nrowsIn, dp) / cellFactor)

    ! if we rounded down, but now we would miss cells, add rows and/or cols
    if ( ncolsOut * rounded_int < ncolsIn ) ncolsOut = ncolsOut + 1_i4
    if ( nrowsOut * rounded_int < nrowsIn ) nrowsOut = nrowsOut + 1_i4

    xllcornerOut = xllcornerIn + real(ncolsIn, dp) * aimingResolution / rounded - real(ncolsOut, dp) * cellsizeOut
    yllcornerOut = yllcornerIn + real(nrowsIn, dp) * aimingResolution / rounded - real(nrowsOut, dp) * cellsizeOut

  end subroutine calculate_grid_properties

end module mo_grid
