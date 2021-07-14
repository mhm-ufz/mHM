!>       \brief Grid-associated routines (File I/O)
!>       \details provides routines for Grid initialization and routines to read Grids from file and dump to file

module mo_grid
  use mo_kind, only : dp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_lowres_level, set_domain_indices, init_advanced_grid_properties, &
          mapCoordinates, geoCoordinates, calculate_grid_properties
  PUBLIC :: write_grid_info    ! write grid to (restart) file
  PUBLIC :: infer_grid_info    ! infer grid from any file
  PUBLIC :: read_grid_info     ! read grid from (restart) file
  PUBLIC :: Grid, GridRemapper     ! basic types

  integer(i4), public :: iFlag_coordinate_sys

  ! -------------------------------------------------------------------
  ! GRID description
  ! -------------------------------------------------------------------
  type Grid
    ! general domain information
    integer(i4) :: ncols     ! Number of columns
    integer(i4) :: nrows     ! Number of rows
    integer(i4) :: nCells     ! Number of rows
    real(dp) :: xllcorner    ! x coordinate of the lowerleft corner
    real(dp) :: yllcorner    ! y coordinate of the lowerleft corner
    real(dp) :: cellsize     ! Cellsize x = cellsize y
    real(dp) :: nodata_value ! Code to define the mask
    real(dp), dimension(:, :), allocatable :: x  ! 2d longitude array (unmasked version is needed for output anyway)
    real(dp), dimension(:, :), allocatable :: y  ! 2d latitude  array (unmasked version is needed for output anyway)
    logical, dimension(:, :), allocatable :: mask  ! the mask for valid cells in the original grid (nrows*ncols)
    ! for referencing values in the nValidCells vector
    integer(i4) :: iStart          ! Starting cell index of a given domain
    integer(i4) :: iEnd            ! Ending cell index of a given domain
    ! dimension(nCells, (x,y) )
    integer(i4), dimension(:, :), allocatable :: CellCoor  ! this is only used for mRM
    real(dp), dimension(:), allocatable :: CellArea  ! area of the cell in sq m
    integer(i4), dimension(:), allocatable :: Id

  end type Grid

  type GridRemapper
    type(Grid), pointer :: high_res_grid
    type(Grid), pointer :: low_res_grid

    ! dimension nCells
    integer(i4), dimension(:), allocatable :: lower_bound  ! 1d index of lower side subgrid
    integer(i4), dimension(:), allocatable :: upper_bound  ! 1d index of upper side subgrid
    integer(i4), dimension(:), allocatable :: left_bound  ! 1d index of left side subgrid
    integer(i4), dimension(:), allocatable :: right_bound  ! 1d index of right side subgrid
    integer(i4), dimension(:), allocatable :: n_subcells   ! 1d numberof valid subgrid cells
    integer(i4), dimension(:, :), allocatable :: lowres_id_on_highres   ! 2d index array of lowres id

  end type GridRemapper


contains

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
  subroutine init_lowres_level(highres, target_resolution, lowres, highres_lowres_remap)

    use mo_common_constants, only : nodata_dp, nodata_i4

    implicit none

    type(Grid), target, intent(in) :: highres  !< high resolution source (input) Grid
    real(dp), intent(in) :: target_resolution  !< resolution if target (output) Grid
    type(Grid), target, intent(inout) :: lowres  !< low resolution target (output) Grid
    type(GridRemapper), intent(inout), optional :: highres_lowres_remap  !< GridRemapper containing remapping information between two grids

    real(dp), dimension(:, :), allocatable :: areaCell0_2D

    real(dp) :: cellFactor

    integer(i4) :: iup, idown

    integer(i4) :: jl, jr

    integer(i4) :: i, j, k, ic, jc, iValue

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
      cellFactor = lowres%cellsize / highres%cellsize

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
      lowres%Id = [ (k, k = 1, lowres%nCells) ]

      allocate(lowres%x(lowres%nrows, lowres%ncols))
      allocate(lowres%y(lowres%nrows, lowres%ncols))

      ! for historic reasons this is in 2d
      lowres%x = spread([((real(iValue) - 0.5_dp) * lowres%cellsize + lowres%xllcorner, iValue=1, lowres%nrows)], &
              2, lowres%ncols)
      lowres%y = spread([((real(iValue) - 0.5_dp) * lowres%cellsize + lowres%yllcorner, iValue=1, lowres%ncols)], &
              1, lowres%nrows)

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

  !>       \brief sets indices for 1d-array of Grids based on nCells attribute
  !>       \details sets the iStart and iEnd attributes of an 1d-array of Grids based on their valid cells (nCells)
  !>       e.g. suppose those nCells of 3 Grids: 10, 100, 33
  !>       those iStart/iEnd attributes are set: 1/10, 11/110, 111/143
  !>       optionally, a list of indices can be given to provide an alternative Grid numbering
  subroutine set_domain_indices(grids, indices)

    type(Grid), intent(inout), dimension(:) :: grids  !< array of Grids
    integer(i4),   intent(in), dimension(:), optional :: indices  !< array of indices with same length as grids

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

  !>       \brief initializes advanced Grid properties
  !>       \details computes the properties
  !>       -  nCells
  !>       -  CellCoor
  !>       -  Id
  !>       -  CellArea
  !>       requires the properties:
  !>       -  mask
  !>       -  cellsize
  !>       -  ncols
  !>       -  nrows
  !>       -  (optionally yllcorner if iFlag_coordinate_sys == 1)
  subroutine init_advanced_grid_properties(new_grid)

    use mo_constants, only : RadiusEarth_dp, TWOPI_dp

    implicit none

    type(Grid), intent(inout) :: new_grid  !< Grid where new attributes are set to

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
    if(iFlag_coordinate_sys .eq. 0) then
      new_grid%CellArea(:) = new_grid%cellsize * new_grid%cellsize

      ! regular lat-lon coordinate system
    else if(iFlag_coordinate_sys .eq. 1) then

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

  end subroutine init_advanced_grid_properties


  !>       \brief Generate map coordinates
  !>       \details Generate 1D arrays of coordinates (x and y) based on cellsize, nrows and ncols properties
  subroutine mapCoordinates(level, y, x)

    type(Grid), intent(in) :: level  !< base Grid to infer grid properties from
    real(dp), intent(out), allocatable, dimension(:) :: x  !< 1d x-coordinate values
    real(dp), intent(out), allocatable, dimension(:) :: y  !< 1d y-coordinate values

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

    y(1) = level%yllcorner + 0.5_dp * cellsize
    do ii = 2, ncols
      y(ii) = y(ii - 1) + cellsize
    end do

  end subroutine mapCoordinates

  !>       \brief Generate geographic coordinates
  !>       \details Generate 2D arrays of coordinates (x and y) based on y and x properties
  subroutine geoCoordinates(level, lat, lon)

    type(Grid), intent(in) :: level  !< base Grid to infer grid properties from
    real(dp), intent(out), allocatable, dimension(:, :) :: lat  !< 2d x-coordinate values
    real(dp), intent(out), allocatable, dimension(:, :) :: lon  !< 2d y-coordinate values


    lat = level%y
    lon = level%x

  end subroutine geoCoordinates

  !>       \brief Calculates basic grid properties at a required coarser level using
  !>       information of a given finer level.
  !>       \details Calculates basic grid properties at a required coarser level (e.g., L11) using
  !>       information of a given finer level (e.g., L0). Basic grid properties such as
  !>       nrows, ncols, xllcorner, yllcorner cellsize are estimated in this
  !>       routine.
  subroutine calculate_grid_properties(nrowsIn, ncolsIn, xllcornerIn, yllcornerIn, cellsizeIn, aimingResolution, &
                                      nrowsOut, ncolsOut, xllcornerOut, yllcornerOut, cellsizeOut)

    use mo_message, only : error_message
    use mo_string_utils, only : num2str

    implicit none

    !> no. of rows at an input level
    integer(i4), intent(in) :: nrowsIn
    !> no. of cols at an input level
    integer(i4), intent(in) :: ncolsIn
    !> xllcorner at an input level
    real(dp), intent(in) :: xllcornerIn
    !> yllcorner at an input level
    real(dp), intent(in) :: yllcornerIn
    !> cell size at an input level
    real(dp), intent(in) :: cellsizeIn
    !> resolution of an output level
    real(dp), intent(in) :: aimingResolution
    !> no. of rows at an output level
    integer(i4), intent(out) :: nrowsOut
    !> no. of cols at an output level
    integer(i4), intent(out) :: ncolsOut
    !> xllcorner at an output level
    real(dp), intent(out) :: xllcornerOut
    !> yllcorner at an output level
    real(dp), intent(out) :: yllcornerOut
    !> cell size at an output level
    real(dp), intent(out) :: cellsizeOut

    real(dp) :: cellFactor, rounded

    cellFactor = aimingResolution / cellsizeIn
    rounded = anint(cellFactor)

    ! TODO: sync this with other comparisons
    if (abs(rounded - cellFactor) > 1.e-6_dp) then
      call error_message('***ERROR: Two resolutions size do not confirm: ', &
              num2str(aimingResolution), ' and ', &
              num2str(cellsizeIn))
    end if

    cellsizeOut = cellsizeIn * rounded
    ncolsOut = nint(real(ncolsIn, dp) / cellFactor)
    nrowsOut = nint(real(nrowsIn, dp) / cellFactor)
    xllcornerOut = xllcornerIn + real(ncolsIn, dp) * cellsizeIn - real(ncolsOut, dp) * cellsizeOut
    yllcornerOut = yllcornerIn + real(nrowsIn, dp) * cellsizeIn - real(nrowsOut, dp) * cellsizeOut

  end subroutine calculate_grid_properties

  !>       \brief write restart file
  !>       \details write restart file for given Grid to given NcDataset
  subroutine write_grid_info(grid_in, level_name, nc)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    !> Grid to be written
    type(Grid), intent(in) :: grid_in
    !> name of level (used to label the grid, e.g. "0", "1", "11", "2")
    character(*), intent(in) :: level_name
    !> NcDataset to write information to
    type(NcDataset), intent(inout) :: nc

    ! dummy for gnu73
    integer(i4), allocatable :: dummy(:, :)
    type(NcDimension) :: rows, cols
    type(NcVariable) :: var


    rows = nc%setDimension("nrows" // trim(level_name), grid_in%nrows)
    cols = nc%setDimension("ncols" // trim(level_name), grid_in%ncols)

    ! now set everything related to the grid
    var = nc%setVariable("L" // trim(level_name) // "_domain_mask", "i32", (/rows, cols/))
    call var%setFillValue(nodata_i4)
    ! transform from logical to i32
    ! ST: where statement is used because gnu73 does not properly translate with merge
    allocate(dummy(size(grid_in%mask, 1), size(grid_in%mask, 2)))
    dummy = 0_i4
    where(grid_in%mask) dummy = 1_i4
    call var%setData(dummy)
    deallocate(dummy)
    call var%setAttribute("long_name", "Mask at level " // trim(level_name))

    var = nc%setVariable("L" // trim(level_name) // "_domain_lat", "f64", (/rows, cols/))
    call var%setFillValue(nodata_dp)
    call var%setData(grid_in%y)
    call var%setAttribute("long_name", "Latitude at level " // trim(level_name))

    var = nc%setVariable("L" // trim(level_name) // "_domain_lon", "f64", (/rows, cols/))
    call var%setFillValue(nodata_dp)
    call var%setData(grid_in%x)
    call var%setAttribute("long_name", "Longitude at level " // trim(level_name))

    var = nc%setVariable("L" // trim(level_name) // "_domain_cellarea", "f64", (/rows, cols/))
    call var%setFillValue(nodata_dp)
    call var%setData(unpack(grid_in%CellArea * 1.0E-6_dp, grid_in%mask, nodata_dp))
    call var%setAttribute("long_name", "Cell area at level " // trim(level_name))

    call nc%setAttribute("xllcorner_L" // trim(level_name), grid_in%xllcorner)
    call nc%setAttribute("yllcorner_L" // trim(level_name), grid_in%yllcorner)
    call nc%setAttribute("cellsize_L" // trim(level_name), grid_in%cellsize)
    call nc%setAttribute("nrows_L" // trim(level_name), grid_in%nrows)
    call nc%setAttribute("ncols_L" // trim(level_name), grid_in%ncols)
    call nc%setAttribute("nCells_L" // trim(level_name), grid_in%nCells)

  end subroutine write_grid_info


  !>       \brief reads complete Grid properties from NetCDF file
  !>       \details reads complete Grid properties from NetCDF file
  subroutine read_grid_info(domainID, inputFile, level_name, new_grid)

    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str
    use mo_os, only: path_isfile

    implicit none

    !> number of domainID
    integer(i4), intent(in) :: domainID
    !> complete path to file
    character(256), intent(in) :: inputFile
    !> name of level (used to label the grid, e.g. "0", "1", "11", "2")
    character(*), intent(in) :: level_name
    !> grid to save information to
    type(Grid), intent(inout) :: new_grid

    ! dummy, 2 dimension I4
    integer(i4), dimension(:, :), allocatable :: dummyI2

    ! dummy, 2 dimension DP
    real(dp), dimension(:, :), allocatable :: dummyD2

    character(256) :: Fname

    type(NcDataset) :: nc

    type(NcVariable) :: var

    integer(i4) :: k


    ! read config
    fname = trim(inputFile)
    call message('    Reading config from     ', trim(adjustl(Fname)), ' ...')

    call path_isfile(path = fname, quiet_ = .true., throwError_ = .true.)
    nc = NcDataset(fname, "r")

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read the grid properties
    call nc%getAttribute("xllcorner_L" // trim(level_name), new_grid%xllcorner)
    call nc%getAttribute("yllcorner_L" // trim(level_name), new_grid%yllcorner)
    call nc%getAttribute("nrows_L" // trim(level_name), new_grid%nrows)
    call nc%getAttribute("ncols_L" // trim(level_name), new_grid%ncols)
    call nc%getAttribute("cellsize_L" // trim(level_name), new_grid%cellsize)
    call nc%getAttribute("nCells_L" // trim(level_name), new_grid%nCells)

    if (.not. allocated(new_grid%mask)) allocate(new_grid%mask(new_grid%nrows, new_grid%ncols))
    if (.not. allocated(new_grid%x)) allocate(new_grid%x(new_grid%nrows, new_grid%ncols))
    if (.not. allocated(new_grid%y)) allocate(new_grid%y(new_grid%nrows, new_grid%ncols))
    ! read L1 mask
    var = nc%getVariable("L" // trim(level_name) // "_domain_mask")
    ! read integer
    call var%getData(dummyI2)
    ! transform to logical
    new_grid%mask = (dummyI2 .eq. 1_i4)

    var = nc%getVariable("L" // trim(level_name) // "_domain_lat")
    call var%getData(new_grid%y)

    var = nc%getVariable("L" // trim(level_name) // "_domain_lon")
    call var%getData(new_grid%x)

    var = nc%getVariable("L" // trim(level_name) // "_domain_cellarea")
    call var%getData(dummyD2)
    if (.not. allocated(new_grid%CellArea)) new_grid%CellArea = pack(dummyD2 / 1.0E-6_dp, new_grid%mask)
    ! new_grid%CellArea = pack(dummyD2 / 1.0E-6_dp, new_grid%mask)

    call nc%close()

    new_grid%Id = (/ (k, k = 1, new_grid%nCells) /)

  end subroutine read_grid_info

  !>       \brief infers Grid properties implicitly from file
  !>       \details infers Grid properties implicitly from file including:
  !>       - xllcorner, yllcorner, nRows, nCols, cellsize, nCells, mask, x, y, cellarea
  !>       this is attempted by scanning for known x and y coordinate names and inferred by mask of given variable
  subroutine infer_grid_info(inputFile, xCoordName, yCoordName, maskVar, new_grid)

    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str
    use mo_os, only: path_isfile

    implicit none

    !> complete path to file
    character(*), intent(in) :: inputFile
    !> name of variable to be used for inferring x-coordinate
    character(*), intent(in) :: xCoordName
    !> name of variable to be used for inferring y-coordinate
    character(*), intent(in) :: yCoordName
    !> name of variable to be used for inferring mask
    character(*), intent(in) :: maskVar
    !> grid to save information to
    type(Grid), intent(inout) :: new_grid

    type(NcDataset) :: nc

    ! read config
    call message('    Inferring grid from     ', trim(inputFile), ' ...')

    call path_isfile(path = inputFile, quiet_ = .true., throwError_ = .true.)
    nc = NcDataset(inputFile, "r")

    call get_coordinate(nc, xCoordName, new_grid%xllcorner, new_grid%nRows, new_grid%cellsize)
    call get_coordinate(nc, yCoordName, new_grid%yllcorner, new_grid%nCols, new_grid%cellsize)
    ! calculate nCells, mask, x, y, cellarea
    call infer_advanced_grid_properties(nc, xCoordName, yCoordName, maskVar, new_grid)

    call nc%close()

  end subroutine infer_grid_info

  !>       \brief retrieves lower bound, nCells and cellsize information of 1d NetCDF coordinate variable
  !>       \details properties used are needed to set Grid properties xllcorner/yllcorner, nRows/nCols and cellsize
  subroutine get_coordinate(nc, coordName, lowerBound, n, cellsize)

    use mo_netcdf, only : NcDataset, NcVariable
    use mo_message, only: error_message

    type(NcDataset), intent(inout) :: nc  !< NcDataset
    character(*), intent(in) :: coordName  !< name of 1d coordinate variable in NcDataset
    real(dp), intent(out) :: lowerBound  !< lower bound of coordinate variable
    integer(i4), intent(out) :: n  !< number of values in coordinate variable
    real(dp), intent(out) :: cellsize  !< stepsize of coordinate values
    
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: varShape
    real(dp), dimension(:), allocatable :: tempValues
    character(256) :: boundVariable
    real(dp), dimension(:, :), allocatable :: dummy

    if (nc%hasVariable(trim(coordName))) then
      ! the coordinate is actually a variable in the NcDataset
      ncVar = nc%getVariable(trim(coordName))
      varShape = ncVar%getShape()
      ! infer the dimensions...
      if (size(varShape) /= 1_i4) then
        call error_message("cannot infer grid from non 1d-coordinate ", trim(coordName))
      end if
      ! the variable is dependent on 1 dimension only
      n = varShape(1)
      allocate(tempValues(n))
      call ncVar%getData(tempValues)
      ! infer directly from bounds
      if (ncVar%hasAttribute('bounds')) then
        call ncVar%getAttribute('bounds', boundVariable)
        ncVar = nc%getVariable(trim(boundVariable))
        call ncVar%getData(dummy)
        lowerBound = dummy(1, 1)
        cellsize = abs((dummy(2, size(dummy, 2)) - lowerBound) / real(n, dp))
        deallocate(dummy)
      elseif (n > 1_i4) then
        cellsize = abs(tempValues(2) - tempValues(1))
        lowerBound = tempValues(1) - 0.5_dp * cellsize
      else
        call error_message("cannot infer cellsize of coordinate ", trim(coordName))
      end if
    else
      call error_message("cannot infer Grid properties by non existing coordinate name ", trim(coordName))
    end if

  end subroutine get_coordinate

  !>       \brief implictily retrieves advanced properties from NetCDF dataset
  !>       \details infers additional properties from nc dataset to set remaining properties of Grid
  subroutine infer_advanced_grid_properties(nc, xCoordName, yCoordName, maskVar, new_grid)

    use mo_netcdf, only : NcDataset, NcVariable
    use mo_message, only: message

    type(NcDataset), intent(inout) :: nc  !< NetCDF dataset to infer properties from
    character(*), intent(in) :: xCoordName  !< 1d coordinate variable name to set x
    character(*), intent(in) :: yCoordName  !< 1d coordinate variable name to set y
    character(*), intent(in) :: maskVar  !< variable name whose mask is used to set mask
    type(Grid), intent(inout) :: new_grid  !< grid to save information to

    type(NcVariable) :: ncVar
    real(dp), dimension(:, :, :), allocatable :: dummyD3
    logical, dimension(:, :, :), allocatable :: maskD3
    real(dp), dimension(:, :), allocatable :: dummyD2
    real(dp), dimension(:), allocatable :: dummyD1
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: i

    if (.not. allocated(new_grid%mask)) allocate(new_grid%mask(new_grid%nrows, new_grid%ncols))
    if (.not. allocated(new_grid%x)) allocate(new_grid%x(new_grid%nrows, new_grid%ncols))
    if (.not. allocated(new_grid%y)) allocate(new_grid%y(new_grid%nrows, new_grid%ncols))

    ! read mask
    ncVar = nc%getVariable(maskVar)
    ncVarShape = ncVar%getshape()

    if (size(ncVarShape) == 3_i4) then
      call ncVar%getData(dummyD3, mask=maskD3)
      new_grid%mask = maskD3(:,:,1)
      deallocate(dummyD3, maskD3)
    elseif (size(ncVarShape) == 2_i4) then
      ! read data
      call ncVar%getData(dummyD2, mask=new_grid%mask)
      deallocate(dummyD2)
    else
      print*, 'Expected 2D or 3D field for inferring mask of grid'
      stop 1
    end if

    ncVar = nc%getVariable(xCoordName)
    call ncVar%getData(dummyD1)
    new_grid%x = spread(dummyD1, 2, new_grid%ncols)

    ncVar = nc%getVariable(yCoordName)
    call ncVar%getData(dummyD1)
    new_grid%y = spread(dummyD1, 1, new_grid%nrows)
    deallocate(dummyD1)

    ! init the remaining properties
    call init_advanced_grid_properties(new_grid)

  end subroutine infer_advanced_grid_properties
end module mo_grid

