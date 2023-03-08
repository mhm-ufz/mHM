!> \file mo_common_restart.f90
!> \brief   \copybrief mo_common_read_data
!> \details \copydetails mo_common_read_data

!> \brief common restart tools
!> \details Routines to deal with grid infos for restart files
!> \authors Robert Schweppe
!> \date Jun 2018
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_restart

  use mo_kind, only : i4, dp
  use mo_message, only: message, error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_grid_info
  PUBLIC :: read_grid_info     ! read restart files for configuration from a given path
  PUBLIC :: read_nLAI_and_check_dims


  !> \brief check consistency of two given items
  INTERFACE check_consistency_element
    MODULE PROCEDURE check_consistency_element_i4, check_consistency_element_dp
  end interface check_consistency_element


CONTAINS


  !> \brief write restart files for each domain
  !> \details write restart files for each domain. For each domain
  !! three restart files are written. These are xxx_states.nc,
  !! xxx_L11_config.nc, and xxx_config.nc (xxx being the three digit
  !! domain index). If a variable is added here, it should also be added
  !! in the read restart routines below.
  !> \changelog
  !! - Stephan Thober     Aug  2015
  !!   - moved write of routing states to mRM
  !! - David Schaefer     Nov  2015
  !!   - mo_netcdf
  !! - Stephan Thober     Nov  2016
  !!   - moved processMatrix to common variables
  !! - Zink M. Demirel C. Mar 2017
  !!   - Added Jarvis soil water stress function at SM process(3)
  !! - Robert Schweppe    Feb 2018
  !!   - Removed all L0 references
  !! - Robert Schweppe    Jun 2018
  !!   - refactoring and reformatting
  !! - Stephan Thober     May 2019
  !!   - where statement for gnu73 to translate level0 mask
  !> \authors Stephan Thober
  !> \date Jun 2014
  subroutine write_grid_info(grid_in, level_name, nc)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_types, only: Grid
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    !> level to be written
    type(Grid), intent(in) :: grid_in

    !> level_id
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


  !> \brief reads configuration apart from Level 11 configuration from a restart directory
  !> \details read configuration variables from a given restart
  !> directory and initializes all configuration variables,
  !> that are initialized in the subroutine initialise,
  !> contained in module mo_startup.
  !> \changelog
  !! - David Schaefer     Nov 2015
  !!   - mo_netcdf
  !! - Zink M. Demirel C. Mar 2017
  !!   - Added Jarvis soil water stress function at SM process(3)
  !! - Robert Schweppe    Feb 2018
  !!   - Removed all L0 references
  !! - Robert Schweppe    Jun 2018
  !!   - refactoring and reformatting
  !! - Stephan Thober     May 2019
  !!   - added allocation check for mask and cellArea because cellArea needs to be read by mRM, but mask is created before by mHM
  !> \authors Stephan Thober
  !> \date Apr 2013
  subroutine read_grid_info(InFile, level_name, new_grid)

    use mo_common_types, only: Grid
    use mo_netcdf, only : NcDataset, NcVariable

    implicit none

    !> Input Path including trailing slash
    character(256), intent(in) :: InFile

    !> level_name (id)
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
    fname = trim(InFile)
    call message('    Reading config from     ', trim(adjustl(Fname)), ' ...')

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


  !> \brief read nubmer of LAI time steps and check dimension configurations read from restart file
  !> \author Robert Schweppe
  !> \date Aug 2019
  !> \author Sebastian Mueller
  !> \date Feb 2023
  subroutine read_nLAI_and_check_dims(iDomain, InFile)

    use mo_mpr_global_variables, only: nLAI, LAIBoundaries ! may read from restart
    use mo_netcdf, only : NcDataset, NcVariable, NcDimension
    use mo_common_constants, only : soilHorizonsVarName, landCoverPeriodsVarName, LAIVarName
    use mo_common_mHM_mRM_variables, only: read_old_style_restart_bounds

    implicit none

    !> domain counter (not ID)
    integer(i4), intent(in) :: iDomain
    !> Input Path including trailing slash
    character(256), intent(in) :: InFile

    character(256) :: fname
    type(NcDataset) :: nc
    type(NcVariable) :: var
    type(NcDimension) :: nc_dim

    integer(i4) :: nSoilHorizons_temp, nLAIs_temp, nLandCoverPeriods_temp
    real(dp), dimension(:), allocatable :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
            LAIBoundaries_temp

    ! dummy, 2 dimension
    real(dp), dimension(:, :), allocatable :: dummyD2, dummyD2_tmp

    integer(i4) :: ii

    ! read config
    fname = trim(InFile)
    call message('    Reading and checking LAI, land-cover and soil-horizons from     ', trim(adjustl(Fname)), ' ...')

    nc = NcDataset(fname, "r")

    ! get the dimensions
    var = nc%getVariable(trim(soilHorizonsVarName)//'_bnds')
    call var%getData(dummyD2_tmp)
    if (allocated(dummyD2)) deallocate(dummyD2)
    if ( read_old_style_restart_bounds ) then
      allocate(dummyD2(size(dummyD2_tmp,2), size(dummyD2_tmp,1)))
      dummyD2 = transpose(dummyD2_tmp)
    else
      allocate(dummyD2(size(dummyD2_tmp,1), size(dummyD2_tmp,2)))
      dummyD2 = dummyD2_tmp
    end if
    deallocate(dummyD2_tmp)
    nSoilHorizons_temp = size(dummyD2, 2)
    allocate(soilHorizonBoundaries_temp(nSoilHorizons_temp+1))
    soilHorizonBoundaries_temp(1:nSoilHorizons_temp) = dummyD2(1, :)
    soilHorizonBoundaries_temp(nSoilHorizons_temp+1) = dummyD2(2, nSoilHorizons_temp)

    ! get the landcover dimension
    var = nc%getVariable(trim(landCoverPeriodsVarName)//'_bnds')
    call var%getData(dummyD2_tmp)
    if (allocated(dummyD2)) deallocate(dummyD2)
    if ( read_old_style_restart_bounds ) then
      allocate(dummyD2(size(dummyD2_tmp,2), size(dummyD2_tmp,1)))
      dummyD2 = transpose(dummyD2_tmp)
    else
      allocate(dummyD2(size(dummyD2_tmp,1), size(dummyD2_tmp,2)))
      dummyD2 = dummyD2_tmp
    end if
    deallocate(dummyD2_tmp)
    nLandCoverPeriods_temp = size(dummyD2, 2)
    allocate(landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1))
    landCoverPeriodBoundaries_temp(1:nLandCoverPeriods_temp) = dummyD2(1, :)
    landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1) = dummyD2(2, nLandCoverPeriods_temp)

    ! get the LAI dimension
    if (nc%hasVariable(trim(LAIVarName)//'_bnds')) then
      var = nc%getVariable(trim(LAIVarName)//'_bnds')
      call var%getData(dummyD2_tmp)
      if (allocated(dummyD2)) deallocate(dummyD2)
      if ( read_old_style_restart_bounds ) then
        allocate(dummyD2(size(dummyD2_tmp,2), size(dummyD2_tmp,1)))
        dummyD2 = transpose(dummyD2_tmp)
      else
        allocate(dummyD2(size(dummyD2_tmp,1), size(dummyD2_tmp,2)))
        dummyD2 = dummyD2_tmp
      end if
      deallocate(dummyD2_tmp)
      nLAIs_temp = size(dummyD2, 2)
      allocate(LAIBoundaries_temp(nLAIs_temp+1))
      LAIBoundaries_temp(1:nLAIs_temp) = dummyD2(1, :)
      LAIBoundaries_temp(nLAIs_temp+1) = dummyD2(2, nLAIs_temp)
    else if (nc%hasDimension('L1_LAITimesteps')) then
      nc_dim = nc%getDimension('L1_LAITimesteps')
      nLAIs_temp = nc_dim%getLength()
      allocate(LAIBoundaries_temp(nLAIs_temp+1))
      LAIBoundaries_temp = [(ii, ii=1, nLAIs_temp+1)]
    else
      call error_message('***ERROR: no LAI information in restart file for reading')
    end if

    ! check LAI for consistency on all domains (-1 indicates first reading)
    if (nLAI == -1_i4) then
      nLAI = nLAIs_temp
      allocate(LAIBoundaries(nLAI + 1_i4))
      LAIBoundaries = LAIBoundaries_temp
    end if

    call check_dimension_consistency(iDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp)

    call nc%close()

  end subroutine read_nLAI_and_check_dims

  !> \brief checks dimension configurations read from restart file
  !> \authors Robert Schweppe
  !> \date Aug 2019
  subroutine check_dimension_consistency(iDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
      nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp)

    use mo_mpr_global_variables, only: nSoilHorizons_mHM, HorizonDepth_mHM, nLAI, LAIBoundaries ! may read from restart
    use mo_common_variables, only: nLCoverScene, LC_year_start, LC_year_end ! read from nml
    use mo_string_utils, only: compress, num2str

    integer(i4), intent(in) :: iDomain

    integer(i4), intent(in) :: nSoilHorizons_temp, nLAIs_temp, nLandCoverPeriods_temp
    real(dp), dimension(:), intent(inout) :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
          LAIBoundaries_temp
    character(256) :: errorString

    integer(i4) :: k

    ! compare local to global
    call check_consistency_element(nLCoverScene, nLandCoverPeriods_temp, 'number of land cover periods', iDomain)
    call check_consistency_element(nSoilHorizons_mHM, nSoilHorizons_temp, 'number of soil horizons', iDomain)
    call check_consistency_element(nLAI, nLAIs_temp, 'number of LAI timesteps', iDomain)

    ! now check the boundaries
    do k=1, nLCoverScene
      errorString = compress(trim(num2str(k)))//'th land cover boundary'
      call check_consistency_element(real(LC_year_start(k), dp), landCoverPeriodBoundaries_temp(k), errorString, iDomain)
    end do
    errorString = 'last land cover boundary (with 1 year added due to real/int conversion) '
    call check_consistency_element(real(LC_year_end(nLCoverScene) + 1_i4, dp), &
          landCoverPeriodBoundaries_temp(nLCoverScene+1), errorString, iDomain)

    ! last soil horizon is spatially variable, so this is not checked yet
    ! first soil horizon 0 and not contained in HorizonDepth_mHM, so skip that, too
    do k=2, nSoilHorizons_mHM
      errorString = compress(trim(num2str(k)))//'th soil horizon boundary'
      call check_consistency_element(HorizonDepth_mHM(k-1), soilHorizonBoundaries_temp(k), errorString, iDomain)
    end do

    do k=1, nLAI+1
      errorString = compress(trim(num2str(k)))//'th LAI period boundary'
      call check_consistency_element(LAIBoundaries(k), LAIBoundaries_temp(k), errorString, iDomain)
    end do

  end subroutine check_dimension_consistency

  !> \copydoc check_consistency_element
  subroutine check_consistency_element_dp(item1, item2, name, iDomain)
    use mo_utils, only: ne
    use mo_string_utils, only: compress, num2str

    real(dp), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iDomain

    if (ne(item1, item2)) then
    call error_message( &
      'The ', trim(name),&
      ' as set in the configuration file (', &
      compress(trim(num2str(item1))), &
      ') does not conform with domain ', &
      compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(item2))), ').')
    end if
  end subroutine check_consistency_element_dp

  !> \copydoc check_consistency_element
  subroutine check_consistency_element_i4(item1, item2, name, iDomain)
    use mo_string_utils, only: compress, num2str

    integer(i4), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iDomain

    if (item1 /= item2) then
    call error_message( &
      'The ', trim(name),&
      ' as set in the configuration file (', &
      compress(trim(num2str(item1))), &
      ') does not conform with domain ', &
      compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(item2))), ').')
    end if
  end subroutine check_consistency_element_i4

end module mo_common_restart
