!>       \file mo_common_read_data.f90

!>       \brief TODO: add description

!>       \details TODO: add description

!>       \authors Robert Schweppe

!>       \date Jun 2018

! Modifications:

module mo_common_read_data
  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_dem, read_lcover

  ! ------------------------------------------------------------------

CONTAINS

  !    NAME
  !        read_dem

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine read_dem(iDomain, level0_iDomain, data_dp_2d)

    use mo_append, only : append
    use mo_common_file, only : varNameDem
    use mo_common_variables, only : dirIn, level0, domainMeta, resolutionHydrology
    use mo_grid, only : infer_grid_info, Grid
    use mo_message, only : error_message, message
    use mo_string_utils, only : num2str
    use mo_netcdf,           only: NcDataset, NcVariable

    implicit none

    ! loop variables
    integer(i4), intent(in) :: iDomain
    type(Grid), intent(inout) :: level0_iDomain
    real(dp), dimension(:, :), allocatable, intent(out) :: data_dp_2d
    integer(i4) :: domainID

    ! file name of file to read
    character(256) :: fName

    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar          ! variables for data form netcdf


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH DOMAIN
    ! ************************************************
    call message('      Reading dem for domain: ', trim(adjustl(num2str(domainMeta%indices(iDomain)))), ' ...')

    fName = trim(dirIn(iDomain)) // trim(varNameDem) // '.nc'
    ! use the dem variable to create the mask
    call infer_grid_info(fName, 'lon', 'lat', trim(varNameDem), level0_iDomain)

    ! check for L0 and L1 scale consistency
    if(resolutionHydrology(iDomain) .LT. level0_iDomain%cellsize) then
      call error_message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)', &
              new_line('a'), &
              '          check set-up (in mhm.nml) for domain: ', trim(adjustl(num2str(domainID))), ' ...')
    end if

    ! read the Dataset
    nc = NcDataset(fname, "r")
    ! get the variable
    ncVar = nc%getVariable(trim(varNameDem))
    ! read data
    call ncVar%getData(data_dp_2d)
    ! close the handler
    call nc%close()

  end subroutine read_dem

  !    NAME
  !        read_lcover

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine read_lcover(iDomain, dataMatrix_i4, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp)

    use mo_append, only : append, paste
    use mo_common_constants, only : nodata_i4
    use mo_common_file, only : varNameLandCover
    use mo_common_variables, only : L0_LCover, dirIn, level0, domainMeta, nLandCoverPeriods
    use mo_grid, only: Grid
    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_netcdf, only: NcDataset, NcVariable

    implicit none

    integer(i4), intent(in) :: iDomain
    integer(i4), dimension(:, :), allocatable, intent(out) :: dataMatrix_i4
    integer(i4), intent(out), optional :: nLandCoverPeriods_temp
    real(dp), dimension(:), allocatable, intent(out), optional :: landCoverPeriodBoundaries_temp

    integer(i4) :: domainID, iVar
    character(256) :: fName

    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    integer(i4), dimension(:, :, :), allocatable :: data_i4_3d
    logical, dimension(:, :, :), allocatable :: mask_3d

    type(Grid), pointer :: level0_iDomain
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar          ! variables for data form netcdf
    real(dp), dimension(:, :), allocatable :: dummyD2

    call message('      Reading lcover for domain: ', trim(adjustl(num2str(domainMeta%indices(iDomain)))), ' ...')
    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    fName = trim(dirIn(iDomain)) // trim(varNameLandCover) // '.nc'
    ! read the Dataset
    nc = NcDataset(fname, "r")
    ! get the variable
    ncVar = nc%getVariable(trim(varNameLandCover))
    call ncVar%getData(data_i4_3d, mask=mask_3d)
    ! LCover read in is realized seperated because of unknown number of scenes
    do iVar = 1, size(data_i4_3d, 3)
      ! put global nodata value into array (probably not all grid cells have values)
      ! this explicit prior allocation is done so that gFortran does not complain with:
      ! "Fortran runtime error: Array bound mismatch for dimension 1 of array 'data_i4_2d' (0/288)"
      allocate(data_i4_2d(size(data_i4_3d, 1), size(data_i4_3d, 2)))
      data_i4_2d = merge(data_i4_3d(:,:,iVar), nodata_i4, mask_3d(:,:,iVar))
      call paste(dataMatrix_i4, pack(data_i4_2d, level0_iDomain%mask), nodata_i4)
      deallocate(data_i4_2d)
    end do

    if (present(nLandCoverPeriods_temp) .and. present(landCoverPeriodBoundaries_temp)) then
      ! get the landcover dimension
      ncVar = nc%getVariable(trim(varNameLandCover)//'_period_bnds')
      call ncVar%getData(dummyD2)
      nLandCoverPeriods_temp = size(dummyD2, 2)
      allocate(landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1))
      landCoverPeriodBoundaries_temp(1:nLandCoverPeriods_temp) = dummyD2(1,:)
      landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1) = dummyD2(2, nLandCoverPeriods_temp)
    end if

    call nc%close()

  end subroutine read_lcover

end module mo_common_read_data
