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

  subroutine read_dem

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_file, only : file_dem, udem
    use mo_common_variables, only : Grid, L0_Basin, L0_elev, dirMorpho, level0, nBasins, nuniquel0Basins, &
                                    resolutionHydrology
    use mo_grid, only : set_basin_indices
    use mo_message, only : message
    use mo_read_spatial_data, only : read_header_ascii, read_spatial_data_ascii
    use mo_string_utils, only : num2str

    implicit none

    ! loop variables
    integer(i4) :: iBasin

    ! file name of file to read
    character(256) :: fName

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    type(Grid), pointer :: level0_iBasin


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH BASIN
    ! ************************************************
    ! allocate necessary variables at Level0
    allocate(level0(nuniquel0Basins))

    do iBasin = 1, nBasins

      level0_iBasin => level0(L0_Basin(iBasin))

      ! check whether L0 data is shared
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          !
          call message('      Using dem of basin ', &
                  trim(adjustl(num2str(L0_Basin(iBasin)))), ' for basin: ',&
                  trim(adjustl(num2str(iBasin))), '...')

          ! DO NOT read L0 data
          cycle

        end if
      end if

      call message('      Reading dem for basin: ', trim(adjustl(num2str(iBasin))), ' ...')

      ! Header (to check consistency)
      fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_dem))
      call read_header_ascii(trim(fName), udem, &
              level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
              level0_iBasin%yllcorner, level0_iBasin%cellsize, level0_iBasin%nodata_value)

      ! check for L0 and L1 scale consistency
      if(resolutionHydrology(iBasin) .LT. level0_iBasin%cellsize) then
        call message()
        call message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)')
        call message('          check set-up (in mhm.nml) for basin: ', trim(adjustl(num2str(iBasin))), ' ...')
        stop
      end if

      ! DEM + overall mask creation
      fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_dem))
      call read_spatial_data_ascii(trim(fName), udem, &
              level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
              level0_iBasin%yllcorner, level0_iBasin%cellsize, data_dp_2d, level0_iBasin%mask)

      ! put global nodata value into array (probably not all grid cells have values)
      data_dp_2d = merge(data_dp_2d, nodata_dp, level0_iBasin%mask)
      ! put data in variable
      call append(L0_elev, pack(data_dp_2d, level0_iBasin%mask))
      ! deallocate arrays
      deallocate(data_dp_2d)

      level0_iBasin%nCells = count(level0_iBasin%mask)

    end do

    call set_basin_indices(level0)

  end subroutine read_dem

  subroutine read_lcover

    USE mo_read_spatial_data, ONLY : &
            read_spatial_data_ascii
    USE mo_append, ONLY : append, paste
    USE mo_string_utils, ONLY : num2str
    USE mo_message, ONLY : message
    !
    USE mo_common_file, ONLY : &
            ulcoverclass ! unit of land cover class map
    use mo_common_variables, ONLY : &
            L0_LCover, & ! classical mHM land cover class (L0)
            level0, & ! grid information (ncols, nrows, ..)
            dirLCover, & ! directories
            L0_Basin, & ! L0_Basin ID
            LCfilename, nLCoverScene, & ! file names and number of land cover scenes
            nBasins, & ! number of basins
            Grid
    USE mo_common_constants, ONLY : nodata_i4                   ! mHM's global nodata vales

    implicit none

    ! local variables
    integer(i4) :: iBasin, iVar                  ! loop variables
    character(256) :: fName                      ! file name of file to read
    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4
    logical, dimension(:, :), allocatable :: mask_2d
    type(Grid), pointer :: level0_iBasin

    do iBasin = 1, nBasins

      level0_iBasin => level0(L0_Basin(iBasin))

      ! check whether L0 data is shared
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          call message('      Using lcover of basin ', &
                  trim(adjustl(num2str(L0_Basin(iBasin)))), ' for basin: ',&
                  trim(adjustl(num2str(iBasin))), '...')
          ! DO NOT read L0 data
          cycle

        end if
      end if

      call message('      Reading lcover for basin: ', trim(adjustl(num2str(iBasin))), ' ...')

      ! LCover read in is realized seperated because of unknown number of scenes
      do iVar = 1, nLCoverScene
        fName = trim(adjustl(dirLCover(iBasin))) // trim(adjustl(LCfilename(iVar)))
        call read_spatial_data_ascii(trim(fName), ulcoverclass, &
                level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
                level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        call paste(dataMatrix_i4, pack(data_i4_2d, level0_iBasin%mask), nodata_i4)
        deallocate(data_i4_2d)
      end do
      call append(L0_LCover, dataMatrix_i4)
      deallocate(dataMatrix_i4)

    end do

  end subroutine read_lcover

end module mo_common_read_data