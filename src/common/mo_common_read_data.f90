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
    use mo_common_variables, only : Grid, L0_Domain, L0_elev, dirMorpho, level0, domainMeta, nuniqueL0Domains, &
                                    resolutionHydrology
    use mo_grid, only : set_basin_indices
    use mo_message, only : message
    use mo_read_spatial_data, only : read_header_ascii, read_spatial_data_ascii
    use mo_string_utils, only : num2str

    implicit none

    ! loop variables
    integer(i4) :: domainID, iDomain

    ! file name of file to read
    character(256) :: fName

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    type(Grid), pointer :: level0_iBasin


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH DOMAIN
    ! ************************************************
    ! allocate necessary variables at Level0
    allocate(level0(nuniqueL0Domains))

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      level0_iBasin => level0(L0_Domain(iDomain))

      ! check whether L0 data is shared
      if (iDomain .gt. 1) then
        ! ToDo: check change
        if (L0_Domain(iDomain) < iDomain) then
          !
          call message('      Using dem of domain ', &
                  trim(adjustl(num2str(domainMeta%indices(L0_Domain(iDomain))))), ' for domain: ',&
                  trim(adjustl(num2str(iDomain))), '...')

          ! DO NOT read L0 data
          cycle

        end if
      end if

      call message('      Reading dem for domain: ', trim(adjustl(num2str(domainID))), ' ...')

      ! Header (to check consistency)
      fName = trim(adjustl(dirMorpho(iDomain))) // trim(adjustl(file_dem))
      call read_header_ascii(trim(fName), udem, &
              level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
              level0_iBasin%yllcorner, level0_iBasin%cellsize, level0_iBasin%nodata_value)

      ! check for L0 and L1 scale consistency
      if(resolutionHydrology(iDomain) .LT. level0_iBasin%cellsize) then
        call message()
        call message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)')
        call message('          check set-up (in mhm.nml) for domain: ', trim(adjustl(num2str(domainID))), ' ...')
        stop
      end if

      ! DEM + overall mask creation
      fName = trim(adjustl(dirMorpho(iDomain))) // trim(adjustl(file_dem))
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

  !    NAME
  !        read_lcover

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine read_lcover

    use mo_append, only : append, paste
    use mo_common_constants, only : nodata_i4
    use mo_common_file, only : ulcoverclass
    use mo_common_variables, only : Grid, L0_Domain, L0_LCover, LCfilename, dirLCover, level0, domainMeta, nLCoverScene
    use mo_message, only : message
    use mo_read_spatial_data, only : read_spatial_data_ascii
    use mo_string_utils, only : num2str

    implicit none

    ! loop variables
    integer(i4) :: domainID, iDomain, iVar

    ! file name of file to read
    character(256) :: fName

    integer(i4), dimension(:, :), allocatable :: data_i4_2d

    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4

    logical, dimension(:, :), allocatable :: mask_2d

    type(Grid), pointer :: level0_iBasin


    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      level0_iBasin => level0(L0_Domain(iDomain))

      ! check whether L0 data is shared
      if (iDomain .gt. 1) then
        if (L0_Domain(iDomain) .eq. L0_Domain(iDomain - 1)) then
          call message('      Using lcover of domain ', &
                  trim(adjustl(num2str(L0_Domain(iDomain)))), ' for domain: ',&
                  trim(adjustl(num2str(domainID))), '...')
          ! DO NOT read L0 data
          cycle

        end if
      end if

      call message('      Reading lcover for domain: ', trim(adjustl(num2str(domainID))), ' ...')

      ! LCover read in is realized seperated because of unknown number of scenes
      do iVar = 1, nLCoverScene
        fName = trim(adjustl(dirLCover(iDomain))) // trim(adjustl(LCfilename(iVar)))
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
