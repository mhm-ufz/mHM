!> \file mo_common_read_data.f90
!> \brief   \copybrief mo_common_read_data
!> \details \copydetails mo_common_read_data

!> \brief Common reading routines
!> \details Routines to read the DEM and landcover files.
!> \authors Robert Schweppe
!> \date Jun 2018
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_read_data
  USE mo_kind, ONLY : i4, dp
  use mo_message, only: message, error_message

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
    use mo_common_types, only: Grid
    use mo_common_variables, only : L0_elev, dirMorpho, level0, domainMeta, &
                                    resolutionHydrology
    use mo_grid, only : set_domain_indices
    use mo_read_spatial_data, only : read_header_ascii, read_spatial_data_ascii
    use mo_string_utils, only : num2str

    implicit none

    ! loop variables
    integer(i4) :: domainID, iDomain

    ! file name of file to read
    character(256) :: fName

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    type(Grid), pointer :: level0_iDomain


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH DOMAIN
    ! ************************************************
    ! allocate necessary variables at Level0
    allocate(level0(domainMeta%nDomains))

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

      ! check whether L0 data is shared
      if (iDomain .gt. 1) then
        if (domainMeta%L0DataFrom(iDomain) < iDomain) then
          !
          call message('      Using dem of domain ', &
                  trim(adjustl(num2str(domainMeta%indices(domainMeta%L0DataFrom(iDomain))))), ' for domain: ',&
                  trim(adjustl(num2str(iDomain))), '...')

          ! DO NOT read L0 data
          cycle

        end if
      end if

      call message('      Reading dem for domain: ', trim(adjustl(num2str(domainID))), ' ...')

      ! Header (to check consistency)
      fName = trim(adjustl(dirMorpho(iDomain))) // trim(adjustl(file_dem))
      call read_header_ascii(trim(fName), udem, &
              level0_iDomain%nrows, level0_iDomain%ncols, level0_iDomain%xllcorner, &
              level0_iDomain%yllcorner, level0_iDomain%cellsize, level0_iDomain%nodata_value)

      ! check for L0 and L1 scale consistency
      if(resolutionHydrology(iDomain) .LT. level0_iDomain%cellsize) then
        call error_message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)', raise=.false.)
        call error_message('          check set-up (in mhm.nml) for domain: ', trim(adjustl(num2str(domainID))), ' ...')
      end if

      ! DEM + overall mask creation
      fName = trim(adjustl(dirMorpho(iDomain))) // trim(adjustl(file_dem))
      call read_spatial_data_ascii(trim(fName), udem, &
              level0_iDomain%nrows, level0_iDomain%ncols, level0_iDomain%xllcorner, &
              level0_iDomain%yllcorner, level0_iDomain%cellsize, data_dp_2d, level0_iDomain%mask)

      ! put global nodata value into array (probably not all grid cells have values)
      data_dp_2d = merge(data_dp_2d, nodata_dp, level0_iDomain%mask)
      ! put data in variable
      call append(L0_elev, pack(data_dp_2d, level0_iDomain%mask))
      ! deallocate arrays
      deallocate(data_dp_2d)

      level0_iDomain%nCells = count(level0_iDomain%mask)

    end do

    call set_domain_indices(level0, indices=domainMeta%L0DataFrom)

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
    use mo_common_types, only: Grid
    use mo_common_variables, only : L0_LCover, LCfilename, dirLCover, level0, domainMeta, nLCoverScene
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

    type(Grid), pointer :: level0_iDomain


    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

      ! check whether L0 data is shared
      ! ToDo: check change
      if (domainMeta%L0DataFrom(iDomain) < iDomain) then
        call message('      Using lcover of domain ', &
                trim(adjustl(num2str(domainMeta%L0DataFrom(iDomain)))), ' for domain: ',&
                trim(adjustl(num2str(domainID))), '...')
        ! DO NOT read L0 data
        cycle

      end if

      call message('      Reading lcover for domain: ', trim(adjustl(num2str(domainID))), ' ...')

      ! LCover read in is realized seperated because of unknown number of scenes
      do iVar = 1, nLCoverScene
        fName = trim(adjustl(dirLCover(iDomain))) // trim(adjustl(LCfilename(iVar)))
        call read_spatial_data_ascii(trim(fName), ulcoverclass, &
                level0_iDomain%nrows, level0_iDomain%ncols, level0_iDomain%xllcorner, &
                level0_iDomain%yllcorner, level0_iDomain%cellsize, data_i4_2d, mask_2d)
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        call paste(dataMatrix_i4, pack(data_i4_2d, level0_iDomain%mask), nodata_i4)
        deallocate(data_i4_2d)
      end do
      call append(L0_LCover, dataMatrix_i4)
      deallocate(dataMatrix_i4)

    end do

  end subroutine read_lcover

end module mo_common_read_data
