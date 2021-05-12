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
    use mo_common_file, only : varNameDem
    use mo_common_variables, only : Grid,  L0_elev, dirMorpho, level0, domainMeta, &
                                    resolutionHydrology
    use mo_grid, only : set_domain_indices, infer_grid_info
    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_netcdf,           only: NcDataset, NcVariable

    implicit none

    ! loop variables
    integer(i4) :: domainID, iDomain

    ! file name of file to read
    character(256) :: fName

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    type(Grid), pointer :: level0_iDomain
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar          ! variables for data form netcdf


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

      fName = trim(dirMorpho(iDomain)) // trim(varNameDem) // '.nc'
      ! use the dem variable to create the mask
      call infer_grid_info(fName, 'lon', 'lat', trim(varNameDem), level0_iDomain)

      ! check for L0 and L1 scale consistency
      if(resolutionHydrology(iDomain) .LT. level0_iDomain%cellsize) then
        call message()
        call message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)')
        call message('          check set-up (in mhm.nml) for domain: ', trim(adjustl(num2str(domainID))), ' ...')
        stop
      end if

      ! read the Dataset
      nc = NcDataset(fname, "r")
      ! get the variable
      ncVar = nc%getVariable(trim(varNameDem))
      ! read data
      call ncVar%getData(data_dp_2d)
      ! put data in variable
      call append(L0_elev, pack(data_dp_2d, level0_iDomain%mask))
      ! deallocate arrays
      deallocate(data_dp_2d)

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
    use mo_common_file, only : varNameLandCover
    use mo_common_variables, only : Grid, L0_LCover, LCfilename, dirLCover, level0, domainMeta, nLCoverScene
    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_netcdf, only: NcDataset, NcVariable

    implicit none

    integer(i4) :: domainID, iDomain, iVar
    character(256) :: fName

    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    integer(i4), dimension(:, :, :), allocatable :: data_i4_3d
    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4
    logical, dimension(:, :, :), allocatable :: mask_3d

    type(Grid), pointer :: level0_iDomain
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar          ! variables for data form netcdf


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

      fName = trim(dirLCover(iDomain)) // trim(varNameLandCover) // '.nc'
      ! read the Dataset
      nc = NcDataset(fname, "r")
      ! get the variable
      ncVar = nc%getVariable(trim(varNameLandCover))
      print*, 'got ', trim(varNameLandCover)
      call ncVar%getData(data_i4_3d, mask=mask_3d)
      print*, 'got data ', shape(data_i4_3d)
      ! LCover read in is realized seperated because of unknown number of scenes
      do iVar = 1, nLCoverScene
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_3d(:,:,iVar), nodata_i4, mask_3d(:,:,iVar))
        call paste(dataMatrix_i4, pack(data_i4_2d, level0_iDomain%mask), nodata_i4)
        deallocate(data_i4_2d)
      end do
      call append(L0_LCover, dataMatrix_i4)
      deallocate(dataMatrix_i4)

      call nc%close()
    end do

  end subroutine read_lcover

end module mo_common_read_data
