!> \file mo_prepare_gridded_lai.f90
!> \brief \copybrief mo_prepare_gridded_lai
!> \details \copydetails mo_prepare_gridded_lai

!> \brief Prepare daily LAI fields (e.g., MODIS data) for mHM
!> \details Prepare daily LAI fields(e.g., MODIS data) for mHM
!> \authors John Craven & Rohini Kumar
!> \date Aug 2013
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_prepare_gridded_LAI

  ! This module provides routines to read daily gridded LAI data.

  ! Written  John Craven & Rohini Kumar, August 2013
  ! Modified from mo_meteo_forcings

  USE mo_kind, ONLY : i4, dp
  use mo_message, only: error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prepare_gridded_daily_LAI_data
  PUBLIC :: prepare_gridded_mean_monthly_LAI_data

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        prepare_gridded_daily_LAI_data

  !    PURPOSE
  !>       \brief Prepare gridded daily LAI data

  !>       \details Prepare gridded daily LAI data at Level-0 (e.g., using MODIS datasets)

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "logical, dimension(:, :) :: mask"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "type(period), optional :: LAIPer_iDomain"

  !    HISTORY
  !>       \authors John Craven & Rohini Kumar

  !>       \date Aug 2013

  ! Modifications:
  ! Matthias Cuntz & Juliane Mai Nov 2014 - use meteo reading routines
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine prepare_gridded_daily_LAI_data(iDomain, nrows, ncols, mask, LAIPer_iDomain)

    use mo_append, only : append
    use mo_common_types, only: period
    use mo_mpr_global_variables, only : L0_gridded_LAI, dirgridded_LAI, inputFormat_gridded_LAI, &
            nLAI, LAIBoundaries, timeStep_LAI_input
    use mo_read_nc, only : read_nc

    implicit none

    ! domain Id
    integer(i4), intent(in) :: iDomain, nrows, ncols

    logical, dimension(:, :), intent(in) :: mask

    type(period), intent(in), optional :: LAIPer_iDomain

    integer(i4) :: ncells, iLAI

    ! data at level-0 [nRow X nCols X nTimeSteps]
    real(dp), dimension(:, :, :), allocatable :: LAI0_3D

    ! data at level-0 [nCells X nTimeSteps]
    real(dp), dimension(:, :), allocatable :: LAI0_2D


    ! select case depending on input data format
    SELECT CASE(trim(inputFormat_gridded_LAI))

    ! netcdf file input option
    CASE('nc')
      CALL read_nc(dirgridded_LAI(iDomain), nRows, nCols, &
              'lai', mask, LAI0_3D, target_period = LAIPer_iDomain, &
              lower = 1.00E-10_dp, upper = 30.0_dp, nctimestep = timeStep_LAI_input)
    CASE DEFAULT
      call error_message('***ERROR: No recognized input format')

    END SELECT

    ! pack variables
    nCells = count(mask)
    ! only set if not yet allocated (e.g. domain 1)
    if (.not. allocated(LAIBoundaries)) then
      nLAI = size(LAI0_3D, 3)
      allocate(LAIBoundaries(nLAI+1))
      LAIBoundaries = [(iLAI, iLAI=1, nLAI+1)]
    end if
    allocate(LAI0_2D(nCells, nLAI))

    do iLAI = 1, nLAI
      LAI0_2D(:, iLAI) = pack(LAI0_3D(:, :, iLAI), MASK = mask(:, :))
    end do

    ! append to Global variable
    call append(L0_gridded_LAI, LAI0_2D(:, :))

    !free space
    deallocate(LAI0_2D, LAI0_3D)

  end subroutine prepare_gridded_daily_LAI_data

  ! ------------------------------------------------------------------

  !    NAME
  !        prepare_gridded_mean_monthly_LAI_data

  !    PURPOSE
  !>       \brief prepare_gridded_mean_monthly_LAI_data

  !>       \details Long term mean monthly gridded LAI data at Level-0 (e.g., using MODIS datasets)
  !>       The netcdf file should contain 12 (calender months) gridded fields of climatological
  !>       LAI data at the input L0 data resolution.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "integer(i4) :: iDomain, nrows, ncols" domain Id
  !>       \param[in] "logical, dimension(:, :) :: mask"

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Dec 2016

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine prepare_gridded_mean_monthly_LAI_data(iDomain, nrows, ncols, mask)

    use mo_append, only : append
    use mo_mpr_global_variables, only : L0_gridded_LAI, dirgridded_LAI, nLAI, LAIBoundaries
    use mo_ncread, only : Get_NcDim, Get_NcVar, Get_NcVarAtt
    use mo_string_utils, only : num2str
    use mo_utils, only : eq

    implicit none

    ! domain Id
    integer(i4), intent(in) :: iDomain, nrows, ncols

    logical, dimension(:, :), intent(in) :: mask

    integer(i4) :: ncells, iLAI

    ! data at level-0 [nRow X nCols X nTimeSteps]
    real(dp), dimension(:, :, :), allocatable :: LAI0_3D

    ! data at level-0 [nCells X nTimeSteps]
    real(dp), dimension(:, :), allocatable :: LAI0_2D

    integer(i4) :: t

    ! name of NetCDF file
    character(256) :: fName

    ! netcdf attribute values
    character(256) :: AttValues

    ! datatype of attribute
    integer(i4) :: datatype

    ! dimension for NetCDF file
    integer(i4), dimension(5) :: dimen

    ! data nodata value
    real(dp) :: nodata_value


    fName = trim(dirgridded_LAI(iDomain)) // trim('lai.nc')

    ! get dimensions
    dimen = Get_NcDim(trim(fName), 'lai')
    if ((dimen(1) .ne. nRows) .or. (dimen(2) .ne. nCols)) then
       call error_message('***ERROR: read_nc: mHM generated x and y are not matching NetCDF dimensions')
    end if
    if (dimen(3) .ne. 12) then
       call error_message('***ERROR: read_nc: the time dimenion of LAI NetCDF file under the option-1 is not 12')
    end if

    ! determine no data value
    call Get_NcVarAtt(trim(fName), 'lai', '_FillValue', AttValues, dtype = datatype)
    ! convert to number
    read(AttValues, *) nodata_value

    call Get_NcVar(trim(fName), 'lai', LAI0_3D)

    ! start checking values
    do t = 1, dimen(3)
      ! checking for nodata values if optional nocheck is given
      if (any(eq(LAI0_3D(:, :, t), nodata_value) .and. (mask))) then
        call error_message('***ERROR: read_nc: nodata value within domain ', raise=.false.)
        call error_message('          boundary in variable: ', 'lai', raise=.false.)
        call error_message('          at timestep         : ', trim(num2str(t)))
      end if
      ! optional check
      if (any((LAI0_3D(:, :, t) .lt. 0.0_dp) .AND. mask(:, :))) then
        call error_message('***ERROR: read_nc: values in variable lai are lower than ', trim(num2str(0, '(F7.2)')), raise=.false.)
        call error_message('          at timestep  : ', trim(num2str(t)), raise=.false.)
        call error_message('File: ', trim(fName), raise=.false.)
        call error_message('Minval at timestep: ', trim(num2str(minval(LAI0_3D(:, :, t)), '(F7.2)')), raise=.false.)
        call error_message('Total minval: ', trim(num2str(minval(LAI0_3D(:, :, :)), '(F7.2)')))
      end if

      if (any((LAI0_3D(:, :, t) .gt. 30.0_dp) .AND. mask(:, :))) then
        call error_message('***ERROR: read_nc: values in variable lai are greater than ', trim(num2str(30, '(F7.2)')), raise=.false.)
        call error_message('          at timestep  : ', trim(num2str(t)), raise=.false.)
        call error_message('File: ', trim(fName), raise=.false.)
        call error_message('Maxval at timestep: ', trim(num2str(maxval(LAI0_3D(:, :, t)), '(F7.2)')), raise=.false.)
        call error_message('Total maxval: ', trim(num2str(maxval(LAI0_3D(:, :, :)), '(F7.2)')))
      end if
    end do

    ! pack variables
    nCells = count(mask)
    ! only set if not yet allocated (e.g. domain 1)
    if (.not. allocated(LAIBoundaries)) then
      nLAI = size(LAI0_3D, 3)
      allocate(LAIBoundaries(nLAI+1))
      LAIBoundaries = [(iLAI, iLAI=1, nLAI+1)]
    end if
    allocate(LAI0_2D(nCells, nLAI))
    do iLAI = 1, nLAI
      LAI0_2D(:, iLAI) = pack(LAI0_3D(:, :, iLAI), MASK = mask(:, :))
    end do

    ! append to Global variable
    call append(L0_gridded_LAI, LAI0_2D(:, :))

    !free space
    deallocate(LAI0_2D, LAI0_3D)

  end subroutine prepare_gridded_mean_monthly_LAI_data


END MODULE mo_prepare_gridded_LAI
