!>       \file mo_prepare_gridded_lai.f90

!>       \brief Prepare daily LAI fields (e.g., MODIS data) for mHM

!>       \details Prepare daily LAI fields(e.g., MODIS data) for mHM

!>       \authors John Craven & Rohini Kumar

!>       \date Aug 2013

! Modifications:

MODULE mo_prepare_gridded_LAI

  ! This module provides routines to read daily gridded LAI data.

  ! Written  John Craven & Rohini Kumar, August 2013
  ! Modified from mo_meteo_forcings

  USE mo_kind, ONLY : i4, dp

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

  subroutine prepare_gridded_daily_LAI_data(iDomain, nrows, ncols, mask, LAIPer_iDomain, nLAIs_temp, LAIBoundaries_temp)

    use mo_append, only : append
    use mo_common_datetime_type, only : period, timeStep_LAI_input
    use mo_message, only : message
    use mo_mpr_global_variables, only : L0_gridded_LAI, dirgridded_LAI, inputFormat_gridded_LAI
    use mo_read_nc, only : read_nc

    implicit none

    ! domain Id
    integer(i4), intent(in) :: iDomain, nrows, ncols

    logical, dimension(:, :), intent(in) :: mask

    type(period), intent(in), optional :: LAIPer_iDomain
    integer(i4), intent(out) :: nLAIs_temp
    real(dp), dimension(:), allocatable, intent(out) :: LAIBoundaries_temp

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
      call message()
      call message('***ERROR: No recognized input format')
      stop 1

    END SELECT

    ! pack variables
    nCells = count(mask)

    nLAIs_temp = size(LAI0_3D, 3)
    allocate(LAIBoundaries_temp(nLAIs_temp+1))
    LAIBoundaries_temp = [(iLAI, iLAI=1, nLAIs_temp+1)]
    allocate(LAI0_2D(nCells, nLAIs_temp))

    do iLAI = 1, nLAIs_temp
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

  subroutine prepare_gridded_mean_monthly_LAI_data(iDomain, nrows, ncols, mask, nLAIs_temp, LAIBoundaries_temp)

    use mo_append, only : append
    use mo_message, only : message
    use mo_mpr_global_variables, only : L0_gridded_LAI, dirgridded_LAI
    use mo_ncread, only : Get_NcDim, Get_NcVar, Get_NcVarAtt
    use mo_string_utils, only : num2str
    use mo_utils, only : eq

    implicit none

    ! domain Id
    integer(i4), intent(in) :: iDomain, nrows, ncols

    logical, dimension(:, :), intent(in) :: mask
    integer(i4), intent(out) :: nLAIs_temp
    real(dp), dimension(:), allocatable, intent(out) :: LAIBoundaries_temp

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
      stop '***ERROR: read_nc: mHM generated x and y are not matching NetCDF dimensions'
    end if
    if (dimen(3) .ne. 12) then
      stop '***ERROR: read_nc: the time dimenion of LAI NetCDF file under the option-1 is not 12'
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
        call message('***ERROR: read_nc: nodata value within domain ')
        call message('          boundary in variable: ', 'lai')
        call message('          at timestep         : ', trim(num2str(t)))
        stop
      end if
      ! optional check
      if (any((LAI0_3D(:, :, t) .lt. 0.0_dp) .AND. mask(:, :))) then
        call message('***ERROR: read_nc: values in variable lai are lower than ', trim(num2str(0, '(F7.2)')))
        call message('          at timestep  : ', trim(num2str(t)))
        call message('File: ', trim(fName))
        call message('Minval at timestep: ', trim(num2str(minval(LAI0_3D(:, :, t)), '(F7.2)')))
        call message('Total minval: ', trim(num2str(minval(LAI0_3D(:, :, :)), '(F7.2)')))
        stop
      end if

      if (any((LAI0_3D(:, :, t) .gt. 30.0_dp) .AND. mask(:, :))) then
        call message('***ERROR: read_nc: values in variable lai are greater than ', trim(num2str(30, '(F7.2)')))
        call message('          at timestep  : ', trim(num2str(t)))
        call message('File: ', trim(fName))
        call message('Maxval at timestep: ', trim(num2str(maxval(LAI0_3D(:, :, t)), '(F7.2)')))
        call message('Total maxval: ', trim(num2str(maxval(LAI0_3D(:, :, :)), '(F7.2)')))
        stop
      end if
    end do

    ! pack variables
    nCells = count(mask)

    nLAIs_temp = size(LAI0_3D, 3)
    allocate(LAIBoundaries_temp(nLAIs_temp+1))
    LAIBoundaries_temp = [(iLAI, iLAI=1, nLAIs_temp+1)]

    allocate(LAI0_2D(nCells, nLAIs_temp))
    do iLAI = 1, nLAIs_temp
      LAI0_2D(:, iLAI) = pack(LAI0_3D(:, :, iLAI), MASK = mask(:, :))
    end do

    ! append to Global variable
    call append(L0_gridded_LAI, LAI0_2D(:, :))

    !free space
    deallocate(LAI0_2D, LAI0_3D)

  end subroutine prepare_gridded_mean_monthly_LAI_data


END MODULE mo_prepare_gridded_LAI
