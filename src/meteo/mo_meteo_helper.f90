!> \file mo_meteo_helper.f90
!> \brief \copybrief mo_meteo_helper
!> \details \copydetails mo_meteo_helper

!> \brief Prepare meteorological forcings data for mHM.
!> \details Prepare meteorological forcings data for mHM.
!> \authors Rohini Kumar
!> \date Jan 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_meteo
MODULE mo_meteo_helper

  ! This module provides routines to read meteorological data.

  ! Written  Rohini Kumar, Jan 2013
  ! Modified

  USE mo_kind, ONLY : i4, dp
  use mo_message, only: message, error_message

  IMPLICIT NONE

  PRIVATE

  public :: meteo_forcings_wrapper
  public :: meteo_weights_wrapper
  public :: chunk_config

CONTAINS

  !> \brief Prepare meteorological forcings data for mHM at Level-1
  !> \details Prepare meteorological forcings data for mHM, which include
  !! 1) Reading meteo. datasets at their native resolution for every Domain
  !! 2) Perform aggregation or disaggregation of meteo. datasets from their
  !! native resolution (level-2) to the required hydrologic resolution (level-1)
  !! 3) Pad the above datasets of every Domain to their respective global ones
  !> \changelog
  !! - Stephan Thober Jun 2014
  !!   - changed to readPer
  !! - Stephan Thober Feb 2016
  !!   - refactored deallocate statements
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - used by meteo-handler
  !!   - now independent of global variables
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine meteo_forcings_wrapper( &
    iDomain, dataPath, inputFormat, dataOut1, readPer, nTstepForcingDay, level1, level2, lower, upper, ncvarName, bound_error)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_types, only : period, grid
    use mo_read_nc, only : read_nc
    use mo_meteo_spatial_tools, only : spatial_aggregation, spatial_disaggregation

    implicit none

    !> Domain Id
    integer(i4), intent(in) :: iDomain
    !> Data path where a given meteo. variable is stored
    character(len = *), intent(in) :: dataPath
    !> only 'nc' possible at the moment
    character(len = *), intent(in) :: inputFormat
    !> Packed meterological variable for the whole simulation period
    real(dp), dimension(:, :), allocatable, intent(inout) :: dataOut1
    !> start and end dates of read period
    type(period), intent(in) :: readPer
    !> Number of forcing intervals per day
    integer(i4), intent(inout) :: nTstepForcingDay
    !> grid information at hydrologic level
    type(Grid), dimension(:), intent(in) :: level1
    !> Reference of the metereological variables
    type(Grid), dimension(:), intent(in) :: level2
    !> Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower
    !> Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper
    !> name of the variable (for .nc files)
    character(len = *), optional, intent(in) :: ncvarName
    !> .FALSE. to only warn about bound (lower, upper) violations, default = .TRUE. - raise an error
    logical, optional, intent(in) :: bound_error

    logical, dimension(:, :), allocatable :: mask1
    integer(i4) :: ncells1
    integer(i4) :: nrows2, ncols2
    logical, dimension(:, :), allocatable :: mask2
    ! meteo data at level-2
    real(dp), dimension(:, :, :), allocatable :: L2_data
    ! meteo data at level-1
    real(dp), dimension(:, :, :), allocatable :: L1_data
    ! packed meteo data at level-1 from 3D to 2D
    real(dp), dimension(:, :), allocatable :: L1_data_packed
    integer(i4) :: nTimeSteps
    ! level-1_resolution/level-2_resolution
    real(dp) :: cellFactorHbyM
    integer(i4) :: t

    ! get basic Domain information at level-1
    nCells1 = level1(iDomain)%nCells
    mask1 = level1(iDomain)%mask

    ! make  basic Domain information at level-2
    nrows2 = level2(iDomain)%nrows
    ncols2 = level2(iDomain)%ncols
    mask2 = level2(iDomain)%mask


    select case (trim(inputFormat))
      case('nc')
        ! Fortran 2008 allowes to pass optional dummy arguments to another routine (when they are optional there too)
        CALL read_nc(dataPath, nRows2, nCols2, ncvarName, mask2, L2_data, target_period=readPer, &
                     lower=lower, upper=upper, is_meteo=.True., bound_error=bound_error, nTstepForcingDay=nTstepForcingDay)
      case DEFAULT
        call error_message('***ERROR: meteo_forcings_wrapper: Not recognized input format')
    end select
    ! cellfactor to decide on the upscaling or downscaling of meteo. fields
    cellFactorHbyM = level1(iDomain)%cellsize / level2(iDomain)%cellsize

    ! upscaling & packing
    if(cellFactorHbyM .gt. 1.0_dp) then
      call spatial_aggregation(L2_data, level2(iDomain)%cellsize, level1(iDomain)%cellsize, mask1, mask2, L1_data)
      ! downscaling
    elseif(cellFactorHbyM .lt. 1.0_dp) then
      call spatial_disaggregation(L2_data, level2(iDomain)%cellsize, level1(iDomain)%cellsize, mask1, mask2, L1_data)
      ! nothing
    else
      allocate(L1_data(size(L2_data, 1), size(L2_data, 2), size(L2_data, 3)))
      L1_data(:, :, :) = L2_data(:, :, :)
    end if
    ! free memory immediately
    deallocate(L2_data)

    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate(L1_data_packed(nCells1, nTimeSteps))
    do t = 1, nTimeSteps

      L1_data_packed(:, t) = pack(L1_data(:, :, t), MASK = mask1(:, :))

    end do

    ! free memory immediately
    deallocate(L1_data)

    ! append
    call append(dataOut1, L1_data_packed(:, :), nodata_dp)

    !free space
    deallocate(L1_data_packed)

  end subroutine meteo_forcings_wrapper


  !> \brief Prepare weights for meteorological forcings data for mHM at Level-1
  !> \details Prepare meteorological weights data for mHM, which include
  !! 1) Reading meteo. weights datasets at their native resolution for every Domain
  !! 2) Perform aggregation or disaggregation of meteo. weights datasets from their
  !! native resolution (level-2) to the required hydrologic resolution (level-1)
  !! 3) Pad the above datasets of every Domain to their respective global ones
  !> \changelog
  !! - Stephan Thober May 2017
  !!   - updated documentation
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Oldrich Rakovec Aug 2018
  !!   - adding message about reading meteo_weights
  !! - Sebastian Müller Mar 2023
  !!   - used by meteo-handler
  !!   - now independent of global variables
  !> \authors Stephan Thober & Rohini Kumar
  !> \date Jan 2017
  subroutine meteo_weights_wrapper(iDomain, read_meteo_weights, dataPath, dataOut1, level1, level2, lower, upper, ncvarName)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_types, only : grid
    use mo_read_nc, only : read_weights_nc
    use mo_meteo_spatial_tools, only : spatial_aggregation, spatial_disaggregation

    implicit none

    !> Domain Id
    integer(i4), intent(in) :: iDomain
    !> Flag for reading meteo weights
    logical, intent(in) :: read_meteo_weights
    !> Data path where a given meteo. variable is stored
    character(len = *), intent(in) :: dataPath
    !> Packed meterological variable for the whole simulation period
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: dataOut1
    !> grid information at hydrologic level
    type(Grid), dimension(:), intent(in) :: level1
    !> Reference of the metereological variables
    type(Grid), dimension(:), intent(in) :: level2
    !> Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower
    !> Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper
    !> name of the variable (for .nc files)
    character(len = *), optional, intent(in) :: ncvarName

    logical, dimension(:, :), allocatable :: mask1
    integer(i4) :: ncells1
    integer(i4) :: nrows2, ncols2
    logical, dimension(:, :), allocatable :: mask2
    ! meteo weights data at level-2
    real(dp), dimension(:, :, :, :), allocatable :: L2_data
    ! meteo weights data at level-1
    real(dp), dimension(:, :, :, :), allocatable :: L1_data
    ! packed meteo weights data at level-1 from 4D to 3D
    real(dp), dimension(:, :, :), allocatable :: L1_data_packed
    integer(i4) :: nMonths, nHours
    ! level-1_resolution/level-2_resolution
    real(dp) :: cellFactorHbyM
    integer(i4) :: t, j

    ! get basic Domain information at level-1
    nCells1 = level1(iDomain)%nCells
    mask1 = level1(iDomain)%mask

    ! make  basic Domain information at level-2
    nrows2 = level2(iDomain)%nrows
    ncols2 = level2(iDomain)%ncols
    mask2 = level2(iDomain)%mask

    if (read_meteo_weights) then
      call message('  read_meteo_weights = .TRUE. ... Reading meteo weights ... ')
      CALL read_weights_nc(dataPath, nRows2, nCols2, ncvarName, L2_data, mask2, lower=lower, upper=upper)

      ! cellfactor to decide on the upscaling or downscaling of meteo. fields
      cellFactorHbyM = level1(iDomain)%cellsize / level2(iDomain)%cellsize

      ! upscaling & packing
      if(cellFactorHbyM .gt. 1.0_dp) then
        call spatial_aggregation(L2_data, level2(iDomain)%cellsize, level1(iDomain)%cellsize, mask1, mask2, L1_data)
        ! downscaling
      elseif(cellFactorHbyM .lt. 1.0_dp) then
        call spatial_disaggregation(L2_data, level2(iDomain)%cellsize, level1(iDomain)%cellsize, mask1, mask2, L1_data)
        ! nothing
      else
        L1_data = L2_data
      end if
      ! free memory immediately
      deallocate(L2_data)

      ! pack variables
      nMonths = size(L1_data, dim = 3)
      nHours = size(L1_data, dim = 4)
      allocate(L1_data_packed(nCells1, nMonths, nHours))
      do t = 1, nMonths
        do j = 1, nHours
          L1_data_packed(:, t, j) = pack(L1_data(:, :, t, j), MASK = mask1(:, :))
        end do
      end do
      ! free memory immediately
      deallocate(L1_data)
    else
      ! dummy allocation
      allocate(L1_data_packed(nCells1, 12, 24))
      L1_data_packed = nodata_dp
    end if

    ! append
    call append(dataOut1, L1_data_packed)

    !free space
    deallocate(L1_data_packed)

  end subroutine meteo_weights_wrapper


  !> \brief determines the start date, end date, and read_flag given Domain id and current timestep
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - used by meteo-handler
  !!   - now independent of global variables
  !> \authors Stephan Thober
  !> \date Jun 2014
  subroutine chunk_config(iDomain, tt, nTstepDay, simPer, timestep, timeStep_model_inputs, read_flag, readPer)

    use mo_common_constants, only : nodata_dp
    use mo_common_types, only: period

    implicit none

    !> current Domain
    integer(i4), intent(in) :: iDomain
    !> current timestep
    integer(i4), intent(in) :: tt
    !> Number of time intervals per day
    integer(i4), intent(in) :: nTstepDay
    !> warmPer + evalPer
    type(period), dimension(:), intent(in) :: simPer
    !> [h] simulation time step (= TS) in [h]
    integer(i4), intent(in) :: timeStep
    !> frequency for reading meteo input
    integer(i4), dimension(:), intent(in) :: timeStep_model_inputs
    !> indicate whether reading data should be read
    logical, intent(out) :: read_flag
    ! start and end dates of reading Period
    type(period), intent(out) :: readPer

    ! initialize
    read_flag = .false.
    if (tt .eq. 1_i4) then
      readPer%julStart = int(nodata_dp, i4)
      readPer%julEnd = int(nodata_dp, i4)
      readPer%dstart = int(nodata_dp, i4)
      readPer%mstart = int(nodata_dp, i4)
      readPer%ystart = int(nodata_dp, i4)
      readPer%dend = int(nodata_dp, i4)
      readper%mend = int(nodata_dp, i4)
      readper%yend = int(nodata_dp, i4)
      readPer%Nobs = int(nodata_dp, i4)
    end if

    ! evaluate date and timeStep_model_inputs to get read_flag -------
    read_flag = is_read(iDomain, tt, nTstepDay, simPer, timestep, timeStep_model_inputs)
    !
    ! determine start and end date of chunk to read
    if (read_flag) call chunk_size(iDomain, tt, nTstepDay, simPer, timeStep_model_inputs, readPer)
    !
  end subroutine chunk_config


  !> \brief evaluate whether new chunk should be read at this timestep
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - used by meteo-handler
  !!   - now independent of global variables
  !> \authors Stephan Thober
  !> \date Jun 2014
  function is_read(iDomain, tt, nTstepDay, simPer, timestep, timeStep_model_inputs)

    use mo_common_types, only : period
    use mo_julian, only : caldat

    implicit none

    !> current Domain
    integer(i4), intent(in) :: iDomain
    !> current time step
    integer(i4), intent(in) :: tt
    !> Number of time intervals per day
    integer(i4), intent(in) :: nTstepDay
    !> warmPer + evalPer
    type(period), dimension(:), intent(in) :: simPer
    !> [h] simulation time step (= TS) in [h]
    integer(i4), intent(in) :: timeStep
    !> frequency for reading meteo input
    integer(i4), dimension(:), intent(in) :: timeStep_model_inputs

    logical :: is_read

    ! number of simulated days
    integer(i4) :: Ndays
    ! day
    integer(i4) :: day
    ! months
    integer(i4) :: month
    ! years
    integer(i4) :: year
    ! number of simulated days one timestep before
    integer(i4) :: Ndays_before
    ! day one simulated timestep before
    integer(i4) :: day_before
    ! month one simulated timestep before
    integer(i4) :: month_before
    ! year one simulated timestep before
    integer(i4) :: year_before

    ! initialize
    is_read = .false.

    ! special case for first timestep
    if (tt .eq. 1_i4) then
      is_read = .true.
    else
      ! check if a new day started by comparing the day of the current time step (Ndays)
      ! with the one before (Ndays_before)
      Ndays = ceiling((real(tt, dp)) / real(nTstepDay, dp))
      Ndays_before = ceiling((real(tt, dp) - 1.0_dp) / real(nTstepDay, dp))

      ! evaluate cases of given timeStep_model_inputs
      select case(timeStep_model_inputs(iDomain))
      case(0)  ! only at the beginning of the period
        if (tt .eq. 1_i4) is_read = .true.
      case(1 :) ! every timestep with frequency timeStep_model_inputs
        if (mod((tt - 1) * timestep, timeStep_model_inputs(iDomain) * 24) .eq. 0_i4) is_read = .true.
      case(-1) ! every day
        if (Ndays .ne. Ndays_before) is_read = .true.
      case(-2) ! every month
        if (Ndays .ne. Ndays_before) then
          ! calculate months
          call caldat(simPer(iDomain)%julStart + Ndays - 1, dd = day, mm = month, yy = year)
          call caldat(simPer(iDomain)%julStart + Ndays_before - 1, dd = day_before, mm = month_before, yy = year_before)
          if (month .ne. month_before) is_read = .true.
        end if
      case(-3) ! every year
        if (Ndays .ne. Ndays_before) then
          ! calculate months
          call caldat(simPer(iDomain)%julStart + Ndays - 1, dd = day, mm = month, yy = year)
          call caldat(simPer(iDomain)%julStart + Ndays_before - 1, dd = day_before, mm = month_before, yy = year_before)
          if (year .ne. year_before) is_read = .true.
        end if
      case default ! not specified correctly
        call error_message('ERROR*** mo_meteo_helper: function is_read: timStep_model_inputs not specified correctly!')
      end select
    end if

  end function is_read


  !> \brief calculate beginning and end of read Period, i.e. that is length of current chunk to read
  !> \changelog
  !! - Stephan Thober  Jan 2015
  !!   - added iDomain
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - used by meteo-handler
  !!   - now independent of global variables
  !> \authors Stephan Thober
  !> \date Jun 2014
  subroutine chunk_size(iDomain, tt, nTstepDay, simPer, timeStep_model_inputs, readPer)

    use mo_common_types, only: period
    use mo_julian, only : caldat, julday

    implicit none

    !> current Domain to process
    integer(i4), intent(in) :: iDomain
    !> current time step
    integer(i4), intent(in) :: tt
    !> Number of time intervals per day
    integer(i4), intent(in) :: nTstepDay
    !> warmPer + evalPer
    type(period), dimension(:), intent(in) :: simPer
    !> frequency for reading meteo input
    integer(i4), dimension(:), intent(in) :: timeStep_model_inputs
    !> start and end dates of read Period
    type(period), intent(out) :: readPer

    ! number of simulated days
    integer(i4) :: Ndays
    ! day
    integer(i4) :: day
    ! months
    integer(i4) :: month
    ! years
    integer(i4) :: year

    ! calculate date of start date
    Ndays = ceiling(real(tt, dp) / real(nTstepDay, dp))

    ! get start date
    readPer%julStart = simPer(iDomain)%julStart + Ndays - 1

    ! calculate end date according to specified frequency
    select case (timeStep_model_inputs(iDomain))
    case(0)  ! length of chunk has to cover whole period
      readPer%julEnd = simPer(iDomain)%julEnd
    case(1 :) ! every timestep with frequency timeStep_model_inputs
      readPer%julEnd = readPer%julStart + timeStep_model_inputs(iDomain) - 1
    case(-1) ! every day
      readPer%julEnd = readPer%julStart
    case(-2) ! every month
      ! calculate date
      call caldat(simPer(iDomain)%julStart + Ndays, dd = day, mm = month, yy = year)
      ! increment month
      if (month .eq. 12) then
        month = 1
        year = year + 1
      else
        month = month + 1
      end if
      readPer%julEnd = julday(dd = 1, mm = month, yy = year) - 1
    case(-3) ! every year
      ! calculate date
      call caldat(simPer(iDomain)%julStart + Ndays, dd = day, mm = month, yy = year)
      readPer%julEnd = julday(dd = 31, mm = 12, yy = year)
    case default ! not specified correctly
      call error_message('ERROR*** mo_meteo_helper: chunk_size: timStep_model_inputs not specified correctly!')
    end select

    ! end date should not be greater than end of simulation period
    readPer%julEnd = min(readPer%julEnd, simPer(iDomain)%julEnd)

    ! calculate the dates of the start and end dates
    call caldat(readPer%julStart, dd = readPer%dstart, mm = readPer%mstart, yy = readPer%ystart)
    call caldat(readPer%julEnd, dd = readPer%dEnd, mm = readPer%mend, yy = readPer%yend)
    readPer%Nobs = readPer%julEnd - readPer%julstart + 1

  end subroutine chunk_size

END MODULE mo_meteo_helper
