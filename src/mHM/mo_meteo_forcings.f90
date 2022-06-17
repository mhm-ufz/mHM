!>       \file mo_meteo_forcings.f90

!>       \brief Prepare meteorological forcings data for mHM.

!>       \details Prepare meteorological forcings data for mHM.

!>       \authors Rohini Kumar

!>       \date Jan 2012

! Modifications:

MODULE mo_meteo_forcings

  ! This module provides routines to read meteorological data.

  ! Written  Rohini Kumar, Jan 2013
  ! Modified

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prepare_meteo_forcings_data

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        prepare_meteo_forcings_data

  !    PURPOSE
  !>       \brief Prepare meteorological forcings data for a given variable

  !>       \details Prepare meteorological forcings data for a given variable.
  !>       Internally this subroutine calls another routine meteo_wrapper
  !>       for different meterological variables

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id
  !>       \param[in] "integer(i4) :: tt"     current timestep

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Jan 2013

  ! Modifications:
  ! Matthias Zink,   Jun 2013 - addded NetCDf reader
  ! Rohini Kumar,    Aug 2013 - name changed "inputFormat" to inputFormat_meteo_forcings
  ! Matthias Zink,   Feb 2014 - added read in for different PET processes (process 5)
  ! Stephan Thober,  Jun 2014 - add chunk_config for chunk read,
  !                             copied L2 initialization to mo_startup
  ! Stephan Thober,  Nov 2016 - moved processMatrix to common variables
  ! Stephan Thober,  Jan 2017 - added subroutine for meteo_weights
  ! Robert Schweppe  Jun 2018 - refactoring and reformatting

  subroutine prepare_meteo_forcings_data(iDomain, domainID, tt)

    use mo_common_mhm_mrm_variables, only : readPer
    use mo_common_variables, only : domainMeta, processMatrix
    use mo_global_variables, only : L1_absvappress, L1_netrad, L1_pet, L1_pet_weights, L1_pre, L1_pre_weights, L1_temp, &
                                    L1_temp_weights, L1_tmax, L1_tmin, L1_windspeed, dirMaxTemperature, &
                                    dirMinTemperature, dirNetRadiation, dirPrecipitation, dirReferenceET, dirTemperature, &
                                    dirabsVapPressure, dirwindspeed, dirRadiation, &
                                    inputFormat_meteo_forcings, read_meteo_weights, &
                                    timeStep_model_inputs, &
                                    L1_ssrd, L1_strd, L1_tann  ! riv-temp related
    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, timer_stop

    implicit none

    ! Domain number
    integer(i4), intent(in) :: iDomain

    ! Domain ID
    integer(i4), intent(in) :: domainID

    ! current timestep
    integer(i4), intent(in) :: tt

    ! indicate whether data should be read
    logical :: read_flag


    ! configuration of chunk_read
    call chunk_config(iDomain, tt, read_flag, readPer)

    ! only read, if read_flag is true
    if (read_flag) then

      ! read weights for hourly disaggregation of temperature
      if (tt .eq. 1) then
        ! TODO-RIV-TEMP: No NC files for weights for radiation at the moment
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for tavg     ...')
        call meteo_weights_wrapper(iDomain, read_meteo_weights, dirTemperature(iDomain), &
                L1_temp_weights, ncvarName = 'tavg_weight')

        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for pet     ...')
        call meteo_weights_wrapper(iDomain, read_meteo_weights, dirReferenceET(iDomain), &
                L1_pet_weights, ncvarName = 'pet_weight')

        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for pre     ...')
        call meteo_weights_wrapper(iDomain, read_meteo_weights, dirPrecipitation(iDomain), &
                L1_pre_weights, ncvarName = 'pre_weight')
      end if

      ! free L1 variables if chunk read is activated
      if (timeStep_model_inputs(iDomain) .ne. 0) then
        if (allocated(L1_pre)) deallocate(L1_pre)
        if (allocated(L1_temp)) deallocate(L1_temp)
        if (allocated(L1_pet)) deallocate(L1_pet)
        if (allocated(L1_tmin)) deallocate(L1_tmin)
        if (allocated(L1_tmax)) deallocate(L1_tmax)
        if (allocated(L1_netrad)) deallocate(L1_netrad)
        if (allocated(L1_absvappress)) deallocate(L1_absvappress)
        if (allocated(L1_windspeed)) deallocate(L1_windspeed)
      end if

      !  Domain characteristics and read meteo header
      if (timeStep_model_inputs(iDomain) .eq. 0) then
        call message('  Reading meteorological forcings for Domain: ', trim(adjustl(num2str(domainID))), ' ...')
        call timer_start(1)
      end if

      ! precipitation
      if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read precipitation        ...')
      call meteo_forcings_wrapper(iDomain, dirPrecipitation(iDomain), inputFormat_meteo_forcings, &
              L1_pre, lower = 0.0_dp, upper = 1000._dp, ncvarName = 'pre')

      ! temperature
      if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read temperature          ...')
      call meteo_forcings_wrapper(iDomain, dirTemperature(iDomain), inputFormat_meteo_forcings, &
              L1_temp, lower = -100._dp, upper = 100._dp, ncvarName = 'tavg')

      ! read input for PET (process 5) depending on specified option
      ! 0 - input, 1 - Hargreaves-Samani, 2 - Priestley-Taylor, 3 - Penman-Monteith
      select case (processMatrix(5, 1))

      case(-1 : 0) ! pet is input
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read pet                  ...')
        call meteo_forcings_wrapper(iDomain, dirReferenceET(iDomain), inputFormat_meteo_forcings, &
                L1_pet, lower = 0.0_dp, upper = 1000._dp, ncvarName = 'pet')
        ! allocate PET and dummies for mhm_call
        if ((iDomain.eq.domainMeta%nDomains) .OR. (timeStep_model_inputs(iDomain) .NE. 0)) then
          allocate(L1_tmin(1, 1)); allocate(L1_tmax(1, 1)); allocate(L1_netrad(1, 1))
          allocate(L1_absvappress(1, 1)); allocate(L1_windspeed(1, 1))
        end if

      case(1) ! Hargreaves-Samani formulation (input: minimum and maximum Temperature)
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read min. temperature     ...')
        call meteo_forcings_wrapper(iDomain, dirMinTemperature(iDomain), inputFormat_meteo_forcings, &
                L1_tmin, lower = -100.0_dp, upper = 100._dp, ncvarName = 'tmin')
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read max. temperature     ...')
        call meteo_forcings_wrapper(iDomain, dirMaxTemperature(iDomain), inputFormat_meteo_forcings, &
                L1_tmax, lower = -100.0_dp, upper = 100._dp, ncvarName = 'tmax')
        ! allocate PET and dummies for mhm_call
        if ((iDomain .eq. domainMeta%nDomains) .OR. (timeStep_model_inputs(iDomain) .NE. 0)) then
          allocate(L1_pet    (size(L1_tmax, dim = 1), size(L1_tmax, dim = 2)))
          allocate(L1_netrad(1, 1)); allocate(L1_absvappress(1, 1)); allocate(L1_windspeed(1, 1))
        end if

      case(2) ! Priestley-Taylor formulation (input: net radiation)
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read net radiation        ...')
        call meteo_forcings_wrapper(iDomain, dirNetRadiation(iDomain), inputFormat_meteo_forcings, &
                L1_netrad, lower = -500.0_dp, upper = 1500._dp, ncvarName = 'net_rad')
        ! allocate PET and dummies for mhm_call
        if ((iDomain .eq. domainMeta%nDomains) .OR. (timeStep_model_inputs(iDomain) .NE. 0)) then
          allocate(L1_pet    (size(L1_netrad, dim = 1), size(L1_netrad, dim = 2)))
          allocate(L1_tmin(1, 1)); allocate(L1_tmax(1, 1))
          allocate(L1_absvappress(1, 1)); allocate(L1_windspeed(1, 1))
        end if

      case(3) ! Penman-Monteith formulation (input: net radiationm absulute vapour pressure, windspeed)
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read net radiation        ...')
        call meteo_forcings_wrapper(iDomain, dirNetRadiation(iDomain), inputFormat_meteo_forcings, &
                L1_netrad, lower = -500.0_dp, upper = 1500._dp, ncvarName = 'net_rad')
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read absolute vapour pressure  ...')
        call meteo_forcings_wrapper(iDomain, dirabsVapPressure(iDomain), inputFormat_meteo_forcings, &
                L1_absvappress, lower = 0.0_dp, upper = 15000.0_dp, ncvarName = 'eabs')
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read windspeed            ...')
        call meteo_forcings_wrapper(iDomain, dirwindspeed(iDomain), inputFormat_meteo_forcings, &
                L1_windspeed, lower = 0.0_dp, upper = 250.0_dp, ncvarName = 'windspeed')
        ! allocate PET and dummies for mhm_call
        if ((iDomain.eq.domainMeta%nDomains) .OR. (timeStep_model_inputs(iDomain) .NE. 0)) then
          allocate(L1_pet    (size(L1_absvappress, dim = 1), size(L1_absvappress, dim = 2)))
          allocate(L1_tmin(1, 1)); allocate(L1_tmax(1, 1))
        end if
      end select

      ! long/short-wave radiation and annual mean temperature for river-temperature routing
      if ( processMatrix(11, 1) .ne. 0 ) then
        ! free L1 variables if chunk read is activated
        if (timeStep_model_inputs(iDomain) .ne. 0) then
          if (allocated(L1_ssrd)) deallocate(L1_ssrd)
          if (allocated(L1_strd)) deallocate(L1_strd)
          if (allocated(L1_tann)) deallocate(L1_tann)
        end if
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read short-wave radiation ...')
          call meteo_forcings_wrapper( &
                iDomain, dirRadiation(iDomain), inputFormat_meteo_forcings, &
                L1_ssrd, lower = 0.0_dp, upper = 1500._dp, ncvarName = 'ssrd')
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read long-wave radiation ...')
          call meteo_forcings_wrapper( &
                iDomain, dirRadiation(iDomain), inputFormat_meteo_forcings, &
                L1_strd, lower = 0.0_dp, upper = 1500._dp, ncvarName = 'strd')
        if (timeStep_model_inputs(iDomain) .eq. 0) call message('    read annual mean temperature ...')
          call meteo_forcings_wrapper( &
                iDomain, dirTemperature(iDomain), inputFormat_meteo_forcings, &
                L1_tann, lower = -100.0_dp, upper = 100._dp, ncvarName = 'tann')
      end if

      if (timeStep_model_inputs(iDomain) .eq. 0) then
        call timer_stop(1)
        call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')
      end if
    end if

  end subroutine prepare_meteo_forcings_data


  ! ------------------------------------------------------------------

  !    NAME
  !        meteo_forcings_wrapper

  !    PURPOSE
  !>       \brief Prepare meteorological forcings data for mHM at Level-1

  !>       \details Prepare meteorological forcings data for mHM, which include
  !>       1) Reading meteo. datasets at their native resolution for every Domain
  !>       2) Perform aggregation or disaggregation of meteo. datasets from their
  !>       native resolution (level-2) to the required hydrologic resolution (level-1)
  !>       3) Pad the above datasets of every Domain to their respective global ones

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"             Domain Id
  !>       \param[in] "character(len = *) :: dataPath"    Data path where a given meteo. variable is stored
  !>       \param[in] "character(len = *) :: inputFormat" only 'nc' possible at the moment

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :) :: dataOut1" Packed meterological variable for the whole simulation
  !>       period

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: lower"               Lower bound for check of validity of data values
  !>       \param[in] "real(dp), optional :: upper"               Upper bound for check of validity of data values
  !>       \param[in] "character(len = *), optional :: ncvarName" name of the variable (for .nc files)

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Jan 2013

  ! Modifications:
  ! Stephan Thober Jun 2014 - changed to readPer
  ! Stephan Thober Feb 2016 - refactored deallocate statements
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine meteo_forcings_wrapper(iDomain, dataPath, inputFormat, dataOut1, lower, upper, ncvarName)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_mhm_mrm_variables, only : readPer
    use mo_common_variables, only : level1
    use mo_global_variables, only : level2
    use mo_read_nc, only : read_nc
    use mo_spatial_agg_disagg_forcing, only : spatial_aggregation, spatial_disaggregation

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    ! Data path where a given meteo. variable is stored
    character(len = *), intent(in) :: dataPath

    ! only 'nc' possible at the moment
    character(len = *), intent(in) :: inputFormat

    ! Packed meterological variable for the whole simulation period
    real(dp), dimension(:, :), allocatable, intent(inout) :: dataOut1

    ! Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower

    ! Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper

    ! name of the variable (for .nc files)
    character(len = *), optional, intent(in) :: ncvarName

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
      if(present(lower) .AND. (.not. present(upper))) then
        CALL read_nc(dataPath, nRows2, nCols2, ncvarName, mask2, L2_data, target_period = readPer, &
                lower = lower, is_meteo=.True.)
      end if
      !
      if(present(upper) .AND. (.not. present(lower))) then
         CALL read_nc(dataPath, nRows2, nCols2, ncvarName, mask2, L2_data, target_period = readPer, &
              upper = upper, is_meteo=.True.)
      end if
      !
      if(present(lower) .AND. present(upper)) then
         CALL read_nc(dataPath, nRows2, nCols2, ncvarName, mask2, L2_data, target_period = readPer, &
              lower = lower, upper = upper, is_meteo=.True.)
      end if
      !
      if((.not. present(lower)) .AND. (.not. present(upper))) then
         CALL read_nc(dataPath, nRows2, nCols2, ncvarName, mask2, L2_data, target_period = readPer, is_meteo=.True.)
      end if
    case DEFAULT
      stop '***ERROR: meteo_forcings_wrapper: Not recognized input format'
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


  ! ------------------------------------------------------------------
  !    NAME
  !        meteo_weights_wrapper

  !    PURPOSE
  !>       \brief Prepare weights for meteorological forcings data for mHM at Level-1

  !>       \details Prepare meteorological weights data for mHM, which include
  !>       1) Reading meteo. weights datasets at their native resolution for every Domain
  !>       2) Perform aggregation or disaggregation of meteo. weights datasets from their
  !>       native resolution (level-2) to the required hydrologic resolution (level-1)
  !>       3) Pad the above datasets of every Domain to their respective global ones

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"          Domain Id
  !>       \param[in] "logical :: read_meteo_weights"  Flag for reading meteo weights
  !>       \param[in] "character(len = *) :: dataPath" Data path where a given meteo. variable is stored

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :, :) :: dataOut1" Packed meterological variable for the whole
  !>       simulation period

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: lower"               Lower bound for check of validity of data values
  !>       \param[in] "real(dp), optional :: upper"               Upper bound for check of validity of data values
  !>       \param[in] "character(len = *), optional :: ncvarName" name of the variable (for .nc files)

  !    HISTORY
  !>       \authors Stephan Thober & Rohini Kumar

  !>       \date Jan 2017

  ! Modifications:
  ! Stephan Thober May 2017 - updated documentation
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Oldrich Rakovec Aug 2018 - adding message about reading meteo_weights

  subroutine meteo_weights_wrapper(iDomain, read_meteo_weights, dataPath, dataOut1, lower, upper, ncvarName)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : level1
    use mo_global_variables, only : level2
    use mo_read_nc, only : read_weights_nc
    use mo_spatial_agg_disagg_forcing, only : spatial_aggregation, spatial_disaggregation
    use mo_message, only : message

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    ! Flag for reading meteo weights
    logical, intent(in) :: read_meteo_weights

    ! Data path where a given meteo. variable is stored
    character(len = *), intent(in) :: dataPath

    ! Packed meterological variable for the whole simulation period
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: dataOut1

    ! Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower

    ! Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper

    ! name of the variable (for .nc files)
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
      if(present(lower) .AND. (.not. present(upper))) then
        CALL read_weights_nc(dataPath, nRows2, nCols2, ncvarName, L2_data, mask2, lower = lower)
      end if
      !
      if(present(upper) .AND. (.not. present(lower))) then
        CALL read_weights_nc(dataPath, nRows2, nCols2, ncvarName, L2_data, mask2, upper = upper)
      end if
      !
      if(present(lower) .AND. present(upper)) then
        CALL read_weights_nc(dataPath, nRows2, nCols2, ncvarName, L2_data, mask2, lower = lower, upper = upper)
      end if

      if((.not. present(lower)) .AND. (.not. present(upper))) then
        CALL read_weights_nc(dataPath, nRows2, nCols2, ncvarName, L2_data, mask2)
      end if

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


  ! ------------------------------------------------------------------
  !    NAME
  !        chunk_config

  !    PURPOSE
  !>       \brief determines the start date, end date, and read_flag given Domain id and current timestep

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" current Domain
  !>       \param[in] "integer(i4) :: tt"     current timestep

  !    INTENT(OUT)
  !>       \param[out] "logical :: read_flag"    indicate whether reading data should be read
  !>       \param[out] "type(period) :: readPer" start and end dates of reading Period

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Jun 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine chunk_config(iDomain, tt, read_flag, readPer)

    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : period
    use mo_kind, only : i4

    implicit none

    ! current Domain
    integer(i4), intent(in) :: iDomain

    ! current timestep
    integer(i4), intent(in) :: tt

    ! indicate whether reading data should be read
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
    read_flag = is_read(iDomain, tt)
    !
    ! determine start and end date of chunk to read
    if (read_flag) call chunk_size(iDomain, tt, readPer)
    !
  end subroutine chunk_config
  ! ------------------------------------------------------------------

  !    NAME
  !        is_read

  !    PURPOSE
  !>       \brief evaluate whether new chunk should be read at this timestep

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" current Domain
  !>       \param[in] "integer(i4) :: tt"     current time step

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Jun 2014

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function is_read(iDomain, tt)

    use mo_common_mhm_mrm_variables, only : nTstepDay, simPer, timestep
    use mo_global_variables, only : timeStep_model_inputs
    use mo_julian, only : caldat
    use mo_kind, only : i4
    use mo_message, only : message

    implicit none

    ! current Domain
    integer(i4), intent(in) :: iDomain

    ! current time step
    integer(i4), intent(in) :: tt

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
        call message('ERROR*** mo_meteo_forcings: function is_read: timStep_model_inputs not specified correctly!')
        stop
      end select
    end if

  end function is_read
  ! ------------------------------------------------------------------

  !    NAME
  !        chunk_size

  !    PURPOSE
  !>       \brief calculate beginning and end of read Period, i.e. that
  !>       is length of current chunk to read

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" current Domain to process
  !>       \param[in] "integer(i4) :: tt"     current time step

  !    INTENT(OUT)
  !>       \param[out] "type(period) :: readPer" start and end dates of read Period

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Jun 2014

  ! Modifications:
  ! Stephan Thober  Jan 2015 - added iDomain
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine chunk_size(iDomain, tt, readPer)

    use mo_common_mhm_mrm_variables, only : nTstepDay, simPer
    use mo_common_variables, only : period
    use mo_global_variables, only : timeStep_model_inputs
    use mo_julian, only : caldat, julday
    use mo_kind, only : i4
    use mo_message, only : message

    implicit none

    ! current time step
    integer(i4), intent(in) :: tt

    ! current Domain to process
    integer(i4), intent(in) :: iDomain

    ! start and end dates of read Period
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
      call message('ERROR*** mo_meteo_forcings: chunk_size: timStep_model_inputs not specified correctly!')
      stop
    end select

    ! end date should not be greater than end of simulation period
    readPer%julEnd = min(readPer%julEnd, simPer(iDomain)%julEnd)

    ! calculate the dates of the start and end dates
    call caldat(readPer%julStart, dd = readPer%dstart, mm = readPer%mstart, yy = readPer%ystart)
    call caldat(readPer%julEnd, dd = readPer%dEnd, mm = readPer%mend, yy = readPer%yend)
    readPer%Nobs = readPer%julEnd - readPer%julstart + 1

  end subroutine chunk_size
  !
END MODULE mo_meteo_forcings
