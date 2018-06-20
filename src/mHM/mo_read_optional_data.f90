!>       \file mo_read_optional_data.f90

!>       \brief Read optional data for mHM calibration.

!>       \details Data have to be provided in resolution of the hydrology.

!>       \authors Matthias Zink

!>       \date Mar 2015

! Modifications:

MODULE mo_read_optional_data

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_soil_moisture, read_basin_avg_TWS, read_neutrons, read_evapotranspiration

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_soil_moisture

  !    PURPOSE
  !>       \brief Read soil moisture data from NetCDF file for calibration

  !>       \details This routine reads oberved soil moisture fields which are used for model
  !>       calibration. The soil moisture file is expected to be called "sm.nc" with
  !>       a variable "sm" inside. The data are read only for the evaluation period
  !>       they are intended to be used for calibration. Soil moisture data are only
  !>       read if one of the corresponding objective functions is chosen.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin" Basin Id

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Mar 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_soil_moisture(iBasin)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : level1
    use mo_global_variables, only : L1_sm, L1_sm_mask, dirSoil_moisture, nTimeSteps_L1_sm, timeStep_sm_input
    use mo_message, only : message
    use mo_read_forcing_nc, only : read_forcing_nc
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, timer_stop

    implicit none

    ! Basin Id
    integer(i4), intent(in) :: iBasin

    ! loop  vars packing L1_data to L1_data_packed
    integer(i4) :: t

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! mask of level 1 for packing
    logical, dimension(:, :), allocatable :: mask1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! data at level-1
    real(dp), dimension(:, :, :), allocatable :: L1_data

    ! packed data at level-1 from 3D to 2D
    real(dp), dimension(:, :), allocatable :: L1_data_packed

    ! mask at level-1
    logical, dimension(:, :, :), allocatable :: L1_mask

    ! packed mask at level-1 from 3D to 2D
    logical, dimension(:, :), allocatable :: L1_mask_packed


    ! get basic basin information at level-1
    nrows1 = level1(iBasin)%nrows
    ncols1 = level1(iBasin)%ncols
    ncells1 = level1(iBasin)%ncells
    mask1 = level1(iBasin)%mask

    !  basin characteristics and read meteo header
    call message('  Reading soil moisture for basin:           ', trim(adjustl(num2str(iBasin))), ' ...')
    call timer_start(1)
    call read_forcing_nc(dirSoil_moisture(iBasin), nRows1, nCols1, 'sm', mask1, L1_data, &
            target_period = evalPer(iBasin), nctimestep = timeStep_sm_input, nocheck = .TRUE., maskout = L1_mask)
    ! pack variables
    nTimeSteps_L1_sm = size(L1_data, 3)
    allocate(L1_data_packed(nCells1, nTimeSteps_L1_sm))
    allocate(L1_mask_packed(nCells1, nTimeSteps_L1_sm))
    do t = 1, nTimeSteps_L1_sm
      L1_data_packed(:, t) = pack(L1_data(:, :, t), MASK = mask1(:, :))
      L1_mask_packed(:, t) = pack(L1_mask(:, :, t), MASK = mask1(:, :))
    end do

    ! append
    call append(L1_sm, L1_data_packed(:, :), fill_value = nodata_dp)
    call append(L1_sm_mask, L1_mask_packed(:, :), fill_value = .FALSE.)

    ! for multi basin calibration number of time steps may vary for different basins
    if (iBasin .GT. 1) nTimeSteps_L1_sm = size(L1_sm, 2)

    !free space
    deallocate(L1_data, L1_data_packed)

    call timer_stop(1)
    call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')

  end subroutine read_soil_moisture

  ! ---------------------------------------------------------------------------

  !    NAME
  !        read_basin_avg_TWS

  !    PURPOSE
  !>       \brief Read basin average TWS timeseries from file, the same way runoff is read

  !>       \details Read basin average TWS timeseries
  !>       Allocate global basin_avg_TWS variable that contains the simulated values after the simulation.

  !    HISTORY
  !>       \authors Oldrich Rakovec

  !>       \date Oct 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_basin_avg_TWS

    use mo_append, only : paste
    use mo_common_constants, only : nodata_dp
    use mo_common_mHM_mRM_variables, only : evalPer, nTstepDay, opti_function, optimize, simPer
    use mo_common_variables, only : nBasins
    use mo_file, only : utws
    use mo_global_variables, only : basin_avg_TWS_obs, basin_avg_TWS_sim, nMeasPerDay_TWS
    use mo_message, only : message
    use mo_read_timeseries, only : read_timeseries
    use mo_string_utils, only : num2str

    implicit none

    integer(i4) :: iBasin

    integer(i4) :: maxTimeSteps

    ! file name of file to read
    character(256) :: fName

    integer(i4), dimension(3) :: start_tmp, end_tmp

    real(dp), dimension(:), allocatable :: data_dp_1d

    logical, dimension(:), allocatable :: mask_1d


    ! ************************************************
    ! INITIALIZE TWS
    ! ************************************************
    maxTimeSteps = maxval(simPer(1 : nBasins)%julEnd - simPer(1 : nBasins)%julStart + 1) * nTstepDay
    allocate(basin_avg_TWS_sim(maxTimeSteps, nBasins))
    basin_avg_TWS_sim = nodata_dp

    ! ************************************************
    ! READ basin average TWS TIME SERIES
    ! ************************************************
    !
    do iBasin = 1, nBasins
      call message('  Reading basin average TWS for basin:     ', trim(adjustl(num2str(iBasin))), ' ...')

      ! get start and end dates
      start_tmp = (/evalPer(iBasin)%yStart, evalPer(iBasin)%mStart, evalPer(iBasin)%dStart/)
      end_tmp = (/evalPer(iBasin)%yEnd, evalPer(iBasin)%mEnd, evalPer(iBasin)%dEnd  /)
      fName = trim(adjustl(basin_avg_TWS_obs%fname(iBasin)))
      call read_timeseries(trim(fName), utws, &
              start_tmp, end_tmp, optimize, opti_function, &
              data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay_TWS)
      data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
      call paste(basin_avg_TWS_obs%TWS, data_dp_1d, nodata_dp)
      deallocate (data_dp_1d)
    end do

  end subroutine read_basin_avg_TWS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_neutrons

  !    PURPOSE
  !>       \brief Read neutrons data from NetCDF file for calibration

  !>       \details This routine reads oberved neutron fields which are used for model
  !>       calibration. The neutrons file is expected to be called "neutrons.nc" with
  !>       a variable "neutrons" inside. The data are read only for the evaluation period
  !>       they are intended to be used for calibration. Neutrons data are only
  !>       read if one of the corresponding objective functions is chosen.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin" Basin Id

  !    HISTORY
  !>       \authors Martin Schroen

  !>       \date Jul 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_neutrons(iBasin)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : level1
    use mo_global_variables, only : L1_neutronsdata, L1_neutronsdata_mask, dirNeutrons, nTimeSteps_L1_neutrons, timeStep_neutrons_input
    use mo_message, only : message
    use mo_read_forcing_nc, only : read_forcing_nc
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, &
                         timer_start, timer_stop

    implicit none

    ! Basin Id
    integer(i4), intent(in) :: iBasin

    ! loop  vars packing L1_data to L1_data_packed
    integer(i4) :: t

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! mask of level 1 for packing
    logical, dimension(:, :), allocatable :: mask1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! data at level-1
    real(dp), dimension(:, :, :), allocatable :: L1_data

    ! packed data at level-1 from 3D to 2D
    real(dp), dimension(:, :), allocatable :: L1_data_packed

    ! mask at level-1
    logical, dimension(:, :, :), allocatable :: L1_mask

    ! packed mask at level-1 from 3D to 2D
    logical, dimension(:, :), allocatable :: L1_mask_packed


    ! get basic basin information at level-1
    nrows1 = level1(iBasin)%nrows
    ncols1 = level1(iBasin)%ncols
    ncells1 = level1(iBasin)%ncells
    mask1 = level1(iBasin)%mask

    !  basin characteristics and read meteo header
    call message('  Reading neutrons for basin:           ', trim(adjustl(num2str(iBasin))), ' ...')
    call timer_start(1)
    call read_forcing_nc(dirNeutrons(iBasin), nRows1, nCols1, 'neutrons', mask1, L1_data, &
            target_period = evalPer(iBasin), nctimestep = timeStep_neutrons_input, nocheck = .TRUE., maskout = L1_mask)
    ! pack variables
    nTimeSteps_L1_neutrons = size(L1_data, 3)
    allocate(L1_data_packed(nCells1, nTimeSteps_L1_neutrons))
    allocate(L1_mask_packed(nCells1, nTimeSteps_L1_neutrons))
    do t = 1, nTimeSteps_L1_neutrons
      L1_data_packed(:, t) = pack(L1_data(:, :, t), MASK = mask1(:, :))
      L1_mask_packed(:, t) = pack(L1_mask(:, :, t), MASK = mask1(:, :))
    end do

    ! append
    call append(L1_neutronsdata, L1_data_packed(:, :), fill_value = nodata_dp)
    call append(L1_neutronsdata_mask, L1_mask_packed(:, :), fill_value = .FALSE.)

    ! for multi basin calibration number of time steps may vary for different basins
    if (iBasin .GT. 1) nTimeSteps_L1_neutrons = size(L1_neutronsdata, 2)

    !free space
    deallocate(L1_data, L1_data_packed)

    call timer_stop(1)
    call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')

  end subroutine read_neutrons

  ! ------------------------------------------------------------------

  !    NAME
  !        read_evapotranspiration

  !    PURPOSE
  !>       \brief Read evapotranspiration data from NetCDF file for calibration

  !>       \details This routine reads oberved evapotranspiration fields which are used for model
  !>       calibration. The evapotranspiration file is expected to be called "et.nc" with
  !>       a variable "et" inside. The data are read only for the evaluation period
  !>       they are intended to be used for calibration. Evapotranspiration data are only
  !>       read if one of the corresponding objective functions is chosen.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin" Basin Id

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date Feb 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_evapotranspiration(iBasin)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : level1
    use mo_global_variables, only : L1_et, L1_et_mask, dirEvapotranspiration, nTimeSteps_L1_et, timeStep_et_input
    use mo_message, only : message
    use mo_read_forcing_nc, only : read_forcing_nc
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, &
                         timer_stop

    implicit none

    ! Basin Id
    integer(i4), intent(in) :: iBasin

    ! loop  vars packing L1_data to L1_data_packed
    integer(i4) :: t

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! mask of level 1 for packing
    logical, dimension(:, :), allocatable :: mask1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! data at level-1
    real(dp), dimension(:, :, :), allocatable :: L1_data

    ! packed data at level-1 from 3D to 2D
    real(dp), dimension(:, :), allocatable :: L1_data_packed

    ! mask at level-1
    logical, dimension(:, :, :), allocatable :: L1_mask

    ! packed mask at level-1 from 3D to 2D
    logical, dimension(:, :), allocatable :: L1_mask_packed


    ! get basic basin information at level-1
    nrows1 = level1(iBasin)%nrows
    ncols1 = level1(iBasin)%ncols
    ncells1 = level1(iBasin)%ncells
    mask1 = level1(iBasin)%mask

    !  basin characteristics and read meteo header
    call message('  Reading evapotranspiration for basin:           ', trim(adjustl(num2str(iBasin))), ' ...')
    call timer_start(1)
    call read_forcing_nc(dirEvapotranspiration(iBasin), nRows1, nCols1, 'et', mask1, L1_data, &
            target_period = evalPer(iBasin), nctimestep = timeStep_et_input, nocheck = .TRUE., maskout = L1_mask)

    ! pack variables
    nTimeSteps_L1_et = size(L1_data, 3)
    allocate(L1_data_packed(nCells1, nTimeSteps_L1_et))
    allocate(L1_mask_packed(nCells1, nTimeSteps_L1_et))
    do t = 1, nTimeSteps_L1_et
      L1_data_packed(:, t) = pack(L1_data(:, :, t), MASK = mask1(:, :))
      L1_mask_packed(:, t) = pack(L1_mask(:, :, t), MASK = mask1(:, :))
    end do

    ! append
    call append(L1_et, L1_data_packed(:, :), fill_value = nodata_dp)
    call append(L1_et_mask, L1_mask_packed(:, :), fill_value = .FALSE.)

    ! for multi basin calibration number of time steps may vary for different basins
    if (iBasin .GT. 1) nTimeSteps_L1_et = size(L1_et, 2)

    !free space
    deallocate(L1_data, L1_data_packed)

    call timer_stop(1)
    call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')

  end subroutine read_evapotranspiration

END MODULE mo_read_optional_data
