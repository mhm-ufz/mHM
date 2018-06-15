!> \file mo_mrm_read_data.f90

!> \brief This module contains all routines to read mRM data from file.

!> \details 

!> \authors Stephan Thober
!> \date Aug 2015

module mo_mrm_read_data
  use mo_kind, only : i4, dp
  implicit none
  public :: mrm_read_L0_data
  public :: mrm_read_discharge
  public :: mrm_read_total_runoff
  private
contains

  ! ------------------------------------------------------------------

  !     NAME
  !         mrm_read_L0_data

  !     PURPOSE
  !>        \brief read L0 data from file
  !
  !>        \details With the exception of L0_mask, L0_elev, and L0_LCover, all
  !>        L0 variables are read from file. The former three are only read if they
  !>        are not provided as variables.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(:), target, optional :: L0_mask - L0 mask"
  !>        \param[in] "real(dp), dimension(:), target, optional :: L0_elev - L0 elevation"
  !>        \param[in] "integer(i4), dimension(:,:), target, optional :: L0_LCover - L0 land cover"
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>        None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Juliane Mai, Matthias Zink, and Stephan Thober
  !>        \date Aug 2015
  !         Modified, Sep 2015 - Stephan Thober, added L0_mask, L0_elev, and L0_LCover
  subroutine mrm_read_L0_data(do_reinit, do_readlatlon, do_readlcover)
    use mo_common_constants, only : nodata_i4! mRM's global nodata vales
    use mo_append, only : append
    use mo_string_utils, only : num2str
    use mo_message, only : message
    use mo_read_spatial_data, only : read_spatial_data_ascii
    use mo_mrm_file, only : &
            file_facc, ufacc, & ! file name and unit of flow acc map
            file_fdir, ufdir, & ! file name and unit of flow dir map
            file_gaugeloc, ugaugeloc ! file name and unit of gauge locations m
    use mo_mrm_global_variables, only : &
            L0_fAcc, & ! flow accumulation on input resolution (L0)
            L0_fDir, & ! flow direction on input resolution (L0)
            L0_gaugeLoc, & ! location of evaluation gauges on input resolution (L0)
            L0_InflowGaugeLoc, & ! location of inflow gauges on input resolution (L0)
            basin_mrm ! basin information for single basins
    use mo_common_variables, only : &
            L0_LCover, &
            Grid, &
            level0, & ! level0 information
            dirMorpho, & ! directories
            L0_Basin, &
            nBasins, &
            processMatrix ! process description
    use mo_common_read_data, only : read_dem, read_lcover
    USE mo_read_latlon, ONLY : read_latlon

    implicit none

    ! optional input variables
    logical, intent(in) :: do_reinit
    logical, intent(in) :: do_readlatlon
    logical, intent(in) :: do_readlcover

    ! local variables
    integer(i4) :: iBasin
    integer(i4) :: iVar
    integer(i4) :: iGauge
    character(256) :: fname
    integer(i4) :: nunit
    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4
    logical, dimension(:, :), allocatable :: mask_2d
    logical, dimension(:, :), allocatable :: mask_global
    type(Grid), pointer :: level0_iBasin


    ! ************************************************
    ! READ SPATIAL DATA FOR EACH BASIN
    ! ************************************************

    if (do_reinit) then
      call read_dem()
    end if

    if (do_readlcover .and. processMatrix(8, 1) .eq. 1) then
      call read_lcover()
    else if (do_readlcover .and. processMatrix(8, 1) .eq. 2) then
      allocate(dataMatrix_i4(count(mask_global), 1))
      dataMatrix_i4 = nodata_i4
      call append(L0_LCover, dataMatrix_i4)
      ! free memory
      deallocate(dataMatrix_i4)
    end if

    do iBasin = 1, nBasins

      level0_iBasin => level0(L0_Basin(iBasin))

      ! check whether L0 data is shared
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          !
          call message('      Using data of basin ', &
                  trim(adjustl(num2str(L0_Basin(iBasin)))), ' for basin: ',&
                  trim(adjustl(num2str(iBasin))), '...')
          cycle
          !
        end if
      end if
      !
      call message('      Reading data for basin: ', trim(adjustl(num2str(iBasin))), ' ...')

      if (do_readlatlon) then
        ! read lat lon coordinates of each basin
        call read_latlon(iBasin, "lon_l0", "lat_l0", "level0", level0_iBasin)
      end if

      ! read fAcc, fDir, gaugeLoc
      do iVar = 1, 3
        select case (iVar)
        case(1) ! flow accumulation
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_facc))
          nunit = ufacc
        case(2) ! flow direction
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_fdir))
          nunit = ufdir
        case(3) ! location of gauging stations
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc))
          nunit = ugaugeloc
        end select
        !
        ! reading and transposing
        call read_spatial_data_ascii(trim(fName), nunit, &
                level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
                level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)

        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        ! put data into global L0 variable
        select case (iVar)
        case(1) ! flow accumulation
          call append(L0_fAcc, pack(data_i4_2d, level0_iBasin%mask))
        case(2) ! flow direction
          ! rotate flow direction and any other variable with directions
          ! NOTE: ONLY when ASCII files are read
          call rotate_fdir_variable(data_i4_2d)
          ! append
          call append(L0_fDir, pack(data_i4_2d, level0_iBasin%mask))
        case(3) ! location of evaluation and inflow gauging stations
          ! evaluation gauges
          ! Input data check
          do iGauge = 1, basin_mrm(iBasin)%nGauges
            ! If gaugeId is found in gauging location file?
            if (.not. any(data_i4_2d .EQ. basin_mrm(iBasin)%gaugeIdList(iGauge))) then
              call message()
              call message('***ERROR: Gauge ID "', trim(adjustl(num2str(basin_mrm(iBasin)%gaugeIdList(iGauge)))), &
                      '" not found in ')
              call message('          Gauge location input file: ', &
                      trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc)))
              stop
            end if
          end do

          call append(L0_gaugeLoc, pack(data_i4_2d, level0_iBasin%mask))

          ! inflow gauges
          ! if no inflow gauge for this subbasin exists still matirx with dim of subbasin has to be paded
          if (basin_mrm(iBasin)%nInflowGauges .GT. 0_i4) then
            ! Input data check
            do iGauge = 1, basin_mrm(iBasin)%nInflowGauges
              ! If InflowGaugeId is found in gauging location file?
              if (.not. any(data_i4_2d .EQ. basin_mrm(iBasin)%InflowGaugeIdList(iGauge))) then
                call message()
                call message('***ERROR: Inflow Gauge ID "', &
                        trim(adjustl(num2str(basin_mrm(iBasin)%InflowGaugeIdList(iGauge)))), &
                        '" not found in ')
                call message('          Gauge location input file: ', &
                        trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_gaugeloc)))
                stop 1
              end if
            end do
          end if

          call append(L0_InflowGaugeLoc, pack(data_i4_2d, level0_iBasin%mask))

        end select
        !
        ! deallocate arrays
        deallocate(data_i4_2d, mask_2d)
        !
      end do
    end do

  end subroutine mrm_read_L0_data
  ! ---------------------------------------------------------------------------

  !      NAME
  !          mrm_read_discharge

  !>        \brief Read discharge timeseries from file

  !>        \details Read Observed discharge at the outlet of a catchment
  !>        and at the inflow of a catchment. Allocate global runoff
  !>        variable that contains the simulated runoff after the simulation.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: iBasin"  basin id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author  Matthias Zink & Stephan Thober
  !         \date    Aug 2015
  subroutine mrm_read_discharge()
    use mo_message, only : message
    use mo_append, only : paste
    use mo_string_utils, only : num2str
    use mo_read_timeseries, only : read_timeseries
    use mo_mrm_file, only : udischarge
    use mo_common_constants, only : nodata_dp
    use mo_mrm_global_variables, only : &
            mRM_runoff, & ! variable storing runoff for each gauge
            nGaugesTotal, gauge, nMeasPerDay, & ! evaluaton gauging station information
            nInflowGaugesTotal, InflowGauge ! inflow stations information
    use mo_common_variables, only : &
            nBasins
    use mo_common_mHM_mRM_variables, only : &
            optimize, & ! optimizeation flag for some error checks
            opti_function, &   ! opti_function that determines to what data to calibrate
            evalPer, & ! model evaluation period (for discharge read in)
            simPer, & ! model simulation period (for inflow read in)
            nTstepDay

    implicit none

    ! local variables
    integer(i4) :: iGauge
    integer(i4) :: iBasin
    integer(i4) :: maxTimeSteps
    character(256) :: fName ! file name of file to read
    integer(i4), dimension(3) :: start_tmp, end_tmp
    real(dp), dimension(:), allocatable :: data_dp_1d
    logical, dimension(:), allocatable :: mask_1d

    !----------------------------------------------------------
    ! INITIALIZE RUNOFF
    !----------------------------------------------------------
    maxTimeSteps = maxval(simPer(1 : nBasins)%julEnd - simPer(1 : nBasins)%julStart + 1) * nTstepDay
    allocate(mRM_runoff(maxTimeSteps, nGaugesTotal))
    mRM_runoff = nodata_dp

    ! READ GAUGE DATA
    do iGauge = 1, nGaugesTotal
      ! get basin id
      iBasin = gauge%basinId(iGauge)
      ! get start and end dates
      start_tmp = (/evalPer(iBasin)%yStart, evalPer(iBasin)%mStart, evalPer(iBasin)%dStart/)
      end_tmp = (/evalPer(iBasin)%yEnd, evalPer(iBasin)%mEnd, evalPer(iBasin)%dEnd  /)
      ! evaluation gauge
      fName = trim(adjustl(gauge%fname(iGauge)))
      call read_timeseries(trim(fName), udischarge, &
              start_tmp, end_tmp, optimize, opti_function, &
              data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
      data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
      call paste(gauge%Q, data_dp_1d, nodata_dp)
      deallocate (data_dp_1d)
    end do
    !
    ! inflow gauge
    !
    ! in mhm call InflowGauge%Q has to be initialized -- dummy allocation with period of basin 1 and initialization
    if (nInflowGaugesTotal .EQ. 0) then
      allocate(data_dp_1d(maxval(simPer(:)%julEnd - simPer(:)%julStart + 1)))
      data_dp_1d = nodata_dp
      call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
    else
      do iGauge = 1, nInflowGaugesTotal
        ! get basin id
        iBasin = InflowGauge%basinId(iGauge)
        ! get start and end dates
        start_tmp = (/simPer(iBasin)%yStart, simPer(iBasin)%mStart, simPer(iBasin)%dStart/)
        end_tmp = (/simPer(iBasin)%yEnd, simPer(iBasin)%mEnd, simPer(iBasin)%dEnd  /)
        ! inflow gauge
        fName = trim(adjustl(InflowGauge%fname(iGauge)))
        call read_timeseries(trim(fName), udischarge, &
                start_tmp, end_tmp, optimize, opti_function, &
                data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
        if (.NOT. (all(mask_1d))) then
          call message()
          call message('***ERROR: Nodata values in inflow gauge time series. File: ', trim(fName))
          call message('          During simulation period from ', num2str(simPer(iBasin)%yStart) &
                  , ' to ', num2str(simPer(iBasin)%yEnd))
          stop
        end if
        data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
        call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
        deallocate (data_dp_1d)
      end do
    end if

  end subroutine mrm_read_discharge

  ! ---------------------------------------------------------------------------

  !      NAME
  !          mrm_read_total_runoff

  !>         \brief read simulated runoff that is to be routed

  !>         \details read spatio-temporal field of total runoff that has been
  !>         simulated by a hydrologic model or land surface model. This 
  !>         total runoff will then be aggregated to the level 11 resolution
  !>         and then routed through the stream network.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: iBasin"  basin id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>        \note The file containing total runoff must be named total_runoff.nc.
  !>        This file must contain a double precision float variable with the name
  !>        "total_runoff". There must also be an integer variable time with the
  !>        units hours, days, months or years.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author  Stephan Thober
  !         \date    Sep 2015
  !     MODIFIED, Feb 2016, Stephan Thober - refactored deallocate statements
  !               Sep 2016, Stephan Thober - added ALMA convention
  subroutine mrm_read_total_runoff(iBasin)
    use mo_append, only : append
    use mo_read_forcing_nc, only : read_forcing_nc
    use mo_mrm_global_variables, only : &
            dirTotalRunoff, & ! directory of total_runoff file for each basin
            L1_total_runoff_in, & ! simulated runoff at L1
            filenameTotalRunoff, & ! filename
            varnameTotalRunoff ! varname
    use mo_common_variables, only : ALMA_convention, level1
    use mo_common_mHM_mRM_variables, only : timestep, simPer
    use mo_common_constants, only : nodata_dp, HourSecs

    implicit none

    ! input variables
    integer(i4), intent(in) :: iBasin

    ! local variables
    integer(i4) :: tt
    integer(i4) :: nTimeSteps
    integer(i4) :: nctimestep ! tell nc file to read daily or hourly values
    real(dp), dimension(:, :, :), allocatable :: L1_data ! read data from file
    real(dp), dimension(:, :), allocatable :: L1_data_packed

    if (timestep .eq. 1) nctimestep = -4 ! hourly input
    if (timestep .eq. 24) nctimestep = -1 ! daily input
    call read_forcing_nc(trim(dirTotalRunoff(iBasin)), level1(iBasin)%nrows, level1(iBasin)%ncols, &
            varnameTotalRunoff, level1(iBasin)%mask, L1_data, target_period = simPer(iBasin), &
            nctimestep = nctimestep, filename = filenameTotalRunoff)
    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate(L1_data_packed(level1(iBasin)%nCells, nTimeSteps))
    do tt = 1, nTimeSteps
      L1_data_packed(:, tt) = pack(L1_data(:, :, tt), mask = level1(iBasin)%mask)
    end do
    ! free space immediately
    deallocate(L1_data)

    ! convert if ALMA conventions have been given
    if (ALMA_convention) then
      ! convert from kg m-2 s-1 to mm TST-1
      ! 1 kg m-2 -> 1 mm depth
      ! multiply with time to obtain per timestep
      L1_data_packed = L1_data_packed * timestep * HourSecs
      ! ! dump to file
      ! call dump_netcdf('test.nc', L1_data_packed)
    else
      ! convert from mm hr-1 to mm TST-1
      L1_data_packed = L1_data_packed * timestep
    end if

    ! append
    call append(L1_total_runoff_in, L1_data_packed(:, :), nodata_dp)

    !free space
    deallocate(L1_data_packed)

  end subroutine mrm_read_total_runoff

  ! ------------------------------------------------------------------
  ! Rotate fdir variable to the new coordinate system
  ! L. Samaniego & R. Kumar
  ! ------------------------------------------------------------------
  subroutine rotate_fdir_variable(x)
    USE mo_common_constants, ONLY : nodata_i4    ! mHM's global nodata vales

    implicit none

    integer(i4), dimension(:, :), intent(INOUT) :: x

    ! local
    integer(i4) :: i, j

    !-------------------------------------------------------------------
    ! NOTE:
    !
    !     Since the DEM was transposed from (lat,lon) to (lon,lat), i.e.
    !     new DEM = transpose(DEM_old), then
    !
    !     the flow direction X (which was read) for every i, j cell needs
    !     to be rotated as follows
    !
    !                 X(i,j) = [R] * {uVector}
    !
    !     with
    !     {uVector} denoting the fDir_old (read) in vector form, and
    !               e.g. dir 8 => {-1, -1, 0 }
    !     [R]       denting a full rotation matrix to transform the flow
    !               direction into the new coordinate system (lon,lat).
    !
    !     [R] = [rx][rz]
    !
    !     with
    !           |      1       0      0 |
    !     [rx] =|      0   cos T  sin T | = elemental rotation along x axis
    !           |      0  -sin T  cos T |
    !
    !           |  cos F   sin F      0 |
    !     [rz] =| -sin F   cos F      0 | = elemental rotation along z axis
    !           |      0       0      1 |
    !
    !     and T = pi, F = - pi/2
    !     thus
    !          !  0  -1   0 |
    !     [R] =| -1   0   0 |
    !          |  0   0  -1 |
    !     making all 8 directions the following transformation were
    !     obtained.
    !-------------------------------------------------------------------

    do i = 1, size(x, 1)
      do j = 1, size(x, 2)
        if (x(i, j)  .eq. nodata_i4) cycle
        select case (x(i, j))
        case(1)
          x(i, j) = 4
        case(2)
          x(i, j) = 2
        case(4)
          x(i, j) = 1
        case(8)
          x(i, j) = 128
        case(16)
          x(i, j) = 64
        case(32)
          x(i, j) = 32
        case(64)
          x(i, j) = 16
        case(128)
          x(i, j) = 8
        end select
      end do
    end do

  end subroutine rotate_fdir_variable

end module mo_mrm_read_data
