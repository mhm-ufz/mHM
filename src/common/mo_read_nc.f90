!>       \file mo_read_nc.f90
!>       \brief Reads forcing input data.
!>       \details This module is to read input data contained in netcdf files, e.g. temperature, precipitation,
!>       total_runoff, lai. Timesteps can be hourly, daily, monthly, and annual. The module provides a subroutine
!>       for NetCDF files only.
!>       \authors Juliane Mai
!>       \date Dec 2012

module mo_read_nc
  use mo_kind, only : dp, i4
  use mo_netcdf,           only: NcDataset, NcVariable, NcDimension
  use mo_message,          only: message
  use mo_string_utils,     only: num2str
  use mo_utils,            only: eq, ne, flip

  implicit none
  public :: read_nc
  public :: read_const_nc
  public :: read_weights_nc
  public :: check_sort_order
  public :: common_check_dimension_consistency
  public :: check_dimension_consistency
  interface check_consistency_element
    module procedure check_consistency_element_i4, &
            check_consistency_element_dp
  end interface check_consistency_element


  ! TODO: fyppify this, then generate docstrings for each procedure
  interface check_sort_order
    procedure check_sort_order_2DF64, check_sort_order_3DF64, check_sort_order_4DF64, &
            check_sort_order_2DI4, check_sort_order_3DI4, check_sort_order_4DI4
  end interface

  private
  !
contains


  !>       \brief Reads forcing input in NetCDF file format.
  !>       \details Reads netCDF forcing files.
  !>       First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>       specified period are read from the specified directory.
  !>       If the optional lower and/or upper bound for the data values is given, the read data are checked for
  !>       validity.
  !>       The program is stopped if any value lies out of range.
  !>       If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>       Additionally in this case an mask of vild data points can be received from the routine in maskout.
  subroutine read_nc(folder, nRows, nCols, varName, mask, data, target_period, lower, upper, nctimestep, &
                            fileName, nocheck, maskout)

    use mo_common_datetime_type, only: period

    implicit none

    !> Name of the folder where data are stored
    character(len = *), intent(in) :: folder
    !> Number of datapoints in longitudinal direction
    integer(i4), intent(in) :: nRows
    !> Number of datapoints in latitudinal  direction
    integer(i4), intent(in) :: nCols
    !> Name of variable name to read
    character(len = *), intent(in) :: varName
    !> mask of valid data fields
    logical, dimension(:, :), intent(in) :: mask
    !> Period the data are needed for
    type(period), optional, intent(in) :: target_period
    !> Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower
    !> Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper
    !> timestep in netcdf file
    integer(i4), optional, intent(in) :: nctimestep
    !> name of file, defaults to varName
    character(256), optional, intent(in) :: fileName
    !> .TRUE. if check for nodata values deactivateddefault = .FALSE. - check is done
    logical, optional, intent(in) :: nocheck
    !> Data matrix dim_1 = longitude, dim_2 = latitude, dim_3 = time
    real(dp), dimension(:, :, :), allocatable, intent(out) :: data
    !> mask of validdata points
    logical, dimension(:, :, :), allocatable, optional, intent(out) :: maskout

    ! netcdf file
    type(NcDataset) :: nc

    ! variables for data and time form netcdf
    type(NcVariable) :: var, time_var

    ! shape of NetCDF variable
    integer(i4), allocatable, dimension(:) :: var_shape

    ! index for selecting time vector
    integer(i4) :: time_start

    ! length of vector of selected time values
    integer(i4) :: time_cnt

    ! name of NetCDF file
    character(256) :: fName

    ! loop variable
    integer(i4) :: i, j, k

    ! data nodata value
    real(dp) :: nodata_value

    ! check if model domain is covered by data
    logical :: checking

    ! check if model domain is covered by data
    integer(i4) :: inctimestep


    ! check optional nctimestep
    inctimestep = -1
    if (present(nctimestep)) inctimestep = nctimestep

    ! default value for performing checks on read input
    checking = .TRUE.
    if (present(nocheck)) checking = .NOT. nocheck

    ! default: fName = varname + '.nc'
    ! optional: fName = filename + '.nc'
    fName = varName
    if (present(fileName)) then
      fName = trim(folder) // trim(fileName) // '.nc'
    else
      fName = trim(folder) // trim(fName) // '.nc'
    end if

    ! read the Dataset
    nc = NcDataset(fname, "r")
    ! get the variable
    var = nc%getVariable(trim(varName))

    ! get dimensions and check if plane is correct
    var_shape = var%getShape()
    if ((var_shape(1) /= nRows) .or. (var_shape(2) /= nCols)) then
      print*, '***ERROR: read_nc: mHM generated x and y: ', nRows, nCols , &
              'are not matching NetCDF dimensions: ', var_shape(1), var_shape(2)
      stop 1
    end if

    ! determine no data value, use _FillValue first, fall back to missing_value
    if (var%hasAttribute("_FillValue")) then
      call var%getAttribute('_FillValue', nodata_value)
    else if (var%hasAttribute("missing_value")) then
      call var%getAttribute('missing_value', nodata_value)
    else
      stop '***ERROR: read_nc: there must be either the attribute "missing_value" or "_FillValue"'
    end if

    ! get time variable
    time_var = nc%getVariable('time')
    ! read the time vector and get start index and count of selection
    call get_time_vector_and_select(time_var, fname, inctimestep, time_start, time_cnt, target_period)
    ! extract data and select time slice
    call var%getData(data, start = (/1, 1, time_start/), cnt = (/nRows, nCols, time_cnt/))

    ! flip the data if any dimension is not sorted correctly
    call check_sort_order(data, var)

    ! save output mask if optional maskout is given
    if (present(maskout)) then
      allocate(maskout(var_shape(1), var_shape(2), var_shape(3)))
      maskout = ne(data(:, :, :), nodata_value)
    end if

    ! start checking values
    do i = 1, size(data, dim = 3)
      ! neglect checking for nodata values if optional nocheck is given
      if (checking) then
        if (any(eq(data(:, :, i), nodata_value) .and. (mask))) then
          call message('***ERROR: read_nc: nodata value within domain ')
          call message('          boundary in variable: ', trim(varName))
          call message('          at timestep         : ', trim(num2str(i)))
          do j = 1, size(data, dim = 2)
            do k = 1, size(data, dim = 1)
              if (eq(data(k, j, i), nodata_value) .and. (mask(k, j))) then
                print*, 'at index: ', k, j, ' data is ', data(k, j, i), &
                        ' with nodata_value ', nodata_value, ' and mask ', mask(k, j)
                stop 1
              end if
            end do
          end do
        end if
      end if
      ! optional check
      if (present(lower)) then
        if (any((data(:, :, i) .lt. lower) .AND. mask(:, :))) then
          call message('***ERROR: read_nc: values in variable "', &
                  trim(varName), &
                  '" are lower than ', trim(num2str(lower, '(F7.2)')))
          call message('          at timestep  : ', trim(num2str(i)))
          call message('File: ', trim(fName))
          call message('Minval at timestep: ', trim(num2str(minval(data(:, :, i)))))
          call message('Total minval: ', trim(num2str(minval(data(:, :, :)))))
          stop
        end if
      end if

      if (present(upper)) then
        if (any((data(:, :, i) .gt. upper) .AND. mask(:, :))) then
          call message('***ERROR: read_nc: values in variable "', &
                  trim(varName), &
                  '" are greater than ', trim(num2str(upper, '(F7.2)')))
          call message('          at timestep  : ', trim(num2str(i)))
          call message('File: ', trim(fName))
          call message('Maxval at timestep: ', trim(num2str(maxval(data(:, :, i)))))
          call message('Total maxval: ', trim(num2str(maxval(data(:, :, :)))))
          print*, data(:, :, i)
          stop
        end if
      end if

    end do

    call nc%close()

  end subroutine read_nc

  !>        \brief Reads time independent forcing input in NetCDF file format.
  !>        \details Reads time independent netCDF (forcing) files.  \n
  !>        First, the dimensions given are cross-checked. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.\n
  !>        If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>        Additionally in this case an mask of valid data points can be received from the routine in maskout.
  !>        \note Files have to be called like defined in mo_files. Furthermore the variable names have to be called
  !>              like they are defined in the declaration of this subroutine. The NetCDF file has to have 2 dimensions:
  !>              1. x, 2. y, It is expected that the variables (especially) within the NetCDF files contain an
  !>              unit attribute. The timestep has to be equidistant.
  subroutine read_const_nc(folder, varName, data, fileName, nRows, nCols, maskout)

    character(len=*),                      intent(in)  :: folder  !< folder where data are stored
    character(len=*),                      intent(in)  :: varName !< name of NetCDF variable
    real(dp), dimension(:,:), allocatable, intent(out) :: data    !< data read in
    ! name of file, defaults to varName
    character(256), optional, intent(in) :: fileName
    integer(i4),                           intent(in), optional  :: nRows   !< number of rows of data fields
    integer(i4),                           intent(in), optional  :: nCols   !< number of columns of data fields
    logical, dimension(:, :), allocatable, optional, intent(out) :: maskout  !< mask of valid data points

    ! local variables
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: var          ! variables for data form netcdf
    integer(i4), allocatable, dimension(:) :: var_shape    ! shape of NetCDF variable

    character(256)                         :: fName        ! name of NetCDF file
    real(dp)                               :: nodata_value ! data nodata value

    fName = varName
    if (present(fileName)) then
      fName = trim(folder) // trim(fileName) // '.nc'
    else
      fName = trim(folder) // trim(fName) // '.nc'
    end if

    ! read the Dataset
    nc = NcDataset(fname, "r")
    ! get the variable
    var = nc%getVariable(trim(varName))

    ! determine no data value, use _FillValue first, fall back to missing_value
    if (var%hasAttribute("_FillValue")) then
      call var%getAttribute('_FillValue', nodata_value)
    else if (var%hasAttribute("missing_value")) then
      call var%getAttribute('missing_value', nodata_value)
    else
      stop '***ERROR: read_const_nc: there must be either the attribute "missing_value" or "_FillValue"'
    end if

    ! get dimensions and check if plane is correct
    var_shape = var%getShape()
    if (present(nRows) .and. present(nCols)) then
      if ( var_shape(1) /= nRows .or. var_shape(2) /= nCols) then
        print*, '***ERROR: read_forcing_nc: mHM generated x and y: ', nRows, nCols , &
              'are not matching NetCDF dimensions: ', var_shape(1), var_shape(2)
        stop 1
      end if
      ! extract data and select time slice
      call var%getData(data, start=(/1,1/), cnt=(/nRows,nCols/))
    else
      call var%getData(data)
    end if

    ! save output mask if optional maskout is given
    if (present(maskout)) then
      allocate(maskout(var_shape(1), var_shape(2)))
      maskout = ne(data(:, :), nodata_value)
    end if

    call nc%close()

  end subroutine read_const_nc

  !>       \brief Reads weights for meteo forcings input in NetCDF file format.
  !>       \details Reads netCDF weight files.
  !>       First, the dimensions given are cross-checked. If the optional lower
  !>       and/or upper bound for the data values is given, the read data are checked for validity.
  !>       The program is stopped if any value lies out of range.
  !>       If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>       Additionally in this case an mask of valid data points can be received from the routine in maskout.
  subroutine read_weights_nc(folder, nRows, nCols, varName, data, mask, lower, upper, nocheck, maskout, fileName)

    !> Name of the folder where data are stored
    character(len = *), intent(in) :: folder
    !> Number of datapoints in longitudinal direction
    integer(i4), intent(in) :: nRows
    !> Number of datapoints in latitudinal  direction
    integer(i4), intent(in) :: nCols
    !> Name of variable name to read
    character(len = *), intent(in) :: varName
    !> mask of valid data fields
    logical, dimension(:, :), intent(in) :: mask
    !> Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: lower
    !> Upper bound for check of validity of data values
    real(dp), optional, intent(in) :: upper
    !> .TRUE. if check for nodata values deactivateddefault = .FALSE. - check is done
    logical, optional, intent(in) :: nocheck
    !> name of variable, defaults to fileName
    character(256), optional, intent(in) :: fileName
    !> Data matrixdim_1 = longitude, dim_2 = latitude, dim_3 = months, dim_4 = hours
    real(dp), dimension(:, :, :, :), allocatable, intent(out) :: data
    !> !> mask of validdata points
    logical, dimension(:, :, :, :), allocatable, optional, intent(out) :: maskout

    ! name of NetCDF file
    character(256) :: fName

    ! loop variable
    integer(i4) :: i

    ! loop variable
    integer(i4) :: j, k, l

    ! data nodata value
    real(dp) :: nodata_value

    ! check if model domain is covered by data
    logical :: checking

    ! container for Netcdf data
    type(NcDataset) :: nc

    ! container for Netcdf variable
    type(NcVariable) :: var

    ! shape of NetCDF variable
    integer(i4), allocatable, dimension(:) :: var_shape


    checking = .TRUE.
    if (present(nocheck)) checking = .NOT. nocheck

    fName = varName
    if (present(fileName)) then
      fName = trim(folder) // fileName
    else
      fName = trim(folder) // trim(fName) // '.nc'
    end if

    nc = NcDataset(fname, "r")
    var = nc%getVariable(trim(varName))

    ! get dimensions
    var_shape = var%getShape()
    if ((var_shape(1) /= nRows) .or. (var_shape(2) /= nCols)) then
      print*, '***ERROR: read_nc: mHM generated x and y: ', nRows, nCols , &
              'are not matching NetCDF dimensions: ', var_shape(1), var_shape(2)
      stop 1
    end if

    ! determine no data value
    call var%getAttribute('missing_value', nodata_value)

    ! extract data
    call var%getData(data)

    ! flip the data if any dimension is not sorted correctly
    call check_sort_order(data, var)

    ! save output mask if optional maskout is given
    if (present(maskout)) then
      allocate(maskout(var_shape(1), var_shape(2), var_shape(3), var_shape(4)))
      maskout = ne(data(:, :, :, :), nodata_value)
    end if

    ! start checking values
    do i = 1, var_shape(3)
      do j = 1, var_shape(4)
        ! neglect checking for naodata values if optional nocheck is given
        if (checking) then
          if (any(eq(data(:, :, i, j), nodata_value) .and. (mask))) then
            call message('***ERROR: read_nc: nodata value within domain ')
            call message('          boundary in variable: ', trim(varName))
            call message('          at hour         : ', trim(num2str(i)))
            do k = 1, size(data, dim = 2)
              do l = 1, size(data, dim = 1)
                if (eq(data(l, k, j, i), nodata_value) .and. (mask(l, k))) then
                  print*, 'at index: ', l, k, ' data is ', data(l, k, j, i), &
                          ' with nodata_value ', nodata_value, ' and mask ', mask(l, k)
                  stop 1
                end if
              end do
            end do
          end if
        end if
        ! optional check
        if (present(lower)) then
          if (any((data(:, :, i, j) .lt. lower) .AND. mask(:, :))) then
            call message('***ERROR: read_nc: values in variable "', &
                    trim(varName), &
                    '" are lower than ', trim(num2str(lower, '(F7.2)')))
            call message('          at hour  : ', trim(num2str(i)))
            call message('File: ', trim(fName))
            call message('Minval at hour: ', trim(num2str(minval(data(:, :, i, j)), '(F7.2)')))
            call message('Total minval: ', trim(num2str(minval(data(:, :, :, :)), '(F7.2)')))
            stop
          end if
        end if

        if (present(upper)) then
          if (any((data(:, :, i, j) .gt. upper) .AND. mask(:, :))) then
            call message('***ERROR: read_nc: values in variable "', &
                    trim(varName), &
                    '" are greater than ', trim(num2str(upper, '(F7.2)')))
            call message('          at hour  : ', trim(num2str(i)))
            call message('File: ', trim(fName))
            call message('Maxval at hour: ', trim(num2str(maxval(data(:, :, i, j)), '(F7.2)')))
            call message('Total maxval: ', trim(num2str(maxval(data(:, :, :, :)), '(F7.2)')))
            stop
          end if
        end if

      end do
    end do

    call nc%close()

  end subroutine read_weights_nc

  !>       \brief Extract time vector
  !>       \details Extract time vector in unit julian hours and get supposed time step in hours
  subroutine get_time_vector_and_select(var, fname, inctimestep, time_start, time_cnt, target_period)

    use mo_common_datetime_type, only: period
    use mo_constants, only : DayHours, DaySecs, YearDays
    use mo_julian, only : caldat, dec2date, julday
    use mo_kind, only : i8
    use mo_string_utils, only : DIVIDE_STRING

    implicit none

    !> variable of interest
    type(NcVariable), intent(in) :: var
    !> fname of ncfile for error message
    character(256), intent(in) :: fname
    !> flag for requested time step
    integer(i4), intent(in) :: inctimestep
    !> time_start index of time selection
    integer(i4), intent(out) :: time_start
    !> time_count of indexes of time selection
    integer(i4), intent(out) :: time_cnt
    !> reference period
    type(period), intent(in), optional :: target_period

    ! reference time of NetCDF
    integer(i4) :: yRef, dRef, mRef, hRef, jRef

    ! netcdf attribute values
    character(256) :: AttValues

    ! dummies for netcdf attribute handling
    character(256), dimension(:), allocatable :: strArr, date, time

    ! native time step converter in ncfile
    integer(i8) :: time_step_seconds

    ! time vector
    integer(i8), allocatable, dimension(:) :: time_data

    ! period of ncfile, for clipping
    type(period) :: nc_period, clip_period

    integer(i4) :: ncJulSta1, dd, n_time

    integer(i4) :: mmcalstart, mmcalend, yycalstart, yycalend

    integer(i4) :: mmncstart, yyncstart

    ! helper variable for error output
    integer(i4) :: hstart_int, hend_int

    ! helper variable for error output
    character(256) :: error_msg


    call var%getAttribute('units', AttValues)
    ! AttValues looks like "<unit> since YYYY-MM-DD[ HH:MM:SS]"
    ! split at space
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)

    ! determine reference time at '-' and convert to integer
    call DIVIDE_STRING(trim(strArr(3)), '-', date)
    read(date(1), *) yRef
    read(date(2), *) mRef
    read(date(3), *) dRef

    jRef = julday(dd = dRef, mm = mRef, yy = yRef)

    ! if existing also read in the time (only hour so far)
    hRef = 0
    if(size(strArr) .gt. 3) then
      call DIVIDE_STRING(trim(strArr(4)), ':', time)
      read(time(1), *) hRef
    end if

    ! determine the step_size
    if (strArr(1) .EQ. 'days') then
      time_step_seconds = int(DaySecs)
    else if (strArr(1) .eq. 'hours') then
      time_step_seconds = int(DaySecs / DayHours)
    else if (strArr(1) .eq. 'minutes') then
      time_step_seconds = int(DaySecs / DayHours / 60._dp)
    else if (strArr(1) .eq. 'seconds') then
      time_step_seconds = 1_i8
    else
      call message('***ERROR: Please provide the input data in (days, hours, minutes, seconds) ', &
              'since YYYY-MM-DD[ HH:MM:SS] in the netcdf file. Found: ', trim(AttValues))
      stop 1
    end if

    ! get the time vector
    call var%getData(time_data)
    ! convert array from units since to seconds
    time_data = time_data * time_step_seconds

    ! check for length of time vector, needs to be at least of length 2, otherwise step width check fails
    if (size(time_data) .le. 1) then
      call message('***ERROR: length of time dimension needs to be at least 2 in file: ' // trim(fname))
      stop 1
    end if

    ! check for equal timesteps and timestep must not be multiple of native timestep
    error_msg = '***ERROR: time_steps are not equal over all times in file and/or do not conform to' // &
            ' requested timestep in file (' // trim(fname) // ') : '

    ! compare the read period from ncfile to the period required
    ! convert julian second information back to date via conversion to float
    ! the 0.5_dp is for the different reference of fractional julian days, hours are truncated
    n_time = size(time_data)
    call dec2date(time_data(1) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dStart, nc_period%mStart, &
            nc_period%yStart, hstart_int)
    nc_period%julStart = int(time_data(1) / DaySecs + jRef + hRef / 24._dp)
    call dec2date(time_data(n_time) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dEnd, nc_period%mEnd, &
            nc_period%yEnd, hend_int)
    nc_period%julEnd = int(time_data(n_time) / DaySecs + jRef + hRef / 24._dp)

    ! if no target period is present, use the whole time period
    if (present(target_period)) then
      clip_period = target_period
    else
      clip_period = nc_period
    end if

    ! prepare the selection and check for required time_step
    select case(inctimestep)
    case(-1) ! daily
      ! difference must be 1 day
      if (.not. all(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs - 1._dp) .lt. 1.e-6)) then
        call message(error_msg // trim('daily'))
        stop 1
      end if
      ncJulSta1 = nc_period%julStart
      time_start = clip_period%julStart - ncJulSta1 + 1_i4
      time_cnt = clip_period%julEnd - clip_period%julStart + 1_i4
    case(-2) ! monthly
      ! difference must be between 28 and 31 days
      if (any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .gt. 31._dp) .or. &
              any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .lt. 28._dp)) then
        call message(error_msg // trim('monthly'))
        stop 1
      end if

      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! monthly timesteps are usually set by month end, so for beginning, we need 1st of month
      ncJulSta1 = julday(1, mmncstart, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = (yycalstart * 12 + mmcalstart) - (yyncstart * 12 + mmncstart) + 1_i4
      time_cnt = (yycalend * 12 + mmcalend) - (yycalstart * 12 + mmcalstart) + 1_i4
    case(-3) ! yearly
      ! difference must be between 365 and 366 days
      if (any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .gt. (YearDays + 1._dp)) .or. &
              any(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / DaySecs) .lt. YearDays)) then
        call message(error_msg // 'yearly')
        stop 1
      end if
      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! yearly timesteps are usually set by year end, so for beginning, we need 1st of year
      ncJulSta1 = julday(1, 1, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = yycalstart - yyncstart + 1_i4
      time_cnt = yycalend - yycalstart + 1_i4
    case(-4) ! hourly
      ! difference must be 1 hour
      if (.not. all(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / 3600._dp - 1._dp) .lt. 1.e-6)) then
        call message(error_msg // 'hourly')
        stop 1
      end if
      ncJulSta1 = nc_period%julStart
      time_start = (clip_period%julStart - ncJulSta1) * 24_i4 + 1_i4 ! convert to hours; always starts at one
      time_cnt = (clip_period%julEnd - clip_period%julStart + 1_i4) * 24_i4 ! convert to hours
    case default ! no output at all
      call message('***ERROR: read_nc: unknown nctimestep switch.')
      stop 1
    end select

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 .LE. clip_period%julStart) .AND. (nc_period%julEnd .GE. clip_period%julEnd))) then
      call message('***ERROR: read_nc: time period of input data: ', trim(fname), &
              '          is not matching modelling period.')
      stop 1
    end if

  end subroutine get_time_vector_and_select

  subroutine check_sort_order_2DF64(data, var)
    real(dp), dimension(:, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_2DF64

  subroutine check_sort_order_3DF64(data, var)
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_3DF64

  subroutine check_sort_order_4DF64(data, var)
    real(dp), dimension(:, :, :, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_4DF64

  subroutine check_sort_order_3DI4(data, var)
    integer(i4), dimension(:, :, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_3DI4

  subroutine check_sort_order_2DI4(data, var)
    integer(i4), dimension(:, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_2DI4

  subroutine check_sort_order_4DI4(data, var)
    integer(i4), dimension(:, :, :, :), allocatable, intent(inout) :: data
    type(NcVariable), intent(in) :: var

    type(NcDimension), dimension(:), allocatable :: ncDims
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: ncVarShape
    integer(i4) :: iDim
    real(dp), dimension(:), allocatable :: data1D
    real(dp), dimension(:, :), allocatable :: data2D
    character(256) :: dimName

    ! get all dimensions of variable
    ncDims = var%getDimensions()
    do iDim=1, size(ncDims)
      ! get the dimension name
      dimName = ncDims(iDim)%getName()
      ! is it also a variable, then we need to get the values
      if (var%parent%hasVariable(trim(dimName))) then
        ncVar = var%parent%getVariable(trim(dimName))
        ncVarShape = ncVar%getShape()
        ! now check the number of dimensions for each
        if (size(ncVarShape) == 1_i4) then
          call ncVar%getData(data1D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data1D(1) > data1D(size(data1D))) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data1D)
        else if (size(ncVarShape) == 2_i4) then
          call ncVar%getData(data2D)
          ! if it is sorted in a descending way, sort along the dimension
          if (data2D(1, 1) > data2D(size(data2D, 1), 1)) then
            call flip(data, iDim)
            print*, 'sorted coordinate: ', trim(dimName), ' for variable: ', trim(var%getName())
          end if
          deallocate(data2D)
        end if
      end if
      ! if it is not Variable it is always integer and sorted in an ascending order
    end do

  end subroutine check_sort_order_4DI4

  subroutine common_check_dimension_consistency(iDomain, uniqueIDomain, boundaries, select_indices)
    use mo_common_variables, only: nLandCoverPeriods, landCoverPeriodBoundaries
    use mo_string_utils, only: compress
    use mo_common_datetime_type, only: simPer, LCyearId

    integer(i4), intent(in) :: iDomain, uniqueIDomain

    real(dp), dimension(:), intent(inout) :: boundaries
    integer(i4), dimension(:), intent(out), allocatable :: select_indices
    logical, dimension(size(boundaries) - 1) :: select_indices_mask
    logical, dimension(size(boundaries)) :: select_indices_temp

    integer(i4) :: select_index, iBoundary, LCyearStart, LCyearEnd

    allocate(select_indices(size(boundaries) - 1))
    select_index = 0_i4
    select_indices_mask = .false.
    LCyearStart = lbound(LCyearId, dim=1)
    LCyearEnd = ubound(LCyearId, dim=1)

    ! set the correct indices to use
    do iBoundary=1, size(boundaries) - 1
      ! check for overlap ((StartA <= EndB) and (EndA >= StartB))
      ! https://stackoverflow.com/questions/325933/
      ! TODO: MPR reinstate, if all landcover periods are to be set
      !if ((boundaries(iBoundary) <= simPer(iDomain)%yend) .and. &
      !        ((boundaries(iBoundary+1)-1) >= simPer(iDomain)%ystart)) then
        ! advance counter
        select_index = select_index + 1_i4
        ! select this iBoundary from dimension
        select_indices_mask(iBoundary) = .true.
        ! set the correct LCyearId
        LCyearId(&
                maxval([int(boundaries(iBoundary)), LCyearStart]):&
                minval([int(boundaries(iBoundary+1)), LCyearEnd]), uniqueIDomain) = select_index
        ! set the boundaries as parsed
        landCoverPeriodBoundaries(select_index, uniqueIDomain) = int(boundaries(iBoundary))
        landCoverPeriodBoundaries(select_index + 1, uniqueIDomain) = int(boundaries(iBoundary + 1))
      !end if
    end do
    select_indices = pack([(iBoundary, iBoundary=1, size(boundaries) - 1)], select_indices_mask)

    ! check if number of periods in namelist is enough
    if (nLandCoverPeriods < select_index) then
      call message('The number of selected land cover periods for domain ', compress(num2str(iDomain)), &
              ' is bigger than allowed (', &
              compress(num2str(nLandCoverPeriods)), '). Please set nLandCoverPeriods in namelist.')
      stop 1
    end if

    ! check if both start and end are covered
    select_indices_temp = [select_indices_mask, .false.]
    if (minval(boundaries, mask=select_indices_temp) > simPer(iDomain)%ystart) then
      call message('The selected land cover periods for domain ', compress(trim(num2str(iDomain))), &
              ' (', compress(trim(num2str(minval(boundaries, mask=select_indices_temp)))), &
              ') do not cover the beginning of the simulation period (', &
              compress(trim(num2str(simPer(iDomain)%ystart))), ').')
      stop 1
    end if
    select_indices_temp = [.false., select_indices_mask]
    if (maxval(boundaries, mask=select_indices_temp) < simPer(iDomain)%yend) then
      call message('The selected land cover periods for domain ', compress(trim(num2str(iDomain))), &
              ' (', compress(trim(num2str(maxval(boundaries, mask=select_indices_temp)))), &
              ' ) do not cover the end of the simulation period (', &
              compress(trim(num2str(simPer(iDomain)%yend))), ').')
      stop 1
    end if
  end subroutine common_check_dimension_consistency

  subroutine check_dimension_consistency(iDomain, uniqueIDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp, &
          landCoverSelect, check_all_arg)
    use mo_global_variables, only: nSoilHorizons, soilHorizonBoundaries, nLAIs, LAIBoundaries
    use mo_common_variables, only: nLandCoverPeriods
    use mo_string_utils, only: compress, num2str
    use mo_utils, only: ne
    use mo_message, only: message

    integer(i4), intent(in) :: iDomain, uniqueIDomain

    integer(i4), intent(in) :: nSoilHorizons_temp, nLAIs_temp, nLandCoverPeriods_temp
    real(dp), dimension(:), intent(inout) :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
            LAIBoundaries_temp
    integer(i4), dimension(:), allocatable, intent(out) :: landCoverSelect
    logical, intent(in), optional :: check_all_arg

    integer(i4) :: k
    logical :: check_all

    check_all = .true.
    if (present(check_all_arg)) check_all = check_all_arg
    if (iDomain == 1 .and. check_all) then
      ! set local to global
      nSoilHorizons = nSoilHorizons_temp
      nLAIs = nLAIs_temp
      nLandCoverPeriods = nLandCoverPeriods_temp
      ! TODO: MPR remove if clause here
      if (.not. allocated(soilHorizonBoundaries)) allocate(soilHorizonBoundaries(nSoilHorizons))
      soilHorizonBoundaries = soilHorizonBoundaries_temp
      allocate(LAIBoundaries(nLAIs+1))
      LAIBoundaries = LAIBoundaries_temp
    else
      ! check if it conforms with global
      if (nSoilHorizons /= nSoilHorizons_temp) then
        call message('The number of soil horizons for domain 1 (', compress(trim(num2str(nSoilHorizons))), &
                ') does not conform with the number of soil horizons for domain ', &
                compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(nSoilHorizons_temp))), ').')
        stop 1
      end if
      if (nLAIs /= nLAIs_temp) then
        call message('The number of soil horizons for domain 1 (', compress(trim(num2str(nLAIs))), &
                ') does not conform with the number of soil horizons for domain ', &
                compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(nLAIs_temp))), ').')
        stop 1
      end if
      do k=1, nSoilHorizons+1
        if (ne(soilHorizonBoundaries(k), soilHorizonBoundaries_temp(k))) then
          call message('The ',compress(trim(num2str(k))),'th soil horizon boundary for domain 1 (', &
                  compress(trim(num2str(soilHorizonBoundaries(k)))), &
                  ') does not conform with domain ', &
                  compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(soilHorizonBoundaries_temp(k)))), ').')
          stop 1
        end if
      end do
      do k=1, nLAIs+1
        if (ne(LAIBoundaries(k), LAIBoundaries_temp(k))) then
          call message('The ',compress(trim(num2str(k))),'th LAI period boundary for domain 1 (', &
                  compress(trim(num2str(LAIBoundaries(k)))), ') does not conform with domain ', &
                  compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(LAIBoundaries_temp(k)))), ').')
          stop 1
        end if
      end do
    end if
    call common_check_dimension_consistency(iDomain, uniqueIDomain, landCoverPeriodBoundaries_temp, landCoverSelect)

  end subroutine check_dimension_consistency

  subroutine check_consistency_element_dp(item1, item2, name, iDomain)
    use mo_utils, only: ne
    use mo_string_utils, only: compress, num2str
    use mo_message, only: message

    real(dp), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iDomain

    if (ne(item1, item2)) then
      call message('The ', trim(name),&
                  ' as set in the configuration file (', &
                  compress(trim(num2str(item1))), &
                  ') does not conform with domain ', &
                  compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(item2))), ').')
      stop 1
    end if
  end subroutine check_consistency_element_dp

  subroutine check_consistency_element_i4(item1, item2, name, iDomain)
    use mo_utils, only: ne
    use mo_string_utils, only: compress, num2str
    use mo_message, only: message

    integer(i4), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iDomain

    if (item1 /= item2) then
      call message('The ', trim(name),&
                  ' as set in the configuration file (', &
                  compress(trim(num2str(item1))), &
                  ') does not conform with domain ', &
                  compress(trim(num2str(iDomain))), ' (', compress(trim(num2str(item2))), ').')
      stop 1
    end if
  end subroutine check_consistency_element_i4

  subroutine check_consistency()
    use mo_global_variables, only : nSoilHorizons
    use mo_common_variables, only : opti_function, optimize
    use mo_global_variables, only : nSoilHorizons_sm_input
    use mo_string_utils, only: num2str
    use mo_message, only: message

    if (optimize) then
      select case (opti_function)
      case(10 : 13, 28)
        ! soil moisture
        if (nSoilHorizons_sm_input > nSoilHorizons) then
          call message()
          call message('***ERROR: Number of soil horizons representative for input soil moisture exceeded')
          call message('          defined number of soil horizions in mHM: ', &
                  adjustl(trim(num2str(nSoilHorizons))), '!')
          stop 1
        end if
      end select
    end if

  end subroutine check_consistency

end module mo_read_nc
