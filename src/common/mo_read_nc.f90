!> \file mo_read_nc.f90
!> \brief \copybrief mo_read_nc
!> \details \copydetails mo_read_nc

!> \brief Reads forcing input data.
!> \details This module is to read forcing input data contained in netcdf files, e.g. temperature, precipitation,
!! total_runoff, lai. Timesteps can be hourly, daily, monthly, and annual. The module provides a subroutine
!! for NetCDF files only. First, the dimensions given are cross-checked with header.txt information. Second,
!! the data of the specified period are read from the specified directory.
!! If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
!! The program is stopped if any value lies out of range.
!> \changelog
!! - Stephan Thober  Sep 2015
!!   - separated routines for netcdf files from routines for binary files
!! - Stephan Thober  Jan 2017
!!   - added reading weights for disaggregation of daily meteorological values to hourly ones
!! - Robert Schweppe Nov 2017
!!   - switched to mo_netcdf library and restuctured routines
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!! - Sebastian Müller Mar 2023
!!   - documentation update
!!   - added bound_error to control if an error is raised if boundaries are violated
!> \authors Juliane Mai
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_read_nc
  use mo_message, only: error_message

  implicit none
  public :: read_nc
  public :: read_const_nc
  public :: read_weights_nc
  private
  !
contains

  !> \brief Reads forcing input in NetCDF file format.
  !> \details Reads netCDF forcing files.
  !! First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !! specified period are read from the specified directory.
  !! If the optional lower and/or upper bound for the data values is given, the read data are checked for
  !! validity.
  !! The program is stopped if any value lies out of range.
  !! If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !! Additionally in this case an mask of vild data points can be received from the routine in maskout.
  !> \changelog
  !! - Stephan Thober     Nov 2013
  !!   - only read required chunk from nc file
  !! - Matthias Cuntz & Juliane Mai Nov 2014 - read daily, monthly or yearly files
  !! - Matthias Zink      Mar 2014
  !!   - added optional nocheck flag and optional maskout
  !! - Stephan Thober     Sep 2015
  !!   - added read for hourly data
  !! - Robert Schweppe    Nov 2017
  !!   - switched to mo_netcdf library and restuctured routines
  !! - Robert Schweppe    Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller   Mar 2023
  !!   - added bound_error to control if an error is raised if boundaries are violated
  !!   - add nTstepForcingDay as argument to be independent of global variables
  !> \authors Matthias Zink
  !> \date May 2013
  subroutine read_nc(folder, nRows, nCols, varName, mask, data, target_period, lower, upper, nctimestep, &
                            fileName, nocheck, maskout, is_meteo, bound_error, nTstepForcingDay)

    use mo_constants, only : nodata_i4
    use mo_common_types, only: period
    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str
    use mo_utils, only : eq, ne

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
    !> Data matrixdim_1 = longitude, dim_2 = latitude, dim_3 = time
    real(dp), dimension(:, :, :), allocatable, intent(out) :: data
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
    !> .TRUE. if check for nodata values deactivated, default = .FALSE. - check is done
    logical, optional, intent(in) :: nocheck
    !> mask of validdata points
    logical, dimension(:, :, :), allocatable, optional, intent(out) :: maskout
    !> logical whether meteorology is currently read
    logical, optional, intent(in) :: is_meteo
    !> .FALSE. to only warn about bound (lower, upper) violations, default = .TRUE. - raise an error
    logical, optional, intent(in) :: bound_error
    !> Number of datapoints in longitudinal direction
    integer(i4), optional, intent(inout) :: nTstepForcingDay

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
    integer(i4) :: i
    ! data nodata value
    real(dp) :: nodata_value
    ! check if model domain is covered by data
    logical :: checking
    ! check if model domain is covered by data
    integer(i4) :: inctimestep
    ! raise error or only warn for bound violations
    logical :: bound_error_

    ! default value for performing checks on read input
    checking = .TRUE.
    if (present(nocheck)) checking = .NOT. nocheck

    ! default value for errors on bound violations
    bound_error_ = .TRUE.
    if (present(bound_error)) bound_error_ = bound_error

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
    if ((var_shape(1) .ne. nRows) .or. (var_shape(2) .ne. nCols)) then
      call error_message('***ERROR: read_nc: mHM generated x and y are not matching NetCDF dimensions')
    end if

    ! determine no data value, use _FillValue first, fall back to missing_value
    if (var%hasAttribute("_FillValue")) then
      call var%getAttribute('_FillValue', nodata_value)
    else if (var%hasAttribute("missing_value")) then
      call var%getAttribute('missing_value', nodata_value)
    else
      call error_message('***ERROR: read_nc: there must be either the attribute "missing_value" or "_FillValue"')
    end if

    ! get time variable
    time_var = nc%getVariable('time')
    ! read the time vector and get start index and count of selection
    call get_time_vector_and_select(time_var, fname, inctimestep, time_start, time_cnt, target_period)

    if (present(is_meteo)) then
      if (.not. present(nTstepForcingDay)) call error_message("read_nc: meteo file given but 'nTstepForcingDay' not passed")
      if (is_meteo) then
        select case(inctimestep)
          case(-1) ! daily
            if (nTstepForcingDay .eq. nodata_i4) then
              nTstepForcingDay = 1_i4
            else if (nTstepForcingDay .ne. 1_i4) then
              call error_message('***ERROR: read_forcing_nc: expected daily input forcing, but read something else. ' // &
                    'Ensure all input time steps have the same units across all input files')
            end if
          case(-4) ! hourly
            if (nTstepForcingDay .eq. nodata_i4) then
              nTstepForcingDay = 24_i4
            else if (nTstepForcingDay .ne. 24_i4) then
              call error_message('***ERROR: read_forcing_nc: expected hourly input forcing, but read something else. ' // &
                    'Ensure all input time steps have the same units across all input files')
            end if
          case default ! no output at all
            call error_message('***ERROR: read_nc: unknown nctimestep switch.')
        end select
      end if
    end if

    ! check optional nctimestep
    if (present(nctimestep)) then
       if (inctimestep .ne. nctimestep) then
          call error_message('***ERROR: provided timestep ' // num2str(nctimestep) //&
                       ' does not match with the one in file ' // num2str(inctimestep), raise=.false.)
          call error_message('File: ' // trim(fname))
       end if
    end if

    ! extract data and select time slice
    call var%getData(data, start = (/1, 1, time_start/), cnt = (/nRows, nCols, time_cnt/))

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
          call error_message('***ERROR: read_nc: nodata value within domain ', raise=.false.)
          call error_message('          boundary in variable: ', trim(varName), raise=.false.)
          call error_message('          at timestep         : ', trim(num2str(i)))
        end if
      end if
      ! optional check
      if (present(lower)) then
        if (any((data(:, :, i) .lt. lower) .AND. mask(:, :))) then
          call error_message('***ERROR: read_nc: values in variable "', &
                  trim(varName), '" are lower than ', trim(num2str(lower, '(F7.2)')), raise=.false.)
          call error_message('          at timestep  : ', trim(num2str(i)), raise=.false.)
          call error_message('File: ', trim(fName), raise=.false.)
          call error_message('Minval at timestep: ', trim(num2str(minval(data(:, :, i)))), raise=.false.)
          call error_message('Total minval: ', trim(num2str(minval(data(:, :, :)))), raise=bound_error_)
        end if
      end if

      if (present(upper)) then
        if (any((data(:, :, i) .gt. upper) .AND. mask(:, :))) then
          call error_message('***ERROR: read_nc: values in variable "', &
                  trim(varName), '" are greater than ', trim(num2str(upper, '(F7.2)')), raise=.false.)
          call error_message('          at timestep  : ', trim(num2str(i)), raise=.false.)
          call error_message('File: ', trim(fName), raise=.false.)
          call error_message('Maxval at timestep: ', trim(num2str(maxval(data(:, :, i)))), raise=.false.)
          call error_message('Total maxval: ', trim(num2str(maxval(data(:, :, :)))), raise=bound_error_)
        end if
      end if

    end do

  end subroutine read_nc


  !> \brief Reads time independent forcing input in NetCDF file format.
  !> \details Reads time independent netCDF forcing files.  \n
  !! First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !! specified period are read from the specified directory.
  !! If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !! The program is stopped if any value lies out of range.\n
  !! If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !! Additionally in this case an mask of vild data points can be received from the routine in maskout.
  !!
  !! \note Files have to be called like defined in mo_files. Furthermore the variable names have to be called
  !!       like they are defined in the declaration of this subroutine. The NetCDF file has to have 2 dimensions:
  !!       1. x, 2. y, It is expected that the variables (especially)within the NetCDF files contain an
  !!       unit attribute. The timestep has to be equidistant.
  !!
  !> \author Lennart Schueler, heavily influenced by read_nc
  !> \date May 2018
  subroutine read_const_nc(folder, nRows, nCols, varName, data, fileName)

    use mo_kind,             only: i4, dp
    use mo_netcdf,           only: NcDataset, NcVariable

    implicit none

    character(len=*),                      intent(in)  :: folder   !< folder where data are stored
    integer(i4),                           intent(in)  :: nRows    !< number of rows of data fields:
    integer(i4),                           intent(in)  :: nCols    !< number of columns of data fields:
    character(len=*),                      intent(in)  :: varName  !< name of NetCDF variable
    real(dp), dimension(:,:), allocatable, intent(out) :: data     !< data read in
    character(256),              optional, intent(in)  :: fileName !< name of file, defaults to varName

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

    ! get dimensions and check if plane is correct
    var_shape = var%getShape()
    if ( (var_shape(1) .ne. nRows) .or. (var_shape(2) .ne. nCols) ) then
       call error_message('***ERROR: read_const_nc: mHM generated x and y are not matching NetCDF dimensions')
    end if

    ! determine no data value, use _FillValue first, fall back to missing_value
    if (var%hasAttribute("_FillValue")) then
      call var%getAttribute('_FillValue', nodata_value)
    else if (var%hasAttribute("missing_value")) then
      call var%getAttribute('missing_value', nodata_value)
    else
       call error_message('***ERROR: read_const_nc: there must be either the attribute "missing_value" or "_FillValue"')
    end if

    ! extract data and select time slice
    call var%getData(data, start=(/1,1/), cnt=(/nRows,nCols/))

  end subroutine read_const_nc


  !> \brief Reads weights for meteo forcings input in NetCDF file format.
  !> \details Reads netCDF weight files.
  !! First, the dimensions given are cross-checked with header.txt information. If the optional lower
  !! and/or upper bound for the data values is given, the read data are checked for validity.
  !! The program is stopped if any value lies out of range.
  !! If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !! Additionally in this case an mask of vild data points can be received from the routine in maskout.
  !> \changelog
  !! - Robert Schweppe    Nov 2017
  !!   - switched to mo_netcdf library and restuctured routine
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - added bound_error to control if an error is raised if boundaries are violated
  !> \authors Stephan Thober & Matthias Zink
  !> \date Jan 2017
  subroutine read_weights_nc(folder, nRows, nCols, varName, data, mask, lower, upper, nocheck, maskout, fileName, bound_error)

    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str
    use mo_utils, only : eq, ne

    implicit none

    !> Name of the folder where data are stored
    character(len = *), intent(in) :: folder
    !> Number of datapoints in longitudinal direction
    integer(i4), intent(in) :: nRows
    !> Number of datapoints in latitudinal  direction
    integer(i4), intent(in) :: nCols
    !> Name of variable name to read
    character(len = *), intent(in) :: varName
    !> Data matrixdim_1 = longitude, dim_2 = latitude, dim_3 = months, dim_4 = hours
    real(dp), dimension(:, :, :, :), allocatable, intent(out) :: data
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
    !> ! mask of validdata points
    logical, dimension(:, :, :, :), allocatable, optional, intent(out) :: maskout
    !> .FALSE. to only warn about bound (lower, upper) violations, default = .TRUE. - raise an error
    logical, optional, intent(in) :: bound_error

    ! name of NetCDF file
    character(256) :: fName
    ! loop variable
    integer(i4) :: i
    ! loop variable
    integer(i4) :: j
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
    ! raise error or only warn for bound violations
    logical :: bound_error_

    checking = .TRUE.
    if (present(nocheck)) checking = .NOT. nocheck

    ! default value for errors on bound violations
    bound_error_ = .TRUE.
    if (present(bound_error)) bound_error_ = bound_error

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
    if ((var_shape(1) .ne. nRows) .or. (var_shape(2) .ne. nCols)) then
       call error_message('***ERROR: read_weights_nc: mHM generated x and y are not matching NetCDF dimensions')
    end if

    ! determine no data value
    call var%getAttribute('missing_value', nodata_value)

    ! extract data
    call var%getData(data)

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
            call error_message('***ERROR: read_weights_nc: nodata value within domain ', raise=.false.)
            call error_message('          boundary in variable: ', trim(varName), raise=.false.)
            call error_message('          at hour         : ', trim(num2str(i)))
          end if
        end if
        ! optional check
        if (present(lower)) then
          if (any((data(:, :, i, j) .lt. lower) .AND. mask(:, :))) then
            call error_message('***ERROR: read_weights_nc: values in variable "', &
                    trim(varName), '" are lower than ', trim(num2str(lower, '(F7.2)')), raise=.false.)
            call error_message('          at hour  : ', trim(num2str(i)), raise=.false.)
            call error_message('File: ', trim(fName), raise=.false.)
            call error_message('Minval at hour: ', trim(num2str(minval(data(:, :, i, j)), '(F7.2)')), raise=.false.)
            call error_message('Total minval: ', trim(num2str(minval(data(:, :, :, :)), '(F7.2)')), raise=bound_error_)
          end if
        end if

        if (present(upper)) then
          if (any((data(:, :, i, j) .gt. upper) .AND. mask(:, :))) then
            call error_message('***ERROR: read_weights_nc: values in variable "', &
                    trim(varName), '" are greater than ', trim(num2str(upper, '(F7.2)')), raise=.false.)
            call error_message('          at hour  : ', trim(num2str(i)), raise=.false.)
            call error_message('File: ', trim(fName), raise=.false.)
            call error_message('Maxval at hour: ', trim(num2str(maxval(data(:, :, i, j)), '(F7.2)')), raise=.false.)
            call error_message('Total maxval: ', trim(num2str(maxval(data(:, :, :, :)), '(F7.2)')), raise=bound_error_)
          end if
        end if

      end do
    end do

  end subroutine read_weights_nc


  !> \brief Extract time vector in unit julian hours and get supposed time step in hours
  !> \changelog
  !! - Matthias Cuntz & Juliane Mai Nov 2014
  !!   - time int or double
  !! - Stephan Thober               Sep 2015
  !!   - added read for hourly data
  !! - Robert Schweppe              Nov 2017
  !!   - restructured routine, reads vector now
  !! - Maren Kaluza                 May 2018
  !!   - fixed bug in time reading
  !! - Stephan Thober               Aug 2020
  !!   - fixed hourly reading
  !! - Stephan Thober               Jan 2022
  !!   - deactivated monthly and annual reading added nTstepForcingDay for hourly reading
  !> \authors Matthias Zink
  !> \date Oct 2012
  subroutine get_time_vector_and_select(var, fname, inctimestep, time_start, time_cnt, target_period)

    use mo_common_types, only: period
    use mo_constants, only : DayHours, DaySecs, YearDays
    use mo_julian, only : caldat, dec2date, julday
    use mo_kind, only : dp, i4, i8
    use mo_netcdf, only : NcVariable
    use mo_string_utils, only : DIVIDE_STRING

    implicit none

    !> variable of interest
    type(NcVariable), intent(in) :: var
    !> fname of ncfile for error message
    character(256), intent(in) :: fname
    !> flag for requested time step
    integer(i4), intent(out) :: inctimestep
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
    integer(i8), allocatable, dimension(:) :: time_diff
    ! period of ncfile, for clipping
    type(period) :: nc_period, clip_period
    integer(i4) :: ncJulSta1, dd, n_time
    integer(i4) :: mmcalstart, mmcalend, yycalstart, yycalend
    integer(i4) :: mmncstart, yyncstart
    ! helper variable for error output
    integer(i4) :: hstart_int, hend_int

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
      call error_message('***ERROR: Please provide the input data in (days, hours, minutes, seconds) ', &
              'since YYYY-MM-DD[ HH:MM:SS] in the netcdf file. Found: ', trim(AttValues))
    end if

    ! get the time vector
    call var%getData(time_data)
    ! convert array from units since to seconds
    time_data = time_data * time_step_seconds

    ! check for length of time vector, needs to be at least of length 2, otherwise step width check fails
    if (size(time_data) .le. 1) then
      call error_message('***ERROR: length of time dimension needs to be at least 2 in file: ' // trim(fname))
    end if

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

    ! calculate input resolution
    allocate(time_diff(n_time - 1))
    time_diff = (time_data(2 : n_time) - time_data(1 : n_time - 1)) / int(DaySecs, i8)
    ! difference must be 1 day
    if (all(abs(time_diff - 1._dp) .lt. 1._dp)) then
       inctimestep = -1 ! daily
    ! difference must be between 28 and 31 days
    else if (all(abs(time_diff) .lt. 32._dp) .and. all(abs(time_diff) .gt. 27._dp)) then
       inctimestep = -2 ! monthly
    ! difference must be between 365 and 366 days
    else if ((all(abs(time_diff) .lt. YearDays + 2)) .and. (all(abs(time_diff) .gt. YearDays - 1._dp))) then
       inctimestep = -3 ! yearly
    ! difference must be 1 hour
    else if (all(abs((time_data(2 : n_time) - time_data(1 : n_time - 1)) / 3600._dp - 1._dp) .lt. 1.e-6)) then
       inctimestep = -4 ! hourly
    else
       call error_message('***ERROR: time_steps are not equal over all times in file and/or do not conform to' // &
       ' requested timestep in file (', trim(fname), ')')
    end if

    ! prepare the selection and check for required time_step
    select case(inctimestep)
    case(-1) ! daily
      ncJulSta1 = nc_period%julStart
      time_start = clip_period%julStart - ncJulSta1 + 1_i4
      time_cnt = clip_period%julEnd - clip_period%julStart + 1_i4
    case(-2) ! monthly
      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! monthly timesteps are usually set by month end, so for beginning, we need 1st of month
      ncJulSta1 = julday(1, mmncstart, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = (yycalstart * 12 + mmcalstart) - (yyncstart * 12 + mmncstart) + 1_i4
      time_cnt = (yycalend * 12 + mmcalend) - (yycalstart * 12 + mmcalstart) + 1_i4
    case(-3) ! yearly
      call caldat(clip_period%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! yearly timesteps are usually set by year end, so for beginning, we need 1st of year
      ncJulSta1 = julday(1, 1, yyncstart)
      call caldat(clip_period%julEnd, dd, mmcalend, yycalend)
      time_start = yycalstart - yyncstart + 1_i4
      time_cnt = yycalend - yycalstart + 1_i4
    case(-4) ! hourly
      ncJulSta1 = nc_period%julStart
      time_start = (clip_period%julStart - ncJulSta1) * 24_i4 + 1_i4 ! convert to hours; always starts at one
      time_cnt = (clip_period%julEnd - clip_period%julStart + 1_i4) * 24_i4 ! convert to hours
    case default ! no output at all
      call error_message('***ERROR: read_nc: unknown nctimestep switch.')
    end select

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 .LE. clip_period%julStart) .AND. (nc_period%julEnd .GE. clip_period%julEnd))) then
      call error_message('***ERROR: read_nc: time period of input data: ', trim(fname), &
              '          is not matching modelling period.')
    end if

    ! free memory
    deallocate(time_diff)

  end subroutine get_time_vector_and_select

end module mo_read_nc
