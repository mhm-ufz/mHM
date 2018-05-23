!> \file mo_read_forcing_nc.f90

!> \brief Reads forcing input data.

!> \details This module is to read forcing input data contained in netcdf files, e.g. temperature, precipitation,
!> total_runoff, lai. Timesteps can be hourly, daily, monthly, and annual. The module provides a subroutine
!> for NetCDF files only. First, the dimensions given are cross-checked with header.txt information. Second,
!> the data of the specified period are read from the specified directory.\n
!> If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
!> The program is stopped if any value lies out of range.

!> \authors Juliane Mai
!> \date Dec 2012
!  Modified Sep 2015, Stephan Thober  - separated routines for netcdf files from routines for binary files
!           Jan 2017, Stephan Thober  - added reading weights for disaggregation of daily meteorological values to hourly ones
!           Nov 2017, Robert Schweppe - switched to mo_netcdf library and restuctured routines

module mo_read_forcing_nc
  implicit none
  public :: read_forcing_nc
  public :: read_weights_nc
  private
  !
contains


  ! ------------------------------------------------------------------

  !     NAME
  !         read_forcing_nc

  !     PURPOSE
  !>        \brief Reads forcing input in NetCDF file format.

  !>        \details Reads netCDF forcing files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.\n
  !>        If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>        Additionally in this case an mask of vild data points can be received from the routine in maskout.

  !     CALLING SEQUENCE
  !         periode%dStart   = 2_i4                                                     ! day
  !         periode%mStart   = 2_i4                                                     ! month
  !         periode%yStart   = 1972_i4                                                  ! year
  !         periode%dEnd     = 7_i4                                                     ! day
  !         periode%mEnd     = 8_i4                                                     ! month
  !         periode%yEnd     = 1977_i4                                                  ! year
  !         periode%julStart = NDAYS(periode%dStart, periode%mStart, periode%yStart)    ! julian day starting
  !         periode%julEnd   = NDAYS(periode%dEnd,   periode%mEnd,   periode%yEnd  )    ! julian day ending
  !         periode%nObs     = periode%julEnd - periode%julStart + 1_i4                 ! total number of observation

  !         call read_forcing_bin('old_code/sub_00020/input/forcing/pre/' , 21_i4, 28_i4, periode, precipitation,   &
  !                             lower=  0.0_dp)
  !         call read_forcing_bin('old_code/sub_00020/input/forcing/pet/' , 21_i4, 28_i4, periode, pot_evapo_trans, &
  !                             lower=  0.0_dp, upper=20.0_dp)
  !         call read_forcing_bin('old_code/sub_00020/input/forcing/tavg/', 21_i4, 28_i4, periode, temp_average,    &
  !                             lower=-50.0_dp, upper=50.0_dp)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: folder"        Name of the folder where data are stored
  !>        \param[in] "integer(i4)      :: nRows"         Number of datapoints in longitudinal direction
  !>        \param[in] "integer(i4)      :: nCols"         Number of datapoints in latitudinal  direction
  !>        \param[in] "type(period)     :: periode"       Period the data are needed for
  !>        \param[in] "character(len=*) :: varName"       Name of variable name to read
  !>        \param[in] "logical, dimension(:,:) :: mask"   mask of valid data fields

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:,:,:) :: data"     Data matrix
  !>                                                             dim_1 = longitude, dim_2 = latitude, dim_3 = time

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp),   optional, intent(in) :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp),   optional, intent(in) :: upper"    Upper bound for check of validity of data values
  !>        \param[in] "logical,    optional, intent(in) :: nocheck"  .TRUE. if check for nodata values deactivated
  !>                                                                   default = .FALSE. - check is done
  !>        \param[in] "integer(i4) optional, intent(in) :: nctimestep" timestep in netcdf file

  
  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[in] "logical, dimension(:,:,:), allocatable,  optional, intent(out) :: maskout"  ! mask of valid
  !>                                                                                                  data points 

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>        \note Files have to be called like defined in mo_files. Furthermore the variable names have to be called
  !>              like they are defined in the declaration of this subroutine. The NetCDF file has to have 3 dimensions:
  !>              1. x, 2. y, 3. t. It is expected that the variables (especially)within the NetCDF files contain an
  !>              unit attribute. The timestep has to be equidistant.

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date May 2013
  !               modified Stephan Thober     Nov 2013 - only read required
  !                                                      chunk from nc file
  !                        Matthias Cuntz & Juliane Mai Nov 2014 - read daily, monthly or yearly files
  !                        Matthias Zink      Mar 2014 - added optional nocheck flag and optional maskout
  !                        Stephan Thober     Sep 2015 - added read for hourly data
  !                        Robert Schweppe    Nov 2017 - switched to mo_netcdf library and restuctured routines

  subroutine read_forcing_nc(folder, nRows, nCols, periode, varName, data, mask, lower, upper, nctimestep, &
          nocheck, maskout, fileName)

    use mo_kind,             only: i4, dp
    use mo_common_variables, only: period
    use mo_message,          only: message
    use mo_netcdf,           only: NcDataset, NcVariable, NcDimension
    use mo_string_utils,     only: num2str
    use mo_utils,            only: eq, ne

    implicit none

    character(len=*),                                  intent(in)  :: folder     ! folder where data are stored
    integer(i4),                                       intent(in)  :: nRows      ! number of rows of data fields:
    integer(i4),                                       intent(in)  :: nCols      ! number of columns of data fields:
    type(period),                                      intent(in)  :: periode    ! time period
    character(len=*),                                  intent(in)  :: varName    ! name of NetCDF variable
    logical, dimension(:,:),                           intent(in)  :: mask       ! mask of valid data fields
    real(dp),                                optional, intent(in)  :: lower      ! lower bound for data points
    real(dp),                                optional, intent(in)  :: upper      ! upper bound for data points
    integer(i4),                             optional, intent(in)  :: nctimestep ! -1: daily (default);
    !                                                                            ! -2: monthly;
    !                                                                            ! -3: yearly;
    !                                                                            ! -4: hourly; DEPRECATED!!!
    character(256),                          optional, intent(in)  :: fileName   ! name of variable, defaults to fileName
    logical,                                 optional, intent(in)  :: nocheck    ! .TRUE. if check for nodata values deactivated
    !                                                                            ! default = .FALSE. - check is done
    real(dp), dimension(:,:,:), allocatable,           intent(out) :: data       ! data read in
    logical,  dimension(:,:,:), allocatable, optional, intent(out) :: maskout    ! mask of data to read

    !
    ! local variables
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: var, time_var! variables for data and time form netcdf
    integer(i4), allocatable, dimension(:) :: var_shape    ! shape of NetCDF variable
    integer(i4)                            :: time_start   ! index for selecting time vector
    integer(i4)                            :: time_cnt     ! length of vector of selected time values

    character(256)                         :: fName        ! name of NetCDF file
    integer(i4)                            :: i            ! loop variable
    real(dp)                               :: nodata_value ! data nodata value
    logical                                :: checking     ! check if model domain is covered by data
    integer(i4)                            :: inctimestep  ! check if model domain is covered by data

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
    if ( (var_shape(1) .ne. nRows) .or. (var_shape(2) .ne. nCols) ) then
       stop '***ERROR: read_forcing_nc: mHM generated x and y are not matching NetCDF dimensions'
    end if

    ! determine no data value, use _FillValue first, fall back to missing_value
    if (var%hasAttribute("_FillValue")) then
      call var%getAttribute('_FillValue', nodata_value)
    else if (var%hasAttribute("missing_value")) then
      call var%getAttribute('missing_value', nodata_value)
    else
      stop '***ERROR: read_forcing_nc: there must be either the attribute "missing_value" or "_FillValue"'
    end if

    ! get time variable
    time_var = nc%getVariable('time')
    ! read the time vector and get start index and count of selection
    call get_time_vector_and_select(time_var, fname, inctimestep, periode, time_start, time_cnt)
    ! extract data and select time slice
    call var%getData(data, start=(/1,1,time_start/), cnt=(/nRows,nCols,time_cnt/))

    ! save output mask if optional maskout is given
    if (present(maskout)) then
       allocate(maskout(var_shape(1), var_shape(2), var_shape(3)))
       maskout = ne(data(:,:,:), nodata_value)
    end if

    ! start checking values
    do i = 1, size(data, dim=3)
       ! neglect checking for nodata values if optional nocheck is given
       if (checking) then
          if (any(eq(data(:,:,i),nodata_value) .and. (mask))) then
             call message('***ERROR: read_forcing_nc: nodata value within basin ')
             call message('          boundary in variable: ', trim(varName))
             call message('          at timestep         : ', trim(num2str(i)))
             stop
          end if
       end if
       ! optional check
       if (present(lower)) then
          if ( any( (data(:,:,i) .lt. lower) .AND. mask(:,:) )  ) then
             call message('***ERROR: read_forcing_nc: values in variable "', &
                  trim(varName),                                           &
                  '" are lower than ', trim(num2str(lower,'(F7.2)')) )
             call message('          at timestep  : ', trim(num2str(i)))
             call message('File: ', trim(fName))
             call message('Minval at timestep: ', trim(num2str(minval(data(:,:,i)))))
             call message('Total minval: ', trim(num2str(minval(data(:,:,:)))))
             stop
          end if
       end if

       if (present(upper)) then
          if ( any( (data(:,:,i) .gt. upper) .AND. mask(:,:) )  ) then
             call message('***ERROR: read_forcing_nc: values in variable "',  &
                  trim(varName),                                           &
                  '" are greater than ', trim(num2str(upper,'(F7.2)')) )
             call message('          at timestep  : ', trim(num2str(i)))
             call message('File: ', trim(fName))
             call message('Maxval at timestep: ', trim(num2str(maxval(data(:,:,i)))))
             call message('Total maxval: ', trim(num2str(maxval(data(:,:,:)))))
             print*, data(:,:,i)
             stop
          end if
       end if

    end do

  end subroutine read_forcing_nc


  ! ------------------------------------------------------------------

  !     NAME
  !         read_weights_nc

  !     PURPOSE
  !>        \brief Reads weights for meteo forcings input in NetCDF file format.

  !>        \details Reads netCDF weight files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. If the optional lower
  !>        and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.\n
  !>        If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>        Additionally in this case an mask of vild data points can be received from the routine in maskout.

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "character(len=*)        :: folder"        Name of the folder where data are stored
  !>        \param[in] "integer(i4)             :: nRows"         Number of datapoints in longitudinal direction
  !>        \param[in] "integer(i4)             :: nCols"         Number of datapoints in latitudinal  direction
  !>        \param[in] "character(len=*)        :: varName"       Name of variable name to read
  !>        \param[in] "logical, dimension(:,:) :: mask"          mask of valid data fields

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:,:,:,:) :: data"     Data matrix
  !>                                                               dim_1 = longitude, dim_2 = latitude, dim_3 = months, dim_4 = hours

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), optional, intent(in)  :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional, intent(in)  :: upper"    Upper bound for check of validity of data values
  !>        \param[in] "logical,  optional, intent(in)  :: nocheck"  .TRUE. if check for nodata values deactivated
  !>                                                                  default = .FALSE. - check is done
  
  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[in] "logical, dimension(:,:,:,:), allocatable,  optional, intent(out) :: maskout"  ! mask of valid
  !>                                                                                                  data points 

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober & Matthias Zink
  !>        \date Jan 2017
  !               modified Robert Schweppe    Nov 2017 - switched to mo_netcdf library and restuctured routine

  
  subroutine read_weights_nc(folder, nRows, nCols, varName, data, mask, lower, upper, nocheck, maskout, fileName)

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_netcdf,           only: NcDataset, NcVariable

    use mo_string_utils,      only: num2str
    use mo_utils,             only: eq, ne

    implicit none

    character(len=*),                                    intent(in)  :: folder  ! folder where data are stored
    integer(i4),                                         intent(in)  :: nRows   ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4),                                         intent(in)  :: nCols   ! number of columns of data fields:
    ! LATITUDE dimension
    character(len=*),                                    intent(in)  :: varName ! name of NetCDF variable
    logical, dimension(:,:),                             intent(in)  :: mask    ! mask of valid data fields
    real(dp),                                  optional, intent(in)  :: lower   ! lower bound for data points
    real(dp),                                  optional, intent(in)  :: upper   ! upper bound for data points
    logical,                                   optional, intent(in)  :: nocheck ! .TRUE. if check for nodata values deactivated
    !                                                                           ! default = .FALSE. - check is done
    character(256),                            optional, intent(in)  :: fileName! name of variable, defaults to fileName
    real(dp), dimension(:,:,:,:), allocatable,           intent(out) :: data    ! data read in
    !                                                                           ! dim 1 = rows
    !                                                                           ! dim 2 = cols
    !                                                                           ! dim 3 = months
    !                                                                           ! dim 4 = hours
    logical,  dimension(:,:,:,:), allocatable, optional, intent(out) :: maskout ! mask of data to read

    ! local variables
    character(256)                         :: fName        ! name of NetCDF file
    integer(i4)                            :: i            ! loop variable
    integer(i4)                            :: j            ! loop variable
    real(dp)                               :: nodata_value ! data nodata value
    logical                                :: checking     ! check if model domain is covered by data
    type(NcDataset)                        :: nc           ! container for Netcdf data
    type(NcVariable)                       :: var          ! container for Netcdf variable
    integer(i4), allocatable, dimension(:) :: var_shape    ! shape of NetCDF variable

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
    if ( (var_shape(1) .ne. nRows) .or. (var_shape(2) .ne. nCols) ) then
       stop '***ERROR: read_forcing_nc: mHM generated x and y are not matching NetCDF dimensions'
    end if

    ! determine no data value
    call var%getAttribute('missing_value', nodata_value)

    ! extract data
    call var%getData(data)

    ! save output mask if optional maskout is given
    if (present(maskout)) then
       allocate(maskout(var_shape(1), var_shape(2), var_shape(3), var_shape(4)))
       maskout = ne(data(:,:,:,:),nodata_value)
    end if
       
    ! start checking values
    do i = 1, var_shape(3)
       do j = 1, var_shape(4)
          ! neglect checking for naodata values if optional nocheck is given
          if (checking) then
             if (any(eq(data(:,:,i,j),nodata_value) .and. (mask))) then
                call message('***ERROR: read_forcing_nc: nodata value within basin ')
                call message('          boundary in variable: ', trim(varName))
                call message('          at hour         : ', trim(num2str(i)))
                stop
             end if
          end if
          ! optional check
          if (present(lower)) then
             if ( any( (data(:,:,i,j) .lt. lower) .AND. mask(:,:) )  ) then
                call message('***ERROR: read_forcing_nc: values in variable "', &
                     trim(varName),                                           &
                     '" are lower than ', trim(num2str(lower,'(F7.2)')) )
                call message('          at hour  : ', trim(num2str(i)))
                call message('File: ', trim(fName))
                call message('Minval at hour: ', trim(num2str(minval(data(:,:,i,j)),'(F7.2)')))
                call message('Total minval: ', trim(num2str(minval(data(:,:,:,:)),'(F7.2)')))
                stop
             end if
          end if

          if (present(upper)) then
             if ( any( (data(:,:,i,j) .gt. upper) .AND. mask(:,:) )  ) then
                call message('***ERROR: read_forcing_nc: values in variable "',  &
                     trim(varName),                                           &
                     '" are greater than ', trim(num2str(upper,'(F7.2)')) )
                call message('          at hour  : ', trim(num2str(i)))
                call message('File: ', trim(fName))
                call message('Maxval at hour: ', trim(num2str(maxval(data(:,:,i,j)),'(F7.2)')))
                call message('Total maxval: ', trim(num2str(maxval(data(:,:,:,:)),'(F7.2)')))
                stop
             end if
          end if
       
       end do
    end do

  end subroutine read_weights_nc

  ! ------------------------------------------------------------------

  !     NAME
  !         get_time_vector_and_select


  !     PURPOSE
  !         Extract time vector in unit julian hours and get supposed time step in hours

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Oct 2012
  !         Modified  Matthias Cuntz & Juliane Mai Nov 2014 - time int or double
  !                   Stephan Thober               Sep 2015 - added read for hourly data
  !                   Robert Schweppe              Nov 2017 - restructured routine, reads vector now
  !                   Maren Kaluza                 May 2018 - fixed bug in time reading


   subroutine get_time_vector_and_select(var, fname, inctimestep, periode, time_start, time_cnt)
    !
    use mo_kind,         only: i4, i8, dp
    use mo_julian,       only: julday, caldat, dec2date
    use mo_message,      only: message
    use mo_netcdf,       only: NcVariable
    use mo_string_utils, only: DIVIDE_STRING
    use mo_common_variables, only: period
    use mo_constants,    only: DaySecs, DayHours, YearDays
    !
    implicit none
    ! intent(in)
    type(NcVariable),             intent(in)  :: var           ! variable of interest
    character(256),               intent(in)  :: fname         ! fname of ncfile for error message
    integer(i4),                  intent(in)  :: inctimestep   ! flag for requested time step
    type(period),                 intent(in)  :: periode       ! reference period
    ! intent(out)
    integer(i4),                  intent(out) :: time_start    ! time_start index of time selection
    integer(i4),                  intent(out) :: time_cnt      ! time_count of indexes of time selection

    ! local
    integer(i4)                               :: yRef, dRef, mRef, hRef, jRef ! reference time of NetCDF
    character(256)                            :: AttValues                    ! netcdf attribute values
    character(256), dimension(:), allocatable :: strArr, date, time           ! dummies for netcdf attribute handling
    integer(i8)                               :: time_step_seconds            ! native time step converter in ncfile
    integer(i8), allocatable, dimension(:)    :: time_data                    ! time vector
    type(period)                              :: nc_period                    ! period of ncfile

    ! time helpers
    integer(i4)                               :: ncJulSta1, dd, n_time
    integer(i4)                               :: mmcalstart, mmcalend, yycalstart, yycalend
    integer(i4)                               :: mmncstart, yyncstart
    integer(i4)                               :: hstart_int, hend_int ! helper variable for error output
    character(256)                            :: error_msg            ! helper variable for error output

    call var%getAttribute('units', AttValues)
    ! AttValues looks like "<unit> since YYYY-MM-DD[ HH:MM:SS]"
    ! split at space
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)

    ! determine reference time at '-' and convert to integer
    call DIVIDE_STRING(trim(strArr(3)), '-', date)
    read(date(1),*) yRef
    read(date(2),*) mRef
    read(date(3),*) dRef

    jRef = julday(dd=dRef, mm=mRef, yy=yRef)

    ! if existing also read in the time (only hour so far)
    hRef = 0
    if( size(strArr) .gt. 3 ) then
      call DIVIDE_STRING(trim(strArr(4)), ':', time)
      read(time(1),*) hRef
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
       stop
    end if

    ! get the time vector
    call var%getData(time_data)
    ! convert array from units since to seconds
    time_data = time_data * time_step_seconds

    ! check for length of time vector, needs to be at least of length 2, otherwise step width check fails
    if (size(time_data) .le. 1) then
      call message('***ERROR: length of time dimension needs to be at least 2 in file: '//trim(fname))
      stop
    end if

    ! check for equal timesteps and timestep must not be multiple of native timestep
    error_msg = '***ERROR: time_steps are not equal over all times in file and/or do not conform to'// &
                     ' requested timestep in file ('//trim(fname)//') : '

    ! compare the read period from ncfile to the period required
    ! convert julian second information back to date via conversion to float
    ! the 0.5_dp is for the different reference of fractional julian days, hours are truncated
    n_time = size(time_data)
    call dec2date(time_data(1) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dStart, nc_period%mStart, &
            nc_period%yStart, hstart_int)
    nc_period%julStart=int(time_data(1) / DaySecs + jRef + hRef / 24._dp)
    call dec2date(time_data(n_time) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, nc_period%dEnd,   nc_period%mEnd, &
            nc_period%yEnd, hend_int)
    nc_period%julEnd=int(time_data(n_time) / DaySecs + jRef + hRef / 24._dp)

    ! prepare the selection and check for required time_step
    select case(inctimestep)
    case(-1) ! daily
      ! difference must be 1 day
      if (.not. all(abs((time_data(2:n_time) - time_data(1:n_time-1)) / DaySecs - 1._dp) .lt. 1.e-6 )) then
        call message(error_msg//trim('daily'))
        stop
      end if
      ncJulSta1 = nc_period%julStart
      time_start = periode%julStart-ncJulSta1 + 1_i4
      time_cnt = periode%julEnd - periode%julStart + 1_i4
    case(-2) ! monthly
      ! difference must be between 28 and 31 days
      if ( any(abs((time_data(2:n_time) - time_data(1:n_time-1)) / DaySecs) .gt. 31._dp) .or. &
          any(abs((time_data(2:n_time) - time_data(1:n_time-1)) / DaySecs) .lt. 28._dp)) then
        call message(error_msg//trim('monthly'))
        stop
      end if

      call caldat(periode%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! monthly timesteps are usually set by month end, so for beginning, we need 1st of month
      ncJulSta1 = julday(1, mmncstart, yyncstart)
      call caldat(periode%julEnd, dd, mmcalend, yycalend)
      time_start = (yycalstart*12+mmcalstart) - (yyncstart*12+mmncstart) + 1_i4
      time_cnt = (yycalend*12+mmcalend) - (yycalstart*12+mmcalstart) + 1_i4
    case(-3) ! yearly
      ! difference must be between 365 and 366 days
      if ( any(abs((time_data(2:n_time) - time_data(1:n_time-1)) / DaySecs) .gt. (YearDays + 1._dp)) .or. &
      any(abs((time_data(2:n_time) - time_data(1:n_time-1)) / DaySecs) .lt. YearDays)) then
        call message(error_msg//'yearly')
        stop
      end if
      call caldat(periode%julStart, dd, mmcalstart, yycalstart)
      call caldat(nc_period%julStart, dd, mmncstart, yyncstart)
      ! yearly timesteps are usually set by year end, so for beginning, we need 1st of year
      ncJulSta1 = julday(1, 1, yyncstart)
      call caldat(periode%julEnd, dd, mmcalend, yycalend)
      time_start = yycalstart - yyncstart + 1_i4
      time_cnt = yycalend - yycalstart + 1_i4
    case(-4) ! hourly
      ! difference must be 1 hour
      if (.not. all(abs((time_data(2:n_time) - time_data(1:n_time-1) ) / 3600._dp - 1._dp) .lt. 1.e-6 )) then
        call message(error_msg//'hourly')
        stop
      end if
      ncJulSta1 = nc_period%julStart
      time_start = (periode%julStart-ncJulSta1) * 24_i4 + 1_i4 ! convert to hours; always starts at one
      time_cnt = (periode%julEnd - periode%julStart + 1_i4) * 24_i4 ! convert to hours
    case default ! no output at all
       call message('***ERROR: read_forcing_nc: unknown nctimestep switch.')
       stop
    end select

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 .LE. periode%julStart) .AND. (nc_period%julEnd .GE. periode%julEnd))) then
       call message('***ERROR: read_forcing_nc: time period of input data: ', trim(fname), &
            '          is not matching modelling period.')
       stop
    end if

   end subroutine get_time_vector_and_select

end module mo_read_forcing_nc

