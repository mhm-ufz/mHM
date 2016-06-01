!> \file mo_read_forcing_nc.f90

!> \brief Reads forcing input data.

!> \details This module is to read forcing input data contained in netcdf files, e.g. temperature, precipitation,
!> total_runoff, lai. Timesteps can be hourly, daily, monthly, and annual. The module provides a subroutine
!> for NetCDF files only. First, the dimensions given are cross-checked with header.txt information. Second,
!> the data of the specified period are read from the specified directory. The names of files in this directory
!> have to be always "<var_name.nc".\n
!> If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
!> The program is stopped if any value lies out of range.

!> \authors Juliane Mai
!> \date Dec 2012
!  Modified Sep 2015, Stephan Thober - separated routines for netcdf files from routines for binary files

module mo_read_forcing_nc
  implicit none
  public :: read_forcing_nc
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

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:,:,:) :: data"     Data matrix
  !>                                                             dim_1 = longitude, dim_2 = latitude, dim_3 = time

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), optional, intent(in)  :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional, intent(in)  :: upper"    Upper bound for check of validity of data values
  !>        \param[in] "logical,  optional, intent(in)  :: nocheck"  .TRUE. if check for nodata values deactivated
  !>                                                                  default = .FALSE. - check is done
  
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
  !>              unit attribute. The timestep has to be equidistant

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

  subroutine read_forcing_nc(folder, nRows, nCols, periode, varName, data, mask, lower, upper, nctimestep, nocheck, maskout)

    use mo_kind,             only: i4, dp
    use mo_common_variables, only: period
    use mo_julian,           only: caldat, julday
    use mo_message,          only: message
    use mo_ncread,           only: Get_NcDim, Get_NcVar, Get_NcVarAtt

    use mo_string_utils,      only: num2str
    use mo_utils,             only: eq, ne

    implicit none

    character(len=*),                                  intent(in)  :: folder    ! folder where data are stored
    integer(i4),                                       intent(in)  :: nRows     ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4),                                       intent(in)  :: nCols     ! number of columns of data fields:
    ! LATITUDE dimension
    type(period),                                      intent(in)  :: periode   ! time period
    character(len=*),                                  intent(in)  :: varName   ! name of NetCDF variable
    logical, dimension(:,:),                           intent(in)  :: mask      ! mask of valid data fields
    real(dp),                                optional, intent(in)  :: lower     ! lower bound for data points
    real(dp),                                optional, intent(in)  :: upper     ! upper bound for data points
    integer(i4),                             optional, intent(in)  :: nctimestep ! -1: daily (default); -2:monthly; -3:yearly
    logical,                                 optional, intent(in)  :: nocheck    ! .TRUE. if check for nodata values deactivated
    !                                                                            ! default = .FALSE. - check is done
    real(dp), dimension(:,:,:), allocatable,           intent(out) :: data      ! data read in
    logical, dimension(:,:,:), allocatable,  optional, intent(out) :: maskout    ! mask of data to read

    !
    ! local variables
    character(256)            :: fName        ! name of NetCDF file
    character(256)            :: AttValues    ! netcdf attribute values
    integer(i4)               :: i            ! loop variable
    integer(i4)               :: ncJulSta     ! start time of nc dataset
    integer(i4)               :: ncJulEnd     ! end time of nc dataset
    integer(i4)               :: datatype     ! datatype of attribute
    integer(i4), dimension(5) :: dimen        ! dimension for NetCDF file
    real(dp)                  :: nodata_value ! data nodata value
    integer(i4)               :: dim3         ! time dim of input data
    integer(i4)               :: ncdim3start  ! start of reading in nc file
    integer(i4)               :: inctimestep  ! local nctimestep
    logical                   :: checking     ! check if model domain is covered by data
    ! time helpers
    integer(i4)               :: julStart1, julEnd1, ncJulSta1, ncJulEnd1, dd
    integer(i4)               :: mmcalstart, mmcalend, yycalstart, yycalend
    integer(i4)               :: mmncstart, mmncend, yyncstart, yyncend

    ! check optional nctimestep
    inctimestep = -1
    if (present(nctimestep)) inctimestep = nctimestep

    checking = .TRUE.
    if (present(nocheck)) checking = .NOT. nocheck
    
    fName = trim(folder) // trim(varName) // '.nc'
    ! get dimensions
    dimen = Get_NcDim(trim(fName), trim(varName))

    if ( (dimen(1) .ne. nRows) .or. (dimen(2) .ne. nCols) ) then
       stop '***ERROR: read_forcing_nc: mHM generated x and y are not matching NetCDF dimensions'
    end if

    ! determine no data value
    call Get_NcVarAtt(fName, varName, '_FillValue', AttValues, dtype=datatype)
    ! convert to number
    read(AttValues, *) nodata_value

    ! get time intervall & check time steps
    call get_time(fName, varName, ncJulSta, ncJulEnd, nctimestep=inctimestep)

    !
    select case(inctimestep)
    case(-1) ! daily
       julStart1 = periode%julStart
       ncJulSta1 = ncJulSta
       julEnd1   = periode%julEnd
       ncJulEnd1 = ncJulEnd
       dim3 = periode%julEnd - periode%julStart + 1_i4
       ncdim3start = periode%julStart-ncJulSta+1
    case(-2) ! monthly
       call caldat(periode%julStart, dd, mmcalstart, yycalstart)
       julStart1 = julday(1, mmcalstart, yycalstart)
       call caldat(ncJulSta, dd, mmncstart, yyncstart)
       ncJulSta1 = julday(1, mmncstart, yyncstart)
       call caldat(periode%julEnd, dd, mmcalend, yycalend)
       julEnd1 = julday(1, mmcalend, yycalend)
       call caldat(ncJulEnd, dd, mmncend, yyncend)
       ncJulEnd1 = julday(1, mmncend, yyncend)
       dim3 = (yycalend*12+mmcalend) - (yycalstart*12+mmcalstart) + 1_i4
       ncdim3start = (yycalstart*12+mmcalstart) - (yyncstart*12+mmncstart) + 1_i4
    case(-3) ! yearly
       call caldat(periode%julStart, dd, mmcalstart, yycalstart)
       julStart1 = julday(1, 1, yycalstart)
       call caldat(ncJulSta, dd, mmncstart, yyncstart)
       ncJulSta1 = julday(1, 1, yyncstart)
       call caldat(periode%julEnd, dd, mmcalend, yycalend)
       julEnd1 = julday(1, 1, yycalend)
       call caldat(ncJulEnd, dd, mmncend, yyncend)
       ncJulEnd1 = julday(1, 1, yyncend)
       dim3 = yycalend - yycalstart + 1_i4
       ncdim3start = yycalstart - yyncstart + 1_i4
    case(-4) ! hourly
       julStart1 = periode%julStart
       ncJulSta1 = ncJulSta
       julEnd1   = periode%julEnd
       ncJulEnd1 = ncJulEnd
       dim3 = (periode%julEnd - periode%julStart + 1_i4) * 24_i4 ! convert to hours
       ncdim3start = (periode%julStart-ncJulSta) * 24_i4 + 1_i4 ! convert to hours; always starts at one
    case default ! no output at all
       call message('***ERROR: read_forcing_nc: unknown nctimestep switch.')
       stop
    end select

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 .LE. julStart1) .AND. (ncJulEnd1 .GE. julEnd1))) then
       call message('***ERROR: read_forcing_nc: time period of input data: ', trim(varName), &
            '          is not matching modelling period.')
       stop
    end if

    !
    ! alloc and read
    allocate(data(dimen(1), dimen(2), dim3))
    call Get_NcVar(trim(fName), trim(varName), data, &
              start = (/ 1_i4, 1_i4, ncdim3start /), &
            a_count = (/ dimen(1), dimen(2), dim3 /) )

    ! save output mask if optional maskout is given
    if (present(maskout)) then
       allocate(maskout(dimen(1), dimen(2), dim3))
       maskout = ne(data(:,:,:),nodata_value)
    end if
       
    ! start checking values
    do i = 1, dim3
       ! neglect checking for naodata values if optional nocheck is given
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
             call message('Minval at timestep: ', trim(num2str(minval(data(:,:,i)),'(F7.2)')))
             call message('Total minval: ', trim(num2str(minval(data(:,:,:)),'(F7.2)')))
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
             call message('Maxval at timestep: ', trim(num2str(maxval(data(:,:,i)),'(F7.2)')))
             call message('Total maxval: ', trim(num2str(maxval(data(:,:,:)),'(F7.2)')))
             stop
          end if
       end if

    end do

  end subroutine read_forcing_nc

  !     PORPOSE
  !         Determine data time interval & check timesteps

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Oct 2012
  !         Modified  Matthias Cuntz & Juliane Mai Nov 2014 - time int or double
  !                   Stephan Thober Sep 2015 - added read for hourly data

  subroutine get_time(fName, vName, julStart, julEnd, nctimestep)
    !
    use mo_kind,         only: i4, dp
    use mo_julian,       only: date2dec
    use mo_message,      only: message
    use mo_NcRead,       only: Get_NcVar, Get_NcDim, Get_NcVarAtt
    use mo_utils,        only: ne
    use mo_string_utils, only: DIVIDE_STRING
    use netcdf,          only: NF90_INT, NF90_DOUBLE, NF90_NOWRITE
    use netcdf,          only: nf90_open, nf90_close, nf90_inq_varid, nf90_inquire_variable
    !
    implicit none
    !
    character(len=*)            , intent(in)  :: fName               ! name of NetCDF file
    character(len=*)            , intent(in)  :: vName               ! name of variable
    integer(i4)                 , intent(out) :: julStart
    integer(i4)                 , intent(out) :: julEnd
    integer(i4),        optional, intent(in)  :: nctimestep          ! -1: daily (default); -2:monthly; -3:yearly
    !
    integer(i4)                               :: i
    integer(i4)                               :: yRef, dRef, mRef    ! reference time of NetCDF (unit attribute of
    integer(i4)                               :: datatype            ! datatype of attribute
    integer(i4),    dimension(5)              :: dimen
    !
    integer(i4),   dimension(:), allocatable  :: timesteps           ! time variable of NetCDF
    real(dp),      dimension(:), allocatable  :: itimesteps          ! time variable of NetCDF
    !
    character(256)                            :: AttValues           ! netcdf attribute values
    character(256), dimension(:), allocatable :: strArr              ! dummy for netcdf attribute handling
    character(256), dimension(:), allocatable :: date                ! dummy for netcdf attribute handling date
    real(dp)                                  :: jday_frac           ! julian day from dec2date
    integer(i4)                               :: inctimestep         ! local nctimestep
    integer(i4) :: ncid    ! id of input stream
    integer(i4) :: varid   ! id of variable to be read
    integer(i4) :: status  ! netcdf inquire return
    real(dp) :: deltaT ! diff between single time steps in NetCDF in days
    real(dp) :: nTStepDay ! number of timesteps per day
    !
    ! check optional nctimestep
    inctimestep = -1 !ST: to do inctimestep should be determined from file
    nTStepDay = -1.
    if (present(nctimestep)) inctimestep = nctimestep

    dimen = Get_NCDim(fName, trim(vName))
    ! get unit attribute of variable 'time'
    call Get_NcVarAtt(fName, 'time', 'units', AttValues, dtype=datatype)
    ! AttValues looks like "<unit> since YYYY-MM-DD HH:MM:SS"
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)
    !
    ! determine reference time and convert to integer
    call DIVIDE_STRING(trim(strArr(3)), '-', date)
    read(date(1),*) yRef
    read(date(2),*) mRef
    read(date(3),*) dRef
    jday_frac = date2dec(dd=dRef, mm=mRef, yy=yRef)
    !
    status = nf90_open(fName, NF90_NOWRITE, ncid)
    status = nf90_inq_varid(ncid, 'time', varid)
    status = nf90_inquire_variable(ncid, varid, xtype=datatype)
    status = nf90_close(ncid)
    status = status + 0  ! this is only to make the variable used
    if ((datatype .eq. NF90_INT) .or. (datatype .eq. NF90_DOUBLE)) then
       if (datatype .eq. NF90_INT) then
          allocate(timesteps(dimen(3)))
          call Get_NcVar(fName, 'time', timesteps)
       else
          allocate(itimesteps(dimen(3)))
          allocate(timesteps(dimen(3)))
          call Get_NcVar(fName, 'time', itimesteps)
          timesteps = nint(itimesteps, i4)
       endif
    else
       call message('***ERROR: data type of time must be NF90_INT or NF90_DOUBLE in netcdf file.')
       stop
    end if
    !
    ! strArr(1) is <unit>
    if (strArr(1) .EQ. 'days') then
       nTStepDay = 1
    else if (strArr(1) .eq. 'hours') then
       nTStepDay = 24
    else
       call message('***ERROR: Please provide the input data in (days or hours) since YYYY-MM-DD HH:MM:SS in ', trim(vName))
       stop
    end if
    !
    ! check consistency of timesteps
    select case(inctimestep)
    case(-1) ! daily
       do i = 2, dimen(3)
          deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
          ! deltaT has to be one but with integer conversion
          if ( ne(deltaT, 1._dp) ) then
             call message('***ERROR: ', trim(vName),' must have daily time steps.')
             stop
          end if
       end do
    case(-2) ! monthly
       do i = 2, dimen(3)
          deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
          ! deltaT has to be one but with integer conversion
          if (( deltaT .lt. 28._dp) .or. ( deltaT .gt. 31._dp)) then
             call message('***ERROR: ', trim(vName),' must have monthly time steps.')
             stop
          end if
       end do
    case(-3) ! yearly
       do i = 2, dimen(3)
          deltaT = real(timesteps(i) - timesteps(i-1), dp) / nTStepDay
          ! deltaT has to be one but with integer conversion
          if ( deltaT .lt. 360._dp) then
             call message('***ERROR: ', trim(vName),' must have yearly time steps.')
             stop
          end if
       end do
    case(-4) ! hourly
       do i = 2, dimen(3)
          deltaT = real(timesteps(i) - timesteps(i-1), dp)
          ! deltaT has to be one but with integer conversion
          if ( ne(deltaT, 1._dp) ) then
             call message('***ERROR: ', trim(vName),' must have hourly time steps.')
             stop
          end if
       end do
    case default
       call message('***ERROR: get_time: unknown nctimestep switch.')
       stop
    end select
    !
    ! determine starting and ending julian day of the dataset
    i  = nint(jday_frac, i4 )
    ! calculate days
    julStart = i + timesteps(1) / int(nTStepDay, i4)
    julEnd   = i + timesteps(dimen(3)) / int(nTStepDay, i4)
    !
  end subroutine get_time

end module mo_read_forcing_nc

