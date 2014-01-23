!> \file mo_read_lai.f90

!> \brief MODIS LAI as a field of data.

!> \details This module is to read MODIS LAI input data.
!> The module provides a subroutine for NetCDF files.\n
!> First, the dimensions given are cross-checked with header.txt information. Second, the data of the specified period are
!> read from the specified directory. The names of files in this directory have to be always "YYYY.bin" or "YYYY.nc".\n
!> If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
!> The program is stopped if any value lies out of range.

!> \authors John Craven
!> \date August 2013

MODULE mo_read_lai

  ! This module provides routines to read MODIS data.

  ! Written  John Craven, Aug 2013

  USE mo_kind,   ONLY: i4, sp, dp
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_lai_bin   ! Read bin files
  PUBLIC :: read_lai_nc    ! Read netCDF files

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_lai_bin

  !     PURPOSE
  !>        \brief Reads binary gridded fields of LAI.

  !>        \details Reads binary LAI gridded field files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        The names of files in this directory have to be always "YYYY.bin" (file ending .bin can be set in mo_file).\n
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.

  !     CALLING SEQUENCE
  !         periode%dStart  = 2_i4                                                     !           day
  !         periode%mStart  = 2_i4                                                     !           month
  !         periode%yStart  = 1972_i4                                                  !           year
  !         periode%dEnd    = 7_i4                                                     !           day
  !         periode%mEnd    = 8_i4                                                     !           month
  !         periode%yEnd    = 1977_i4                                                  !           year
  !         periode%julStart  = NDAYS(periode%dStart, periode%mStart, periode%yStart)  !           julian day starting
  !         periode%julEnd    = NDAYS(periode%dEnd,   periode%mEnd,   periode%yEnd  )  !           julian day ending
  !         periode%nObs    = periode%julEnd - periode%julStart + 1_i4                 !           total number of observations

  !         call read_meteo_bin('old_code/sub_00020/input/lai/' , 21_i4, 28_i4, periode, lai,   &
  !                             lower=  0.0_dp)

  !     INTENT(IN)
  !>        \param[in] "character(len=*)        :: folder"   Name of the folder where data are stored
  !>        \param[in] "integer(i4)             :: nRows"    Number of datapoints in longitudinal direction
  !>        \param[in] "integer(i4)             :: nCols"    Number of datapoints in latitudinal  direction
  !>        \param[in] "type(period)            :: periode"  Period the data are needed for
  !>        \param[in] "logical, dimension(:,:) :: mask"     Mask of valid field for checking lower and upper bounds

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:,:,:) :: data"     Data matrix
  !>                                                             dim_1 = longitude, dim_2 = latitude, dim_3 = time

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), optional)           :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional)           :: upper"    Upper bound for check of validity of data values

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>        \note The values stored in YYYY.bin are of single precision, i.e. 4 bytes per value.

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author J. Craven & R. Kumar
  !>        \date Aug. 2013

  subroutine read_lai_bin(folder, nRows, nCols, periode, data, mask, lower, upper)

    use mo_global_variables, only: period
    use mo_julian,           only: date2dec
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_file,             only: file_lai_header, ulai_header, &
                                   file_lai_binary_end, ulai
    implicit none

    character(len=*),                                  intent(in)  :: folder    ! folder where data are stored
    integer(i4),                                       intent(in)  :: nRows     ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4),                                       intent(in)  :: nCols     ! number of columns of data fields:
    ! LATITUDE dimension
    type(period),                                      intent(in)  :: periode   ! time period
    logical, dimension(:,:),                           intent(in)  :: mask      ! mask of valid data fields
    real(dp), dimension(:,:,:), allocatable,           intent(out) :: data      ! data read in
    real(dp),                                optional, intent(in)  :: lower     ! lower bound for data points
    real(dp),                                optional, intent(in)  :: upper     ! upper bound for data points

    ! local variables
    integer(i4)                           :: i
    character(10)                         :: dummy
    integer(i4)                           :: nc, nr
    character(4)                          :: yearStr
    integer(i4)                           :: julStartYear, julEndYear       ! start and end in julianDay of year
    integer(i4)                           :: julStartPeriod, julEndPeriod   ! start and end in julianDay of Period
    integer(i4)                           :: cummulativeDays      ! which day will be read in
    integer(i4)                           :: day                  ! counter: actual day
    integer(i4)                           :: year                 ! counter: actual year
    integer(i4)                           :: record_number        ! counter: record number in file
    real(sp), dimension(:,:), allocatable :: tmp                  ! data of a single day
    real(sp)                              :: nodata_value
    character(256)                        :: fName
    real(dp)                              :: jday_frac

    ! checking input nCols and nRows with header file of data files stored under folder/header.txt
    fName = trim(adjustl(folder) )//trim(adjustl(file_lai_header))
    open(unit=ulai_header, file=trim(fName), status='old')
    read(ulai_header, *) dummy, nc
    read(ulai_header, *) dummy, nr
    read(ulai_header, *) dummy
    read(ulai_header, *) dummy
    read(ulai_header, *) dummy
    read(ulai_header, *) dummy, nodata_value
    dummy = dummy//''   ! only to avoid warning

    close(ulai_header)
    if ( (nc .ne. nRows) .or. (nr .ne. nCols) ) then
       stop 'read_lai_bin: mHM generated nRows and nCols do not match header.txt information'
    end if

    allocate( data(nRows, nCols, periode%julEnd - periode%julStart + 1_i4) )

    cummulativeDays = 0_i4
    yearLoop: do year=periode%yStart, periode%yEnd

       write(yearStr, '(I4)') year
       fName = trim(folder) // trim(yearStr) // trim(file_lai_binary_end)
       open(unit=ulai, file=trim(fName), form='unformatted', access='direct', recl=4*nRows*nCols)

       ! julian day of strating and end year
       jday_frac    = date2dec(dd=01, mm=01, yy=year)
       julStartYear = nint(jday_frac)

       jday_frac   = date2dec(dd=31, mm=12, yy=year)
       julEndYear  = nint(jday_frac)


       ! Julian Days for the read-in period within this year
       ! - same as above in intermediate years
       ! - different if start or end of period is not first and last day of the year
       julStartPeriod   = Max(julStartYear, periode%julStart)
       julEndPeriod     = Min(julEndYear,   periode%julEnd)
       ! DaysOfYear = julEndYear - julStartYear + 1_i4

       dayLoop: do day = julStartPeriod, julEndPeriod
          cummulativeDays = cummulativeDays + 1_i4
          record_number   = day - julStartYear + 1_i4
          allocate(tmp(nCols,nRows))
          read(ulai,rec=record_number) tmp
          data(:,:,cummulativeDays) = real(transpose(tmp(:,:)), dp)
          deallocate(tmp)
       end do dayLoop

       close(ulai)

    end do yearLoop

    ! start checking values
    do i = 1, (periode%julEnd - periode%julStart + 1_i4)
       ! for no data value
       if(  any(                                                                  &
            ( abs(data(:,:,i) - nodata_value) .LT. tiny(1.0_dp) )  .AND.       &
            mask(:,:)                                                             &
            )                                                                     &
            ) then
          call message('***ERROR: nodata value within basin boundary in folder ', &
               folder, trim( num2str(year,'(I4.4)') )          )
       end if


       ! optinal check
       if (present(lower)) then
          if(  any( (data(:,:,i) .lt. lower) .AND. mask(:,:) )  ) then
             call message('read_lai_bin: ERROR occured: values in file "',     &
                  folder, trim(num2str(year,'(I4.4)')), file_lai_binary_end,   &
                  '" are lower than ', trim(num2str(lower,'(F7.2)')) )
             stop
          end if
       end if

       if (present(upper)) then
          if(  any( (data(:,:,i) .gt. upper) .AND. mask(:,:) )  ) then
             call message('read_lai_bin: ERROR occured: values in file "', &
                  folder, trim(num2str(year,'(I4.4)')), file_lai_binary_end, &
                  '" are greater than ', trim(num2str(upper,'(F7.2)')) )
             stop
          end if
       end if

    end do

  end subroutine read_lai_bin
  
  ! ------------------------------------------------------------------

  !     NAME
  !         read_lai_nc

  !     PURPOSE
  !>        \brief Reads MODIS LAI input in NetCDF file format.

  !>        \details Reads netCDF MODIS LAI files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        The names of files in this directory have to be always "YYYY.nc".\n
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.

  !     CALLING SEQUENCE
  !         periode%dStart   = 2_i4                                                     ! day
  !         periode%mStart   = 2_i4                                                     ! month
  !         periode%yStart   = 1972_i4                                                  ! year
  !         periode%dEnd     = 7_i4                                                     ! day
  !         periode%mEnd     = 8_i4                                                     ! month
  !         periode%yEnd     = 1977_i4                                                  ! year
  !         periode%julStart = NDAYS(periode%dStart, periode%mStart, periode%yStart)    ! julian day starting
  !         periode%julEnd   = NDAYS(periode%dEnd,   periode%mEnd,   periode%yEnd  )    ! julian day ending
  !         periode%nObs     = periode%julEnd - periode%julStart + 1_i4                 ! total number of observations

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
  !>        \param[in] "logical, dimension(:,:),optional :: mask"     mask of valid field for checking bounds
  !>        \param[in] "real(dp), optional               :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional               :: upper"    Upper bound for check of validity of data values

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>        \note Files have to be called as defined in mo_files. Furthermore the variable names have to be called
  !>              like they are defined in the declaration of this subroutine. The NetCDF file has to have 3 dimensions:
  !>              1. x, 2. y, 3. t. It is expected that the variables (especially)within the NetCDF files contain an 
  !>              unit attribute. The timestep has to be equidistant 

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY 
  !>        \author John Craven
  !>        \date August 2013

  subroutine read_lai_nc(folder, nRows, nCols, periode, varName, laidata, mask, lower, upper)

    use mo_global_variables, only: period
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_ncread,           only: Get_NcDim, Get_NcVar, Get_NcVarAtt

    implicit none

    character(len=*),                                  intent(in)  :: folder    ! folder where data are stored
   
    integer(i4),                                       intent(in)  :: nRows     ! number of rows of data fields:
    integer(i4),                                       intent(in)  :: nCols     ! number of columns of data fields:

    type(period),                                      intent(in)  :: periode   ! time period
    character(len=*),                                  intent(in)  :: varName   ! name of NetCDF variable
    real(dp), dimension(:,:,:), allocatable,           intent(out) :: laidata      ! data read in
    logical, dimension(:,:),                           intent(in)  :: mask      ! mask of valid data fields
    real(dp),                                optional, intent(in)  :: lower     ! lower bound for data points
    real(dp),                                optional, intent(in)  :: upper     ! upper bound for data points
    !
    ! local variables
    character(256)                                                 :: fName     ! name of NetCDF file
    character(256)                                                 :: AttValues ! netcdf attribute values
    integer(i4)                                                    :: i         ! loop variable
    integer(i4)                                                    :: ncJulSta  ! start time of nc dataset
    integer(i4)                                                    :: ncJulEnd  ! end time of nc dataset
    integer(i4)                                                    :: datatype  ! datatype of attribute
    integer(i4), dimension(5)                                      :: dimen     ! dimension for NetCDF file
    real(dp)                                                       :: nodata_value ! data nodata value
    real(dp), dimension(:,:,:), allocatable                        :: tmp          ! data read in for all nctimesteps

    fName = trim(folder) // trim(varName) // '.nc'
    
    ! get dimensions
    dimen = Get_NcDim( trim(fName), trim(varName) )

    if ( (dimen(1) .ne. nRows) .or. (dimen(2) .ne. nCols) ) then
       stop '***ERROR: read_lai_nc: mHM generated x and y do not match NetCDF dimensions'
    end if
    !
    allocate(     tmp(dimen(1), dimen(2),  dimen(3)                               ) )
    allocate( laidata(dimen(1), dimen(2), periode%julEnd - periode%julStart + 1_i4) )

    ! determine no data value           
    call Get_NcVarAtt( trim(fName), trim(varName), '_FillValue', AttValues, dtype=datatype)

    ! convert to number
    read(AttValues, *) nodata_value 

    call Get_NcVar(trim(fName), trim(varName), tmp)

    ! get time interval & check time steps
    call get_time(fName, varName, ncJulSta, ncJulEnd)

    ! check if time periods overlay, put data of relevant time period 
    if ((ncJulSta .LE. periode%julStart) .AND. (ncJulEnd .GT. periode%julEnd)) then
       laidata(:,:,:) = tmp(:,:,periode%julStart-ncJulSta+1:periode%julEnd-ncJulSta+1)
    else
       call message('***ERROR: read_lai_nc: time period of input data: ', trim(varName), &
                    ' is not matching modelling period.')
       stop
    end if

    ! start checking values
    do i = 1, (periode%julEnd - periode%julStart + 1_i4)
       ! for no data value
       if(  any(                                                              &
            ( abs(laidata(:,:,i) - nodata_value) .LT. tiny(1.0_dp) )  .AND.   &
            mask(:,:)                                                         &
            )                                                                 &
            ) then
          call message('***ERROR: read_lai_nc: nodata value within basin ')
          call message('          boundary in variable: ', trim(varName))
          stop
       end if

       ! optional check
       if (present(lower)) then
          if(  any( (laidata(:,:,i) .lt. lower) .AND. mask(:,:) )  ) then
             call message('***ERROR: read_lai_nc: values in variable "', &
                  trim(varName),                                           &
                  '" are lower than ', trim(num2str(lower,'(F7.2)')) )
             stop
          end if
       end if

       if (present(upper)) then
          if(  any( (laidata(:,:,i) .gt. upper) .AND. mask(:,:) )  ) then
             call message('***ERROR: read_lai_nc: values in variable"',  &
                  trim(varName),                                           &
                  '" are greater than ', trim(num2str(upper,'(F7.2)')) )
             stop
          end if
       end if

    end do

    deallocate(tmp)

  end subroutine read_lai_nc

 !
  !     PORPOSE
  !         Determine data time interval & check timesteps
  
  !     CALLING SEQUENCE
  !         BIAS(NetCDF_filename, VariableName, startJulianDay, endJulianDay)
  
  !     RESTRICTIONS
  !         Input values must be floating points.
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !         Written,  Matthias Zink, Oct 2012

  subroutine get_time(fName, vName, julStart, julEnd)
    !
    use mo_julian,       only: date2dec
    use mo_message,      only: message
    use mo_NcRead,       only: Get_NcVar, Get_NcDim, Get_NcVarAtt
    use mo_string_utils, only: DIVIDE_STRING
    !
    implicit none
    !
    character(len=*)            , intent(in)  :: fName               ! name of NetCDF file
    character(len=*)            , intent(in)  :: vName               ! name of variable
    integer(i4)                 , intent(out) :: julStart 
    integer(i4)                 , intent(out) :: julEnd
    !
    integer(i4)                               :: i
    integer(i4)                               :: deltaT              ! diff between single time steps in NetCDF
    integer(i4)                               :: yRef, dRef, mRef    ! reference time of NetCDF (unit attribute of
    integer(i4)                               :: datatype            ! datatype of attribute
    integer(i4),    dimension(5)              :: dimen 
    !
    integer(i4),   dimension(:), allocatable  :: timesteps           ! time variable of NetCDF
    !
    character(256)                            :: AttValues           ! netcdf attribute values
    character(256), dimension(:), allocatable :: strArr              ! dummy for netcdf attribute handling
    real(dp)                                 :: jday_frac            ! julian day from dec2date
    !
    dimen = Get_NCDim(fName, trim(vName))
    ! get unit attribute of variable 'time'
    call Get_NcVarAtt(fName, 'time', 'units', AttValues, dtype=datatype)
    ! AttValues looks like "<unit> since YYYY-MM-DD HH:MM:SS"
    call DIVIDE_STRING(trim(AttValues), ' ', strArr)
    !
    ! strArr(1) is <unit>
    if (strArr(1) .EQ. 'days') then
       ! determine reference time and convert to integer
       call DIVIDE_STRING(trim(strArr(3)), '-', strArr) 
       read(strArr(1),*) yRef
       read(strArr(2),*) mRef
       read(strArr(3),*) dRef
       !
       allocate(timesteps(dimen(3)))
       call Get_NcVar(fName, 'time', timesteps)
       !
       ! check if timestep is one day
       do i = 2, dimen(3)
          deltaT = timesteps(i) - timesteps(i-1)
          ! deltaT has to be one but with inetger conversion 1000
          if ( deltaT .NE. 1_i4) then
             call message('***ERROR: Timestep has to be equidistant as one day in ', trim(vName))
             stop
          end if
       end do
       !
       ! determine starting and ending julian day of the dataset
       jday_frac = date2dec(dd=dRef, mm=mRef, yy=yRef)
       i  = nint(jday_frac, i4 ) 
       julStart = i + timesteps(1)       
       julEnd   = i + timesteps(dimen(3))
    else 
       call message('***ERROR: Please provide the input data ', trim(vName) , 'in days.')
       stop
    end if
    !
  end subroutine get_time

END MODULE mo_read_lai
