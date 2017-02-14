!> \file mo_read_meteo.f90

!> \brief Reads meteorological input data.

!> \details This module is to read meteorological input data, e.g. temperature, precipitation.\n
!> The module provides a subroutine for binary files.\n
!> First, the dimensions given are cross-checked with header.txt information. Second, the data of the specified period are
!> read from the specified directory. The names of files in this directory have to be always "YYYY.bin".\n
!> If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
!> The program is stopped if any value lies out of range.

!> \authors Juliane Mai
!> \date Dec 2012
!  Modified Sep 2015, Stephan Thober - separated routines for netcdf files from routines for binary files

MODULE mo_read_meteo

  ! This module provides routines to read meteorological data.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Rohini Kumar, Feb 2013
  !          Stephan Thober, Sep 2015 - moved read of meteorological netcdf files to
  !                                     mo_read_forcing_nc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_meteo_bin   ! Read binary files

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_meteo_bin

  !     PURPOSE
  !>        \brief Reads binary meteorological files.

  !>        \details Reads binary meteorological files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        The names of files in this directory have to be always "YYYY.bin" (file ending .bin can be set in mo_file).\n
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.

  !     CALLING SEQUENCE
  !         periode%dStart  = 2_i4                                                    !           day
  !         periode%mStart  = 2_i4                                                    !           month
  !         periode%yStart  = 1972_i4                                                 !           year
  !         periode%dEnd    = 7_i4                                                    !           day
  !         periode%mEnd    = 8_i4                                                    !           month
  !         periode%yEnd    = 1977_i4                                                 !           year
  !         periode%julStart  = NDAYS(periode%dStart, periode%mStart, periode%yStart)   !           julian day starting
  !         periode%julEnd    = NDAYS(periode%dEnd,   periode%mEnd,   periode%yEnd  )   !           julian day ending
  !         periode%nObs    = periode%julEnd - periode%julStart + 1_i4                    !           total number of observations

  !         call read_meteo_bin('old_code/sub_00020/input/meteo/pre/' , 21_i4, 28_i4, periode, precipitation,   &
  !                             lower=  0.0_dp)
  !         call read_meteo_bin('old_code/sub_00020/input/meteo/pet/' , 21_i4, 28_i4, periode, pot_evapo_trans, &
  !                             lower=  0.0_dp, upper=20.0_dp)
  !         call read_meteo_bin('old_code/sub_00020/input/meteo/tavg/', 21_i4, 28_i4, periode, temp_average,    &
  !                             lower=-50.0_dp, upper=50.0_dp)

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
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, Stephan Thober, Jun 2014 -- added julstart and julend

  subroutine read_meteo_bin(folder, nRows, nCols, periode, data, mask, &
       lower, upper )

    use mo_kind,             only: i4, dp, sp
    use mo_common_variables, only: period
    use mo_julian,           only: date2dec
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_file,             only: file_meteo_header, umeteo_header, &
         file_meteo_binary_end, umeteo

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
    fName = trim(adjustl(folder) )//trim(adjustl(file_meteo_header))
    open (unit=umeteo_header, file=trim(fName), status='old')
    read (umeteo_header, *) dummy, nc
    read (umeteo_header, *) dummy, nr
    read (umeteo_header, *) dummy
    read (umeteo_header, *) dummy
    read (umeteo_header, *) dummy
    read (umeteo_header, *) dummy, nodata_value
    dummy = dummy//''   ! only to avoid warning

    close(umeteo_header)
    if ( (nc .ne. nRows) .or. (nr .ne. nCols) ) then
       stop 'read_meteo_bin: mHM generated nRows and nCols are not matching header.txt information'
    end if

    allocate( data(nRows, nCols, periode%julEnd - periode%julStart + 1_i4) )

    cummulativeDays = 0_i4
    yearLoop: do year=periode%yStart, periode%yEnd

       write(yearStr, '(I4)') year
       fName = trim(folder) // trim(yearStr) // trim(file_meteo_binary_end)
       open(unit=umeteo, file=trim(fName), &
            form='unformatted', access='direct', recl=4*nRows*nCols)

       ! julian day of starting and end year
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
          read(umeteo,rec=record_number) tmp
          data(:,:,cummulativeDays) = real(transpose(tmp(:,:)), dp)
          deallocate(tmp)
       end do dayLoop

       close(umeteo)

    end do yearLoop

    ! start checking values
    do i = 1, size( data, dim = 3 )
       ! for no data value
       if(  any(                                                                  &
            ( abs(data(:,:,i) - nodata_value) .LT. tiny(1.0_dp) )  .AND.       &
            mask(:,:)                                                             &
            )                                                                     &
            ) then
          call message('***ERROR: nodata value within basin boundary in timestep: ', &
               trim( num2str(i,'(I5.5)') )          )
       end if


       ! optinal check
       if (present(lower)) then
          if(  any( (data(:,:,i) .lt. lower) .AND. mask(:,:) )  ) then
             call message('read_meteo_bin: ERROR occured: values at timestep: ',     &
                  trim(num2str(i,'(I5.5)')), &
                  ' are lower than ', trim(num2str(lower,'(F7.2)')) )
             stop
          end if
       end if

       if (present(upper)) then
          if(  any( (data(:,:,i) .gt. upper) .AND. mask(:,:) )  ) then
             call message('read_meteo_bin: ERROR occured: values at timestep ', &
                  trim(num2str(i,'(I5.5)')), &
                  ' are greater than ', trim(num2str(upper,'(F7.2)')) )
             stop
          end if
       end if

    end do

  end subroutine read_meteo_bin

END MODULE mo_read_meteo
