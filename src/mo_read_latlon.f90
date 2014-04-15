!> \file mo_read_latlon.f90

!> \brief reading latitude and longitude coordinates for each basin

!> \authors Stephan Thober
!> \date Nov 2013

MODULE mo_read_latlon

  ! This module provides routines for reading latitude and longitude coordinates
  ! from file.

  ! Written  Stephan Thober, Nov 2013

  USE mo_kind, ONLY: i4, dp

  ! Of course
  IMPLICIT NONE

  PUBLIC :: read_latlon

  PRIVATE

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_latlon

  !     PURPOSE
  !>        \brief reads latitude and longitude coordinates

  !>        \details reads latitude and longitude coordinates from
  !>        netcdf file for each basin and appends it to the global
  !>        variables latitude and longitude.

  !     CALLING SEQUENCE
  !         call read_latlon(ii)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: ii"        basin index

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
  !>        File name of the basins must be xxx_latlon.nc, where
  !>        xxx is the basin id. Variable names in the netcdf file
  !>        have to be 'lat' for latitude and 'lon' for longitude.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date   Nov 2013

  subroutine read_latlon(ii)
    
    USE mo_global_variables, ONLY: dirLatLon, latitude, longitude, level1
    USE mo_append,           ONLY: append
    USE mo_message,          ONLY: message
    USE mo_ncread,           ONLY: get_NcVar, get_NcDim

    implicit none

    integer(i4), intent(in) :: ii ! basin index

    ! local variables
    character(256)                        :: fname ! file name
    integer(i4), dimension(5)             :: dl    ! dimension lengths
    real(dp), dimension(:,:), allocatable :: dummy ! dummy variable

    ! construct filename
    fname = trim( dirLatLon(ii) ) 

    ! read dimension length of variable in netcdf File
    dl = get_NcDim( trim(fname), 'lat' )
    
    ! consistency check
    if ( (dl(1) /= level1%nrows(ii) ) .or. &
         (dl(2) /= level1%ncols(ii) ) ) then
         call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
         stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!'
    end if
    
    allocate(dummy(dl(1), dl(2)))
    
    ! read dummy variable
    call get_NcVar( trim(fname), 'lat', dummy )

    ! append it to global variables
    call append( latitude, pack( dummy, .true. ))
    deallocate(dummy)

    ! read dimension length of variable in netcdf File
    dl = get_NcDim( trim(fname), 'lon' )
    
    ! consistency check
    if ( (dl(1) /= level1%nrows(ii) ) .or. &
         (dl(2) /= level1%ncols(ii) ) ) then
         call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
         stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!'
    end if
    
    allocate(dummy(dl(1), dl(2)))
    
    ! read dummy variable
    call get_NcVar( trim(fname), 'lon', dummy )

    ! append it to global variables
    call append( longitude, pack( dummy, .true. ))
    deallocate(dummy)

  end subroutine read_latlon

END MODULE mo_read_latlon
