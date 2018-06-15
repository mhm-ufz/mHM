!> \file mo_read_latlon.f90

!> \brief reading latitude and longitude coordinates for each basin

!> \authors Stephan Thober
!> \date Nov 2013

MODULE mo_read_latlon

  ! This module provides routines for reading latitude and longitude coordinates
  ! from file.

  ! Written  Stephan Thober, Nov 2013

  USE mo_kind, ONLY : i4, dp

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
  !         modified, Stephan Thober, Sep 2015 - added latitude and longitude for level 0
  !                   Stephan Thober, Oct 2015 - added L1_rect_latitude and L1_rect_longitude
  !                   David Schaefer, May 2016 - removed ncread dependency
  !                   Robert Schweppe, Mar 2018 - major rewrite

  subroutine read_latlon(ii, lon_var_name, lat_var_name, level_name, level)

    use mo_common_variables, only : fileLatLon, Grid
    USE mo_message, ONLY : message
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    integer(i4), intent(in) :: ii ! basin index
    character(*), intent(in) :: lon_var_name
    character(*), intent(in) :: lat_var_name
    character(*), intent(in) :: level_name
    type(Grid), intent(inout) :: level

    ! local variables
    character(256) :: fname ! file name
    real(dp), dimension(:, :), allocatable :: dummy ! dummy variable
    type(NcDataset) :: nc
    type(NcVariable) :: var

    ! construct filename
    fname = trim(fileLatLon(ii))

    nc = NcDataset(fname, "r")

    ! -------------------------------------------------------------------------
    ! READ LEVEL LATITUDE / LONGITUDE
    ! -------------------------------------------------------------------------
    var = nc%getVariable(trim(lat_var_name))
    call var%getData(dummy)
    ! consistency check
    if ((size(dummy, dim = 1) .NE. level%nrows) .or. &
            (size(dummy, dim = 2) .NE. level%ncols)) then
      call message('   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file for ', trim(level_name), &
              ' in basin ', trim(adjustl(num2str(ii))), '!')
      call message('  Latitude expected to have following dimensions ... rows:', &
              trim(adjustl(num2str(level%nrows))), ', cols:', trim(adjustl(num2str(level%ncols))))
      call message('  Latitude provided ... rows:', &
              trim(adjustl(num2str(size(dummy, dim = 1)))), ', cols:', trim(adjustl(num2str(size(dummy, dim = 2)))))
      stop 1
    end if
    level%y = dummy

    var = nc%getVariable(trim(lon_var_name))
    call var%getData(dummy)
    ! consistency check
    if ((size(dummy, dim = 1) .NE. level%nrows) .or. &
            (size(dummy, dim = 2) .NE. level%ncols)) then
      call message('   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file for ', trim(level_name), &
              ' in basin ', trim(adjustl(num2str(ii))), '!')
      call message('  Longitude expected to have following dimensions ... rows:', &
              trim(adjustl(num2str(level%nrows))), ', cols:', trim(adjustl(num2str(level%ncols))))
      call message('  Longitude provided ... rows:', &
              trim(adjustl(num2str(size(dummy, dim = 1)))), ', cols:', trim(adjustl(num2str(size(dummy, dim = 2)))))
      stop 1
    end if
    level%x = dummy

    call nc%close()

  end subroutine read_latlon

END MODULE mo_read_latlon
