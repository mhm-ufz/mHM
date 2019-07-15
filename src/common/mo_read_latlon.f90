!>       \file mo_read_latlon.f90

!>       \brief reading latitude and longitude coordinates for each domain

!>       \details TODO: add description

!>       \authors Stephan Thober

!>       \date Nov 2013

! Modifications:

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

  !    NAME
  !        read_latlon

  !    PURPOSE
  !>       \brief reads latitude and longitude coordinates

  !>       \details reads latitude and longitude coordinates from
  !>       netcdf file for each domain and appends it to the global
  !>       variables latitude and longitude.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: ii"            domain indexFile name of the domains must be xxx_latlon.nc, wherexxx
  !>       is the domain id. Variable names in the netcdf filehave to be 'lat' for latitude and 'lon' for longitude.
  !>       \param[in] "character(*) :: lon_var_name"
  !>       \param[in] "character(*) :: lat_var_name"
  !>       \param[in] "character(*) :: level_name"

  !    INTENT(INOUT)
  !>       \param[inout] "type(Grid) :: level"

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Nov 2013

  ! Modifications:
  ! Stephan Thober, Sep 2015 - added latitude and longitude for level 0
  ! Stephan Thober, Oct 2015 - added L1_rect_latitude and L1_rect_longitude
  ! David Schaefer, May 2016 - removed ncread dependency
  ! Robert Schweppe, Mar 2018 - major rewrite
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_latlon(ii, lon_var_name, lat_var_name, level_name, level)

    use mo_common_variables, only : Grid, fileLatLon
    use mo_message, only : message
    use mo_netcdf, only : NcDataset, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    ! domain indexFile name of the domains must be xxx_latlon.nc, wherexxx is the domain id. Variable names in the netcdf
    ! filehave to be 'lat' for latitude and 'lon' for longitude.
    integer(i4), intent(in) :: ii

    character(*), intent(in) :: lon_var_name

    character(*), intent(in) :: lat_var_name

    character(*), intent(in) :: level_name

    type(Grid), intent(inout) :: level

    ! file name
    character(256) :: fname

    ! dummy variable
    real(dp), dimension(:, :), allocatable :: dummy

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
              ' in domain ', trim(adjustl(num2str(ii))), '!')
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
              ' in domain ', trim(adjustl(num2str(ii))), '!')
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
