!> \file mo_read_latlon.f90

!> \brief reading latitude and longitude coordinates for each basin

!> \authors Stephan Thober
!> \date Nov 2013

MODULE mo_mrm_read_latlon

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
  !         modified, Stephan Thober, Sep 2015 - added latitude and longitude for level 0
  !                   Stephan Thober, Oct 2015 - added L1_rect_latitude and L1_rect_longitude
  !                   Stephan Thober, Oct 2015 - adapted to mRM

  subroutine read_latlon(ii)
    
    USE mo_mrm_global_variables, ONLY: fileLatLon, level11, &
         L0_latitude, L0_longitude, level0, basin_mrm, L0_Basin, &
         L11_rect_latitude, L11_rect_longitude
    USE mo_append,           ONLY: append
    USE mo_message,          ONLY: message
    USE mo_ncread,           ONLY: get_NcVar, get_NcDim
    use mo_string_utils,     only: num2str

    implicit none

    integer(i4), intent(in) :: ii ! basin index

    ! local variables
    character(256)                        :: fname ! file name
    integer(i4), dimension(5)             :: dl    ! dimension lengths
    real(dp), dimension(:,:), allocatable :: dummy ! dummy variable
    logical, dimension(:,:), allocatable  :: mask

    ! construct filename
    fname = trim( fileLatLon(ii) ) 

    ! -------------------------------------------------------------------------
    ! READ LEVEL 0 LATITUDE / LONGITUDE
    ! -------------------------------------------------------------------------
    if (ii .eq. 1) then
       ! create mask for level 0
       mask = reshape(basin_mrm%L0_mask(basin_mrm%L0_iStartMask(ii):basin_mrm%L0_iEndMask(ii)), &
            (/level0%nrows(ii), level0%ncols(ii)/))
       ! read dimension length of variable in netcdf File
       dl = get_NcDim( trim(fname), 'lat_l0' )
    
       ! consistency check
       if ( (dl(1) .NE. level0%nrows(ii) ) .or. &
            (dl(2) .NE. level0%ncols(ii) ) ) then
          call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
          call message( '  Latlon expected to have following dimensions ... level0%nrows(ii):',               &
               trim(adjustl(num2str(level0%nrows(ii)))),', level0%ncols(ii):', trim(adjustl(num2str(level0%ncols(ii)))))
          call message( '  Latlon provided ... level0%nrows(ii):',               &
               trim(adjustl(num2str(dl(1)))),', level0%ncols(ii):', trim(adjustl(num2str(dl(2)))))
          stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!' 
       end if
       
       allocate(dummy(dl(1), dl(2)))
       
       ! read dummy variable
       call get_NcVar( trim(fname), 'lat_l0', dummy )
       
       ! append it to global variables
       call append( L0_latitude, pack( dummy, mask ))
       deallocate(dummy)
       
       ! read dimension length of variable in netcdf File
       dl = get_NcDim( trim(fname), 'lon_l0' )
       
       ! consistency check
       if ( (dl(1) .NE. level0%nrows(ii) ) .or. &
            (dl(2) .NE. level0%ncols(ii) ) ) then
          call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
          call message( '  Latlon expected to have following dimensions ... level0%nrows(ii):',               &
               trim(adjustl(num2str(level0%nrows(ii)))),', level0%ncols(ii):', trim(adjustl(num2str(level0%ncols(ii)))))
          call message( '  Latlon provided ... level0%nrows(ii):',               &
               trim(adjustl(num2str(dl(1)))),', level0%ncols(ii):', trim(adjustl(num2str(dl(2)))))
          stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!'
       end if
       
       allocate(dummy(dl(1), dl(2)))

       ! read dummy variable
       call get_NcVar( trim(fname), 'lon_l0', dummy )

       ! append it to global variables
       call append( L0_longitude, pack( dummy, .true. ))
       deallocate(dummy)
    else if (L0_Basin(ii) .ne. L0_Basin(ii - 1)) then
       ! create mask for level 0
       mask = reshape(basin_mrm%L0_mask(basin_mrm%L0_iStartMask(ii):basin_mrm%L0_iEndMask(ii)), &
            (/level0%nrows(ii), level0%ncols(ii)/))
       ! read dimension length of variable in netcdf File
       dl = get_NcDim( trim(fname), 'lat_l0' )
    
       ! consistency check
       if ( (dl(1) .NE. level0%nrows(ii) ) .or. &
            (dl(2) .NE. level0%ncols(ii) ) ) then
          call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
          call message( '  Latlon expected to have following dimensions ... level0%nrows(ii):',               &
               trim(adjustl(num2str(level0%nrows(ii)))),', level0%ncols(ii):', trim(adjustl(num2str(level0%ncols(ii)))))
          call message( '  Latlon provided ... level0%nrows(ii):',               &
               trim(adjustl(num2str(dl(1)))),', level0%ncols(ii):', trim(adjustl(num2str(dl(2)))))
          stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!' 
       end if
       
       allocate(dummy(dl(1), dl(2)))
       
       ! read dummy variable
       call get_NcVar( trim(fname), 'lat_l0', dummy )
       
       ! append it to global variables
       call append( L0_latitude, pack( dummy, mask ))
       deallocate(dummy)
       
       ! read dimension length of variable in netcdf File
       dl = get_NcDim( trim(fname), 'lon_l0' )
       
       ! consistency check
       if ( (dl(1) .NE. level0%nrows(ii) ) .or. &
            (dl(2) .NE. level0%ncols(ii) ) ) then
          call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
          call message( '  Latlon expected to have following dimensions ... level0%nrows(ii):',               &
               trim(adjustl(num2str(level0%nrows(ii)))),', level0%ncols(ii):', trim(adjustl(num2str(level0%ncols(ii)))))
          call message( '  Latlon provided ... level0%nrows(ii):',               &
               trim(adjustl(num2str(dl(1)))),', level0%ncols(ii):', trim(adjustl(num2str(dl(2)))))
          stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!'
       end if
       
       allocate(dummy(dl(1), dl(2)))

       ! read dummy variable
       call get_NcVar( trim(fname), 'lon_l0', dummy )

       ! append it to global variables
       call append( L0_longitude, pack( dummy, mask ))
       deallocate(dummy)
    end if
    ! clean up
    if (allocated(mask)) deallocate(mask)
    mask = reshape(basin_mrm%L11_mask(basin_mrm%L11_iStartMask(ii):basin_mrm%L11_iEndMask(ii)), &
         (/level11%nrows(ii), level11%ncols(ii)/))

    ! -------------------------------------------------------------------------
    ! READ LEVEL 11 LATITUDE / LONGITUDE
    ! -------------------------------------------------------------------------
    ! read dimension length of variable in netcdf File
    dl = get_NcDim( trim(fname), 'lat_l11' )
    
    ! consistency check
    if ( (dl(1) .NE. level11%nrows(ii) ) .or. &
         (dl(2) .NE. level11%ncols(ii) ) ) then
       call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
       call message( '  Latlon expected to have following dimensions ... level11%nrows(ii):',               &
            trim(adjustl(num2str(level11%nrows(ii)))),', level11%ncols(ii):', trim(adjustl(num2str(level11%ncols(ii)))))
       call message( '  Latlon provided ... level11%nrows(ii):',               &
            trim(adjustl(num2str(dl(1)))),', level11%ncols(ii):', trim(adjustl(num2str(dl(2)))))
       stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!' 
    end if
    
    allocate(dummy(dl(1), dl(2)))
    
    ! read dummy variable
    call get_NcVar( trim(fname), 'lat_l11', dummy )

    ! append it to global variables
    call append( L11_rect_latitude, pack(dummy, .true.))
    deallocate(dummy)

    ! read dimension length of variable in netcdf File
    dl = get_NcDim( trim(fname), 'lon_l11' )
    
    ! consistency check
    if ( (dl(1) .NE. level11%nrows(ii) ) .or. &
         (dl(2) .NE. level11%ncols(ii) ) ) then
         call message( '   ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!')
         call message( '  Latlon expected to have following dimensions ... level11%nrows(ii):',               &
              trim(adjustl(num2str(level11%nrows(ii)))),', level11%ncols(ii):', trim(adjustl(num2str(level11%ncols(ii)))))
         call message( '  Latlon provided ... level11%nrows(ii):',               &
              trim(adjustl(num2str(dl(1)))),', level11%ncols(ii):', trim(adjustl(num2str(dl(2)))))
         stop '    ***ERROR: subroutine mo_read_latlon: size mismatch in latlon file!'
    end if
    
    allocate(dummy(dl(1), dl(2)))
    
    ! read dummy variable
    call get_NcVar( trim(fname), 'lon_l11', dummy )

    ! append it to global variables
    call append( L11_rect_longitude, pack(dummy, .true.))
    deallocate(dummy, mask)

  end subroutine read_latlon

END MODULE mo_mrm_read_latlon
