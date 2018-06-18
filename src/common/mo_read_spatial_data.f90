!>       \file mo_read_spatial_data.f90

!>       \brief Reads spatial input data.

!>       \details This module is to read spatial input data, e.g. dem, aspect, flow direction.
!>       The module provides a subroutine for ASCII files. 
!>       (Subroutine for NetCDF files will come with release 5.1).
!>       The data are read from the specified directory.

!>       \authors s Juliane Mai

!>       \date Dec 2012

! Modifications:

MODULE mo_read_spatial_data

  ! This module provides routines to read spatial data.

  ! Written  Juliane Mai, Jan 2013
  ! Modified 

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PUBLIC :: read_header_ascii           ! Reads header of ASCII files
  PUBLIC :: read_spatial_data_ascii     ! Read ASCII  files
  ! PUBLIC :: read_spatial_data_nc      ! Read netCDF files -> will be implemented in release 5.1

  ! ------------------------------------------------------------------

  INTERFACE  read_spatial_data_ascii
    MODULE PROCEDURE read_spatial_data_ascii_i4, read_spatial_data_ascii_dp
  END INTERFACE read_spatial_data_ascii

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_spatial_data_ascii_dp

  !    PURPOSE
  !>       \brief Reads spatial data files of ASCII format.

  !>       \details Reads spatial input data, e.g. dem, aspect, flow direction.

  !    INTENT(IN)
  !>       \param[in] "character(len = *) :: filename" Name of file and its location
  !>       \param[in] "integer(i4) :: fileunit"        File unit for open file
  !>       \param[in] "integer(i4) :: header_nCols"    Reference number of columns
  !>       \param[in] "integer(i4) :: header_nRows"    Reference number of rows
  !>       \param[in] "real(dp) :: header_xllcorner"   Reference lower left corner (x)
  !>       \param[in] "real(dp) :: header_yllcorner"   Reference lower left corner (y)
  !>       \param[in] "real(dp) :: header_cellsize"    Reference cell size [m]

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: data" Data matrixdim_1 = longitude, dim_2 = latitude
  !>       \param[out] "logical, dimension(:, :) :: mask"  Mask of data matrixdim_1 = longitude, dim_2 = latitude

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Jan 2013

  ! Modifications:
  ! Matthias Zink  Feb 2013 - , added interface and routine for datatype i4
  ! David Schaefer Mar 2015 - , removed double allocation of temporary data

  subroutine read_spatial_data_ascii_dp(filename, fileunit, header_ncols, header_nrows, header_xllcorner, &
                                       header_yllcorner, header_cellsize, data, mask)
    implicit none

    ! Name of file and its location
    character(len = *), intent(in) :: filename

    ! File unit for open file
    integer(i4), intent(in) :: fileunit

    ! Reference number of rows
    integer(i4), intent(in) :: header_nRows

    ! Reference number of columns
    integer(i4), intent(in) :: header_nCols

    ! Reference lower left corner (x)
    real(dp), intent(in) :: header_xllcorner

    ! Reference lower left corner (y)
    real(dp), intent(in) :: header_yllcorner

    ! Reference cell size [m]
    real(dp), intent(in) :: header_cellsize

    ! Data matrixdim_1 = longitude, dim_2 = latitude
    real(dp), dimension(:, :), allocatable, intent(out) :: data

    ! Mask of data matrixdim_1 = longitude, dim_2 = latitude
    logical, dimension(:, :), allocatable, intent(out) :: mask

    ! number of rows of data fields:
    integer(i4) :: file_nRows

    ! number of columns of data fields:
    integer(i4) :: file_nCols

    ! file read in lower left corner
    real(dp) :: file_xllcorner

    ! file read in lower left corner
    real(dp) :: file_yllcorner

    ! file read in cellsize
    real(dp) :: file_cellsize

    ! file read in nodata value
    real(dp) :: file_nodata

    integer(i4) :: i, j

    ! data
    real(dp), dimension(:, :), allocatable :: tmp_data

    ! mask
    logical, dimension(:, :), allocatable :: tmp_mask


    ! compare headers always with reference header (intent in)
    call read_header_ascii(filename, fileunit, &
            file_ncols, file_nrows, &
            file_xllcorner, file_yllcorner, file_cellsize, file_nodata)
    if ((file_ncols .ne. header_ncols)) &
            stop 'read_spatial_data_ascii: header not matching with reference header: ncols'
    if ((file_nrows .ne. header_nrows)) &
            stop 'read_spatial_data_ascii: header not matching with reference header: nrows'
    if ((abs(file_xllcorner - header_xllcorner) .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: xllcorner'
    if ((abs(file_yllcorner - header_yllcorner) .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: yllcorner'
    if ((abs(file_cellsize - header_cellsize)   .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: cellsize'

    ! allocation and initialization of matrices
    allocate(tmp_data(file_nrows, file_ncols))
    tmp_data = file_nodata
    allocate(tmp_mask(file_nrows, file_ncols))
    tmp_mask = .true.

    ! read in
    ! recl is only a rough estimate on bytes per line in the ascii
    ! default for nag: recl=1024(byte) which is not enough for 100s of columns
    open (unit = fileunit, file = filename, action = 'read', status = 'old', recl = 48 * file_ncols)
    ! (a) skip header
    do i = 1, 6
      read(fileunit, *)
    end do
    ! (b) read data
    do i = 1, file_nrows
      read(fileunit, *) (tmp_data(i, j), j = 1, file_ncols)
    end do
    close(fileunit)

    ! set mask .false. if nodata value appeared
    where (abs(tmp_data - file_nodata) .lt. tiny(1.0_dp))
      tmp_mask = .false.
    end where

    ! transpose of data due to longitude-latitude ordering
    allocate(data(file_ncols, file_nrows))
    data = transpose(tmp_data)
    deallocate(tmp_data)

    allocate(mask(file_ncols, file_nrows))
    mask = transpose(tmp_mask)
    deallocate(tmp_mask)

  end subroutine read_spatial_data_ascii_dp

  subroutine read_spatial_data_ascii_i4(filename, fileunit, &
          header_ncols, header_nrows, &
          header_xllcorner, header_yllcorner, &
          header_cellsize, &
          data, mask)

    implicit none

    character(len = *), intent(in) :: filename         ! filename with location
    integer(i4), intent(in) :: fileunit         ! unit for opening the file
    integer(i4), intent(in) :: header_nRows     ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4), intent(in) :: header_nCols     ! number of columns of data fields:
    ! LATITUDE dimension
    real(dp), intent(in) :: header_xllcorner ! header read in lower left corner
    real(dp), intent(in) :: header_yllcorner ! header read in lower left corner
    real(dp), intent(in) :: header_cellsize  ! header read in cellsize
    integer(i4), dimension(:, :), allocatable, intent(out) :: data             ! data
    logical, dimension(:, :), allocatable, intent(out) :: mask             ! mask

    ! local variables
    integer(i4) :: file_nRows     ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4) :: file_nCols     ! number of columns of data fields:
    ! LATITUDE dimension
    real(dp) :: file_xllcorner ! file read in lower left corner
    real(dp) :: file_yllcorner ! file read in lower left corner
    real(dp) :: file_cellsize  ! file read in cellsize
    real(dp) :: file_nodata    ! file read in nodata value
    integer(i4) :: i, j
    integer(i4), dimension(:, :), allocatable :: tmp_data       ! data
    logical, dimension(:, :), allocatable :: tmp_mask       ! mask

    ! compare headers always with reference header (intent in)
    call read_header_ascii(filename, fileunit, &
            file_ncols, file_nrows, &
            file_xllcorner, file_yllcorner, file_cellsize, file_nodata)
    if ((file_ncols .ne. header_ncols)) &
            stop 'read_spatial_data_ascii: header not matching with reference header: ncols'
    if ((file_nrows .ne. header_nrows)) &
            stop 'read_spatial_data_ascii: header not matching with reference header: nrows'
    if ((abs(file_xllcorner - header_xllcorner) .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: xllcorner'
    if ((abs(file_yllcorner - header_yllcorner) .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: yllcorner'
    if ((abs(file_cellsize - header_cellsize)   .gt. tiny(1.0_dp))) &
            stop 'read_spatial_data_ascii: header not matching with reference header: cellsize'

    ! allocation and initialization of matrices
    allocate(tmp_data(file_nrows, file_ncols))
    tmp_data = int(file_nodata, i4)
    allocate(tmp_mask(file_nrows, file_ncols))
    tmp_mask = .true.

    ! read in
    ! recl is only a rough estimate on bytes per line in the ascii
    ! default for nag: recl=1024(byte) which is not enough for 100s of columns
    open (unit = fileunit, file = filename, action = 'read', status = 'old', recl = 48 * file_ncols)
    ! (a) skip header
    do i = 1, 6
      read(fileunit, *)
    end do
    ! (b) read data
    do i = 1, file_nrows
      read(fileunit, *) (tmp_data(i, j), j = 1, file_ncols)
    end do
    close(fileunit)

    ! set mask .false. if nodata value appeared
    where (tmp_data .EQ. int(file_nodata, i4))
      tmp_mask = .false.
    end where

    ! transpose of data due to longitude-latitude ordering
    allocate(data(file_ncols, file_nrows))
    data = transpose(tmp_data)
    deallocate(tmp_data)

    allocate(mask(file_ncols, file_nrows))
    mask = transpose(tmp_mask)
    deallocate(tmp_mask)

  end subroutine read_spatial_data_ascii_i4

  ! ------------------------------------------------------------------

  !     NAME
  !         read_header_ascii

  !     PURPOSE
  !>        \brief Reads header lines of ASCII files.

  !>        \details Reads header lines of ASCII files, e.g. dem, aspect, flow direction.

  !     CALLING SEQUENCE
  !         call read_header_ascii(trim(adjustl(dirMorpho))//file_dem, udem, &
  !                                header_nRows, header_nCols, header_xllcorner, header_yllcorner, &
  !                                header_cellsize, header_nodata)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: filename"          Name of file and its location
  !>        \param[in] "integer(i4)      :: fileunit"          File unit for open file 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)      :: header_ncols"      Reference number of columns 
  !>        \param[out] "integer(i4)      :: header_nrows"      Reference number of rows
  !>        \param[out] "real(dp)         :: header_xllcorner"  Reference lower left corner (x)
  !>        \param[out] "real(dp)         :: header_yllcorner"  Reference lower left corner (y)
  !>        \param[out] "integer(i4)      :: header_cellsize"   Reference cell size [m]
  !>        \param[out] "real(dp)         :: header_nodata"     Reference nodata value

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jan 2013

  subroutine read_header_ascii(filename, fileunit, &
          header_ncols, header_nrows, header_xllcorner, header_yllcorner, header_cellsize, header_nodata)

    implicit none

    character(len = *), intent(in) :: filename         ! filename with location
    integer(i4), intent(in) :: fileunit         ! unit for opening the file
    integer(i4), intent(out) :: header_nRows     ! number of rows of data fields:
    ! LONGITUDE dimension
    integer(i4), intent(out) :: header_nCols     ! number of columns of data fields:
    ! LATITUDE dimension
    real(dp), intent(out) :: header_xllcorner ! header read in lower left corner
    real(dp), intent(out) :: header_yllcorner ! header read in lower left corner
    real(dp), intent(out) :: header_cellsize  ! header read in cellsize
    real(dp), intent(out) :: header_nodata    ! header read in nodata value

    ! local variables
    character(5) :: dummy

    ! reading header from a file
    open (unit = fileunit, file = filename, status = 'old')
    read (fileunit, *) dummy, header_nCols
    read (fileunit, *) dummy, header_nRows
    read (fileunit, *) dummy, header_xllcorner
    read (fileunit, *) dummy, header_yllcorner
    read (fileunit, *) dummy, header_cellsize
    read (fileunit, *) dummy, header_nodata
    close(fileunit)
    dummy = dummy // ''   ! only to avoid warning

  end subroutine read_header_ascii

END MODULE mo_read_spatial_data
