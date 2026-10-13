!> \file mo_meteo_weights.f90
!> \brief \copybrief mo_meteo_weights
!> \details \copydetails mo_meteo_weights

!> \brief Read NetCDF weights for meteorological disaggregation.
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
!> \ingroup f_meteo
module mo_meteo_weights
  use mo_message, only: error_message
  implicit none
  private
  public :: read_weights_nc
contains
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

    character(len = *), intent(in) :: folder !< Name of the folder where data are stored
    integer(i4), intent(in) :: nRows !< Number of datapoints in longitudinal direction
    integer(i4), intent(in) :: nCols !< Number of datapoints in latitudinal  direction
    character(len = *), intent(in) :: varName !< Name of variable name to read
    real(dp), dimension(:, :, :, :), allocatable, intent(out) :: data !< Data matrixdim_1 = longitude, dim_2 = latitude, dim_3 = months, dim_4 = hours
    logical, dimension(:, :), intent(in) :: mask !< mask of valid data fields
    real(dp), optional, intent(in) :: lower !< Lower bound for check of validity of data values
    real(dp), optional, intent(in) :: upper !< Upper bound for check of validity of data values
    logical, optional, intent(in) :: nocheck !< .TRUE. if check for nodata values deactivateddefault = .FALSE. - check is done
    character(256), optional, intent(in) :: fileName !< name of variable, defaults to fileName
    logical, dimension(:, :, :, :), allocatable, optional, intent(out) :: maskout !< ! mask of validdata points
    logical, optional, intent(in) :: bound_error !< .FALSE. to only warn about bound (lower, upper) violations, default = .TRUE. - raise an error

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
end module mo_meteo_weights
