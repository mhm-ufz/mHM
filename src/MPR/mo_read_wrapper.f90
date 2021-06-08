!>       \file mo_read_wrapper.f90

!>       \brief Wrapper for all reading routines.

!>       \details This module is to wrap up all reading routines.
!>       The general written reading routines are used to store now the read data into global variables.

!>       \authors Juliane Mai, Matthias Zink

!>       \date Jan 2013

! Modifications:

MODULE mo_read_wrapper

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PUBLIC :: read_data            ! reads all available data

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_data

  !    PURPOSE
  !>       \brief Reads data.

  !>       \details The namelists are already read by read_config call.
  !>       All LUTs are read from their respective directory and information within those
  !>       files are shared across all domains to be modeled.

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "type(period), dimension(:), optional :: LAIPer"

  !    HISTORY
  !>       \authors Juliane Mai & Matthias Zink

  !>       \date Feb 2013

  ! Modifications:
  ! Luis Samaniego               Feb 2013 - rotate fdir variable to the new coordinate system
  ! Rohini Kumar                 Aug 2013 - name changed from "L0_LAI" to "L0_LCover_LAI"
  ! Rohini Kumar                 Aug 2013 - added dirSoil_LUT and dirGeology_LUT, and changed to
  !                                         read datapaths and variables made accordingly
  ! Rohini Kumar                 Aug 2013 - added iFlag_LAI_data_format to handle LAI options, and changed
  !                                         within the code made accordingly
  ! Rohini  Kumar                Sep 2013 - read input data for routing processes according to process_matrix flag
  ! Matthias Zink                Mar 2014 - added inflow gauge
  ! Kumar & Schroen              Apr 2014 - added check for consistency of L0 and L1 spatial resolution
  ! Stephan Thober               Jun 2014 - added perform_mpr for omitting L0 read
  ! Matthias Cuntz & Juliane Mai Nov 2014 - LAI input from daily, monthly or yearly files
  ! Stephan Thober               Aug 2015 - moved routing related variables and routines to mRM
  ! Rohini Kumar                 Mar 2016 - options to handle different soil databases
  ! Matthias Zink                Mar 2014 - added subroutine for consistency check
  ! Stephan Thober               Nov 2016 - moved processMatrix to common variables
  ! Rohini Kumar                 Dec 2016 - option to handle monthly mean gridded fields of LAI
  ! Robert Schweppe              Jun 2018 - refactoring and reformatting

  subroutine read_data(LAIPer)

    use mo_append, only : append, paste
    use mo_constants, only : YearMonths
    use mo_file, only : varNameLAIclass, varNameSlope, varNameAspect, varNameSoilClass, varNameGeoClass
    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_read_data, only : read_dem, read_lcover
    use mo_common_variables, only : Grid, dirCommonFiles, dirMorpho, &
                                    global_parameters, level0, domainMeta, period, processMatrix
    use mo_message, only : message
    use mo_mpr_file, only : file_geolut, file_hydrogeoclass, &
                            file_laiclass, file_lailut, file_soil_database, file_soil_database_1, &
                            ugeolut, ulailut
    use mo_mpr_global_variables, only : GeoUnitKar, &
                                        GeoUnitList, L0_asp, L0_geoUnit, L0_gridded_LAI, L0_slope, L0_soilId, LAILUT, &
                                        LAIUnitList, iFlag_soilDB, nGeoUnits, nLAI, nLAIclass, nSoilHorizons_mHM, soilDB, &
                                        LAIBoundaries
    use mo_common_datetime_type, only: timeStep_LAI_input
    use mo_prepare_gridded_lai, only : prepare_gridded_daily_LAI_data, prepare_gridded_mean_monthly_LAI_data
    use mo_read_latlon, only : read_latlon
    use mo_read_lut, only : read_geoformation_lut, read_lai_lut
    use mo_soil_database, only : read_soil_LUT
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, &
                         timer_stop
    use mo_netcdf, only: NcDataset, NcVariable
    use mo_read_nc, only: read_const_nc

    implicit none

    type(period), dimension(:), intent(in), optional :: LAIPer

    ! loop variables
    integer(i4) :: domainID, iDomain, iVar, iHorizon, iMon, itimer, ll

    ! dummy variable
    integer(i4) :: nH

    ! file name of file to read
    character(256) :: fName, varName

    real(dp), dimension(:, :), allocatable :: data_dp_2d

    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    !! real(dp), dimension(:, :, :), allocatable :: data_dp_3d_nc

    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4

    logical, dimension(:, :), allocatable :: mask_2d
    !! logical, dimension(:, :, :), allocatable :: mask_3d_nc

    integer(i4), dimension(:), allocatable :: dummy_i4

    type(Grid), pointer :: level0_iDomain

    real(dp), parameter :: slope_minVal = 0.01_dp

    real(dp), parameter :: aspect_minVal = 1.00_dp
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar          ! variables for data form netcdf


    call message('  Reading data ...')
    itimer = 1
    call timer_start(itimer)

    call message('    Reading dem and lcover ...')
    call read_dem()
    call read_lcover()

    ! ************************************************
    ! READ LOOKUP TABLES
    ! ************************************************
    !
    ! Soil LUT
    if(iFlag_soilDB .eq. 0) then
      fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_soil_database))
    else if(iFlag_soilDB .eq. 1) then
      fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_soil_database_1))
    end if
    call read_soil_LUT(trim(fName))

    ! Geological formation LUT
    fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_geolut))
    call read_geoformation_lut(trim(fName), ugeolut, nGeoUnits, GeoUnitList, GeoUnitKar)

    ! LAI LUT
    if(timeStep_LAI_input .EQ. 0) then
      fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_lailut))
      call read_lai_lut(trim(fName), ulailut, nLAIclass, LAIUnitList, LAILUT)
    end if

    domains: do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

      call message('    Reading data for domain: ', trim(adjustl(num2str(domainID))), ' ...')
      ! check whether L0 data is shared
      ! ToDo: check change
      if (domainMeta%L0DataFrom(iDomain) < iDomain) then
        !
        call message('      Using data of domain ', &
                trim(adjustl(num2str(domainMeta%indices(domainMeta%L0DataFrom(iDomain))))), ' for domain: ',&
                trim(adjustl(num2str(domainID))), '...')
        ! DO NOT read L0 data
        cycle

      end if

      itimer = 2
      call timer_start(itimer)

      ! read slope and aspect - datatype real
      call message('      Reading slope ...')
      call read_const_nc(&
              folder=dirMorpho(iDomain),&
              varName=varNameSlope, &
              data=data_dp_2d, &
              nRows=level0_iDomain%nrows, &
              nCols=level0_iDomain%ncols, &
              maskout=mask_2d &
      )
      ! put global nodata value into array (probably not all grid cells have values)
      data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d)
      call append(L0_slope, pack(data_dp_2d, level0_iDomain%mask))
      call message('      Reading aspect ...')
      call read_const_nc(&
              folder=dirMorpho(iDomain),&
              varName=varNameAspect, &
              data=data_dp_2d, &
              nRows=level0_iDomain%nrows, &
              nCols=level0_iDomain%ncols, &
              maskout=mask_2d &
      )
      ! put global nodata value into array (probably not all grid cells have values)
      data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d)
      call append(L0_asp, pack(data_dp_2d, level0_iDomain%mask))

      deallocate(data_dp_2d, mask_2d)

      ! read datatype integer
      ! ***** CHANGE IS MADE HERE TO ACCOMODATE READ SEVERAL TYPES OF SOIL DATABASES
      ! change from the earlier code where everything was done in the DO loop (**see below **)
      ! here everything is more explicit and no do loop appears as was the case in previous version

      ! read soilID in both options
      call message('      Reading soil ids ...')

      nH = 1 !> by default; when iFlag_soilDB = 0
      if(iFlag_soilDB .eq. 1) nH = nSoilHorizons_mHM
      ! modified way to read multiple horizons specific soil class
      do iHorizon = 1, nH
        if(iFlag_soilDB .eq. 0) then
          fName = trim(adjustl(dirMorpho(iDomain))) // trim(varNameSoilClass) // '.nc'
          varName = trim(varNameSoilClass)
        else if(iFlag_soilDB .eq. 1) then
          write(varName, 172) iHorizon
          172             format('soil_class_horizon_', i2.2)
          fName = trim(adjustl(dirMorpho(iDomain))) // trim(adjustl(varName)) // '.nc'
        end if
        ! read the Dataset
        nc = NcDataset(fname, "r")
        ! get the variable
        ncVar = nc%getVariable(trim(varName))
        ! read data
        call ncVar%getData(data_i4_2d, mask=mask_2d)
        if ( size(data_i4_2d, 1) /= level0_iDomain%nrows .or. size(data_i4_2d, 2) /= level0_iDomain%ncols) then
          print*, '***ERROR: read_forcing_nc: mHM generated x and y: ', level0_iDomain%nrows, level0_iDomain%ncols , &
                'are not matching NetCDF dimensions: ', shape(data_i4_2d)
          stop 1
        end if
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        ! put data in variable
        call paste(dataMatrix_i4, pack(data_i4_2d, level0_iDomain%mask), nodata_i4)
        ! deallocate arrays
        deallocate(data_i4_2d)
        ! close the handler
        call nc%close()
      end do
      call append(L0_soilId, dataMatrix_i4)
      deallocate(dataMatrix_i4)

      ! read geoUnit
      fName = trim(adjustl(dirMorpho(iDomain))) // trim(varNameGeoClass) // '.nc'
      ! reading and transposing
      ! read the Dataset
      nc = NcDataset(fname, "r")
      ! get the variable
      ncVar = nc%getVariable(trim(varNameGeoClass))
      ! read data
      call ncVar%getData(data_i4_2d, mask=mask_2d)
      if ( size(data_i4_2d, 1) /= level0_iDomain%nrows .or. size(data_i4_2d, 2) /= level0_iDomain%ncols) then
        print*, '***ERROR: read_forcing_nc: mHM generated x and y: ', level0_iDomain%nrows, level0_iDomain%ncols , &
              'are not matching NetCDF dimensions: ', shape(data_i4_2d)
        stop 1
      end if
      ! put global nodata value into array (probably not all grid cells have values)
      data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
      ! put global nodata value into array (probably not all grid cells have values)
      call append(L0_geoUnit, pack(data_i4_2d, level0_iDomain%mask))
      ! deallocate arrays
      deallocate(data_i4_2d, mask_2d)
      ! close the handler
      call nc%close()

      ! LAI values
      call message('      Reading LAI ...')
      select case (timeStep_LAI_input)
      case(1) ! long term mean monthly gridded fields
        call prepare_gridded_mean_monthly_LAI_data(iDomain, level0_iDomain%nrows, level0_iDomain%ncols, level0_iDomain%mask)

      case(0) ! long term mean monthly values per class with LUT
        ! only set if not yet allocated (e.g. domain 1)
        if (.not. allocated(LAIBoundaries)) then
          nLAI = nint(YearMonths, i4)
          allocate(LAIBoundaries(nLAI+1))
          LAIBoundaries = [(iMon, iMon=1, nLAI+1)]
        end if

        fName = trim(dirMorpho(iDomain)) // trim(varNameLAIclass) // '.nc'
        ! read the Dataset
        nc = NcDataset(fname, "r")
        ! get the variable
        ncVar = nc%getVariable(trim(varNameLAIclass))
        ! read data
        call ncVar%getData(data_i4_2d, mask=mask_2d)
        if ( size(data_i4_2d, 1) /= level0_iDomain%nrows .or. size(data_i4_2d, 2) /= level0_iDomain%ncols) then
          print*, '***ERROR: read_forcing_nc: mHM generated x and y: ', level0_iDomain%nrows, level0_iDomain%ncols , &
                'are not matching NetCDF dimensions: ', shape(data_i4_2d)
          stop 1
        end if
        call nc%close()

        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        allocate(dummy_i4(count(level0_iDomain%mask)))
        dummy_i4 = pack(data_i4_2d, level0_iDomain%mask)
        deallocate(data_i4_2d, mask_2d)

        call check_consistency_lut_map(dummy_i4, LAIUnitList, file_laiclass)

        allocate(data_dp_2d(count(level0_iDomain%mask), nLAI))
        do iMon = 1, nLAI
          ! determine LAIs per month
          do ll = 1, size(LAILUT, dim = 1)
            data_dp_2d(:, iMon) = merge(LAILUT(ll, iMon), data_dp_2d(:, iMon), dummy_i4(:) .EQ. LAIUnitList(ll))
          end do
        end do
        call append(L0_gridded_LAI, data_dp_2d(:, :))
        deallocate(dummy_i4, data_dp_2d)
        ! correction for 0 LAI values to avoid numerical instabilities
        L0_gridded_LAI(:, :) = merge(1.00E-10_dp, L0_gridded_LAI(:, :), L0_gridded_LAI(:, :) .LT. 1.00E-10_dp)
        L0_gridded_LAI(:, :) = merge(30.0_dp, L0_gridded_LAI(:, :), L0_gridded_LAI(:, :) .GT. 30.0_dp)
      case(-3 : -1) ! daily, monthly or yearly gridded fields (time-series)
        call prepare_gridded_daily_LAI_data(iDomain, level0_iDomain%nrows, level0_iDomain%ncols, level0_iDomain%mask, &
        LAIPer(iDomain))

      end select

      ! read lat lon coordinates of each domain
      call message('      Reading latitude/logitude ...')
      call read_latlon(iDomain, "lon_l0", "lat_l0", "level0", level0_iDomain)

      call timer_stop(itimer)
      call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

    end do domains
    itimer = 1
    call timer_stop(itimer)
    call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

    !Soil
    ! determine name od soil class definition file based on input option iFlag_soilDB
    if(iFlag_soilDB .eq. 0) then
      fName = file_soil_database
    else if(iFlag_soilDB .eq. 1) then
      fName = file_soil_database_1
    end if

    ! If you are getting segmentation fault with intel on large file (e.g. here on soils), then
    ! make sure to include following into your Marefile: INTEL_EXCLUDE  := mo_read_wrapper.f90
    call check_consistency_lut_map(reshape(L0_soilId, (/ size(L0_soilId, 1) * size(L0_soilId, 2) /)), &
            soilDB%id(:), fName)
    ! Geology
    call check_consistency_lut_map(L0_geoUnit, GeoUnitList, file_hydrogeoclass, dummy_i4)

    ! deactivate parameters of non existing geological classes in study domain for optimization
    ! loop over geological units in look up list
    do iVar = 1, size(GeoUnitList, 1)
      ! check if unit appears in geological map (dummy_i4 is unique number in L0_geoUnit)
      if (.not. ANY(dummy_i4 .EQ. GeoUnitList(iVar))) then
        ! deactivate optimization flag (dim=4 from global_parameters)
        global_parameters(processMatrix(9, 3) - processMatrix(9, 2) + iVar, 4) = 0
        call message('***WARNING: Geological unit ', trim(adjustl(num2str(GeoUnitList(iVar)))))
        call message('            is not appearing in study domain.')
      end if
    end do

    deallocate(dummy_i4) ! is allocated in subroutine check_consistency_lut_map - geology

    !----------------------------------------------------------------
    ! Correction for slope and aspect -- min value set above
    !----------------------------------------------------------------
    ! keep the colons (:) in the statements because of Intel's reallocation lhs problem
    L0_slope(:) = merge(slope_minVal, L0_slope(:), (L0_slope(:)  .lt.  slope_minVal))
    L0_asp(:) = merge(aspect_minVal, L0_asp(:), (L0_asp(:)    .lt. aspect_minVal))

  end subroutine read_data

  ! ------------------------------------------------------------------

  !    NAME
  !        check_consistency_lut_map

  !    PURPOSE
  !>       \brief Checks if classes in input maps appear in look up tables.

  !>       \details Determines wether a class appearing in the morphological input
  !>       is occuring in the respective look up table. mHM breaks if inconsistencies
  !>       are discovered.

  !    INTENT(IN)
  !>       \param[in] "integer(i4), dimension(:) :: data"        map of study domain
  !>       \param[in] "integer(i4), dimension(:) :: lookuptable" look up table corresponding to map
  !>       \param[in] "character(*) :: filename"                 name of the lut file - ERORR warn

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "integer(i4), dimension(:), optional :: unique_values" array of unique values in dataone

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Nov 2016

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine check_consistency_lut_map(data, lookuptable, filename, unique_values)

    use mo_message, only : message
    use mo_orderpack, only : unista
    use mo_string_utils, only : num2str

    implicit none

    ! map of study domain
    integer(i4), dimension(:), intent(in) :: data

    ! look up table corresponding to map
    integer(i4), dimension(:), intent(in) :: lookuptable

    ! name of the lut file - ERORR warn
    character(*), intent(in) :: filename

    ! array of unique values in dataone
    integer(i4), dimension(:), allocatable, intent(out), optional :: unique_values

    integer(i4) :: n_unique_elements

    integer(i4) :: ielement

    integer(i4), dimension(:), allocatable :: temp


    allocate(temp(size(data, 1)))
    ! copy L0_geoUnit because subroutine unista overwrites the nunit entries of the
    ! input array with the unique array values
    temp = data
    ! retrieve unique values of data
    call unista(temp, n_unique_elements)
    ! check if unit exists in look up table
    do ielement = 1, n_unique_elements
      if (.not. ANY(lookuptable .EQ. temp(ielement))) then
        call message()
        call message('***ERROR: Class ', trim(adjustl(num2str(temp(ielement)))), ' is missing')
        call message('          in input file ', trim(adjustl(filename)), ' ...')
        stop
      end if
    end do

    ! pass unique values if optional argument unique_values is given
    if (present(unique_values)) then
      allocate(unique_values(n_unique_elements))
      unique_values(:) = temp(1 : n_unique_elements)
    end if
    deallocate(temp)

  end subroutine check_consistency_lut_map


END MODULE mo_read_wrapper
