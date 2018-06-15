!> \file mo_read_wrapper.f90

!> \brief Wrapper for all reading routines.

!> \details This module is to wrap up all reading routines.\n
!> The general written reading routines are used to store now the read data into global variables.

!> \authors Juliane Mai, Matthias Zink
!> \date Jan 2013

MODULE mo_read_wrapper

  ! Written  Juliane Mai & Matthias Zink, Jan 2013
  ! Modified
  !          Luis Samaniego, Feb 2013  ! rotate fdir variable to the new coordinate system

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PUBLIC :: read_data            ! reads all available data

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !         read_data

  !     PURPOSE
  !>        \brief Reads data.

  !>        \details The namelists are already read by read_config call./n
  !>                 All LUTs are read from their respective directory and information within those
  !>                 files are shared across all basins to be modeled.
  !     INTENT(IN)
  !         None

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
  !>       \note read_config has to be called before

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai & Matthias Zink
  !>        \date Feb 2013
  !          Modified,
  !                    Luis Samaniego, Feb 2013  - rotate fdir variable to the new coordinate system
  !                    Rohini Kumar,   Aug 2013  - name changed from "L0_LAI" to "L0_LCover_LAI"
  !                    Rohini Kumar,   Aug 2013  - added dirSoil_LUT and dirGeology_LUT, and changed to
  !                                                read datapaths and variables made accordingly
  !                    Rohini Kumar,   Aug 2013  - added iFlag_LAI_data_format to handle LAI options,
  !                                                and changed within the code made accordingly
  !                    Rohini  Kumar,  Sep 2013  - read input data for routing processes according
  !                    Stephan Thober             to process_matrix flag
  !                    Matthias Zink   Mar 2014   added inflow gauge
  !                    Kumar & Schroen Apr 2014  - added check for consistency of L0 and L1 spatial resolution
  !                    Stephan Thober  Jun 2014  - added perform_mpr for omitting L0 read
  !                    Matthias Cuntz &
  !                    Juliane Mai     Nov 2014  - LAI input from daily, monthly or yearly files
  !                    Stephan Thober  Aug 2015  - moved routing related variables and routines to mRM
  !                    Rohini Kumar,   Mar 2016  - options to handle different soil databases
  !                    Matthias Zink   Mar 2014  - added subroutine for consistency check
  !                    Stephan Thober, Nov 2016  - moved processMatrix to common variables
  !                    Rohini Kuamr,   Dec  2016 - option to handle monthly mean gridded fields of LAI
  ! ------------------------------------------------------------------

  subroutine read_data(LAIPer)
    !
    use mo_common_read_data, only : read_dem, read_lcover
    USE mo_read_lut, ONLY : read_lai_lut, &
            read_geoformation_lut
    USE mo_soil_database, ONLY : read_soil_LUT
    USE mo_read_spatial_data, ONLY : read_spatial_data_ascii
    USE mo_read_latlon, ONLY : read_latlon
    USE mo_append, ONLY : append, paste
    USE mo_string_utils, ONLY : num2str
    USE mo_message, ONLY : message
    use mo_timer, only : timer_get, timer_start, timer_stop
    !
    USE mo_mpr_file, ONLY : file_geolut, ugeolut, & ! file name and unit of hydrogeology LuT
            file_lailut, ulailut, & ! file name and unit of LAI LuT
            file_slope, uslope, & ! file name and unit of slope map
            file_aspect, uaspect, & ! file name and unit of aspect map
            file_soilclass, usoilclass, & ! file name and unit of soil class map
            file_hydrogeoclass, uhydrogeoclass, & ! file name and unit of hydrogeo class map
            file_laiclass, ulaiclass, & ! file name and unit of lai class map
            file_soil_database, & ! file name of soil class map (iFlag_soilDB = 0)
            file_soil_database_1 ! file name of soil class map (iFlag_soilDB = 1)
    USE mo_mpr_global_variables, ONLY : nGeoUnits, GeoUnitList, GeoUnitKar, & ! geological class information
            L0_slope, & ! slope on input resolution (L0)
            L0_asp, & ! aspect on input resolution (L0)
            L0_soilId, & ! soil ID on L0 resolution
            L0_geoUnit, & ! hydrogeological class ID on input resolution (L0)
            timeStep_LAI_input, & ! flag on how LAI data has to be read
            iFlag_soilDB, & ! options to handle different types of soil databases
            nLAIclass, LAIUnitList, LAILUT, soilDB, L0_gridded_LAI, &
            nSoilHorizons_mHM, & ! soil horizons info for mHM
            nLAI
    use mo_common_variables, ONLY : processMatrix, & ! Info about which process runs in which option
            level0, & ! grid information (ncols, nrows, ..)
            dirMorpho, & ! directories
            dirCommonFiles, & ! directory of common files
            L0_Basin, & ! L0_Basin ID
            nBasins, & ! number of basins
            Grid, &
            period, &
            global_parameters                   ! global parameters
    USE mo_common_constants, ONLY : nodata_i4, nodata_dp, YearMonths_i4                   ! mHM's global nodata vales
    use mo_prepare_gridded_lai, only : prepare_gridded_daily_LAI_data, prepare_gridded_mean_monthly_LAI_data

    implicit none

    type(period), dimension(:), intent(in), optional :: LAIPer
    ! local variables
    integer(i4) :: iBasin, iVar, iHorizon, iMon, itimer, ll     ! loop variables
    integer(i4) :: nH                         ! dummy variable
    integer(i4) :: nunit                      ! file unit of file to read
    character(256) :: fName                      ! file name of file to read
    real(dp), dimension(:, :), allocatable :: data_dp_2d
    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    integer(i4), dimension(:, :), allocatable :: dataMatrix_i4
    logical, dimension(:, :), allocatable :: mask_2d
    integer(i4), dimension(:), allocatable :: dummy_i4
    type(Grid), pointer :: level0_iBasin

    ! min. value of slope and aspect
    real(dp), parameter :: slope_minVal = 0.01_dp
    real(dp), parameter :: aspect_minVal = 1.00_dp

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

    basins : do iBasin = 1, nBasins

      level0_iBasin => level0(L0_Basin(iBasin))

      call message('    Reading data for basin: ', trim(adjustl(num2str(iBasin))), ' ...')
      ! check whether L0 data is shared
      if (iBasin .gt. 1) then
        if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
          !
          call message('      Using data of basin ', &
                  trim(adjustl(num2str(L0_Basin(iBasin)))), ' for basin: ',&
                  trim(adjustl(num2str(iBasin))), '...')
          ! DO NOT read L0 data
          cycle

        end if
      end if

      itimer = 2
      call timer_start(itimer)

      ! read slope and aspect - datatype real
      nVars_real : do iVar = 1, 2
        select case (iVar)
        case(1) ! slope
          call message('      Reading slope ...')
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_slope))
          nunit = uslope
        case(2) ! aspect
          call message('      Reading aspect ...')
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_aspect))
          nunit = uaspect
        end select

        ! reading
        call read_spatial_data_ascii(trim(fName), nunit, &
                level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
                level0_iBasin%yllcorner, level0_iBasin%cellsize, data_dp_2d, mask_2d)
        ! put global nodata value into array (probably not all grid cells have values)
        data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d)
        ! put data in variable
        select case (iVar)
        case(1) ! slope
          call append(L0_slope, pack(data_dp_2d, level0_iBasin%mask))

        case(2) ! aspect
          call append(L0_asp, pack(data_dp_2d, level0_iBasin%mask))
        end select
        ! deallocate arrays
        deallocate(data_dp_2d, mask_2d)

      end do nVars_real

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
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_soilclass))
        else if(iFlag_soilDB .eq. 1) then
          write(fName, 172) iHorizon
          172             format('soil_class_horizon_', i2.2, '.asc')
          fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(fName))
        end if
        call read_spatial_data_ascii(trim(fName), usoilclass, &
                level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
                level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        call paste(dataMatrix_i4, pack(data_i4_2d, level0_iBasin%mask), nodata_i4)
        deallocate(data_i4_2d)
      end do
      call append(L0_soilId, dataMatrix_i4)
      deallocate(dataMatrix_i4)

      ! read geoUnit
      fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_hydrogeoclass))
      ! reading and transposing
      call read_spatial_data_ascii(trim(fName), uhydrogeoclass, &
              level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
              level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)
      ! put global nodata value into array (probably not all grid cells have values)
      data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
      call append(L0_geoUnit, pack(data_i4_2d, level0_iBasin%mask))
      deallocate(data_i4_2d, mask_2d)

      ! LAI values
      call message('      Reading LAI ...')
      select case (timeStep_LAI_input)
      case(1) ! long term mean monthly gridded fields
        call prepare_gridded_mean_monthly_LAI_data(iBasin, level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%mask)

      case(0) ! long term mean monthly values per class with LUT
        nLAI = YearMonths_i4
        fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_laiclass))
        ! reading and transposing
        call read_spatial_data_ascii(trim(fName), ulaiclass, &
                level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%xllcorner, &
                level0_iBasin%yllcorner, level0_iBasin%cellsize, data_i4_2d, mask_2d)
        ! put global nodata value into array (probably not all grid cells have values)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
        allocate(dummy_i4(count(level0_iBasin%mask)))
        dummy_i4 = pack(data_i4_2d, level0_iBasin%mask)
        deallocate(data_i4_2d, mask_2d)

        call check_consistency_lut_map(dummy_i4, LAIUnitList, file_laiclass)

        allocate(data_dp_2d(count(level0_iBasin%mask), nLAI))
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
        call prepare_gridded_daily_LAI_data(iBasin, level0_iBasin%nrows, level0_iBasin%ncols, level0_iBasin%mask, &
        LAIPer(iBasin))

      end select

      ! read lat lon coordinates of each basin
      call message('      Reading latitude/logitude ...')
      call read_latlon(iBasin, "lon_l0", "lat_l0", "level0", level0_iBasin)

      call timer_stop(itimer)
      call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

    end do basins
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

  !      NAME
  !         check_consistency_lut_map

  !     PURPOSE
  !>        \brief Checks if classes in input maps appear in look up tables.

  !>        \details Determines wether a class appearing in the morphological input
  !>                 is occuring in the respective look up table. mHM breaks if inconsistencies
  !>                 are discovered.
  !     INTENT(IN)
  !>        \param[in] " integer(i4), dimension(:) :: data"          map of study domain
  !>        \param[in] "integer(i4), dimension(:)  :: lookuptable"   look up table corresponding to map
  !>        \param[in] "character(*)               :: filename"      name of the lut file - ERORR warn 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), dimension(:), allocatable :: unique_values" array of unique values in data

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \none

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Nov 2016
  ! ------------------------------------------------------------------

  subroutine check_consistency_lut_map(data, lookuptable, filename, unique_values)

    USE mo_orderpack, ONLY : unista
    USE mo_string_utils, ONLY : num2str
    USE mo_message, ONLY : message

    implicit none

    integer(i4), dimension(:), intent(in) :: data          ! map of study domain
    integer(i4), dimension(:), intent(in) :: lookuptable   ! look up table corresponding to map
    character(*), intent(in) :: filename      ! name of the lut file - ERORR warn
    integer(i4), dimension(:), allocatable, intent(out), optional :: unique_values ! array of unique values in data

    ! local variables
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
