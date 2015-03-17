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

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PUBLIC  :: read_data            ! reads all available data
  PRIVATE :: rotate_fdir_variable ! rotate flow direction to the new coordinate system

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
  !                    Matthias Cuntz & Juliane Mai Nov 2014 - LAI input from daily, monthly or yearly files
  ! ------------------------------------------------------------------

  subroutine read_data
    !
    USE mo_read_lut,           ONLY: read_lai_lut,                        &
                                     read_geoformation_lut
    USE mo_soil_database,      ONLY: read_soil_LUT
    USE mo_read_spatial_data,  ONLY: read_header_ascii,                   &
                                     read_spatial_data_ascii
    USE mo_read_timeseries ,   ONLY: read_timeseries
    USE mo_append,             ONLY: append, paste
    USE mo_string_utils,       ONLY: num2str
    USE mo_message,            ONLY: message
    !
    USE mo_file,               ONLY: file_geolut        , ugeolut,        & ! file name and unit of hydrogeology LuT
                                     file_lailut        , ulailut,        & ! file name and unit of LAI LuT
                                     file_dem           , udem,           & ! file name and unit of elevation map
                                     file_slope         , uslope,         & ! file name and unit of slope map
                                     file_aspect        , uaspect,        & ! file name and unit of aspect map
                                     file_facc          , ufacc,          & ! file name and unit of flow acc map
                                     file_fdir          , ufdir,          & ! file name and unit of flow dir map
                                     file_soilclass     , usoilclass,     & ! file name and unit of soil class map
                                     file_hydrogeoclass , uhydrogeoclass, & ! file name and unit of hydrogeo class map
                                     file_gaugeloc      , ugaugeloc,      & ! file name and unit of gauge locations map
                                     file_laiclass      , ulaiclass,      & ! file name and unit of lai class map
                                     file_soil_database ,                 & ! file name and unit of soil class map
                                     udischarge         ,                 & ! unit of discharge time series
                                     ulcoverclass                           ! unit of land cover class map
    USE mo_global_variables,   ONLY: nGeoUnits, GeoUnitList, GeoUnitKar,  & ! geological class information
                                     L0_Basin,                            & ! L0_Basin ID
                                     L0_elev,                             & ! elevation on input resolution (L0)
                                     L0_slope,                            & ! slope on input resolution (L0)
                                     L0_asp,                              & ! aspect on input resolution (L0)
                                     L0_fAcc,                             & ! flow accumulation on input resolution (L0)
                                     L0_fDir,                             & ! flow direction on input resolution (L0)
                                     L0_soilId,                           & ! soil class ID on input resolution (L0)
                                     L0_geoUnit,                          & ! hydro class ID on input resolution (L0)
                                     L0_gaugeLoc,                         & ! location of evaluation gauges on input resolution (L0)
                                     L0_InflowGaugeLoc,                   & ! location of inflow gauges on input resolution (L0)
                                     L0_LCover_LAI,                       & ! LAI class ID on input resolution (L0)
                                     L0_LCover,                           & ! Normal land cover class ID on input resolution (L0)
                                     dirMorpho, dirLCover,                & ! directories
                                     dirCommonFiles,                      & ! directory of common files
                                     LCfilename, nLCover_scene,           & ! file names and number of land cover scenes
                                     level0,                              & ! grid information (ncols, nrows, ..)
                                     optimize,                            & ! optimizeation flag for some error checks
                                     nGaugesTotal, gauge, nMeasPerDay,    & ! evaluaton gauging station information
                                     nInflowGaugesTotal, InflowGauge,     & ! inflow stations information
                                     nBasins,                             & ! number of basins
                                     basin,                               & ! basin information for single basins
                                     evalPer,                             & ! model evaluation period (for discharge read in)
                                     simPer,                              & ! model simulation period (for inflow read in)
                                     processMatrix,                       & ! identify activated processes
                                     perform_mpr,                         & ! flag indicating whether L0 is read
                                     !timeStep_LAI_input,                  & ! flag on how LAI data has to be read
                                     resolutionHydrology                    ! hydrology resolution (L1 scale)
    USE mo_global_variables,   ONLY: nLAIclass, LAIUnitList, LAILUT,soilDB
    USE mo_mhm_constants,      ONLY: nodata_i4, nodata_dp                   ! mHM's global nodata vales
    
    implicit none

    ! local variables
    integer(i4)                               :: iGauge
    integer(i4)                               :: iBasin, iVar     ! loop variables
    integer(i4)                               :: nunit            ! file unit of file to read
    integer(i4)                               :: nCells           ! number of cells in global_mask
    character(256)                            :: fName            ! file name of file to read
    real(dp), dimension(:),   allocatable     :: data_dp_1d
    real(dp), dimension(:,:), allocatable     :: data_dp_2d
    integer(i4), dimension(:,:), allocatable  :: data_i4_2d
    integer(i4), dimension(:,:), allocatable  :: dataMatrix_i4
    logical, dimension(:),   allocatable      :: mask_1d
    logical, dimension(:,:), allocatable      :: mask_2d
    logical, dimension(:,:), allocatable      :: mask_global
    integer(i4), dimension(3)                 :: start_tmp, end_tmp
    integer(i4), dimension(:),  allocatable   :: dummy_i4
    real(dp),    dimension(:),  allocatable   :: dummy_dp

    ! min. value of slope and aspect
    real(dp), parameter                       :: slope_minVal  = 0.01_dp
    real(dp), parameter                       :: aspect_minVal = 1.00_dp

    ! ************************************************
    ! READ LOOKUP TABLES
    ! ************************************************
    !
    ! Soil LUT
    fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_soil_database))
    call read_soil_LUT( trim(fName), soilDB )

    ! Geological formation LUT
    fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_geolut))
    call read_geoformation_lut(trim(fName), ugeolut, nGeoUnits, GeoUnitList, GeoUnitKar)

    ! LAI LUT
    fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_lailut))
    call read_lai_lut(trim(fName), ulailut, nLAIclass, LAIUnitList, LAILUT)
    ! end if
    ! ************************************************
    ! READ SPATIAL DATA FOR EACH BASIN
    ! ************************************************

    ! allocate necessary variables at Level0
    allocate(level0%nrows       (nBasins))
    allocate(level0%ncols       (nBasins))
    allocate(level0%xllcorner   (nBasins))
    allocate(level0%yllcorner   (nBasins))
    allocate(level0%cellsize    (nBasins))
    allocate(level0%nodata_value(nBasins))
    !
    allocate(basin%L0_iStart    (nBasins))
    allocate(basin%L0_iEnd      (nBasins))
    allocate(basin%L0_iStartMask(nBasins))
    allocate(basin%L0_iEndMask  (nBasins))
    !
    ! allocate necessary variables at Level110
    allocate(basin%L110_iStart    (nBasins))
    allocate(basin%L110_iEnd      (nBasins))
    !
    basins: do iBasin = 1, nBasins
       !
       ! Header (to check consistency)
       fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_dem))
       call read_header_ascii(trim(fName), udem,   &
            level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
            level0%yllcorner(iBasin), level0%cellsize(iBasin), level0%nodata_value(iBasin))

       ! check for L0 and L1 scale consistency
       if( resolutionHydrology(iBasin) .LT. level0%cellsize(iBasin)) then
          call message()
          call message('***ERROR: resolutionHydrology (L1) should be smaller than the input data resolution (L0)')
          call message('          check set-up (in mhm.nml) for basin: ', trim(adjustl(num2str(iBasin))),' ...')
          stop
       end if

       !
       ! DEM + overall mask creation
       fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_dem))
       call read_spatial_data_ascii(trim(fName), udem, &
            level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin),&
            level0%yllcorner(iBasin), level0%cellsize(iBasin), data_dp_2d, mask_global)
       !
       ! Saving indices at Level110 irrespective of whether L0_data is shared or not
       if (iBasin .eq. 1) then
          basin%L110_iStart(iBasin) = 1
          basin%L110_iEnd  (iBasin) = basin%L110_iStart(iBasin) + count(mask_global) - 1
        else
          basin%L110_iStart(iBasin) = basin%L110_iEnd(iBasin-1) + 1
          basin%L110_iEnd  (iBasin) = basin%L110_iStart(iBasin) + count(mask_global) - 1
       end if
       !
       ! check whether L0 data is shared
       if (iBasin .gt. 1) then
          if (L0_Basin(iBasin) .eq. L0_Basin(iBasin - 1)) then
             !
             call message('    Using data of previous basin: ', trim(adjustl(num2str(iBasin))),' ...')
             basin%L0_iStart(iBasin) = basin%L0_iStart(iBasin - 1)
             basin%L0_iEnd  (iBasin) = basin%L0_iEnd(iBasin - 1)
             !
             basin%L0_iStartMask(iBasin) = basin%L0_iStartMask(iBasin - 1 )
             basin%L0_iEndMask  (iBasin) = basin%L0_iEndMask(iBasin - 1 )
             !
             ! DONT read L0 data
             cycle
             !
          end if
       end if
       !
       call message('    Reading data for basin: ', trim(adjustl(num2str(iBasin))),' ...')
       !
       ! create overall mHM mask on L0 and save indices
       nCells = size(mask_global, dim=1)*size(mask_global, dim=2)
       call append( basin%L0_mask, reshape(mask_global, (/nCells/)))
       !
       ! Saving indices of mask and packed data
       if(iBasin .eq. 1) then
          basin%L0_iStart(iBasin) = 1
          basin%L0_iEnd  (iBasin) = basin%L0_iStart(iBasin) + count(mask_global) - 1
          !
          basin%L0_iStartMask(iBasin) = 1
          basin%L0_iEndMask  (iBasin) = basin%L0_iStartMask(iBasin) + nCells - 1
       else
          basin%L0_iStart(iBasin) = basin%L0_iEnd(iBasin-1) + 1
          basin%L0_iEnd  (iBasin) = basin%L0_iStart(iBasin) + count(mask_global) - 1
          !
          basin%L0_iStartMask(iBasin) = basin%L0_iEndMask(iBasin-1) + 1
          basin%L0_iEndMask  (iBasin) = basin%L0_iStartMask(iBasin) + nCells - 1
       end if
       ! Read L0 data, if restart is false
       read_L0_data: if ( perform_mpr ) then
          !
          ! put global nodata value into array (probably not all grid cells have values)
          data_dp_2d = merge(data_dp_2d,  nodata_dp, mask_global)
          ! put data in variable
          call append( L0_elev, pack(data_dp_2d, mask_global) )
          ! deallocate arrays
          deallocate(data_dp_2d)
          !
          ! read slope and aspect - datatype real
          nVars_real: do iVar = 1, 2
             select case (iVar)
             case(1) ! slope
                fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl(file_slope))
                nunit = uslope
             case(2) ! aspect
                fName = trim(adjustl(dirMorpho(iBasin)))// trim(adjustl(file_aspect))
                nunit = uaspect
             end select
             !
             ! reading
             call read_spatial_data_ascii(trim(fName), nunit,                                     &
                  level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
                  level0%yllcorner(iBasin), level0%cellsize(iBasin), data_dp_2d, mask_2d)
             !
             ! put global nodata value into array (probably not all grid cells have values)
             data_dp_2d = merge(data_dp_2d,  nodata_dp, mask_2d)
             !
             ! put data in variable
             select case (iVar)
             case(1) ! slope
                call append( L0_slope, pack(data_dp_2d, mask_global) )
             case(2) ! aspect
                call append( L0_asp, pack(data_dp_2d, mask_global) )
             end select
             !
             ! deallocate arrays
             deallocate(data_dp_2d, mask_2d)
             !
          end do nVars_real
          !
          ! read fAcc, fDir, soilID, geoUnit, gaugeLoc, LAI, LCover - datatype integer
          nVars_integer: do iVar = 1, 6

             ! handle routing related input data
             if( ( (iVar .EQ. 1) .or. (iVar .EQ. 2) .or. (iVar .EQ. 5) ) .and. &
                  (processMatrix(8, 1) .EQ. 0)  ) CYCLE

             ! handle LAI options
             ! if( (iVar .EQ. 6)  .AND. (timeStep_LAI_input < 0) ) CYCLE

             select case (iVar)
             case(1) ! flow accumulation
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_facc))
                nunit = ufacc
             case(2) ! flow direction
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_fdir))
                nunit = ufdir
             case(3) ! soil ID
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_soilclass))
                nunit = usoilclass
             case(4) ! geological ID
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_hydrogeoclass))
                nunit = uhydrogeoclass
             case(5) ! location of gauging stations
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_gaugeloc))
                nunit = ugaugeloc
             case(6) ! LAI classes
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_laiclass))
                nunit = ulaiclass
             end select

             !
             ! reading and transposing
             call read_spatial_data_ascii(trim(fName), nunit,                               &
                  level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
                  level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)

             ! put global nodata value into array (probably not all grid cells have values)
             data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)

             ! put data into global L0 variable
             select case (iVar)
             case(1) ! flow accumulation
                call append( L0_fAcc,    pack(data_i4_2d, mask_global) )
             case(2) ! flow direction
                ! rotate flow direction and any other variable with directions
                ! NOTE: ONLY when ASCII files are read
                call rotate_fdir_variable(data_i4_2d)
                ! append
                call append( L0_fDir,    pack(data_i4_2d, mask_global) )
             case(3) ! soil class ID
                call append( L0_soilId,  pack(data_i4_2d, mask_global) )
             case(4) ! hydrogeological class ID
                call append( L0_geoUnit, pack(data_i4_2d, mask_global) )
             case(5) ! location of evaluation and inflow gauging stations
                ! evaluation gauges
                ! Input data check
                do iGauge = 1, basin%nGauges(iBasin)
                   ! If gaugeId is found in gauging location file?
                   if (.not. any(data_i4_2d .EQ. basin%gaugeIdList(iBasin, iGauge))) then
                      call message()
                      call message('***ERROR: Gauge ID "', trim(adjustl(num2str(basin%gaugeIdList(iBasin, iGauge)))), &
                           '" not found in ' )
                      call message('          Gauge location input file: ', &
                           trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_gaugeloc)))
                      stop
                   end if
                end do

                call append( L0_gaugeLoc, pack(data_i4_2d, mask_global) )

                ! inflow gauges
                ! if no inflow gauge for this subbasin exists still matirx with dim of subbasin has to be paded
                if (basin%nInflowGauges(iBasin) .GT. 0_i4) then
                   ! Input data check
                   do iGauge = 1, basin%nInflowGauges(iBasin)
                      ! If InflowGaugeId is found in gauging location file?
                      if (.not. any(data_i4_2d .EQ. basin%InflowGaugeIdList(iBasin, iGauge))) then
                         call message()
                         call message('***ERROR: Inflow Gauge ID "', &
                              trim(adjustl(num2str(basin%InflowGaugeIdList(iBasin, iGauge)))), &
                              '" not found in ' )
                         call message('          Gauge location input file: ', &
                              trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_gaugeloc)))
                         stop
                      end if
                   end do
                end if

                call append( L0_InflowGaugeLoc, pack(data_i4_2d, mask_global) )

             case(6) ! Land cover related to LAI classes
                call append( L0_LCover_LAI, pack(data_i4_2d, mask_global) )
             end select
             !
             ! deallocate arrays
             deallocate(data_i4_2d, mask_2d)
             !
          end do nVars_integer
          !
       else
          ! if restart is switched on, perform dummy allocation of
          allocate( dummy_dp( count(mask_global) ) )
          allocate( dummy_i4( count(mask_global) ) )
          call append( L0_elev,     dummy_dp )
          call append( L0_slope,    dummy_dp )
          call append( L0_asp,      dummy_dp )
          call append( L0_fAcc,     dummy_i4 )
          call append( L0_fDir,     dummy_i4 )
          call append( L0_soilId,   dummy_i4 )
          call append( L0_geoUnit,  dummy_i4 )
          call append( L0_gaugeLoc, dummy_i4 )
          call append( L0_InflowGaugeLoc, dummy_i4 )
          deallocate( dummy_dp, dummy_i4 )
          ! read L0_LCover_LAI
          fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_laiclass))
          nunit = ulaiclass
          call read_spatial_data_ascii(trim(fName), nunit,                               &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)

          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)
          call append( L0_LCover_LAI, pack(data_i4_2d, mask_global) )
          ! end if
       end if read_L0_data

       !
       ! LCover read in is realized seperated because of unknown number of scenes
       do iVar = 1, nLCover_scene
          fName = trim(adjustl(dirLCover(iBasin)))//trim(adjustl(LCfilename(iVar)))
          call read_spatial_data_ascii(trim(fName), ulcoverclass,                        &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)

          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)
          call paste(dataMatrix_i4, pack(data_i4_2d, mask_global), nodata_i4)
          !
          deallocate(data_i4_2d)
       end do
       !
       call append( L0_LCover, dataMatrix_i4 )
       !
       deallocate(mask_global)
       deallocate(dataMatrix_i4)

    end do basins
    !----------------------------------------------------------------
    ! Correction for slope and aspect -- min value set above
    !----------------------------------------------------------------
    if ( perform_mpr ) then
       ! keep the colons (:) in the statements because of Intel's reallocation lhs problem
       L0_slope(:)  = merge(  slope_minVal, L0_slope(:),  (L0_slope(:)  .lt.  slope_minVal)  )
       L0_asp(:)    = merge( aspect_minVal, L0_asp(:),    (L0_asp(:)    .lt. aspect_minVal)  )
    end if

    ! ************************************************
    ! READ DISCHARGE TIME SERIES
    ! ************************************************
    !
    ! processMatrix(8,1) - process(8)=discharge
    if( processMatrix(8,1) .GE. 1 ) then
       !
       do iGauge = 1, nGaugesTotal
          ! get basin id
          iBasin = gauge%basinId(iGauge)
          ! get start and end dates
          start_tmp = (/evalPer(iBasin)%yStart, evalPer(iBasin)%mStart, evalPer(iBasin)%dStart/)
          end_tmp   = (/evalPer(iBasin)%yEnd,   evalPer(iBasin)%mEnd,   evalPer(iBasin)%dEnd  /)
          ! evaluation gauge
          fName = trim(adjustl(gauge%fname(iGauge)))
          call read_timeseries(trim(fName), udischarge, &
               start_tmp, end_tmp, optimize, &
               data_dp_1d, mask=mask_1d, nMeasPerDay=nMeasPerDay)
          data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
          call paste(gauge%Q, data_dp_1d, nodata_dp )
          deallocate (data_dp_1d)
       end do
       !
    end if


    ! inflow gauge
    !
    ! in mhm call InflowGauge%Q has to be initialized -- dummy allocation with period of basin 1 and initialization
    if (nInflowGaugesTotal .EQ. 0) then
       allocate( data_dp_1d( maxval( simPer(:)%julEnd  - simPer(:)%julStart + 1 ) ) )
       data_dp_1d = nodata_dp
       call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
    else

       do iGauge = 1, nInflowGaugesTotal
          ! get basin id
          iBasin = InflowGauge%basinId(iGauge)
          ! get start and end dates
          start_tmp = (/simPer(iBasin)%yStart, simPer(iBasin)%mStart, simPer(iBasin)%dStart/)
          end_tmp   = (/simPer(iBasin)%yEnd,   simPer(iBasin)%mEnd,   simPer(iBasin)%dEnd  /)
          ! inflow gauge
          fName = trim(adjustl(InflowGauge%fname(iGauge)))
          call read_timeseries(trim(fName), udischarge, &
               start_tmp, end_tmp, optimize, &
               data_dp_1d, mask=mask_1d, nMeasPerDay=nMeasPerDay)
          if ( .NOT. (all(mask_1d)) ) then
             call message()
             call message('***ERROR: Nodata values in inflow gauge time series. File: ', trim(fName))
             call message('          During simulation period from ', num2str(simPer(iBasin)%yStart) &
                  ,' to ', num2str(simPer(iBasin)%yEnd))
             stop
          end if
          data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
          call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
          deallocate (data_dp_1d)
       end do
    end if

  end subroutine read_data

  ! ------------------------------------------------------------------
  ! Rotate fdir variable to the new coordinate system
  ! L. Samaniego & R. Kumar
  ! ------------------------------------------------------------------
  subroutine rotate_fdir_variable( x )
    USE mo_mhm_constants,      ONLY: nodata_i4    ! mHM's global nodata vales

    implicit none

    integer(i4), dimension(:,:), intent(INOUT) :: x
    ! local
    integer(i4)                                :: i, j

    !-------------------------------------------------------------------
    ! NOTE:
    !
    !     Since the DEM was transposed from (lat,lon) to (lon,lat), i.e.
    !     new DEM = transpose(DEM_old), then
    !
    !     the flow direction X (which was read) for every i, j cell needs
    !     to be rotated as follows
    !
    !                 X(i,j) = [R] * {uVector}
    !
    !     with
    !     {uVector} denoting the fDir_old (read) in vector form, and
    !               e.g. dir 8 => {-1, -1, 0 }
    !     [R]       denting a full rotation matrix to transform the flow
    !               direction into the new coordinate system (lon,lat).
    !
    !     [R] = [rx][rz]
    !
    !     with
    !           |      1       0      0 |
    !     [rx] =|      0   cos T  sin T | = elemental rotation along x axis
    !           |      0  -sin T  cos T |
    !
    !           |  cos F   sin F      0 |
    !     [rz] =| -sin F   cos F      0 | = elemental rotation along z axis
    !           |      0       0      1 |
    !
    !     and T = pi, F = - pi/2
    !     thus
    !          !  0  -1   0 |
    !     [R] =| -1   0   0 |
    !          |  0   0  -1 |
    !     making all 8 directions the following transformation were
    !     obtained.
    !-------------------------------------------------------------------

    do i = 1, size(x,1)
       do j = 1, size(x,2)
          if ( x(i,j)  .eq. nodata_i4 ) cycle
          select case ( x(i,j) )
          case(1)
             x(i,j) =   4
          case(2)
             x(i,j) =   2
          case(4)
             x(i,j) =   1
          case(8)
             x(i,j) = 128
          case(16)
             x(i,j) =  64
          case(32)
             x(i,j) =  32
          case(64)
             x(i,j) =  16
          case(128)
             x(i,j) =   8
          end select
       end do
    end do

end subroutine rotate_fdir_variable

END MODULE mo_read_wrapper
