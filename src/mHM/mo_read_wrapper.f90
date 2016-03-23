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
  !                     Rohini Kumar,  Mar 2016  - options to handle different soil databases
  ! ------------------------------------------------------------------

  subroutine read_data
    !
    USE mo_read_lut,           ONLY: read_lai_lut,                        &
                                     read_geoformation_lut
    USE mo_soil_database,      ONLY: read_soil_LUT
    USE mo_read_spatial_data,  ONLY: read_header_ascii,                   &
                                     read_spatial_data_ascii
    USE mo_append,             ONLY: append, paste
    USE mo_string_utils,       ONLY: num2str
    USE mo_message,            ONLY: message
    !
    USE mo_file,               ONLY: file_geolut        , ugeolut,        & ! file name and unit of hydrogeology LuT
                                     file_lailut        , ulailut,        & ! file name and unit of LAI LuT
                                     file_dem           , udem,           & ! file name and unit of elevation map
                                     file_slope         , uslope,         & ! file name and unit of slope map
                                     file_aspect        , uaspect,        & ! file name and unit of aspect map
                                     file_soilclass     , usoilclass,     & ! file name and unit of soil class map
                                     file_hydrogeoclass , uhydrogeoclass, & ! file name and unit of hydrogeo class map
                                     file_laiclass      , ulaiclass,      & ! file name and unit of lai class map
                                     file_soil_database ,                 & ! file name of soil class map (iFlag_soilDB = 0)
                                     file_soil_database_1,                & ! file name of soil class map (iFlag_soilDB = 1)
                                     ulcoverclass                           ! unit of land cover class map
    USE mo_global_variables,   ONLY: nGeoUnits, GeoUnitList, GeoUnitKar,  & ! geological class information
                                     L0_Basin,                            & ! L0_Basin ID
                                     L0_mask,                             & ! global mask variable
                                     L0_elev,                             & ! elevation on input resolution (L0)
                                     L0_slope,                            & ! slope on input resolution (L0)
                                     L0_asp,                              & ! aspect on input resolution (L0)
                                     L0_soilId,                           & ! soil ID on L0 resolution
                                     L0_geoUnit,                          & ! hydrogeological class ID on input resolution (L0)
                                     L0_LCover_LAI,                       & ! LAI class ID on input resolution (L0)
                                     L0_LCover,                           & ! classical mHM land cover class (L0)
                                     dirMorpho, dirLCover,                & ! directories
                                     dirCommonFiles,                      & ! directory of common files
                                     LCfilename, nLCoverScene,            & ! file names and number of land cover scenes
                                     level0,                              & ! grid information (ncols, nrows, ..)
                                     nBasins,                             & ! number of basins
                                     basin,                               & ! basin information for single basins
                                     perform_mpr,                         & ! flag indicating whether L0 is read
                                     !timeStep_LAI_input,                 & ! flag on how LAI data has to be read
                                     iFlag_soilDB,                        & ! options to handle different types of soil databases
                                     nSoilHorizons_mHM,                   & ! soil horizons info for mHM
                                     resolutionHydrology,                 & ! hydrology resolution (L1 scale)
                                     nLAIclass, LAIUnitList, LAILUT,soilDB        
    USE mo_mhm_constants,      ONLY: nodata_i4, nodata_dp                   ! mHM's global nodata vales
    
    implicit none

    ! local variables
    integer(i4)                               :: iBasin, iVar, iHorizon     ! loop variables
    integer(i4)                               :: nH                         ! dummy variable
    integer(i4)                               :: nunit                      ! file unit of file to read
    integer(i4)                               :: nCells                     ! number of cells in global_mask
    character(256)                            :: fName                      ! file name of file to read
    real(dp), dimension(:,:), allocatable     :: data_dp_2d
    integer(i4), dimension(:,:), allocatable  :: data_i4_2d
    integer(i4), dimension(:,:), allocatable  :: dataMatrix_i4
    logical, dimension(:,:), allocatable      :: mask_2d
    logical, dimension(:,:), allocatable      :: mask_global
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
    if( iFlag_soilDB .eq. 0 ) then
       fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_soil_database))
    else if( iFlag_soilDB .eq. 1) then
       fName = trim(adjustl(dirCommonFiles)) // trim(adjustl(file_soil_database_1))
    end if
    call read_soil_LUT( trim(fName), iFlag_soilDB, soilDB)

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
             ! DO NOT read L0 data
             cycle
             !
          end if
       end if
       !
       call message('    Reading data for basin: ', trim(adjustl(num2str(iBasin))),' ...')
       !
       ! create overall mHM mask on L0 and save indices
       nCells = size(mask_global, dim=1)*size(mask_global, dim=2)
       call append( L0_mask, reshape(mask_global, (/nCells/)))
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
          
          ! put global nodata value into array (probably not all grid cells have values)
          data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global)
          ! put data in variable
          call append( L0_elev, pack(data_dp_2d, mask_global) )
          ! deallocate arrays
          deallocate(data_dp_2d)
          
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
             
             ! reading
             call read_spatial_data_ascii(trim(fName), nunit,                                     &
                  level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
                  level0%yllcorner(iBasin), level0%cellsize(iBasin), data_dp_2d, mask_2d)
             ! put global nodata value into array (probably not all grid cells have values)
             data_dp_2d = merge(data_dp_2d,  nodata_dp, mask_2d)
             ! put data in variable
             select case (iVar)
             case(1) ! slope
                call append( L0_slope, pack(data_dp_2d, mask_global) )
             case(2) ! aspect
                call append( L0_asp, pack(data_dp_2d, mask_global) )
             end select
             ! deallocate arrays
             deallocate(data_dp_2d, mask_2d)

          end do nVars_real

          ! read datatype integer
          ! ***** CHANGE IS MADE HERE TO ACCOMODATE READ SEVERAL TYPES OF SOIL DATABASES
          ! change from the earlier code where everything was done in the DO loop (**see below **)
          ! here everything is more explicit and no do loop appears as was the case in previous version
          
          ! read soilID in both options
          nH = 1 !> by default; when iFlag_soilDB = 0
          if(iFlag_soilDB .eq. 1) nH = nSoilHorizons_mHM
          ! modified way to read multiple horizons specific soil class
          do iHorizon = 1, nH
              if( iFlag_soilDB .eq. 0 ) then
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_soilclass)) 
             else if( iFlag_soilDB .eq. 1 ) then
                write(fName, 172) iHorizon
172             format('soil_class_horizon_',i2.2,'.asc')
                fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(fName))
             end if
             call read_spatial_data_ascii(trim(fName), usoilclass,                          &
                  level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
                  level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
             ! put global nodata value into array (probably not all grid cells have values)
             data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)               
             call paste(dataMatrix_i4, pack(data_i4_2d, mask_global), nodata_i4)
             deallocate(data_i4_2d)
          end do
          call append( L0_soilId, dataMatrix_i4 )
          deallocate(dataMatrix_i4)

          ! read geoUnit
          fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_hydrogeoclass))
          ! reading and transposing
          call read_spatial_data_ascii(trim(fName), uhydrogeoclass,                      &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
          call append( L0_geoUnit, pack(data_i4_2d, mask_global) )
          deallocate(data_i4_2d, mask_2d)
          
          ! read LAI related land cover class
          fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_laiclass))
          ! reading and transposing
          call read_spatial_data_ascii(trim(fName), ulaiclass,                           &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)
          call append( L0_LCover_LAI, pack(data_i4_2d, mask_global) )
          deallocate(data_i4_2d, mask_2d)


          ! ! read soilID, geoUnit, LAI - datatype integer
          ! nVars_integer: do iVar = 1, 3
          !    ! handle LAI options
          !    ! if( (iVar .EQ. 3)  .AND. (timeStep_LAI_input < 0) ) CYCLE
          !    select case (iVar)
          !    case(1) ! soil ID
          !       fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_soilclass))
          !       nunit = usoilclass
          !    case(2) ! geological ID
          !       fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_hydrogeoclass))
          !       nunit = uhydrogeoclass
          !    case(3) ! LAI classes
          !       fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_laiclass))
          !       nunit = ulaiclass
          !    end select

          !    ! reading and transposing
          !    call read_spatial_data_ascii(trim(fName), nunit,                               &
          !         level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
          !         level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
          !    ! put global nodata value into array (probably not all grid cells have values)
          !    data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d)

          !    ! put data into global L0 variable
          !    select case (iVar)
          !    case(1) ! soil class ID
          !       call append( L0_soilId,  pack(data_i4_2d, mask_global) )
          !    case(2) ! hydrogeological class ID
          !       call append( L0_geoUnit, pack(data_i4_2d, mask_global) )
          !    case(3) ! Land cover related to LAI classes
          !       call append( L0_LCover_LAI, pack(data_i4_2d, mask_global) )
          !    end select
             
          !    ! deallocate arrays
          !    deallocate(data_i4_2d, mask_2d)
          !    !
          ! end do nVars_integer
          !
       else
          ! if restart is switched on, perform dummy allocation of
          allocate( dummy_dp( count(mask_global) ) )
          allocate( dummy_i4( count(mask_global) ) )
          call append( L0_elev,     dummy_dp )
          call append( L0_slope,    dummy_dp )
          call append( L0_asp,      dummy_dp )
          ! for soil class 
          nH = 1 !> by default; when iFlag_soilDB = 0
          if ( iFlag_soilDB .eq. 1 ) nH = nSoilHorizons_mHM
          do iHorizon = 1, nH
             call paste(dataMatrix_i4, dummy_i4)
          end do
          call append( L0_soilId, dataMatrix_i4 )
          deallocate(dataMatrix_i4)
          !
          call append( L0_geoUnit,  dummy_i4 )
          deallocate( dummy_dp, dummy_i4 )
          
          ! read L0_LCover_LAI
          fName = trim(adjustl(dirMorpho(iBasin)))//trim(adjustl(file_laiclass))
          call read_spatial_data_ascii(trim(fName), ulaiclass,                           &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)
          call append( L0_LCover_LAI, pack(data_i4_2d, mask_global) )
          ! end if
          
       end if read_L0_data

       
       ! LCover read in is realized seperated because of unknown number of scenes
       do iVar = 1, nLCoverScene
          fName = trim(adjustl(dirLCover(iBasin)))//trim(adjustl(LCfilename(iVar)))
          call read_spatial_data_ascii(trim(fName), ulcoverclass,                        &
               level0%nrows(iBasin),     level0%ncols(iBasin), level0%xllcorner(iBasin), &
               level0%yllcorner(iBasin), level0%cellsize(iBasin), data_i4_2d, mask_2d)
          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)
          call paste(dataMatrix_i4, pack(data_i4_2d, mask_global), nodata_i4)
          deallocate(data_i4_2d)
       end do
       call append( L0_LCover, dataMatrix_i4 )
       deallocate(dataMatrix_i4)

       ! deallocate mask
       deallocate(mask_global)

    end do basins
    !----------------------------------------------------------------
    ! assign L0_mask to basin
    !----------------------------------------------------------------
    basin%L0_mask => L0_mask
    
    !----------------------------------------------------------------
    ! Correction for slope and aspect -- min value set above
    !----------------------------------------------------------------
    if ( perform_mpr ) then
       ! keep the colons (:) in the statements because of Intel's reallocation lhs problem
       L0_slope(:)  = merge(  slope_minVal, L0_slope(:),  (L0_slope(:)  .lt.  slope_minVal)  )
       L0_asp(:)    = merge( aspect_minVal, L0_asp(:),    (L0_asp(:)    .lt. aspect_minVal)  )
    end if

  end subroutine read_data
  
END MODULE mo_read_wrapper
