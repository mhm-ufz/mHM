!>       \file mo_mhm_mpr_interface.f90
!>       \brief calls the MPR library for parameter estimation of mHM parameters
!>       \details - MPR is called with global parameters and path to its namelist
!>                - MPR reads all the land surface properties on the level 0 grid(s)
!>                - it then estimates the parameters and scaled them to level 1 (mHM grid (hydrology_resolution))
!>                - parameter names have to conform with mHM parameter names (e.g. start with "L1_")
!>                - level1 grid is inferred from MPR coordinates 'lon_out' and 'lat_out'
!>                - level1 grid mask is inferred from field 'L1_latitude'
!>                - land cover scenes are inferred from MPR coordinate 'land_cover_period_out'
!>                - soil horizons are inferred from MPR coordinate 'horizon_out'
!>                - lai periods are inferred from 3rd coordinate of field 'L1_Max_Canopy_Intercept' (name can be flexible)
!>                - mHM parameters required are:
!>                   - L1_SealedFraction            (rows, cols, lcscenes)
!>                   - L1_Alpha                     (rows, cols, lcscenes)
!>                   - L1_DegDayInc                 (rows, cols, lcscenes)
!>                   - L1_DegDayMax                 (rows, cols, lcscenes)
!>                   - L1_DegDayNoPre               (rows, cols, lcscenes)
!>                   - L1_KarstLoss                 (rows, cols)
!>                   - L1_Max_Canopy_Intercept      (rows, cols, lais)
!>                   - L1_FastFlow                  (rows, cols, lcscenes)
!>                   - L1_SlowFlow                  (rows, cols, lcscenes)
!>                   - L1_kBaseFlow                 (rows, cols, lcscenes)
!>                   - L1_Kperco                    (rows, cols, lcscenes)
!>                   - L1_FieldCap                  (rows, cols, horizons, lcscenes)
!>                   - L1_SatSoilMoisture           (rows, cols, horizons, lcscenes)
!>                   - L1_SoilMoistureExponent      (rows, cols, horizons, lcscenes)
!>                   - L1_PermWiltPoint             (rows, cols, horizons, lcscenes)
!>                   - L1_TempThresh                (rows, cols, lcscenes)
!>                   - L1_UnsatThreshold            (rows, cols)
!>                   - L1_SealedThresh              (rows, cols)
!>                   - L1_Jarvis_Threshold          (rows, cols)
!>                   - L1_PET_LAI_correction_factor (rows, cols, lais, lcscenes)
!>                   - L1_HarSamCoeff               (rows, cols)
!>                   - L1_PrieTayAlpha              (rows, cols, lais)
!>                   - L1_Aerodyn_resist            (rows, cols, lais, lcscenes)
!>                   - L1_Bulk_Surface_Resist       (rows, cols, lais)
!>                   - L1_fAsp                      (rows, cols)
!>                   - L1_latitude                  (rows, cols)
!>                   - L1_fRoots_FC (case 3 == 3)   (rows, cols, horizons, lcscenes)
!>                   - L1_fRoots (case 3 != 3)      (rows, cols, horizons, lcscenes)

module mo_mhm_mpr_interface

  ! some imports relevant for all routines in the module
  use mo_mpr_data_array, only : MPR_DATA_ARRAYS
  use mo_mpr_constants, only : maxNameLength
  use mo_kind, only : i4, dp
  implicit none

  private

  public :: call_mpr

contains

  subroutine call_mpr(parameterValues, parameterNames, grids, doInitArg)
    !< this calls the external mpr library
    !< the parameters stored in MPR_DATA_ARRAYS are explicitly read and iterated over
    !< coordinate information are inferred implicitly through attributes of MPR_DATA_ARRAYS

    use mo_mpr_global_variables, only : out_filename
    use mo_mpr_read_config, only : mpr_read_config
    use mo_mpr_constants, only : maxNoDataArrays
    use mo_mpr_reset, only: reset

    use mo_netcdf, only : NcDataset
    use mo_common_variables, only : write_restart, dummy_global_parameters, dummy_global_parameters_name, domainMeta
    use mo_append, only: append
    use mo_file, only : unamelist_mpr
    use mo_global_variables, only: pathMprNml, are_parameter_initialized
    use mo_grid, only: Grid

    !> the vector of parameter values passed from mHM
    real(dp), dimension(:), intent(in) :: parameterValues
    !> the vector of parameter names passed from mHM
    character(*), dimension(:), intent(in) :: parameterNames
    !> the Grid type that stores the lat-lon grid of the parameters
    type(Grid), dimension(:), intent(inout) :: grids
    !> boolean flag indicating whether initalize grid, land_cover and LAI period coordinates
    logical, intent(in), optional :: doInitArg

    real(dp), dimension(:), allocatable :: parameterValuesConcat
    character(maxNameLength), dimension(:), allocatable :: parameterNamesConcat
    integer(i4), dimension(:), allocatable :: landCoverSelect, laiSelect
    integer(i4) :: iDomain, previousDomain, iDA, nLandCoverPeriodsTemp, nLaisTemp
    type(NcDataset) :: nc

    ! as the mHM configuration might not provide parameters for all data arrays in the mpr.nml -
    ! dummy values for the global parameters missing are inserted here
    parameterValuesConcat = parameterValues
    parameterNamesConcat = parameterNames
    call append(parameterValuesConcat, dummy_global_parameters)
    call append(parameterNamesConcat, dummy_global_parameters_name)

    ! loop over all domains
    previousDomain = 0
    do iDomain = 1, domainMeta%nDomains
      
      if (allocated(landCoverSelect)) deallocate(landCoverSelect)
      if (allocated(laiSelect)) deallocate(laiSelect)
      ! if L0 data is shared ...
      if (domainMeta%L0DataFrom(iDomain) == previousDomain) then
        ! use grid from previousDomain
        grids(iDomain) = grids(previousDomain)
        ! get new landCoverperiod indices
        call init_grid(grids(iDomain), iDomain, landCoverSelect, nLandCoverPeriodsTemp, &
            laiSelect, nLaisTemp, .false.)
        ! use parameters from previousDomain but use iDomains landCoverSeelct (based also on simulation period)
        call init_eff_params_from_data_arrays(grids(iDomain), nLandCoverPeriodsTemp, landCoverSelect, &
                nLaisTemp, laiSelect, doInsertArg=doInitArg)
        cycle
      end if
      
      previousDomain = domainMeta%L0DataFrom(iDomain)
      ! delete MPR_DATA_ARRAYS
      call reset()
      ! read the configuration file
      call mpr_read_config(pathMprNml(iDomain), unamelist_mpr, &
              parameterValues=parameterValuesConcat, parameterNames=parameterNamesConcat)
    
      ! reorder data array to minimize storage requirement and assure that all required input fields have been read before
      ! TODO: add feature to mpr code, that reorder_data_arrays also considers the init of coordinates
      !   (this is the reason that this call is not active, as coordinates are not initialized properly)
      ! call reorder_data_arrays(MPR_DATA_ARRAYS)
    
      if (write_restart) then
        ! open the output file
        nc = NcDataset(out_filename, "w")
      end if
    
      do iDA = 1, maxNoDataArrays
        ! stop if no name is initialized by the configuration file
        if (.not. MPR_DATA_ARRAYS(iDA)%is_initialized) then
          cycle
        end if
        ! execute MPR on each DataArray
        call MPR_DATA_ARRAYS(iDA)%execute()
        if (write_restart .and. MPR_DATA_ARRAYS(iDA)%toFile) then
          ! write the current parameter to the output file
          call MPR_DATA_ARRAYS(iDA)%write(nc)
        end if
      end do
    
      if (write_restart) then
        ! close the output file
        call nc%close()
      end if
    
      ! TODO: MPR we need to make sure there is no data gaps in data arrays... (former L0_check_input)
      ! - are all parameters initialized?
      ! - same grids etc.
      call init_grid(grids(iDomain), iDomain, landCoverSelect, nLandCoverPeriodsTemp, &
              laiSelect, nLaisTemp, doInitArg)
      ! set parameters (moves values from MPR_DATA_ARRAYS to mo_global_variables data structures)
      call init_eff_params_from_data_arrays(grids(iDomain), nLandCoverPeriodsTemp, landCoverSelect, &
              nLaisTemp, laiSelect, doInsertArg=doInitArg)
    end do
    ! delete MPR_DATA_ARRAYS
    call reset()
    ! some sanity checks
    call check_consistency()
    are_parameter_initialized = .true.

  end subroutine call_mpr

  subroutine init_eff_params_from_data_arrays(grid, nLandCoverPeriodsRaw, landCoverSelect, &
          nlaiPeriodsRaw, laiSelect, doInsertArg)
    !< fills and appends the effective parameters (mo_global_variables) based on MPR_DATA_ARRAYS
    !< those parameters are stored in a lat-lon packed mode and are appended on this packed dimension (cells)
    !< however, as land cover periods and LAI periods can vary over domains, empty slices might have to be inserted
    use mo_global_variables, only: nSoilHorizons, L1_fSealed, L1_alpha, L1_degDayInc, L1_degDayMax, &
            L1_degDayNoPre, L1_karstLoss, L1_maxInter, L1_kFastFlow, L1_kSlowFlow, &
            L1_kBaseFlow, L1_kPerco, L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp, L1_jarvis_thresh_c1, &
            L1_petLAIcorFactor, L1_tempThresh, L1_unsatThresh, L1_sealedThresh, L1_wiltingPoint, L1_fAsp, &
            L1_latitude, L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_fRoots, L1_surfResist
    use mo_common_variables, only: processMatrix
    use mo_grid, only: grid_type => Grid
    use mo_common_datetime_type, only: nLandCoverPeriods, nlaiPeriods
    use mo_mpr_constants, only: maxNoDataArrays
    use mo_append, only: append, add_nodata_slice
    use mo_constants, only: nodata_dp
    use mo_mrm_global_variables, only: L0_slope

    !> grid that contains critical information for appending values for domains >=2
    type(grid_type), intent(in) :: grid
    !> how many land cover periods are present in MPR_DATA_ARRAYS
    integer(i4), intent(in) :: nLandCoverPeriodsRaw
    !> indices which land cover periods to select from MPR_DATA_ARRAYS
    integer(i4), dimension(:), allocatable, intent(in) :: landCoverSelect
    !> how many LAI periods are present in MPR_DATA_ARRAYS
    integer(i4), intent(in) :: nlaiPeriodsRaw
    !> indices which LAI periods to select from MPR_DATA_ARRAYS
    integer(i4), dimension(:), allocatable, intent(in) :: laiSelect
    !> optional flag indicating whether
    logical, intent(in), optional :: doInsertArg
    
    integer(i4) :: nCells, iDA, s1, e1, newLandCoverSlices, newlaiSlices, iLC
    logical :: doInsert
    real(dp), dimension(:, :, :), allocatable :: dummy3D, dummy3D_temp
    real(dp), dimension(:, :), allocatable :: dummy2D
  
    doInsert = .true.
    if (present(doInsertArg)) doInsert = doInsertArg
  
    nCells = grid%nCells
    newLandCoverSlices = size(landCoverSelect)
    newlaiSlices = size(laiSelect)
    ! init the parameters
    do iDA = 1, maxNoDataArrays
      ! stop if no name is initialized
      if (MPR_DATA_ARRAYS(iDA)%name(1:3) == 'L1_') then
        if (.not. doInsert) then
          s1 = grid%iStart
          e1 = grid%iEnd
          ! The data arrays are selected based on name
          ! they only exist in a packed 1D-version and thus first need to be reshaped
          ! we know from landCoverSelect that we only need a certain slices of the landCoverPeriods
          ! we select the slices after the reshape by using newLandCoverSlices
          ! then, we need to have nLandCoverPeriods in the arrays last dimensions, newLandCoverSlices may be different
          ! to fulfill that, we add dummy slices in the last dimension
          select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_SealedFraction')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_fSealed(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_Alpha')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_alpha(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayInc')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_degDayInc(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayMax')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_degDayMax(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayNoPre')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_degDayNoPre(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_KarstLoss')
              ! rows, cols
              L1_karstLoss(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_Max_Canopy_Intercept')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              L1_maxInter(s1:e1, 1:newlaiSlices) = dummy2D(:, laiSelect)
            case('L1_FastFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_kFastFlow(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_SlowFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_kSlowFlow(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_kBaseFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_kBaseFlow(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_Kperco')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_kPerco(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_FieldCap')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_soilMoistFC(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_SatSoilMoisture')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_soilMoistSat(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_SoilMoistureExponent')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_soilMoistExp(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_PermWiltPoint')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_wiltingPoint(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_TempThresh')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_tempThresh(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_UnsatThreshold')
              ! rows, cols
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              L1_unsatThresh(s1:e1, 1:newLandCoverSlices) = dummy2D(:, landCoverSelect)
            case('L1_SealedThresh')
              ! rows, cols
              L1_sealedThresh(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_Jarvis_Threshold')
              ! rows, cols
              L1_jarvis_thresh_c1(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_PET_LAI_correction_factor')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw, nLandCoverPeriodsRaw])
              L1_petLAIcorFactor(s1:e1, 1:newlaiSlices, 1:newLandCoverSlices) = dummy3D(:, laiSelect, landCoverSelect)
            case('L1_HarSamCoeff')
              ! rows, cols
              L1_HarSamCoeff(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_PrieTayAlpha')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              L1_PrieTayAlpha(s1:e1, 1:newlaiSlices) = dummy2D(:, laiSelect)
            case('L1_Aerodyn_resist')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw, nLandCoverPeriodsRaw])
              L1_aeroResist(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, laiSelect, landCoverSelect)
            case('L1_Bulk_Surface_Resist')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              L1_surfResist(s1:e1, 1:newlaiSlices) = dummy2D(:, laiSelect)
            case('L1_fAsp')
              ! rows, cols
              L1_fAsp(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_latitude')
              ! rows, cols
              L1_latitude(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
          end select
          ! process-sensitive types of calculation for parameters
          if (processMatrix(3,1) == 3) then
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots_FC')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_fRoots(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            end select
          else
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              L1_fRoots(s1:e1, :, 1:newLandCoverSlices) = dummy3D(:, :, landCoverSelect)
            end select
          end if
        else
          select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_SealedFraction')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_fSealed, dummy2D)
            case('L1_Alpha')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_alpha, dummy2D)
            case('L1_DegDayInc')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_degDayInc, dummy2D)
            case('L1_DegDayMax')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_degDayMax, dummy2D)
            case('L1_DegDayNoPre')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_degDayNoPre, dummy2D)
            case('L1_KarstLoss')
              ! rows, cols
              call append(L1_karstLoss, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_Max_Canopy_Intercept')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              dummy2D = dummy2D(:, laiSelect)
              call add_nodata_slice(dummy2D, nlaiPeriods - newlaiSlices, nodata_dp)
              call append(L1_maxInter, dummy2D)
            case('L1_FastFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_kFastFlow, dummy2D)
            case('L1_SlowFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_kSlowFlow, dummy2D)
            case('L1_kBaseFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_kBaseFlow, dummy2D)
            case('L1_Kperco')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_kPerco, dummy2D)
            case('L1_FieldCap')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_soilMoistFC, dummy3D, idim=1)
            case('L1_SatSoilMoisture')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_soilMoistSat, dummy3D, idim=1)
            case('L1_SoilMoistureExponent')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_soilMoistExp, dummy3D, idim=1)
            case('L1_PermWiltPoint')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_wiltingPoint, dummy3D, idim=1)
            case('L1_TempThresh')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_tempThresh, dummy2D)
            case('L1_UnsatThreshold')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLandCoverPeriodsRaw])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_unsatThresh, dummy2D)
            case('L1_SealedThresh')
              ! rows, cols
              call append(L1_sealedThresh, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_Jarvis_Threshold')
              ! rows, cols
              call append(L1_jarvis_thresh_c1, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_PET_LAI_correction_factor')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_petLAIcorFactor, dummy3D)
            case('L1_HarSamCoeff')
              ! rows, cols
              call append(L1_HarSamCoeff, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_PrieTayAlpha')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              dummy2D = dummy2D(:, laiSelect)
              call add_nodata_slice(dummy2D, nlaiPeriods - newlaiSlices, nodata_dp)
              call append(L1_PrieTayAlpha, dummy2D)
            case('L1_Aerodyn_resist')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, laiSelect, landCoverSelect)
              allocate(dummy3D_temp(nCells, nlaiPeriods, nLandCoverPeriodsRaw))
              do iLC=1, size(landCoverSelect)
                dummy2D = dummy3D(:,:,iLC)
                call add_nodata_slice(dummy2D, nlaiPeriods - newlaiSlices, nodata_dp)
                dummy3D_temp(:, :, iLC) = dummy2D
              end do
              call add_nodata_slice(dummy3D_temp, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_aeroResist, dummy3D_temp)
            case('L1_Bulk_Surface_Resist')
              ! rows, cols, lais
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nlaiPeriodsRaw])
              dummy2D = dummy2D(:, laiSelect)
              call add_nodata_slice(dummy2D, nlaiPeriods - newlaiSlices, nodata_dp)
              call append(L1_surfResist, dummy2D)
            case('L1_fAsp')
              ! rows, cols
              call append(L1_fAsp, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_latitude')
              ! rows, cols
              call append(L1_latitude, MPR_DATA_ARRAYS(iDA)%data)
          end select
          ! process-sensitive types of calculation for parameters
          if (processMatrix(3,1) == 3) then
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots_FC')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_fRoots, dummy3D, idim=1)
            end select
          else
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nLandCoverPeriodsRaw])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newLandCoverSlices, nodata_dp)
              call append(L1_fRoots, dummy3D, idim=1)
            end select
          end if
        end if
      else if (trim(MPR_DATA_ARRAYS(iDA)%name) == 'slope') then
        call append(L0_slope, MPR_DATA_ARRAYS(iDA)%data)
      end if
    end do
  
  end subroutine init_eff_params_from_data_arrays

  subroutine init_grid(newGrid, iDomain, landCoverSelect, nLandCoverPeriodsTemp, laiSelect, nLaisTemp, doInitArg)
    !< initialize the lat-lon grid, the land cover period and lai period objects and return the information
    !< for those periods to select (n and ids)
    use mo_common_variables, only : Grid
    use mo_grid, only : init_advanced_grid_properties
    use mo_mpr_coordinate, only: get_index_in_coordinate, MPR_COORDINATES
    use mo_mpr_utils, only: get_index_in_vector
    use mo_mpr_constants, only : maxNameLength
    use mo_read_nc, only: check_soil_dimension_consistency
    use mo_common_datetime_type, only: simPer, laiPeriods, landCoverPeriods
    use mo_common_constants, only: maxNLais, maxNLcovers, keepUnneededPeriodsLAI, keepUnneededPeriodsLandCover

    !> grid to save information to
    type(Grid), intent(inout) :: newGrid
    !> domain index (if 1, set soil information as reference for all other domains)
    integer(i4), intent(in) :: iDomain
    !> ids for land cover periods in nc file to select
    integer(i4), dimension(:), allocatable, intent(out) :: landCoverSelect
    !> total number of land cover periods in nc file
    integer(i4), intent(out) :: nLandCoverPeriodsTemp
    !> ids for lai periods in nc file to select
    integer(i4), dimension(:), allocatable, intent(out) :: laiSelect
    !> total number of lai periods in nc file
    integer(i4), intent(out) :: nLaisTemp
    !> optional flag whether to initialize the latlon grid (can be shared among domains)
    logical, intent(in), optional :: doInitArg
  
    logical :: doInit
    integer(i4) :: nSoilHorizonsTemp
    real(dp), dimension(:), allocatable :: soilHorizonBoundariesTemp
    integer(i4) :: iDA
    integer(i4), dimension(:), allocatable :: iDimX, iDimY, iDim
    character(maxNameLength) :: laiDimName
    character(2048) :: units

    if (present(doInitArg)) then
      doInit = doInitArg
    else
      doInit = .false.
    end if

    if (doInitArg) then
      ! TODO: make all that more flexible (coordinate and mask detection)
      ! get the x dimension
      iDimX = get_index_in_coordinate('lon_out')
      ! get the y dimension
      iDimY = get_index_in_coordinate('lat_out')

      ! init some main properties
      newGrid%xllcorner = MPR_COORDINATES(iDimX(1))%bounds(1)
      newGrid%yllcorner = MPR_COORDINATES(iDimY(1))%bounds(1)
      newGrid%nrows = int(MPR_COORDINATES(iDimX(1))%count, kind=i4)
      newGrid%ncols = int(MPR_COORDINATES(iDimY(1))%count, kind=i4)
      newGrid%cellsize = MPR_COORDINATES(iDimX(1))%step

      ! allocate all 2d properties, now that we know the dimensionality
      allocate(newGrid%x(newGrid%nrows, newGrid%ncols))
      allocate(newGrid%y(newGrid%nrows, newGrid%ncols))
      allocate(newGrid%mask(newGrid%nrows, newGrid%ncols))

      ! for historic reasons this is in 2d
      newGrid%x = spread(MPR_COORDINATES(iDimX(1))%values(1:), 2, newGrid%ncols)
      newGrid%y = spread(MPR_COORDINATES(iDimY(1))%values(1:), 1, newGrid%nrows)

      ! use a 2d field for inferring the mask (cannot be derived from dimensions)
      iDA = get_index_in_vector('L1_latitude', MPR_DATA_ARRAYS)
      newGrid%mask = reshape(MPR_DATA_ARRAYS(iDA)%reshapedMask, [newGrid%nrows, newGrid%ncols])
      call init_advanced_grid_properties(newGrid)

      ! get the z dimension
      iDim = get_index_in_coordinate('horizon_out')
      nSoilHorizonsTemp = int(MPR_COORDINATES(iDim(1))%count, kind=i4)
      allocate(soilHorizonBoundariesTemp(0: nSoilHorizonsTemp))
      soilHorizonBoundariesTemp = MPR_COORDINATES(iDim(1))%values(0: nSoilHorizonsTemp)

      ! check if soil and LAI are consistent with other domains, set to global if domain == 1
      call check_soil_dimension_consistency(iDomain, nSoilHorizonsTemp, soilHorizonBoundariesTemp)

    end if

    ! use a field for inferring that dimension (name depends on the input data and is flexible)
    iDA = get_index_in_vector('L1_Max_Canopy_Intercept', MPR_DATA_ARRAYS)
    laiDimName = MPR_DATA_ARRAYS(iDA)%coords(3)%coord_p%name

    ! get the LAI dimension
    iDim = get_index_in_coordinate(trim(laiDimName))
    nLaisTemp = int(MPR_COORDINATES(iDim(1))%count, kind=i4)
    call MPR_COORDINATES(iDim(1))%get_attribute('units', units)

    ! converting real values to integer
    call laiPeriods(iDomain)%init(n=nLaisTemp, nMax=maxNLais, name='LAI', simPerArg=simPer(iDomain), &
            units=units, &
            periodValues=nint(MPR_COORDINATES(iDim(1))%values), &
            keepUnneededPeriods=keepUnneededPeriodsLAI, selectIndices=laiSelect)

    ! get the landcover dimension
    iDim = get_index_in_coordinate('land_cover_period_out')
    nLandCoverPeriodsTemp = int(MPR_COORDINATES(iDim(1))%count, kind=i4)
    call MPR_COORDINATES(iDim(1))%get_attribute('units', units)

    ! converting real values to integer
    call landCoverPeriods(iDomain)%init(n=nLandCoverPeriodsTemp, nMax=maxNLcovers, name='land cover', &
            simPerArg=simPer(iDomain), &
            units=units, &
            periodValues=nint(MPR_COORDINATES(iDim(1))%values), &
            keepUnneededPeriods=keepUnneededPeriodsLandCover, selectIndices=landCoverSelect)

  end subroutine init_grid

  subroutine check_consistency()
    !< routine collecting all checks on the parameter fields read from MPR and its coordinates
    use mo_global_variables, only : nSoilHorizons
    use mo_common_variables, only : opti_function, optimize
    use mo_global_variables, only : nSoilHorizons_sm_input
    use mo_string_utils, only: num2str
    use mo_message, only: error_message

    if (optimize) then
      select case (opti_function)
      case(10 : 13, 28)
        ! soil moisture
        if (nSoilHorizons_sm_input > nSoilHorizons) then
          call error_message('***ERROR: Number of soil horizons representative for input soil moisture exceeded', &
                   new_line('a'), '          defined number of soil horizions in mHM: ', &
                  adjustl(trim(num2str(nSoilHorizons))), '!')
        end if
      end select
    end if

  end subroutine check_consistency

end module mo_mhm_mpr_interface
