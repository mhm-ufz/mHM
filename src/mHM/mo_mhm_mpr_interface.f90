!>       \file mo_mhm_mpr_interface.f90
!>       \brief calls the MPR library for parameter estimation of mHM parameters
!>       \details - MPR is called with global parameters and path to its namelist
!>                - MPR reads all the land surface properties on the level 0 grid(s)
!>                - it then estimates the parameters and scaled them to level 1 (mHM grid (hydrology_resolution))
!>                - parameter names have to conform with mHM parameter names
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

  use mo_mpr_data_array, only : MPR_DATA_ARRAYS
  use mo_mpr_constants, only : maxNameLength
  use mo_kind, only : i4, dp
  implicit none

  private

  public :: call_mpr
  ! ------------------------------------------------------------------

contains

  subroutine call_mpr(parameterValues, parameterNames, grids, do_init_arg, opti_domain_indices)
    use mo_mpr_global_variables, only : out_filename
    use mo_mpr_read_config, only : mpr_read_config
    use mo_mpr_constants, only : maxNoDataArrays
    use mo_mpr_reset, only: reset
    ! use mo_mpr_reorder_data_array, only: reorder_data_arrays
    use mo_netcdf, only : NcDataset
    use mo_common_variables, only : write_restart, dummy_global_parameters, dummy_global_parameters_name, domainMeta, &
            nLandCoverPeriods
    use mo_append, only: append
    use mo_file, only : unamelist_mpr
    use mo_global_variables, only: pathMprNml, are_parameter_initialized
    use mo_grid, only: Grid
    use mo_read_nc, only: get_land_cover_period_indices

    real(dp), dimension(:), intent(in) :: parameterValues
    character(*), dimension(:), intent(in) :: parameterNames
    type(Grid), dimension(:), intent(inout) :: grids
    logical, intent(in), optional :: do_init_arg
    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices

    real(dp), dimension(:), allocatable :: parameterValuesConcat
    character(maxNameLength), dimension(:), allocatable :: parameterNamesConcat
    real(dp), dimension(:), allocatable :: landCoverPeriodBoundaries_temp
    integer(i4), dimension(:), allocatable :: landCoverSelect
    integer(i4) :: iDomain, previousDomain, iDA, nSlices
    type(NcDataset) :: nc

    parameterValuesConcat = parameterValues
    parameterNamesConcat = parameterNames
    call append(parameterValuesConcat, dummy_global_parameters)
    call append(parameterNamesConcat, dummy_global_parameters_name)
    
    previousDomain = 0
    do iDomain = 1, domainMeta%nDomains
      
      if (allocated(landCoverSelect)) deallocate(landCoverSelect)
      ! read the config
      if (domainMeta%L0DataFrom(iDomain) == previousDomain) then
        ! use grid from previousDomain
        grids(iDomain) = grids(previousDomain)
        ! get new landCoverperiod indices
        call get_land_cover_period_indices(iDomain, landCoverPeriodBoundaries_temp, landCoverSelect)
        ! use parameters from previousDomain but use iDomains landCoverSeelct (based also on simulation period)
        call init_eff_params_from_data_arrays(iDomain, grids(iDomain), nSlices, landCoverSelect, do_init_arg)
        cycle
      end if
      
      previousDomain = domainMeta%L0DataFrom(iDomain)
      
      ! delete MPR_DATA_ARRAYS
      call reset()

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
        ! stop if no name is initialized
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
    
      ! TODO: we need to make sure there is no data gaps in data arrays... (former L0_check_input)
      call init_grid(grids(iDomain), iDomain, landCoverPeriodBoundaries_temp, nSlices, do_init_arg)
      if (iDomain == 1) then
        nLandCoverPeriods = nSlices
      end if
      ! set new landCoverperiod
      call get_land_cover_period_indices(iDomain, landCoverPeriodBoundaries_temp, landCoverSelect)
      ! set parameters
      call init_eff_params_from_data_arrays(iDomain, grids(iDomain), nSlices, landCoverSelect, do_init_arg)
    end do
    ! delete MPR_DATA_ARRAYS
    call reset()
    
    call check_consistency()
    are_parameter_initialized = .true.

  end subroutine call_mpr

  subroutine init_eff_params_from_data_arrays(iDomain, grid, nSlices, landCoverSelect, do_init_arg)
    use mo_global_variables, only: nSoilHorizons, nLAIs, L1_fSealed, L1_alpha, L1_degDayInc, L1_degDayMax, &
            L1_degDayNoPre, L1_karstLoss, L1_maxInter, L1_kFastFlow, L1_kSlowFlow, &
            L1_kBaseFlow, L1_kPerco, L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp, L1_jarvis_thresh_c1, &
            L1_petLAIcorFactor, L1_tempThresh, L1_unsatThresh, L1_sealedThresh, L1_wiltingPoint, L1_fAsp, &
            L1_latitude, L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_fRoots, L1_surfResist
    use mo_common_variables, only: domainMeta, processMatrix
    use mo_grid, only: grid_type => Grid
    use mo_common_variables, only: nLandCoverPeriods
    use mo_mpr_constants, only: maxNoDataArrays
    use mo_append, only: append, add_nodata_slice
    use mo_constants, only: nodata_dp
    use mo_mrm_global_variables, only: L0_slope

    integer(i4), intent(in) :: iDomain
    type(grid_type), intent(inout) :: grid
    integer(i4), intent(in) :: nSlices
    integer(i4), dimension(:), allocatable, intent(in) :: landCoverSelect
    logical, intent(in), optional :: do_init_arg
    
    integer(i4) :: nCells, iDA, s1, e1, newSlices
    logical :: do_init
    real(dp), dimension(:, :, :), allocatable :: dummy3D
    real(dp), dimension(:, :), allocatable :: dummy2D
  
    do_init = .true.
    if (present(do_init_arg)) do_init = do_init_arg
  
    nCells = grid%nCells
    newSlices = size(landCoverSelect)
    ! init the parameters
    do iDA = 1, maxNoDataArrays
      ! stop if no name is initialized
      if (MPR_DATA_ARRAYS(iDA)%name(1:3) == 'L1_') then
        if (.not. do_init) then
          s1 = grid%iStart
          e1 = grid%iEnd
          ! The data arrays are selected based on name
          ! they only exist in a packed 1D-version and thus first need to be reshaped
          ! we know from landCoverSelect that we only need a certain slices of the landCoverPeriods
          ! we select the slices after the reshape by using newSlices
          ! then, we need to have nLandCoverPeriods in the arrays last dimensions, newSlices may be different
          ! to fulfill that, we add dummy slices in the last dimension
          select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_SealedFraction')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_fSealed(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_Alpha')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_alpha(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayInc')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_degDayInc(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayMax')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_degDayMax(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_DegDayNoPre')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_degDayNoPre(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_KarstLoss')
              ! rows, cols
              L1_karstLoss(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_Max_Canopy_Intercept')
              ! rows, cols, lais
              L1_maxInter(s1:e1, :) = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs])
            case('L1_FastFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_kFastFlow(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_SlowFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_kSlowFlow(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_kBaseFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_kBaseFlow(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_Kperco')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_kPerco(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_FieldCap')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_soilMoistFC(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_SatSoilMoisture')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_soilMoistSat(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_SoilMoistureExponent')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_soilMoistExp(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_PermWiltPoint')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_wiltingPoint(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_TempThresh')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_tempThresh(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_UnsatThreshold')
              ! rows, cols
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              L1_unsatThresh(s1:e1, 1:newSlices) = dummy2D(:, landCoverSelect)
            case('L1_SealedThresh')
              ! rows, cols
              L1_sealedThresh(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_Jarvis_Threshold')
              ! rows, cols
              L1_jarvis_thresh_c1(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_PET_LAI_correction_factor')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs, nSlices])
              L1_petLAIcorFactor(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_HarSamCoeff')
              ! rows, cols
              L1_HarSamCoeff(s1:e1) = MPR_DATA_ARRAYS(iDA)%data
            case('L1_PrieTayAlpha')
              ! rows, cols, lais
              L1_PrieTayAlpha(s1:e1, :) = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs])
            case('L1_Aerodyn_resist')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs, nSlices])
              L1_aeroResist(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            case('L1_Bulk_Surface_Resist')
              ! rows, cols, lais
              L1_surfResist(s1:e1, :) = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs])
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
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_fRoots(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            end select
          else
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              L1_fRoots(s1:e1, :, 1:newSlices) = dummy3D(:, :, landCoverSelect)
            end select
          end if
        else
          select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_SealedFraction')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_fSealed, dummy2D)
            case('L1_Alpha')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_alpha, dummy2D)
            case('L1_DegDayInc')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_degDayInc, dummy2D)
            case('L1_DegDayMax')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_degDayMax, dummy2D)
            case('L1_DegDayNoPre')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_degDayNoPre, dummy2D)
            case('L1_KarstLoss')
              ! rows, cols
              call append(L1_karstLoss, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_Max_Canopy_Intercept')
              ! rows, cols, lais
              call append(L1_maxInter, reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs]))
            case('L1_FastFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_kFastFlow, dummy2D)
            case('L1_SlowFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_kSlowFlow, dummy2D)
            case('L1_kBaseFlow')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_kBaseFlow, dummy2D)
            case('L1_Kperco')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_kPerco, dummy2D)
            case('L1_FieldCap')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_soilMoistFC, dummy3D, idim=1)
            case('L1_SatSoilMoisture')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_soilMoistSat, dummy3D, idim=1)
            case('L1_SoilMoistureExponent')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_soilMoistExp, dummy3D, idim=1)
            case('L1_PermWiltPoint')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_wiltingPoint, dummy3D, idim=1)
            case('L1_TempThresh')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_tempThresh, dummy2D)
            case('L1_UnsatThreshold')
              ! rows, cols, lcscenes
              dummy2D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSlices])
              dummy2D = dummy2D(:, landCoverSelect)
              call add_nodata_slice(dummy2D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_unsatThresh, dummy2D)
            case('L1_SealedThresh')
              ! rows, cols
              call append(L1_sealedThresh, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_Jarvis_Threshold')
              ! rows, cols
              call append(L1_jarvis_thresh_c1, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_PET_LAI_correction_factor')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_petLAIcorFactor, dummy3D)
            case('L1_HarSamCoeff')
              ! rows, cols
              call append(L1_HarSamCoeff, MPR_DATA_ARRAYS(iDA)%data)
            case('L1_PrieTayAlpha')
              ! rows, cols, lais
              call append(L1_PrieTayAlpha, reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs]))
            case('L1_Aerodyn_resist')
              ! rows, cols, lais, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_aeroResist, dummy3D)
            case('L1_Bulk_Surface_Resist')
              ! rows, cols, lais
              call append(L1_surfResist, reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nLAIs]))
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
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_fRoots, dummy3D, idim=1)
            end select
          else
            select case(trim(MPR_DATA_ARRAYS(iDA)%name))
            case('L1_fRoots')
              ! rows, cols, soil, lcscenes
              dummy3D = reshape(MPR_DATA_ARRAYS(iDA)%data, [nCells, nSoilHorizons, nSlices])
              dummy3D = dummy3D(:, :, landCoverSelect)
              call add_nodata_slice(dummy3D, nLandCoverPeriods - newSlices, nodata_dp)
              call append(L1_fRoots, dummy3D, idim=1)
            end select
          end if
        end if
      else if (trim(MPR_DATA_ARRAYS(iDA)%name) == 'slope') then
        call append(L0_slope, MPR_DATA_ARRAYS(iDA)%data)
      end if
    end do
  
  
  end subroutine init_eff_params_from_data_arrays

  subroutine init_grid(new_grid, iDomain, landCoverPeriodBoundaries_temp, nLandCoverPeriods_temp, do_init_arg)
    use mo_common_variables, only : Grid, domainMeta
    use mo_grid, only : init_advanced_grid_properties
    use mo_mpr_coordinate, only: get_index_in_coordinate, MPR_COORDINATES
    use mo_mpr_utils, only: get_index_in_vector
    use mo_mpr_constants, only : maxNoDataArrays, maxNameLength
    use mo_read_nc, only: check_dimension_consistency
  
    ! grid to save information to
    type(Grid), intent(inout) :: new_grid
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), allocatable, intent(out) :: landCoverPeriodBoundaries_temp
    integer(i4), intent(out) :: nLandCoverPeriods_temp
    logical, intent(in), optional :: do_init_arg
  
    logical :: do_init
    integer(i4) :: nSoilHorizons_temp, nLAIs_temp
    real(dp), dimension(:), allocatable :: soilHorizonBoundaries_temp, &
            LAIBoundaries_temp
    integer(i4) :: iDA, uniqueIDomain
    integer(i4), dimension(:), allocatable :: iDim_x, iDim_y, iDim
    character(maxNameLength) :: LAI_dim_name
  
    ! TODO: make all that more flexible (coordinate and mask detection)
    ! get the x dimension
    iDim_x = get_index_in_coordinate('lon_out')
    ! get the y dimension
    iDim_y = get_index_in_coordinate('lat_out')

    ! init some main properties
    new_grid%xllcorner = MPR_COORDINATES(iDim_x(1))%bounds(1)
    new_grid%yllcorner = MPR_COORDINATES(iDim_y(1))%bounds(1)
    new_grid%nrows = MPR_COORDINATES(iDim_x(1))%count
    new_grid%ncols = MPR_COORDINATES(iDim_y(1))%count
    new_grid%cellsize = MPR_COORDINATES(iDim_x(1))%step

    ! allocate all 2d properties, now that we know the dimensionality
    allocate(new_grid%x(new_grid%nrows, new_grid%ncols))
    allocate(new_grid%y(new_grid%nrows, new_grid%ncols))
    allocate(new_grid%mask(new_grid%nrows, new_grid%ncols))

    ! for historic reasons this is in 2d
    new_grid%x = spread(MPR_COORDINATES(iDim_x(1))%values(1:), 2, new_grid%ncols)
    new_grid%y = spread(MPR_COORDINATES(iDim_y(1))%values(1:), 1, new_grid%nrows)

    ! use a 2d field for inferring the mask (cannot be derived from dimensions)
    iDA = get_index_in_vector('L1_latitude', MPR_DATA_ARRAYS)
    new_grid%mask = reshape(MPR_DATA_ARRAYS(iDA)%reshapedMask, [new_grid%nrows, new_grid%ncols])
    call init_advanced_grid_properties(new_grid)

    ! get the z dimension
    iDim = get_index_in_coordinate('horizon_out')
    nSoilHorizons_temp = MPR_COORDINATES(iDim(1))%count
    allocate(soilHorizonBoundaries_temp(0: nSoilHorizons_temp))
    soilHorizonBoundaries_temp = MPR_COORDINATES(iDim(1))%values(0: nSoilHorizons_temp)
  
    ! get the landcover dimension
    iDim = get_index_in_coordinate('land_cover_period_out')
    nLandCoverPeriods_temp = MPR_COORDINATES(iDim(1))%count
    allocate(landCoverPeriodBoundaries_temp(0: nLandCoverPeriods_temp))
    landCoverPeriodBoundaries_temp = MPR_COORDINATES(iDim(1))%values(0: nLandCoverPeriods_temp)

    ! use a field for inferring that dimension (name depends on the input data and is flexible)
    iDA = get_index_in_vector('L1_Max_Canopy_Intercept', MPR_DATA_ARRAYS)
    LAI_dim_name = MPR_DATA_ARRAYS(iDA)%coords(3)%coord_p%name
  
    ! get the LAI dimension
    iDim = get_index_in_coordinate(trim(LAI_dim_name))
    nLAIs_temp = MPR_COORDINATES(iDim(1))%count
    allocate(LAIBoundaries_temp(0: nLAIs_temp))
    LAIBoundaries_temp = MPR_COORDINATES(iDim(1))%values(0: nLAIs_temp)

    ! check if soil and LAI are consistent with other domains, set to global if domain == 1
    call check_dimension_consistency(iDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp)
  
  
  end subroutine init_grid

  subroutine check_consistency()
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
