!>       \file mo_mhm_mpr_interface.f90
!>       \brief TODO: add description
!>       \details TODO: add description

module mo_mhm_mpr_interface

  use mo_mpr_data_array, only : MPR_DATA_ARRAYS
  use mo_mpr_constants, only : maxNameLength
  use mo_kind, only : i4, dp
  implicit none

  private

  public :: call_mpr
  ! public :: init_eff_params_from_data_arrays


  ! ------------------------------------------------------------------

contains

    subroutine call_mpr(parameterValues, parameterNames, grids, do_init_arg, opti_domain_indices)
    ! use mo_mpr_global_variables, only : soilDB, L0_asp, L0_geoUnit, L0_gridded_LAI, L0_slope, L0_slope_emp, L0_soilId
    use mo_mpr_global_variables, only : out_filename
    use mo_mpr_read_config, only : mpr_read_config
    use mo_mpr_constants, only : maxNoDataArrays
    use mo_mpr_reset, only: reset
    use mo_mpr_reorder_data_array, only: reorder_data_arrays
    use mo_netcdf, only : NcDataset
    use mo_common_variables, only : nuniqueL0Domains, write_restart, &
            dummy_global_parameters, dummy_global_parameters_name, domainMeta, level0, l0_l1_remap, L0_LCover, L0_elev
    use mo_append, only: append
    use mo_file, only : unamelist_mpr
    use mo_global_variables, only: pathMprNml, are_parameter_initialized, LAIBoundaries
    use mo_grid, only: Grid
    ! use mo_mpr_eval, only: mpr_eval
    ! use mo_mpr_startup, only: mpr_initialize
    ! use mo_read_wrapper, only: read_data
    use mo_common_datetime_type, only: simPer
    USE mo_timer, ONLY : timer_start, timer_stop
    use mo_restart, only: reset_eff_params

    real(dp), dimension(:), intent(in) :: parameterValues
    character(*), dimension(:), intent(in) :: parameterNames
    type(Grid), dimension(:), intent(inout) :: grids
    logical, intent(in), optional :: do_init_arg
    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices

    logical :: do_init

    real(dp), dimension(:), allocatable :: parameterValuesConcat
    character(maxNameLength), DIMENSION(:), allocatable :: parameterNamesConcat
    integer(i4) :: i, iUniqueDomain
    type(NcDataset) :: nc

    ! deallocate everything in case optimization is active
    ! if (allocated(soilDB%id)) deallocate(soilDB%id)
    ! if (allocated(soilDB%nHorizons)) deallocate(soilDB%nHorizons)
    ! if (allocated(soilDB%is_present)) deallocate(soilDB%is_present)
    ! if (allocated(soilDB%UD)) deallocate(soilDB%UD)
    ! if (allocated(soilDB%LD)) deallocate(soilDB%LD)
    ! if (allocated(soilDB%clay)) deallocate(soilDB%clay)
    ! if (allocated(soilDB%sand)) deallocate(soilDB%sand)
    ! if (allocated(soilDB%DbM)) deallocate(soilDB%DbM)
    ! if (allocated(soilDB%depth)) deallocate(soilDB%depth)
    ! if (allocated(soilDB%RZdepth)) deallocate(soilDB%RZdepth)
    ! if (allocated(soilDB%Wd)) deallocate(soilDB%Wd)
    ! if (allocated(soilDB%nTillHorizons)) deallocate(soilDB%nTillHorizons)
    ! if (allocated(soilDB%thetaS_Till)) deallocate(soilDB%thetaS_Till)
    ! if (allocated(soilDB%thetaS)) deallocate(soilDB%thetaS)
    ! if (allocated(soilDB%Db)) deallocate(soilDB%Db)
    ! if (allocated(soilDB%thetaFC_Till)) deallocate(soilDB%thetaFC_Till)
    ! if (allocated(soilDB%thetaFC)) deallocate(soilDB%thetaFC)
    ! if (allocated(soilDB%thetaPW_Till)) deallocate(soilDB%thetaPW_Till)
    ! if (allocated(soilDB%thetaPW)) deallocate(soilDB%thetaPW)
    ! if (allocated(soilDB%Ks)) deallocate(soilDB%Ks)
    ! if (allocated(L0_asp)) deallocate(L0_asp)
    if (allocated(L0_LCover)) deallocate(L0_LCover)
    if (allocated(L0_elev)) deallocate(L0_elev)
    ! if (allocated(L0_geoUnit)) deallocate(L0_geoUnit)
    ! if (allocated(L0_gridded_LAI)) deallocate(L0_gridded_LAI)
    ! if (allocated(L0_slope)) deallocate(L0_slope)
    ! if (allocated(L0_slope_emp)) deallocate(L0_slope_emp)
    ! if (allocated(L0_soilId)) deallocate(L0_soilId)
    if (allocated(level0)) deallocate(level0)
    if (allocated(LAIBoundaries)) deallocate(LAIBoundaries)
    if (allocated(l0_l1_remap)) deallocate(l0_l1_remap)

    ! do_init = .true.
    ! if (present(do_init_arg)) do_init = do_init_arg

    ! if (do_init) then
    !   call read_data(simPer)
    !   call mpr_initialize()
    ! else
    !   call reset_eff_params()
    ! end if

    ! call mpr_eval(parameterValues, opti_domain_indices)

    parameterValuesConcat = parameterValues
    parameterNamesConcat = parameterNames
    call append(parameterValuesConcat, dummy_global_parameters)
    call append(parameterNamesConcat, dummy_global_parameters_name)
    
    do iUniqueDomain=1, nuniqueL0Domains
      ! read the config
      call mpr_read_config(pathMprNml(iUniqueDomain), unamelist_mpr, &
              parameterValues=parameterValuesConcat, parameterNames=parameterNamesConcat)
    
      ! reorder data array to minimize storage requirement and assure that all required input fields have been read before
      ! call reorder_data_arrays(MPR_DATA_ARRAYS)
    
      if (write_restart) then
        ! open the output file
        nc = NcDataset(out_filename, "w")
      end if
    
      do i = 1, maxNoDataArrays
        ! stop if no name is initialized
        if (.not. MPR_DATA_ARRAYS(i)%is_initialized) then
          cycle
        end if
    
        ! execute MPR on each DataArray
        call MPR_DATA_ARRAYS(i)%execute()
        if (write_restart .and. MPR_DATA_ARRAYS(i)%toFile) then
          ! write the current parameter to the output file
          call MPR_DATA_ARRAYS(i)%write(nc)
        end if
    
      end do
    
      if (write_restart) then
        ! close the output file
        call nc%close()
      end if
    
      ! TODO: we need to make sure there is no data gaps in data arrays... (former L0_check_input)
      call init_eff_params_from_data_arrays(iUniqueDomain, grids(iUniqueDomain), do_init_arg)
      ! delete MPR_DATA_ARRAYS
      call reset()
    end do
    call check_consistency()
    are_parameter_initialized = .true.

  end subroutine call_mpr

  subroutine init_eff_params_from_data_arrays(iUniqueDomain, grid, do_init_arg)
    use mo_global_variables, only: nSoilHorizons, nLAIs, L1_fSealed, L1_alpha, L1_degDayInc, L1_degDayMax, &
            L1_degDayNoPre, L1_karstLoss, L1_maxInter, L1_kFastFlow, L1_kSlowFlow, &
            L1_kBaseFlow, L1_kPerco, L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp, L1_jarvis_thresh_c1, &
            L1_petLAIcorFactor, L1_tempThresh, L1_unsatThresh, L1_sealedThresh, L1_wiltingPoint, L1_fAsp, &
            L1_latitude, L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_fRoots, L1_surfResist
    use mo_common_variables, only: domainMeta, L0_Domain, processMatrix
    use mo_grid, only: grid_type => Grid
    use mo_common_variables, only: nLandCoverPeriods
    use mo_mpr_constants, only: maxNoDataArrays
    use mo_append, only: append, add_nodata_slice
    use mo_constants, only: nodata_dp
  
    integer(i4), intent(in) :: iUniqueDomain
    type(grid_type), intent(inout) :: grid
    logical, intent(in), optional :: do_init_arg
    integer(i4) :: iDomain, nCells, iDA, s1, e1, nSlices, newSlices
    logical :: share_L0, do_init
    integer(i4), dimension(:), allocatable :: landCoverSelect
    real(dp), dimension(:, :, :), allocatable :: dummy3D
    real(dp), dimension(:, :), allocatable :: dummy2D
  
    do_init = .true.
    if (present(do_init_arg)) do_init = do_init_arg
  
    if (allocated(landCoverSelect)) deallocate(landCoverSelect)
    ! init the grid
    call init_grid(grid, iUniqueDomain, landCoverSelect, nSlices, do_init)
  
    nCells = grid%nCells
    newSlices = size(landCoverSelect)
    ! init the parameters
    do iDA = 1, maxNoDataArrays
      print*, MPR_DATA_ARRAYS(iDA)%name !TODO: remove
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
            ! rows, cols, lais, lcscenes
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
            print*, "test", landCoverSelect !TODO: remove
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
      end if
    end do
  
  
  end subroutine init_eff_params_from_data_arrays

  subroutine init_grid(new_grid, iDomain, landCoverSelect, nLandCoverPeriods_temp, do_init_arg)
    use mo_common_variables, only : Grid, L0_Domain, domainMeta
    use mo_grid, only : init_advanced_grid_properties
    use mo_mpr_coordinate, only: get_index_in_coordinate, MPR_COORDINATES
    use mo_mpr_utils, only: get_index_in_vector
    use mo_mpr_constants, only : maxNoDataArrays, maxNameLength
    use mo_read_nc, only: check_dimension_consistency
  
    ! grid to save information to
    type(Grid), intent(inout) :: new_grid
    integer(i4), intent(in) :: iDomain
    integer(i4), dimension(:), intent(out), allocatable :: landCoverSelect
    integer(i4), intent(out) :: nLandCoverPeriods_temp
    logical, intent(in), optional :: do_init_arg
  
    logical :: do_init
    integer(i4) :: nSoilHorizons_temp, nLAIs_temp
    real(dp), dimension(:), allocatable :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
            LAIBoundaries_temp
    integer(i4) :: iDA, uniqueIDomain
    integer(i4), dimension(:), allocatable :: iDim_x, iDim_y, iDim
    character(maxNameLength) :: LAI_dim_name
  
    do_init = .true.
    if (present(do_init_arg)) do_init = do_init_arg
  
    if (do_init) then
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
  
    end if
  
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

    uniqueIDomain = domainMeta%L0DataFrom(iDomain)
    call check_dimension_consistency(iDomain, uniqueIDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp, &
            landCoverSelect, do_init)
  
  
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
