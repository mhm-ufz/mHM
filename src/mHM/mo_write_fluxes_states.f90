!> \file mo_write_fluxes_states.f90
!> \brief \copybrief mo_write_fluxes_states
!> \details \copydetails mo_write_fluxes_states

!> \brief Creates NetCDF output for different fluxes and state variables of mHM.
!> \details NetCDF is first initialized and later on variables are put to the NetCDF.
!> \changelog
!! - David Schaefer       Aug 2015
!!   - major rewrite
!! - O. Rakovec, R. Kumar Nov 2017
!!   - added project description for the netcdf outputs
!> \authors Matthias Zink
!> \date Apr 2013
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
module mo_write_fluxes_states

  use mo_nc_output, only: OutputDataset, OutputVariable, set_attributes, data_dims, data_dtype
  use mo_kind, only : i4, dp
  use mo_string_utils, only : num2str
  use mo_common_constants, only : nodata_dp
  use mo_netcdf, only : NcDataset
  use mo_global_variables, only : &
    output_deflate_level, output_double_precision, output_time_reference, timeStep_model_outputs, outputFlxState
  use mo_common_mHM_mRM_variables, only: timeStep
  use mo_mpr_global_variables, only : nSoilHorizons_mHM
  use mo_common_variables, only : iFlag_cordinate_sys, level1
  use mo_file, only : file_mhm_output

  implicit none

contains

  !> \brief Initialize mHM OutputDataset
  !> \details Create and initialize the output file. If new a new output
  !!       variable needs to be written, this is the first of two
  !!       procedures to change (second: updateDataset)
  !!
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Pallav Shrestha Mar 2020
  !!   - iFlag_cordinate_sys based dimensions (dims)
  !!
  !> \return type(OutputDataset)
  !> \authors Matthias Zink
  !> \date Apr 2013
  function mHM_OutputDataset(iDomain, mask1) result(out)
    implicit none

    !> domain id
    integer(i4), intent(in) :: iDomain
    !> L1 mask to reconstruct the data
    logical, pointer, intent(in), dimension(:, :) :: mask1

    type(OutputDataset) :: out

    integer(i4) :: ii, nn, nCells
    character(3) :: dtype
    character(16), dimension(3) :: dims
    character(16) :: unit
    type(OutputVariable), dimension(size(outputFlxState) * nSoilHorizons_mHM) :: tmpvars

    nCells = level1(iDomain)%nCells

    out = OutputDataset( &
      iDomain=iDomain, &
      level=level1, &
      file_name=file_mhm_output, &
      double_precision=output_double_precision, &
      outputs_frequence=timeStep_model_outputs, &
      time_reference=output_time_reference &
    )

    dtype = data_dtype(output_double_precision)
    unit = fluxesUnit(iDomain)
    dims = data_dims()

    ii = 0

    if (outputFlxState(1)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "interception", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "canopy interception storage", "mm", output_double_precision)
    end if

    if (outputFlxState(2)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "snowpack", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(tmpvars(ii)%nc, "depth of snowpack", "mm", output_double_precision)
    end if

    if (outputFlxState(3)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(out%nc, "SWC_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims, nCells, mask1, output_deflate_level, .true.)
        call set_attributes(tmpvars(ii)%nc, &
                'soil water content of soil layer' // trim(num2str(nn)), "mm", output_double_precision)
      end do
    end if

    if (outputFlxState(4)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(out%nc, "SM_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims, nCells, mask1, output_deflate_level, .true.)
        call set_attributes(tmpvars(ii)%nc, &
                'volumetric soil moisture of soil layer' // trim(num2str(nn)), "mm mm-1", output_double_precision)
      end do
    end if

    if (outputFlxState(5)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "SM_Lall", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "average soil moisture over all layers", "mm mm-1", output_double_precision)
    end if

    if (outputFlxState(6)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "sealedSTW", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "reservoir of sealed areas (sealedSTW)", "mm", output_double_precision)
    end if

    if (outputFlxState(7)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "unsatSTW", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "reservoir of unsaturated zone", "mm", output_double_precision)
    end if

    if (outputFlxState(8)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "satSTW", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "water level in groundwater reservoir", "mm", output_double_precision)
    end if

    if (outputFlxState(18)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "neutrons", dtype, dims, nCells, mask1, output_deflate_level, .true.)
      call set_attributes(&
              tmpvars(ii)%nc, "ground albedo neutrons", "cph", output_double_precision)
    end if

    if (outputFlxState(9)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "PET", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "potential Evapotranspiration", trim(unit), output_double_precision)
    end if

    if (outputFlxState(10)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "aET", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "actual Evapotranspiration", trim(unit), output_double_precision)
    end if

    if (outputFlxState(11)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "Q", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "total runoff generated by every cell", trim(unit), output_double_precision)
    end if

    if (outputFlxState(12)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "QD", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(tmpvars(ii)%nc, &
              "direct runoff generated by every cell (runoffSeal)", trim(unit), output_double_precision)
    end if

    if (outputFlxState(13)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "QIf", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(tmpvars(ii)%nc, &
              "fast interflow generated by every cell (fastRunoff)", trim(unit), output_double_precision)
    end if

    if (outputFlxState(14)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "QIs", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(tmpvars(ii)%nc, &
              "slow interflow generated by every cell (slowRunoff)", trim(unit), output_double_precision)
    end if

    if (outputFlxState(15)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "QB", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "baseflow generated by every cell", trim(unit), output_double_precision)
    end if

    if (outputFlxState(16)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "recharge", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "groundwater recharge", trim(unit), output_double_precision)
    end if

    if (outputFlxState(17)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(&
                out%nc, "soil_infil_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims, nCells, mask1, output_deflate_level)
        call set_attributes(tmpvars(ii)%nc, &
                "infiltration flux from soil layer" // trim(num2str(nn)), unit, output_double_precision)
      end do
    end if

    if (outputFlxState(19)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(&
                out%nc, "aET_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims, nCells, mask1, output_deflate_level)
        call set_attributes(tmpvars(ii)%nc, &
                'actual Evapotranspiration from soil layer' // trim(num2str(nn)), &
                "mm " // trim(unit), output_double_precision)
      end do
    end if

    if (outputFlxState(20)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "preEffect", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "effective precipitation", trim(unit), output_double_precision)
    end if

    if (outputFlxState(21)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              out%nc, "Qsm", dtype, dims, nCells, mask1, output_deflate_level)
      call set_attributes(&
              tmpvars(ii)%nc, "Average liquid water generated from solid to liquid phase change in the snow", &
              trim(unit), output_double_precision)
    end if

    allocate(out%vars(ii))
    out%vars = tmpvars(1 : ii)

  end function mHM_OutputDataset

  !> \brief Update all variables.
  !> \details Call the type bound procedure updateVariable for
  !> all output variables. If a new output
  !> variable needs to be written, this is the second
  !> of two procedures to change (first: newOutputDataset)
  !!
  !> \changelog
  !! - L. Samaniego et al. Dec  2013
  !!   - nullify pointer Matthias Zink,        Feb. 2014
  !!   - added aditional output: pet V. Prykhodk, J. Mai,  Nov. 2014
  !!   - adding new variable infilSoil
  !!   - case 16 David Schaefer      , Jun. 2015
  !!   - major rewrite
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Matthias Zink
  !> \date Apr 2013
  subroutine mHM_updateDataset(nc_mhm, sidx, eidx, L1_fSealed, L1_fNotSealed, L1_inter, L1_snowPack, L1_soilMoist, &
                          L1_soilMoistSat, L1_sealSTW, L1_unsatSTW, L1_satSTW, L1_neutrons, L1_pet, L1_aETSoil, &
                          L1_aETCanopy, L1_aETSealed, L1_total_runoff, L1_runoffSeal, L1_fastRunoff, L1_slowRunoff, &
                          L1_baseflow, L1_percol, L1_infilSoil, L1_preEffect, L1_melt)

    implicit none

    class(OutputDataset), intent(inout), target :: nc_mhm

    integer(i4), intent(in) :: sidx !< start index
    integer(i4), intent(in) :: eidx !< end index
    real(dp), intent(in), dimension(:) :: L1_fSealed
    real(dp), intent(in), dimension(:) :: L1_fNotSealed
    real(dp), intent(in), dimension(:) :: L1_inter
    real(dp), intent(in), dimension(:) :: L1_snowPack
    real(dp), intent(in), dimension(:, :) :: L1_soilMoist
    real(dp), intent(in), dimension(:, :) :: L1_soilMoistSat
    real(dp), intent(in), dimension(:) :: L1_sealSTW
    real(dp), intent(in), dimension(:) :: L1_unsatSTW
    real(dp), intent(in), dimension(:) :: L1_satSTW
    real(dp), intent(in), dimension(:) :: L1_neutrons
    real(dp), intent(in), dimension(:) :: L1_pet
    real(dp), intent(in), dimension(:, :) :: L1_aETSoil
    real(dp), intent(in), dimension(:) :: L1_aETCanopy
    real(dp), intent(in), dimension(:) :: L1_aETSealed
    real(dp), intent(in), dimension(:) :: L1_total_runoff
    real(dp), intent(in), dimension(:) :: L1_runoffSeal
    real(dp), intent(in), dimension(:) :: L1_fastRunoff
    real(dp), intent(in), dimension(:) :: L1_slowRunoff
    real(dp), intent(in), dimension(:) :: L1_baseflow
    real(dp), intent(in), dimension(:) :: L1_percol
    real(dp), intent(in), dimension(:, :) :: L1_infilSoil
    real(dp), intent(in), dimension(:) :: L1_preEffect
    real(dp), intent(in), dimension(:) :: L1_melt

    type(OutputVariable), pointer, dimension(:) :: vars
    integer(i4) :: ii, nn

    ii = 0
    vars => nc_mhm%vars

    if (outputFlxState(1)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_inter(sidx : eidx))
    end if

    if (outputFlxState(2)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_snowPack(sidx : eidx))
    end if

    if (outputFlxState(3)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        call vars(ii)%updateVariable(L1_soilMoist(sidx : eidx, nn))
      end do
    end if

    if (outputFlxState(4)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        call vars(ii)%updateVariable(L1_soilMoist(sidx : eidx, nn) / L1_soilMoistSat(sidx : eidx, nn))
      end do
    end if

    if (outputFlxState(5)) then
      ii = ii + 1
      call vars(ii)%updateVariable(sum(L1_soilMoist(sidx : eidx, :), dim = 2) / sum(L1_soilMoistSat(sidx : eidx, :), dim = 2))
    end if

    if (outputFlxState(6)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_sealSTW(sidx : eidx))
    end if

    if (outputFlxState(7)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_unsatSTW(sidx : eidx))
    end if

    if (outputFlxState(8)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_satSTW(sidx : eidx))
    end if

    if (outputFlxState(18)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_neutrons(sidx : eidx))
    end if

    if (outputFlxState(9)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_pet(sidx : eidx))
    end if

    if (outputFlxState(10)) then
      ii = ii + 1
      call vars(ii)%updateVariable(sum(L1_aETSoil(sidx : eidx, :), dim = 2) * L1_fNotSealed(sidx : eidx) &
              + L1_aETCanopy(sidx : eidx) + L1_aETSealed(sidx : eidx) * L1_fSealed(sidx : eidx))
    end if

    if (outputFlxState(11)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_total_runoff(sidx : eidx))
    end if

    if (outputFlxState(12)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_runoffSeal(sidx : eidx) * L1_fSealed(sidx : eidx))
    end if

    if (outputFlxState(13)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_fastRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
    end if

    if (outputFlxState(14)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_slowRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
    end if

    if (outputFlxState(15)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_baseflow(sidx : eidx) * L1_fNotSealed(sidx : eidx))
    end if

    if (outputFlxState(16)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_percol(sidx : eidx) * L1_fNotSealed(sidx : eidx))
    end if

    if (outputFlxState(17)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        call vars(ii)%updateVariable(L1_infilSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
      end do
    end if

    if (outputFlxState(19)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        call vars(ii)%updateVariable(L1_aETSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
      end do
    end if

    if (outputFlxState(20)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_preEffect(sidx : eidx))
    end if

    if (outputFlxState(21)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L1_melt(sidx : eidx))
    end if

  end subroutine mHM_updateDataset

  !> \brief Generate a unit string
  !> \details Generate the unit string for the output variable netcdf attribute based on modeling timestep
  !!
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !!
  !> \return character(16)
  !> \authors David Schaefer
  !> \date June 2015
  function fluxesUnit(iDomain)

    use mo_common_mhm_mrm_variables, only : nTstepDay, simPer

    implicit none

    integer(i4), intent(in) :: iDomain !< selected domain

    character(16) :: fluxesUnit
    real(dp) :: ntsteps

    if (timestep * timestep_model_outputs .eq. 1) then
      fluxesUnit = 'mm h-1'
    else if (timestep_model_outputs > 1) then
      fluxesUnit = 'mm ' // trim(adjustl(num2str(timestep))) // 'h-1'
    else if (timestep_model_outputs .eq. 0) then
      ntsteps = (simPer(iDomain)%julEnd - simPer(iDomain)%julStart + 1) * nTstepDay
      fluxesUnit = 'mm ' // trim(adjustl(num2str(nint(ntsteps)))) // 'h-1'
    else if (timestep_model_outputs .eq. -1) then
      fluxesUnit = 'mm d-1'
    else if (timestep_model_outputs .eq. -2) then
      fluxesUnit = 'mm month-1'
    else if (timestep_model_outputs .eq. -3) then
      fluxesUnit = 'mm a-1'
    else
      fluxesUnit = ''
    end if

  end function fluxesUnit

end module mo_write_fluxes_states
