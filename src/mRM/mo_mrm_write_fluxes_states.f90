!> \file mo_mrm_write_fluxes_states.f90
!> \brief   \copybrief mo_mrm_write_fluxes_states
!> \details \copydetails mo_mrm_write_fluxes_states

!> \brief Creates NetCDF output for different fluxes and state variables of mHM.
!> \details NetCDF is first initialized and later on variables are put to the NetCDF.
!> \changelog
!! - David Schaefer       Aug 2015
!!   - major rewrite
!! - Stephan Thober       Oct 2015
!!   - adapted to mRM
!! - O. Rakovec, R. Kumar Nov 2017
!!   - added project description for the netcdf outputs
!> \authors Matthias Zink
!> \date Apr 2013
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_write_fluxes_states

  use mo_nc_output, only: OutputDataset, OutputVariable, set_attributes, data_dims, data_dtype
  use mo_kind, only : i4, dp
  use mo_common_constants, only : nodata_dp
  use mo_netcdf, only : NcDataset, NcVariable
  use mo_mrm_global_variables, only : output_deflate_level_mrm, output_double_precision_mrm, &
    output_time_reference_mrm, timeStep_model_outputs_mrm, outputFlxState_mrm, riv_temp_pcs, level11
  use mo_common_mHM_mRM_variables, only: timeStep
  use mo_mrm_file, only : file_mrm_output
  use mo_common_variables, only : iFlag_cordinate_sys, level0
  use mo_String_utils, only : num2str

  implicit none

contains

  !> \brief Initialize mRM OutputDataset
  !> \details Create and initialize the output file. If new a new output
  !! variable needs to be written, this is the first of two
  !! procedures to change (second: updateDataset)
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Mueller Jul 2020
  !!   - added output for river temperature
  !> \return type(OutputDataset)
  !> \authors Matthias Zink
  !> \date Apr 2013
  function mRM_OutputDataset(iDomain, mask) result(out)
    implicit none

    integer(i4), intent(in) :: iDomain !< domain id
    logical, intent(in), pointer, dimension(:, :) :: mask !< L11 mask

    type(OutputDataset) :: out

    integer(i4) :: ii, nCells
    character(3) :: dtype
    character(16), dimension(3) :: dims
    type(OutputVariable), dimension(size(outputFlxState_mrm)) :: tmpvars

    out = OutputDataset( &
      iDomain=iDomain, &
      level=level11, &
      file_name=file_mrm_output, &
      double_precision=output_double_precision_mrm, &
      outputs_frequence=timeStep_model_outputs_mrm, &
      time_reference=output_time_reference_mrm &
    )
    dtype = data_dtype(output_double_precision_mrm)
    dims = data_dims()
    nCells = level11(iDomain)%nCells

    ii = 0

    if (outputFlxState_mrm(1)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(out%nc, "Qrouted", dtype, dims, nCells, mask, output_deflate_level_mrm, .true.)
      call set_attributes(tmpvars(ii)%nc, "routed streamflow", "m3 s-1", output_double_precision_mrm)
    end if

    if (outputFlxState_mrm(2) .AND. riv_temp_pcs%active) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(out%nc, "RivTemp", dtype, dims, nCells, mask, output_deflate_level_mrm, .true.)
      call set_attributes(tmpvars(ii)%nc, "routed river temperature", "degC", output_double_precision_mrm)
    end if

    allocate(out%vars(ii))
    out%vars = tmpvars(1:ii)

  end function mRM_OutputDataset

  !> \brief Update all variables.
  !> \details Call the type bound procedure updateVariable for
  !! all output variables. If a new output
  !! variable needs to be written, this is the second
  !! of two procedures to change (first: newOutputDataset)
  !!
  !> \changelog
  !! - L. Samaniego et al. Dec  2013
  !!   - nullify pointer Matthias Zink,        Feb. 2014
  !!   - added aditional output: pet V. Prykhodk, J. Mai,  Nov. 2014
  !!   - adding new variable infilSoil
  !!   - case 16 David Schaefer      , Jun. 2015
  !!   - major rewrite
  !! - Stephan Thober      Oct  2015
  !!   - adapted to mRM
  !! - Robert Schweppe     Jun  2018
  !!   - refactoring and reformatting
  !! - Sebastian Mueller   Jul  2020
  !!   - add river temperature output (optional)
    !> \authors Matthias Zink
  !> \date Apr 2013
  subroutine mRM_updateDataset(nc_mrm, L11_Qmod, L11_riv_temp)

    use mo_mrm_global_variables, only : outputFlxState_mrm

    implicit none

    class(OutputDataset), intent(inout), target :: nc_mrm
    real(dp), intent(in), dimension(:) :: L11_Qmod
    real(dp), intent(in), dimension(:), optional :: L11_riv_temp

    type(OutputVariable), pointer, dimension(:) :: vars

    integer(i4) :: ii

    ii = 0
    vars => nc_mrm%vars

    if (outputFlxState_mrm(1)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L11_Qmod)
    end if

    if (outputFlxState_mrm(2) .AND. riv_temp_pcs%active .AND. present(L11_riv_temp)) then
      ii = ii + 1
      call vars(ii)%updateVariable(L11_riv_temp)
    end if

  end subroutine mRM_updateDataset

  !> \brief Initialize groundwater coupling OutputDataset
  !> \details Create and initialize the output file. If new a new output
  !! variable needs to be written, this is the first of two
  !! procedures to change (second: updateDataset)
  !> \return type(OutputDataset)
  !> \authors Sebastian Mueller
  !> \date Dec 2022
  function GW_OutputDataset(iDomain, mask) result(out)
    implicit none

    integer(i4), intent(in) :: iDomain !< domain id
    logical, intent(in), target, dimension(:, :) :: mask !< L11 mask

    type(OutputDataset) :: out

    integer(i4) :: nCells
    character(3) :: dtype
    character(16), dimension(3) :: dims

    out = OutputDataset( &
      iDomain=iDomain, &
      level=level0, &
      file_name='mRM_riverhead_' // trim(num2str(iDomain, '(i3.3)')) // '.nc', &
      double_precision=output_double_precision_mrm, &
      outputs_frequence=timeStep_model_outputs_mrm, &
      time_reference=output_time_reference_mrm &
    )
    dtype = data_dtype(output_double_precision_mrm)
    dims = data_dims()
    nCells = level0(iDomain)%nCells

    allocate(out%vars(1))
    out%vars(1) = OutputVariable(out%nc, "riverhead", dtype, dims, nCells, mask, output_deflate_level_mrm, .true.)
    call set_attributes(out%vars(1)%nc, "simulated riverhead at each node at level 0", "m", output_double_precision_mrm)

  end function GW_OutputDataset

  !> \brief Update riverhead.
  !> \details Call the type bound procedure updateVariable for
  !! all output variables. If a new output
  !! variable needs to be written, this is the second
  !! of two procedures to change (first: newOutputDataset)
  !> \authors Sebastian Mueller
  !> \date Dec 2022
  subroutine GW_updateDataset(nc_gw, L0_river_head)

    implicit none

    class(OutputDataset), intent(inout), target :: nc_gw
    real(dp), intent(in), dimension(:) :: L0_river_head

    type(OutputVariable), pointer, dimension(:) :: vars

    vars => nc_gw%vars
    call vars(1)%updateVariable(L0_river_head)

  end subroutine GW_updateDataset

end module mo_mrm_write_fluxes_states
