!>       \file mo_write_fluxes_states.f90

!>       \brief Creates NetCDF output for different fluxes and state variables of mHM.

!>       \details NetCDF is first initialized and later on variables are put to the NetCDF.

!>       \authors Matthias Zink

!>       \date Apr 2013

! Modifications:
! David Schaefer       Aug 2015 - major rewrite
! O. Rakovec, R. Kumar Nov 2017 - added project description for the netcdf outputs

module mo_write_fluxes_states

  use mo_kind, only : i4, dp
  use mo_string_utils, only : num2str
  use mo_common_variables, only : project_details, setup_description, simulation_type, &
          Conventions, contact, mHM_details, history
  use mo_common_constants, only : nodata_dp
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable
  use mo_global_variables, only : output_deflate_level, output_double_precision

  implicit none

  type OutputVariable
    type(NcVariable) :: nc                 !> NcDataset which contains the variable
    logical :: avg = .false.      !> average data before writing
    logical, pointer :: mask(:, :)          !> mask to reconstruct data
    real(dp), allocatable :: data(:)            !> store the data between writes
    integer(i4) :: counter = 0        !> count the number of updateVariable calls

  contains
    procedure, public :: updateVariable
    procedure, public :: writeVariableTimestep

  end type OutputVariable

  ! constructor interface
  interface OutputVariable
    procedure newOutputVariable
  end interface OutputVariable

  type OutputDataset
    integer(i4) :: iDomain      !> domain id
    type(NcDataset) :: nc          !> NcDataset to write
    type(OutputVariable), allocatable :: vars(:)     !> store all created (dynamic) variables
    integer(i4) :: counter = 0 !> count written time steps

  contains
    procedure, public :: updateDataset
    procedure, public :: writeTimestep
    procedure, public :: close

  end type OutputDataset

  ! constructor interface
  interface OutputDataset
    procedure newOutputDataset
  end interface OutputDataset

  private

  public :: OutputDataset

contains

  !------------------------------------------------------------------
  !    NAME
  !        newOutputVariable

  !    PURPOSE
  !>       \brief Initialize OutputVariable

  !>       \details TODO: add description

  !>       \return type(OutputVariable)

  !    INTENT(IN)
  !>       \param[in] "type(NcDataset) :: nc"               -> NcDataset which contains the variable
  !>       \param[in] "character(*) :: name"
  !>       \param[in] "character(*) :: dtype"
  !>       \param[in] "character(16), dimension(3) :: dims"
  !>       \param[in] "integer(i4) :: ncells"               -> number of cells in domain
  !>       \param[in] "logical, dimension(:, :) :: mask"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: avg" -> average the data before writing

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! David Schaefer Nov 2017 - , added NcVariable initialization
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function newOutputVariable(nc, name, dtype, dims, ncells, mask, avg) result(out)
    implicit none

    ! -> NcDataset which contains the variable
    type(NcDataset), intent(in) :: nc

    character(*), intent(in) :: name

    character(*), intent(in) :: dtype

    character(16), intent(in), dimension(3) :: dims

    ! -> number of cells in domain
    integer(i4), intent(in) :: ncells

    logical, intent(in), target, dimension(:, :) :: mask

    ! -> average the data before writing
    logical, intent(in), optional :: avg

    type(OutputVariable) :: out


    allocate(out%data(ncells))
    out%nc = nc%setVariable(name, dtype, dims, deflate_level = output_deflate_level, shuffle = .true.)
    out%data = 0
    out%mask => mask
    if (present(avg)) out%avg = avg
  end function newOutputVariable

  !------------------------------------------------------------------
  !    NAME
  !        updateVariable

  !    PURPOSE
  !>       \brief Update OutputVariable

  !>       \details Add the array given as actual argument
  !>       to the derived type's component 'data'

  !>       \return type(OutputVariable)

  !    INTENT(INOUT)
  !>       \param[inout] "class(OutputVariable) :: self"

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: data"

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine updateVariable(self, data)
    implicit none

    class(OutputVariable), intent(inout) :: self

    real(dp), intent(in), dimension(:) :: data


    self%data = self%data + data
    self%counter = self%counter + 1

  end subroutine updateVariable

  !------------------------------------------------------------------
  !    NAME
  !        writeVariableTimestep

  !    PURPOSE
  !>       \brief Write timestep to file

  !>       \details Write the content of the derived types's component
  !>       'data' to file, average if necessary

  !    INTENT(INOUT)
  !>       \param[inout] "class(OutputVariable) :: self"

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: timestep" -> index along the time dimension of the netcdf variable

  !    HISTORY
  !>       \authors David Schafer

  !>       \date June 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine writeVariableTimestep(self, timestep)
    implicit none

    class(OutputVariable), intent(inout) :: self

    ! -> index along the time dimension of the netcdf variable
    integer(i4), intent(in) :: timestep


    if (self%avg) then
      self%data = self%data / real(self%counter, dp)
    end if
    call self%nc%setData(unpack(self%data, self%mask, nodata_dp), &
            (/1, 1, timestep/))
    self%data = 0
    self%counter = 0

  end subroutine writeVariableTimestep

  !------------------------------------------------------------------
  !    NAME
  !        newOutputDataset

  !    PURPOSE
  !>       \brief Initialize OutputDataset

  !>       \details Create and initialize the output file. If new a new output
  !>       variable needs to be written, this is the first of two
  !>       procedures to change (second: updateDataset)

  !>       \return type(OutputDataset)

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"             -> domain id
  !>       \param[in] "logical, dimension(:, :) :: mask1" -> L1 mask to reconstruct the data
  !>       \param[in] "integer(i4) :: nCells"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Pallav Shrestha Mar 2020 - iFlag_cordinate_sys based dimensions (dims1)

  function newOutputDataset(iDomain, mask1, nCells) result(out)

    use mo_global_variables, only : outputFlxState
    use mo_mpr_global_variables, only : nSoilHorizons_mHM
    use mo_common_variables, only : iFlag_cordinate_sys

    implicit none

    ! -> domain id
    integer(i4), intent(in) :: iDomain

    ! -> L1 mask to reconstruct the data
    logical, pointer, intent(in), dimension(:, :) :: mask1

    integer(i4), intent(in) :: nCells

    type(OutputDataset) :: out

    integer(i4) :: ii, nn

    character(3) :: dtype

    character(16), dimension(3) :: dims1

    character(16) :: unit

    type(NcDataset) :: nc

    type(OutputVariable), dimension(size(outputFlxState) * nSoilHorizons_mHM) :: tmpvars


    if ( output_double_precision ) then
      dtype = "f64"
    else
      dtype = "f32"
    end if
    unit = fluxesUnit(iDomain)

    if (iFlag_cordinate_sys == 0) then
      dims1 = (/"easting ", "northing", "time    "/) ! X & Y coordinate system
    else
      dims1 = (/"lon ", "lat ", "time"/) ! lat & lon coordinate system
    endif

    nc = createOutputFile(iDomain)

    ii = 0

    if (outputFlxState(1)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "interception", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "canopy interception storage", "mm")
    end if

    if (outputFlxState(2)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "snowpack", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(tmpvars(ii), "depth of snowpack", "mm")
    end if

    if (outputFlxState(3)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(nc, "SWC_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims1, nCells, mask1, .true.)
        call writeVariableAttributes(tmpvars(ii), &
                'soil water content of soil layer' // trim(num2str(nn)), "mm")
      end do
    end if

    if (outputFlxState(4)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(nc, "SM_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims1, nCells, mask1, .true.)
        call writeVariableAttributes(tmpvars(ii), &
                'volumetric soil moisture of soil layer' // trim(num2str(nn)), "mm mm-1")
      end do
    end if

    if (outputFlxState(5)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "SM_Lall", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "average soil moisture over all layers", "mm mm-1")
    end if

    if (outputFlxState(6)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "sealedSTW", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "reservoir of sealed areas (sealedSTW)", "mm")
    end if

    if (outputFlxState(7)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "unsatSTW", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "reservoir of unsaturated zone", "mm")
    end if

    if (outputFlxState(8)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "satSTW", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "water level in groundwater reservoir", "mm")
    end if

    if (outputFlxState(18)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "neutrons", dtype, dims1, nCells, mask1, .true.)
      call writeVariableAttributes(&
              tmpvars(ii), "ground albedo neutrons", "cph")
    end if

    if (outputFlxState(9)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "PET", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "potential Evapotranspiration", trim(unit))
    end if

    if (outputFlxState(10)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "aET", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "actual Evapotranspiration", trim(unit))
    end if

    if (outputFlxState(11)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "Q", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "total runoff generated by every cell", trim(unit))
    end if

    if (outputFlxState(12)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "QD", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(tmpvars(ii), &
              "direct runoff generated by every cell (runoffSeal)", trim(unit))
    end if

    if (outputFlxState(13)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "QIf", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(tmpvars(ii), &
              "fast interflow generated by every cell (fastRunoff)", trim(unit))
    end if

    if (outputFlxState(14)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "QIs", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(tmpvars(ii), &
              "slow interflow generated by every cell (slowRunoff)", trim(unit))
    end if

    if (outputFlxState(15)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "QB", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "baseflow generated by every cell", trim(unit))
    end if

    if (outputFlxState(16)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "recharge", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "groundwater recharge", trim(unit))
    end if

    if (outputFlxState(17)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(&
                nc, "soil_infil_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims1, nCells, mask1)
        call writeVariableAttributes(tmpvars(ii), &
                "infiltration flux from soil layer" // trim(num2str(nn)), unit)
      end do
    end if

    if (outputFlxState(19)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
        tmpvars(ii) = OutputVariable(&
                nc, "aET_L" // trim(num2str(nn, '(i2.2)')), &
                dtype, dims1, nCells, mask1)
        call writeVariableAttributes(tmpvars(ii), &
                'actual Evapotranspiration from soil layer' // trim(num2str(nn)), &
                "mm " // trim(unit))
      end do
    end if

    if (outputFlxState(20)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "preEffect", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "effective precipitation", trim(unit))
    end if

    if (outputFlxState(21)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "Qsm", dtype, dims1, nCells, mask1)
      call writeVariableAttributes(&
              tmpvars(ii), "Average liquid water generated from solid to liquid phase change in the snow", trim(unit))
    end if

    out%vars = tmpvars(1 : ii)
    out%nc = nc
    out%iDomain = iDomain

  end function newOutputDataset

  !------------------------------------------------------------------
  !    NAME
  !        updateDataset

  !    PURPOSE
  !>       \brief Update all variables.

  !>       \details Call the type bound procedure updateVariable for
  !>       all output variables. If a new output
  !>       variable needs to be written, this is the second
  !>       of two procedures to change (first: newOutputDataset)

  !    INTENT(INOUT)
  !>       \param[inout] "class(OutputDataset) :: self"

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: sidx, eidx"
  !>       \param[in] "integer(i4) :: sidx, eidx"
  !>       \param[in] "real(dp), dimension(:) :: L1_fSealed"
  !>       \param[in] "real(dp), dimension(:) :: L1_fNotSealed"
  !>       \param[in] "real(dp), dimension(:) :: L1_inter"
  !>       \param[in] "real(dp), dimension(:) :: L1_snowPack"
  !>       \param[in] "real(dp), dimension(:, :) :: L1_soilMoist"
  !>       \param[in] "real(dp), dimension(:, :) :: L1_soilMoistSat"
  !>       \param[in] "real(dp), dimension(:) :: L1_sealSTW"
  !>       \param[in] "real(dp), dimension(:) :: L1_unsatSTW"
  !>       \param[in] "real(dp), dimension(:) :: L1_satSTW"
  !>       \param[in] "real(dp), dimension(:) :: L1_neutrons"
  !>       \param[in] "real(dp), dimension(:) :: L1_pet"
  !>       \param[in] "real(dp), dimension(:, :) :: L1_aETSoil"
  !>       \param[in] "real(dp), dimension(:) :: L1_aETCanopy"
  !>       \param[in] "real(dp), dimension(:) :: L1_aETSealed"
  !>       \param[in] "real(dp), dimension(:) :: L1_total_runoff"
  !>       \param[in] "real(dp), dimension(:) :: L1_runoffSeal"
  !>       \param[in] "real(dp), dimension(:) :: L1_fastRunoff"
  !>       \param[in] "real(dp), dimension(:) :: L1_slowRunoff"
  !>       \param[in] "real(dp), dimension(:) :: L1_baseflow"
  !>       \param[in] "real(dp), dimension(:) :: L1_percol"
  !>       \param[in] "real(dp), dimension(:, :) :: L1_infilSoil"
  !>       \param[in] "real(dp), dimension(:) :: L1_preEffect"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! L. Samaniego et al. Dec  2013 - nullify pointer Matthias Zink,        Feb. 2014
  !                              - added aditional output: pet V. Prykhodk, J. Mai,  Nov. 2014
  !                              - adding new variable infilSoil
  !                              - case 16 David Schaefer      , Jun. 2015
  !                              - major rewrite
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine updateDataset(self, sidx, eidx, L1_fSealed, L1_fNotSealed, L1_inter, L1_snowPack, L1_soilMoist, &
                          L1_soilMoistSat, L1_sealSTW, L1_unsatSTW, L1_satSTW, L1_neutrons, L1_pet, L1_aETSoil, &
                          L1_aETCanopy, L1_aETSealed, L1_total_runoff, L1_runoffSeal, L1_fastRunoff, L1_slowRunoff, &
                          L1_baseflow, L1_percol, L1_infilSoil, L1_preEffect, L1_melt)

    use mo_global_variables, only : outputFlxState
    use mo_mpr_global_variables, only : nSoilHorizons_mHM

    implicit none

    class(OutputDataset), intent(inout), target :: self

    integer(i4), intent(in) :: sidx, eidx

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
    vars => self%vars

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

  end subroutine updateDataset

  !------------------------------------------------------------------
  !    NAME
  !        writeTimestep

  !    PURPOSE
  !>       \brief Write all accumulated data.

  !>       \details Write all accumulated and potentially averaged
  !>       data to disk.

  !>       \return type(OutputVariable)

  !    INTENT(INOUT)
  !>       \param[inout] "class(OutputDataset) :: self"

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: timestep" The model timestep to write

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine writeTimestep(self, timestep)
    implicit none

    class(OutputDataset), intent(inout), target :: self

    ! The model timestep to write
    integer(i4), intent(in) :: timestep

    integer(i4) :: ii

    type(NcVariable) :: tvar


    self%counter = self%counter + 1

    ! add to time variable
    tvar = self%nc%getVariable("time")
    call tvar%setData(timestep, (/self%counter/))

    do ii = 1, size(self%vars)
      call self%vars(ii)%writeVariableTimestep(self%counter)
    end do

  end subroutine writeTimestep

  !------------------------------------------------------------------
  !    NAME
  !        close

  !    PURPOSE
  !>       \brief Close the file

  !>       \details Close the file associated with variable of
  !>       type(OutputDataset)

  !    HISTORY
  !>       \authors Rohini Kumar & Stephan Thober

  !>       \date August 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine close(self)

    use mo_String_utils, only : num2str
    use mo_common_variables, only : dirOut
    use mo_message, only : message

    implicit none

    class(OutputDataset) :: self


    call self%nc%close()
    call message('  OUTPUT: saved netCDF file for domain', trim(num2str(self%iDomain)))
    call message('    to ', trim(dirOut(self%iDomain)))

  end subroutine close

  !------------------------------------------------------------------
  !    NAME
  !        createOutputFile

  !    PURPOSE
  !>       \brief Create and initialize output file for X & Y coordinate system

  !>       \details Create output file, write all non-dynamic variables
  !>       and global attributes for the given domain for X & Y coordinate system

  !>       \return type(NcDataset)

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" -> domain id

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Stephan Thober  Oct 2015 - added actual version of mHM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Pallav Shrestha Mar 2020 - output file lat and lon are 1d or 2d based on coordinate system

  function createOutputFile(iDomain) result(nc)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : dirOut, level1, iFlag_cordinate_sys
    use mo_file, only : version
    use mo_grid, only : geoCoordinates, mapCoordinates
    use mo_julian, only : dec2date

    implicit none

    ! -> domain id
    integer(i4), intent(in) :: iDomain

    type(NcDataset) :: nc

    type(NcDimension), dimension(3) :: dimids1

    type(NcVariable) :: var

    integer(i4) :: day, month, year

    character(1028) :: fname

    character(128) :: unit, date, time, datetime

    real(dp), allocatable, dimension(:) :: easting, northing

    real(dp), allocatable, dimension(:) :: lat1d, lon1d    ! 1D lat lon vectors. Used if coordinate system is lat & lon

    real(dp), allocatable, dimension(:, :) :: lat2d, lon2d ! temporary storage of mHM's 2D latlon array.
                                                           ! Used as 2d lat lon arrays if coordinate system is X & Y
    character(3) :: dtype

    if ( output_double_precision ) then
      dtype = "f64"
    else
      dtype = "f32"
    end if

    fname = trim(dirOut(iDomain)) // 'mHM_Fluxes_States.nc'
    call geoCoordinates(level1(iDomain), lat2d, lon2d)

    nc = NcDataset(trim(fname), "w")

    ! set the horizonal dimensions
    if (iFlag_cordinate_sys == 0) then

      ! X & Y coordinate system; 2D lat lon!
      !============================================================
      call mapCoordinates(level1(iDomain), northing, easting)

      dimids1 = (/ &
        nc%setDimension("easting", size(easting)), &
        nc%setDimension("northing", size(northing)), &
        nc%setDimension("time", 0) &
      /)
      ! northing
      var = nc%setVariable("northing", dtype, (/ dimids1(2) /))
      call var%setData(northing)
      call var%setAttribute("units", "m")
      call var%setAttribute("long_name", "y-coordinate in the given coordinate system")
      ! easting
      var = nc%setVariable("easting", dtype, (/ dimids1(1) /))
      call var%setData(easting)
      call var%setAttribute("units", "m")
      call var%setAttribute("long_name", "x-coordinate in the given coordinate system")
      ! lon
      var = nc%setVariable("lon", dtype, dimids1(1 : 2))
      call var%setFillValue(nodata_dp)
      call var%setData(lon2d)
      call var%setAttribute("units", "degrees_east")
      call var%setAttribute("long_name", "longitude")
      call var%setAttribute("missing_value", nodata_dp)
      ! lat
      var = nc%setVariable("lat", dtype, dimids1(1 : 2))
      call var%setFillValue(nodata_dp)
      call var%setData(lat2d)
      call var%setAttribute("units", "degrees_north")
      call var%setAttribute("long_name", "latitude")
      call var%setAttribute("missing_value", nodata_dp)

    else

      ! lat & lon coordinate system; 1D lat lon!
      !============================================================
      lat1d = lat2d(1, :) ! first row info is sufficient
      lon1d = lon2d(:, 1) ! first column info is sufficient
      dimids1 = (/ &
        nc%setDimension("lon", size(lon1d)), &
        nc%setDimension("lat", size(lat1d)), &
        nc%setDimension("time", 0) &
      /)
      ! lon
      var = nc%setVariable("lon", dtype, (/ dimids1(1) /)) ! sufficient to store lon as vector
      call var%setFillValue(nodata_dp)
      call var%setData(lon1d)
      call var%setAttribute("units", "degrees_east")
      call var%setAttribute("long_name", "longitude")
      call var%setAttribute("missing_value", nodata_dp)
      ! lat
      var = nc%setVariable("lat", dtype, (/ dimids1(2) /)) ! sufficient to store lat as vector
      call var%setFillValue(nodata_dp)
      call var%setData(lat1d)
      call var%setAttribute("units", "degrees_north")
      call var%setAttribute("long_name", "latitude")
      call var%setAttribute("missing_value", nodata_dp)

    endif

    ! set record dimension
    ! time units
    call dec2date(real(evalPer(iDomain)%julStart, dp), dd = day, mm = month, yy = year)
    write(unit, "('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    var = nc%setVariable("time", "i32", (/ dimids1(3) /))
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name", "time")


    ! global attributes
    call date_and_time(date = date, time = time)
    write(datetime, "(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1 : 4), &
            date(5 : 6), date(7 : 8), time(1 : 2), time(3 : 4), time(5 : 6)

    call nc%setAttribute("project", project_details)
    call nc%setAttribute("setup_description", setup_description)
    call nc%setAttribute("simulation_type", simulation_type)
    call nc%setAttribute("Conventions", Conventions)
    call nc%setAttribute("contact", contact)
    call nc%setAttribute("mHM_details", trim(mHM_details) // ", release mHMv" // trim(version))
    call nc%setAttribute("history", trim(datetime) // ", " // history)
    call nc%setAttribute("title", "mHMv"//trim(version)//" "//trim(simulation_type)//" outputs")
    call nc%setAttribute("creation_date", datetime)

  end function createOutputFile

  !------------------------------------------------------------------
  !    NAME
  !        writeVariableAttributes

  !    PURPOSE
  !>       \brief Write output variable attributes

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "type(OutputVariable) :: var"
  !>       \param[in] "character(*) :: long_name, unit" -> variable name
  !>       \param[in] "character(*) :: long_name, unit" -> physical unit

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine writeVariableAttributes(var, long_name, unit)
    implicit none

    type(OutputVariable), intent(inout) :: var

    ! -> physical unit
    character(*), intent(in) :: long_name, unit


    call var%nc%setFillValue(nodata_dp)
    call var%nc%setAttribute("long_name", long_name)
    call var%nc%setAttribute("units", unit)
    call var%nc%setAttribute("scale_factor", 1.0_dp)
    call var%nc%setAttribute("missing_value", nodata_dp)
    call var%nc%setAttribute("coordinates", "lat lon")

  end subroutine writeVariableAttributes


  !------------------------------------------------------------------
  !    NAME
  !        fluxesUnit

  !    PURPOSE
  !>       \brief Generate a unit string

  !>       \details Generate the unit string for the output variable
  !>       netcdf attribute based on modeling timestep

  !>       \return character(16)

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function fluxesUnit(iDomain)

    use mo_common_mhm_mrm_variables, only : nTstepDay, simPer, timestep
    use mo_global_variables, only : timeStep_model_outputs

    implicit none

    integer(i4), intent(in) :: iDomain

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
