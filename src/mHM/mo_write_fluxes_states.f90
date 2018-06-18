!>       \file mo_write_fluxes_states.f90

!>       \brief Creates NetCDF output for different fluxes and state variables of mHM.

!>       \details NetCDF is first initialized and later on variables are put to the NetCDF.

!>       \authors s Matthias Zink

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
    integer(i4) :: ibasin      !> basin id
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
#ifdef pgiFortran154
  public :: newOutputDataset
#endif

contains

  !------------------------------------------------------------------
  !    NAME
  !        newOutputVariable

  !    PURPOSE
  !>       \brief Initialize OutputVariable

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "type(NcDataset) :: nc"               -> NcDataset which contains the variable
  !>       \param[in] "character(*) :: name"                
  !>       \param[in] "character(*) :: dtype"               
  !>       \param[in] "character(16), dimension(3) :: dims" 
  !>       \param[in] "integer(i4) :: ncells"               -> number of cells in basin
  !>       \param[in] "logical, dimension(:, :) :: mask"    

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: avg" -> average the data before writing

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! David Schaefer Nov 2017 - , added NcVariable initialization

  function newOutputVariable(nc, name, dtype, dims, ncells, mask, avg) result(out)
    implicit none

    ! -> NcDataset which contains the variable
    type(NcDataset), intent(in) :: nc

    character(*), intent(in) :: name

    character(*), intent(in) :: dtype

    character(16), intent(in), dimension(3) :: dims

    ! -> number of cells in basin
    integer(i4), intent(in) :: ncells

    logical, intent(in), target, dimension(:, :) :: mask

    ! -> average the data before writing
    logical, intent(in), optional :: avg

    type(OutputVariable) :: out


    allocate(out%data(ncells))
    out%nc = nc%setVariable(name, dtype, dims, deflate_level = 1, shuffle = .true.)
    out%data = 0
    out%mask => mask
    if (present(avg)) out%avg = avg
  end function newOutputVariable

  !------------------------------------------------------------------
  !     NAME
  !         updateVariable
  !
  !     PURPOSE
  !>        \brief Update OutputVariable
  !>        \details Add the array given as actual argument
  !>                 to the derived type's component 'data'
  !
  !     CALLING SEQUENCE
  !         -> with nc of type(OutputVariable):
  !         call var%updateVariable(data)
  !
  !     INTENT(IN)
  !>        \param[in] "type(NcDataset)   :: nc"        -> NcDataset which contains the variable
  !>        \param[in] "integer(i4)       :: ncells"    -> number of cells in basin
  !>        \param[in] "logical, target   :: mask(:,:)" -> mask to reconstruct data
  !>        \param[in] "logical, optional :: avg"       -> average the data before writing
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         \return type(OutputVariable)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine updateVariable(self, data)
    class(OutputVariable), intent(inout) :: self
    real(dp), intent(in) :: data(:)

    self%data = self%data + data
    self%counter = self%counter + 1

  end subroutine updateVariable

  !------------------------------------------------------------------
  !     NAME
  !         writeVariableTimestep
  !
  !     PURPOSE
  !>        \brief Write timestep to file
  !>        \details Write the content of the derived types's component
  !>                 'data' to file, average if necessary
  !
  !     CALLING SEQUENCE
  !         -> with var of type(OutputVariable):
  !         call var%updateVariable(data)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: timestep"
  !>            -> index along the time dimension of the netcdf variable
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schafer
  !>        \date June 2015
  !         Modified:
  !             David Schaefer, Sep. 2015 - bugfix
  subroutine writeVariableTimestep(self, timestep)
    class(OutputVariable), intent(inout) :: self
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
  !     NAME
  !         newOutputDataset
  !
  !     PURPOSE
  !>        \brief Initialize OutputDataset
  !>        \details Create and initialize the output file. If new a new output
  !>                 variable needs to be written, this is the first of two
  !>                 procedures to change (second: updateDataset)
  !
  !     CALLING SEQUENCE
  !         nc = OutputDataset(ibasin, mask1)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: ibasin" -> basin id
  !>        \param[in] "logical     :: mask1"  -> L1 mask to reconstruct the data
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         \return type(OutputDataset)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                               during writing of the netcdf file
  !             Matthias Zink       , Feb. 2014 - added aditional output: pet
  !             V. Prykhodk, J. Mai , Nov. 2014 - adding new variable infilSoil - case 16
  !             David Schaefer      , Jun. 2015 - major rewrite
  !             David Schaefer      , Nov. 2016 - moved NcVariable initialization to newOutputVariable
  function newOutputDataset(ibasin, mask1, nCells) result(out)

    use mo_global_variables, only : outputFlxState
    use mo_mpr_global_variables, only : nSoilHorizons_mHM

    integer(i4), intent(in) :: ibasin
    logical, target, intent(in) :: mask1(:, :)
    integer(i4), intent(in) :: nCells
    type(OutputDataset) :: out
    ! local
    integer(i4) :: ii, nn
    character(3) :: dtype
    character(16) :: dims1(3), unit
    type(NcDataset) :: nc
    type(OutputVariable) :: tmpvars(size(outputFlxState) * nSoilHorizons_mHM)

    dtype = "f64"
    unit = fluxesUnit(ibasin)
    dims1 = (/"easting ", "northing", "time    "/)
    nc = createOutputFile(ibasin)

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

    out%vars = tmpvars(1 : ii)
    out%nc = nc
    out%ibasin = ibasin

  end function newOutputDataset

  !------------------------------------------------------------------
  !     NAME
  !         updateDataset
  !
  !     PURPOSE
  !>        \brief Update all variables.
  !>        \details Call the type bound procedure updateVariable for
  !>                 all output variables. If a new output
  !>                 variable needs to be written, this is the second
  !>                 of two procedures to change (first: newOutputDataset)
  !
  !     CALLING SEQUENCE
  !        with nc of type(OutputDataset):
  !        call nc%updateDataset(&
  !             self         , sidx         , eidx,           ,    &
  !             L1_fSealed   , L1_fNotSealed, L1_inter        ,    &
  !             L1_snowPack  , L1_soilMoist , L1_soilMoistSat ,    &
  !             L1_sealSTW   , L1_unsatSTW  , L1_satSTW       ,    &
  !             L1_neutrons  , L1_pet       , L1_aETSoil      ,    &
  !             L1_aETCanopy , L1_aETSealed , L1_total_runoff ,    &
  !             L1_runoffSeal, L1_fastRunoff, L1_slowRunoff   ,    &
  !             L1_baseflow  , L1_percol    , L1_infilSoil    ,    &
  !             L1_preEffect                )
  !
  !
  !     INTENT(IN)
  !>             \param[in] "sidx"        -> start index of the basin related data in L1_* arguments
  !>             \param[in] "eidx"        -> end index of the basin related data in L1_* arguments
  !>             \param[in] "L1_fSealed"
  !>             \param[in] "L1_fNotSealed"
  !>             \param[in] "L1_inter"
  !>             \param[in] "L1_snowPack"
  !>             \param[in] "L1_soilMoist"
  !>             \param[in] "L1_soilMoistSat"
  !>             \param[in] "L1_sealSTW"
  !>             \param[in] "L1_unsatSTW"
  !>             \param[in] "L1_satSTW"
  !>             \param[in] "L1_neutrons"
  !>             \param[in] "L1_pet"
  !>             \param[in] "L1_aETSoil"
  !>             \param[in] "L1_aETCanopy"
  !>             \param[in] "L1_aETSealed"
  !>             \param[in] "L1_total_runoff"
  !>             \param[in] "L1_runoffSeal"
  !>             \param[in] "L1_fastRunoff"
  !>             \param[in] "L1_slowRunoff"
  !>             \param[in] "L1_baseflow"
  !>             \param[in] "L1_percol"
  !>             \param[in] "L1_infilSoil"
  !>             \param[in] "L1_preEffect"
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                               during writing of the netcdf file
  !             L. Samaniego et al.,  Dec  2013 - nullify pointer
  !             Matthias Zink,        Feb. 2014 - added aditional output: pet
  !             V. Prykhodk, J. Mai,  Nov. 2014 - adding new variable infilSoil - case 16
  !             David Schaefer      , Jun. 2015 - major rewrite
  subroutine updateDataset(&
          self, sidx, eidx, &
          L1_fSealed, L1_fNotSealed, L1_inter, &
          L1_snowPack, L1_soilMoist, L1_soilMoistSat, &
          L1_sealSTW, L1_unsatSTW, L1_satSTW, &
          L1_neutrons, L1_pet, L1_aETSoil, &
          L1_aETCanopy, L1_aETSealed, L1_total_runoff, &
          L1_runoffSeal, L1_fastRunoff, L1_slowRunoff, &
          L1_baseflow, L1_percol, L1_infilSoil, &
          L1_preEffect)

    use mo_global_variables, only : outputFlxState
    use mo_mpr_global_variables, only : nSoilHorizons_mHM

    class(OutputDataset), intent(inout), target :: self
    integer(i4), intent(in) :: sidx, eidx
    real(dp), intent(in) :: L1_fSealed(:)
    real(dp), intent(in) :: L1_fNotSealed(:)
    ! states
    real(dp), intent(in) :: L1_inter(:)
    real(dp), intent(in) :: L1_snowPack(:)
    real(dp), intent(in) :: L1_soilMoist(:, :)
    real(dp), intent(in) :: L1_soilMoistSat(:, :)
    real(dp), intent(in) :: L1_sealSTW(:)
    real(dp), intent(in) :: L1_unsatSTW(:)
    real(dp), intent(in) :: L1_satSTW(:)
    real(dp), intent(in) :: L1_neutrons(:)
    ! fluxes,
    real(dp), intent(in) :: L1_pet(:)
    real(dp), intent(in) :: L1_aETSoil(:, :)
    real(dp), intent(in) :: L1_aETCanopy(:)
    real(dp), intent(in) :: L1_aETSealed(:)
    real(dp), intent(in) :: L1_total_runoff(:)
    real(dp), intent(in) :: L1_runoffSeal(:)
    real(dp), intent(in) :: L1_fastRunoff(:)
    real(dp), intent(in) :: L1_slowRunoff(:)
    real(dp), intent(in) :: L1_baseflow(:)
    real(dp), intent(in) :: L1_percol(:)
    real(dp), intent(in) :: L1_infilSoil(:, :)
    real(dp), intent(in) :: L1_preEffect(:)

    ! local
    type(OutputVariable), pointer :: vars(:)
    integer(i4) :: ii, nn

    ii = 0
    vars => self%vars

    if (outputFlxState(1)) then
      ii = ii + 1
#ifdef pgiFortran
       call updateVariable(vars(ii), L1_inter(sidx : eidx))
#else
       call vars(ii)%updateVariable(L1_inter(sidx : eidx))
#endif
    end if

    if (outputFlxState(2)) then
      ii = ii + 1
#ifdef pgiFortran
       call updateVariable(vars(ii), L1_snowPack(sidx : eidx))
#else
       call vars(ii)%updateVariable(L1_snowPack(sidx : eidx))
#endif
    end if

    if (outputFlxState(3)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
#ifdef pgiFortran
        call updateVariable(vars(ii), L1_soilMoist(sidx : eidx, nn))
#else
        call vars(ii)%updateVariable(L1_soilMoist(sidx : eidx, nn))
#endif
       end do
    end if

    if (outputFlxState(4)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
#ifdef pgiFortran
        call updateVariable(vars(ii), L1_soilMoist(sidx : eidx, nn) &
                / L1_soilMoistSat(sidx : eidx, nn))
#else
        call vars(ii)%updateVariable(L1_soilMoist(sidx : eidx, nn) &
                / L1_soilMoistSat(sidx : eidx, nn))
#endif
       end do
    end if

    if (outputFlxState(5)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), sum(L1_soilMoist(sidx : eidx, :), dim = 2) &
              / sum(L1_soilMoistSat(sidx : eidx, :), dim = 2))
#else
      call vars(ii)%updateVariable(sum(L1_soilMoist(sidx : eidx, :), dim = 2) &
              / sum(L1_soilMoistSat(sidx : eidx, :), dim = 2))
#endif
    end if

    if (outputFlxState(6)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_sealSTW(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_sealSTW(sidx : eidx))
#endif
    end if

    if (outputFlxState(7)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_unsatSTW(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_unsatSTW(sidx : eidx))
#endif
    end if

    if (outputFlxState(8)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_satSTW(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_satSTW(sidx : eidx))
#endif
    end if

    if (outputFlxState(18)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_neutrons(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_neutrons(sidx : eidx))
#endif
    end if

    if (outputFlxState(9)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_pet(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_pet(sidx : eidx))
#endif
    end if

    if (outputFlxState(10)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), sum(L1_aETSoil(sidx : eidx, :), dim = 2) * L1_fNotSealed(sidx : eidx) &
              + L1_aETCanopy(sidx : eidx) + L1_aETSealed(sidx : eidx) * L1_fSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(sum(L1_aETSoil(sidx : eidx, :), dim = 2) * L1_fNotSealed(sidx : eidx) &
              + L1_aETCanopy(sidx : eidx) + L1_aETSealed(sidx : eidx) * L1_fSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(11)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_total_runoff(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_total_runoff(sidx : eidx))
#endif
    end if

    if (outputFlxState(12)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_runoffSeal(sidx : eidx) * L1_fSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_runoffSeal(sidx : eidx) * L1_fSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(13)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_fastRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_fastRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(14)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_slowRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_slowRunoff(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(15)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_baseflow(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_baseflow(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(16)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L1_percol(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#else
      call vars(ii)%updateVariable(L1_percol(sidx : eidx) * L1_fNotSealed(sidx : eidx))
#endif
    end if

    if (outputFlxState(17)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
#ifdef pgiFortran
        call updateVariable(vars(ii), L1_infilSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
#else
        call vars(ii)%updateVariable(L1_infilSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
#endif
       end do
    end if

    if (outputFlxState(19)) then
      do nn = 1, nSoilHorizons_mHM
        ii = ii + 1
#ifdef pgiFortran
        call updateVariable(vars(ii), L1_aETSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
#else
        call vars(ii)%updateVariable(L1_aETSoil(sidx : eidx, nn) * L1_fNotSealed(sidx : eidx))
#endif
       end do
    end if

    if (outputFlxState(20)) then
      ii = ii + 1
#ifdef pgiFortran
        call updateVariable(vars(ii), L1_preEffect(sidx : eidx))
#else
        call vars(ii)%updateVariable(L1_preEffect(sidx : eidx))
#endif
    end if

  end subroutine updateDataset

  !------------------------------------------------------------------
  !     NAME
  !         writeTimestep
  !
  !     PURPOSE
  !>        \brief Write all accumulated data.
  !>        \details Write all accumulated and potentially averaged
  !>                 data to disk.
  !
  !     CALLING SEQUENCE
  !         -> with nc of type(OutputDataset)
  !         call nc%writeTimestep(timestep)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: timestep" The model timestep to write
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         \return type(OutputVariable)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine writeTimestep(self, timestep)
    class(OutputDataset), intent(inout), target :: self
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
  !     NAME
  !         close
  !
  !     PURPOSE
  !>        \brief Close the file
  !>        \details Close the file associated with variable of
  !>                 type(OutputDataset)
  !
  !     CALLING SEQUENCE
  !         -> with nc of type(OutputDataset):
  !         call nc%close()
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Rohini Kumar & Stephan Thober
  !>        \date August 2013
  !         Modified:
  !             David Schaefer, June 2015 - adapted to new structure
  subroutine close(self)

    use mo_String_utils, only : num2str
    use mo_message, only : message
    use mo_common_variables, only : dirOut

    class(OutputDataset) :: self
    call self%nc%close()
    call message('  OUTPUT: saved netCDF file for basin', trim(num2str(self%ibasin)))
    call message('    to ', trim(dirOut(self%ibasin)))

  end subroutine close

  !------------------------------------------------------------------
  !     NAME
  !         createOutputFile
  !
  !     PURPOSE
  !>        \brief Create and initialize output file
  !>        \details Create output file, write all non-dynamic variables
  !>                 and global attributes for the given basin.
  !>
  !
  !     CALLING SEQUENCE
  !         nc = createOutputFile(ibasin)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)     :: ibasin"      -> basin id
  !>        \param[in] "logical, target :: mask1(:,:)"  -> level1 mask
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         \return type(NcDataset)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  !         modified,
  !                  Oct 2015, Stephan Thober - added actual version of mHM
  function createOutputFile(ibasin) result(nc)

    use mo_common_variables, only : dirOut, level1
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_grid, only : mapCoordinates, geoCoordinates
    use mo_file, only : version
    use mo_julian, only : dec2date

    integer(i4), intent(in) :: ibasin
    type(NcDataset) :: nc
    type(NcDimension) :: dimids1(3)
    type(NcVariable) :: var
    integer(i4) :: day, month, year
    character(1028) :: fname
    character(128) :: unit, date, time, datetime
    real(dp), allocatable :: northing(:), easting(:), lat(:, :), lon(:, :)

    fname = trim(dirOut(ibasin)) // 'mHM_Fluxes_States.nc'
    call mapCoordinates(level1(ibasin), northing, easting)
    call geoCoordinates(level1(ibasin), lat, lon)

    nc = NcDataset(trim(fname), "w")
    dimids1 = (/&
            nc%setDimension("easting", size(easting)), &
                    nc%setDimension("northing", size(northing)), &
                    nc%setDimension("time", 0) &
            /)

    ! time units
    call dec2date(real(evalPer(ibasin)%julStart, dp), dd = day, mm = month, yy = year)
    write(unit, "('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    var = nc%setVariable("time", "i32", (/ dimids1(3) /))
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name", "time")

    ! northing
    var = nc%setVariable("northing", "f64", (/ dimids1(2) /))
    call var%setData(northing)
    call var%setAttribute("units", "m or deg. dec.")
    call var%setAttribute("long_name", "y-coordinate in the given coordinate system")

    ! easting
    var = nc%setVariable("easting", "f64", (/ dimids1(1) /))
    call var%setData(easting)
    call var%setAttribute("units", "m or deg. dec.")
    call var%setAttribute("long_name", "x-coordinate in the given coordinate system")

    ! lon
    var = nc%setVariable("lon", "f64", dimids1(1 : 2))
    call var%setData(lon)
    call var%setAttribute("units", "deg. dec.")
    call var%setAttribute("long_name", "longitude")
    call var%setAttribute("missing_value", nodata_dp)

    ! lat
    var = nc%setVariable("lat", "f64", dimids1(1 : 2))
    call var%setData(lat)
    call var%setAttribute("units", "deg. dec.")
    call var%setAttribute("long_name", "latitude")
    call var%setAttribute("missing_value", nodata_dp)

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
  !     NAME
  !         writeVariableAttributes
  !
  !     PURPOSE
  !>        \brief Write output variable attributes
  !
  !     CALLING SEQUENCE
  !         call writeVariableAttributes(var, long_name, unit)
  !
  !     INTENT(IN)
  !>        \param[in] "type(OutputVariable) :: var"
  !>        \param[in] "character(*)         :: long_name"    -> variable name
  !>        \param[in] "character(*)         :: unit"         -> physical unit
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !        None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine writeVariableAttributes(var, long_name, unit)
    type(OutputVariable), intent(in) :: var
    character(*), intent(in) :: long_name, unit

    call var%nc%setAttribute("long_name", long_name)
    call var%nc%setAttribute("unit", unit)
    call var%nc%setAttribute("scale_factor", 1.0_dp)
    call var%nc%setAttribute("missing_value", nodata_dp)
    call var%nc%setAttribute("coordinates", "lat lon")

  end subroutine writeVariableAttributes


  !------------------------------------------------------------------
  !     NAME
  !         fluxesUnit
  !
  !     PURPOSE
  !>        \brief Generate a unit string
  !>        \details Generate the unit string for the output variable
  !>                 netcdf attribute based on modeling timestep
  !
  !     CALLING SEQUENCE
  !         unit = fluxesUnit(iBasin)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: iBasin"    -> basin id
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         \return character(16)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  function fluxesUnit(ibasin)

    use mo_common_mhm_mrm_variables, only : simPer, timestep, nTstepDay
    use mo_global_variables, only : timeStep_model_outputs

    integer(i4), intent(in) :: ibasin
    character(16) :: fluxesUnit
    real(dp) :: ntsteps

    if (timestep * timestep_model_outputs .eq. 1) then
      fluxesUnit = 'mm h-1'
    else if (timestep_model_outputs > 1) then
      fluxesUnit = 'mm ' // trim(adjustl(num2str(timestep))) // 'h-1'
    else if (timestep_model_outputs .eq. 0) then
      ntsteps = (simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1) * nTstepDay
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
