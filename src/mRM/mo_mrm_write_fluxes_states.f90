!>       \file mo_mrm_write_fluxes_states.f90

!>       \brief Creates NetCDF output for different fluxes and state variables of mHM.

!>       \details NetCDF is first initialized and later on variables are put to the NetCDF.

!>       \authors Matthias Zink

!>       \date Apr 2013

! Modifications:
! David Schaefer       Aug 2015 - major rewrite
! Stephan Thober       Oct 2015 - adapted to mRM
! O. Rakovec, R. Kumar Nov 2017 - added project description for the netcdf outputs 

module mo_mrm_write_fluxes_states

  use mo_kind, only : i4, dp
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
    out%nc = nc%setVariable(name, dtype, dims, deflate_level = 1, shuffle = .true.)
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
  !>       \param[in] "integer(i4) :: iDomain"            -> domain id
  !>       \param[in] "logical, dimension(:, :) :: mask"
  !>       \param[in] "integer(i4) :: nCells"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function newOutputDataset(iDomain, mask, nCells) result(out)

    use mo_mrm_global_variables, only : outputFlxState_mrm

    implicit none

    ! -> domain id
    integer(i4), intent(in) :: iDomain

    logical, intent(in), pointer, dimension(:, :) :: mask

    integer(i4), intent(in) :: nCells

    type(OutputDataset) :: out

    integer(i4) :: ii

    character(3) :: dtype

    character(16), dimension(3) :: dims1

    type(NcDataset) :: nc

    type(OutputVariable), dimension(size(outputFlxState_mrm)) :: tmpvars


    dtype = "f64"
    dims1 = (/"easting ", "northing", "time    "/)
    nc = createOutputFile(iDomain)

    ii = 0

    if (outputFlxState_mrm(1)) then
      ii = ii + 1
      tmpvars(ii) = OutputVariable(&
              nc, "Qrouted", dtype, dims1, nCells, mask, .true.)
      call writeVariableAttributes(tmpvars(ii), "routed streamflow", "m3 s-1")
    end if

    ! out = OutputDataset(iDomain, nc, tmpvars(1 : ii))
    allocate(out%vars(ii))
    out%vars = tmpvars(1:ii)
    out%nc = nc
    out%iDomain = iDomain
    ! print*, 'Finished OutputDatasetInit'

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
  !>       \param[in] "integer(i4) :: sidx, eidx"          - start index of the domain related data in L1_* arguments
  !>       \param[in] "integer(i4) :: sidx, eidx"          - end index of the domain related data in L1_* arguments
  !>       \param[in] "real(dp), dimension(:) :: L11_Qmod"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! L. Samaniego et al. Dec  2013 - nullify pointer Matthias Zink,        Feb. 2014 
  !                              - added aditional output: pet V. Prykhodk, J. Mai,  Nov. 2014 
  !                              - adding new variable infilSoil 
  !                              - case 16 David Schaefer      , Jun. 2015 
  !                              - major rewrite
  ! Stephan Thober      Oct  2015 - adapted to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine updateDataset(self, sidx, eidx, L11_Qmod)

    use mo_mrm_global_variables, only : outputFlxState_mrm

    implicit none

    class(OutputDataset), intent(inout), target :: self

    ! - end index of the domain related data in L1_* arguments
    integer(i4), intent(in) :: sidx, eidx

    real(dp), intent(in), dimension(:) :: L11_Qmod

    type(OutputVariable), pointer, dimension(:) :: vars

    integer(i4) :: ii


    ii = 0
    vars => self%vars

    if (outputFlxState_mrm(1)) then
      ii = ii + 1
#ifdef pgiFortran
      call updateVariable(vars(ii), L11_Qmod(sidx : eidx))
#else
      call vars(ii)%updateVariable(L11_Qmod(sidx : eidx))
#endif
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
  ! Stephan Thober Oct  2015 - adapted to mRM
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
  !>       \brief Create and initialize output file

  !>       \details Create output file, write all non-dynamic variables
  !>       and global attributes for the given domain.

  !>       \return type(NcDataset)

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" -> domain id

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date June 2015

  ! Modifications:
  ! Stephan Thober  Oct 2015 - adapted to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function createOutputFile(iDomain) result(nc)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : dirOut
    use mo_grid, only : geoCoordinates, mapCoordinates
    use mo_julian, only : dec2date
    use mo_mrm_file, only : file_mrm_output, version
    use mo_mrm_global_variables, only : level11

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

    real(dp), allocatable, dimension(:,:) :: lat, lon


    fname = trim(dirOut(iDomain)) // trim(file_mrm_output)
    call mapCoordinates(level11(iDomain), northing, easting)
    call geoCoordinates(level11(iDomain), lat, lon)

    nc = NcDataset(trim(fname), "w")
    dimids1 = (/&
            nc%setDimension("easting", size(easting)), &
                    nc%setDimension("northing", size(northing)), &
                    nc%setDimension("time", 0) &
            /)

    ! time units
    call dec2date(real(evalPer(iDomain)%julStart, dp), dd = day, mm = month, yy = year)
    write(unit, "('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    var = nc%setVariable("time", "i32", (/ dimids1(3) /))
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name", "time")

    ! northing
    var = nc%setVariable("northing", "f64", (/ dimids1(2) /))
    call var%setData(northing)
    call var%setAttribute("units", "m")
    call var%setAttribute("long_name", "y-coordinate in the given coordinate system")

    ! easting
    var = nc%setVariable("easting", "f64", (/ dimids1(1) /))
    call var%setData(easting)
    call var%setAttribute("units", "m")
    call var%setAttribute("long_name", "x-coordinate in the given coordinate system")

    ! lon
    var = nc%setVariable("lon", "f64", dimids1(1 : 2))
    call var%setData(lon)
    call var%setAttribute("units", "degrees_east")
    call var%setAttribute("long_name", "longitude")
    call var%setAttribute("missing_value", nodata_dp)

    ! lat
    var = nc%setVariable("lat", "f64", dimids1(1 : 2))
    call var%setData(lat)
    call var%setAttribute("units", "degrees_north")
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
    call nc%setAttribute("mRM_details", trim(mHM_details) // ", release mRMv" // trim(version))
    call nc%setAttribute("history", trim(datetime) // ", " // history)
    call nc%setAttribute("title", "mRMv"//trim(version)//" "//trim(simulation_type)//" outputs")
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

    type(OutputVariable), intent(in) :: var

    ! -> physical unit
    character(*), intent(in) :: long_name, unit


    call var%nc%setAttribute("long_name", long_name)
    call var%nc%setAttribute("unit", unit)
    call var%nc%setAttribute("scale_factor", 1.0_dp)
    call var%nc%setAttribute("missing_value", nodata_dp)
    call var%nc%setAttribute("coordinates", "lat lon")

  end subroutine writeVariableAttributes

end module mo_mrm_write_fluxes_states
