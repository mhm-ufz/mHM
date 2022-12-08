!> \file mo_nc_output.f90
!> \brief \copybrief mo_nc_output
!> \details \copydetails mo_nc_output

!> \brief Creates NetCDF output for different fluxes and state variables of mHM.
!> \details NetCDF is first initialized and later on variables are put to the NetCDF.
!!
!! Modifications:
!! - David Schaefer       Aug 2015 - major rewrite
!! - Stephan Thober       Oct 2015 - adapted to mRM
!! - O. Rakovec, R. Kumar Nov 2017 - added project description for the netcdf outputs
!! - S. Mueller,          Dec 2022 - unified module for mHM and mRM
!!
!> \authors Matthias Zink
!> \date Apr 2013
!> \ingroup f_common
module mo_nc_output

  use mo_kind, only : i4, dp
  use mo_common_variables, only : project_details, setup_description, simulation_type, &
          Conventions, contact, mHM_details, history, dirOut, iFlag_cordinate_sys, Grid
  use mo_file, only : version
  use mo_common_constants, only : nodata_dp
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable
  use mo_common_mHM_mRM_variables, only: timeStep

  implicit none

  public :: OutputDataset, OutputVariable, writeVariableAttributes, data_dims, data_dtype

  private

  type OutputVariable
    type(NcVariable) :: nc                 !< NcDataset which contains the variable
    logical :: avg = .false.      !< average data before writing
    logical, pointer :: mask(:, :)          !< mask to reconstruct data
    real(dp), allocatable :: data(:)            !< store the data between writes
    integer(i4) :: counter = 0        !< count the number of updateVariable calls

  contains
    procedure, public :: updateVariable
    procedure, public :: writeVariableTimestep

  end type OutputVariable

  ! constructor interface
  interface OutputVariable
    procedure newOutputVariable
  end interface OutputVariable

  type OutputDataset
    integer(i4) :: iDomain      !< domain id
    type(NcDataset) :: nc          !< NcDataset to write
    type(OutputVariable), allocatable :: vars(:)     !< store all created (dynamic) variables
    integer(i4) :: counter !< count written time steps
    integer(i4) :: previous_time !< previous time steps for bounds
    integer(i4) :: time_unit_factor !< possible factor to convert hours to minutes when using center as time reference
    integer(i4) :: outputs_frequence !< write out frequence (-3: yearly, -2: monthly, -1: daily, 0: end of run, >0: after n steps)
    integer(i4) :: time_reference !< time stamp reference (0: begin, 1: center, 2: end of time interval)
    logical :: double_precision !< output precision switch for nc files

  contains
    procedure, public :: writeTimestep
    procedure, public :: close

  end type OutputDataset

  ! constructor interface
  interface OutputDataset
    procedure newOutputDataset
  end interface OutputDataset

contains

  !> \brief Output variable dtype for single of double precision.
  !> \return "f64" or "f32"
  character(3) function data_dtype(double_precision)
    implicit none
    logical, intent(in) :: double_precision !< flag to use double precision
    if ( double_precision ) then
      data_dtype = "f64"
    else
      data_dtype = "f32"
    end if
  end function data_dtype

  !> \brief Output variable dimension names.
  !> \return (X, Y, T) names tuple
  function data_dims()
    implicit none
    character(16), dimension(3) :: data_dims
    if (iFlag_cordinate_sys == 0) then
      data_dims = (/"easting ", "northing", "time    "/) ! X & Y coordinate system
    else
      data_dims = (/"lon ", "lat ", "time"/) ! lat & lon coordinate system
    endif
  end function data_dims

  !> \brief Initialize OutputVariable
  !> \details Modifications:
  !! - David Schaefer Nov 2017 - added NcVariable initialization
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !> \return type(OutputVariable)
  !> \authors David Schaefer
  !> \date June 2015
  function newOutputVariable(nc, name, dtype, dims, ncells, mask, deflate_level, avg) result(out)
    implicit none

    type(NcDataset), intent(in) :: nc !< NcDataset which contains the variable
    character(*), intent(in) :: name !< name of the variable
    character(*), intent(in) :: dtype !< data type of the variable
    character(16), intent(in), dimension(3) :: dims !< dimensions of the variable (by name)
    integer(i4), intent(in) :: ncells !< number of cells in domain
    logical, intent(in), target, dimension(:, :) :: mask !< mask of the variable
    integer(i4), intent(in) :: deflate_level !< deflate level for compression
    logical, intent(in), optional :: avg !< flag to average the data before writing

    type(OutputVariable) :: out

    allocate(out%data(ncells))
    out%nc = nc%setVariable(name, dtype, dims, deflate_level = deflate_level, shuffle = .true.)
    out%data = 0
    out%mask => mask
    if (present(avg)) out%avg = avg
  end function newOutputVariable

  !> \brief Update OutputVariable
  !> \details Add the array given as actual argument
  !> to the derived type's component 'data'
  !!
  !! Modifications:
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !!
  !> \return type(OutputVariable)
  !> \authors David Schaefer
  !> \date June 2015
  subroutine updateVariable(self, data)
    implicit none

    class(OutputVariable), intent(inout) :: self
    real(dp), intent(in), dimension(:) :: data !< data for current time step

    self%data = self%data + data
    self%counter = self%counter + 1

  end subroutine updateVariable

  !> \brief Write timestep to file
  !> \details Write the content of the derived types's component 'data' to file, average if necessary
  !!
  !! Modifications:
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !!
  !> \authors David Schafer
  !> \date June 2015
  subroutine writeVariableTimestep(self, current_time_step)
    implicit none

    class(OutputVariable), intent(inout) :: self

    !> index along the time dimension of the netcdf variable
    integer(i4), intent(in) :: current_time_step

    if (self%avg) then
      self%data = self%data / real(self%counter, dp)
    end if

    call self%nc%setData(unpack(self%data, self%mask, nodata_dp), &
            (/1, 1, current_time_step/))
    self%data = 0
    self%counter = 0

  end subroutine writeVariableTimestep

  !> \brief Initialize OutputDataset
  !> \details Create and initialize the output file. If new a new output
  !! variable needs to be written, this is the first of two
  !! procedures to change (second: updateDataset)
  !!
  !! Modifications:
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !! - Sebastian Mueller Jul 2020 - added output for river temperature
  !!
  !> \return type(OutputDataset)
  !> \authors Matthias Zink
  !> \date Apr 2013
  function newOutputDataset( iDomain, level, file_name, double_precision, outputs_frequence, time_reference ) result(out)
    implicit none

    integer(i4), intent(in) :: iDomain !< domain id
    type(Grid), dimension(:), allocatable, target, intent(in) :: level !< level definitions for all domains
    character(*), intent(in) :: file_name !< long name of the variable
    logical, intent(in) :: double_precision !< mask on desired level
    integer(i4), intent(in) :: outputs_frequence !< write out frequence (-3, -2, -1, 0, >0)
    integer(i4), intent(in) :: time_reference !< time stamp reference (0: begin, 1: center, 2: end of time interval)

    type(OutputDataset) :: out

    out%nc = createOutputFile(iDomain, level, file_name, double_precision, outputs_frequence, time_reference)
    out%iDomain = iDomain

    ! check if we need minutes instead of hours as time unit and set the time unit factor accordingly
    if ( (outputs_frequence > 0) &
         .and. (mod(timestep * outputs_frequence, 2) == 1) &
         .and. (time_reference == 1) &
    ) then
      out%time_unit_factor = 60
    else
      out%time_unit_factor = 1
    end if
    out%outputs_frequence = outputs_frequence
    out%time_reference = time_reference
    out%double_precision = double_precision
    out%previous_time = 0
    out%counter = 0

  end function newOutputDataset

  !> \brief Write all accumulated data.
  !> \details Write all accumulated and potentially averaged data to disk.
  !!
  !! Modifications:
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !!
  !> \authors David Schaefer
  !> \date June 2015
  subroutine writeTimestep(self, current_time_step)
    implicit none

    class(OutputDataset), intent(inout), target :: self

    !> The model timestep to write
    integer(i4), intent(in) :: current_time_step

    integer(i4) :: ii

    type(NcVariable) :: tvar

    self%counter = self%counter + 1

    ! add to time variable
    tvar = self%nc%getVariable("time")
    select case( self%time_reference )
      case(0)
        call tvar%setData(self%previous_time * self%time_unit_factor, (/self%counter/))
      case(1)
        call tvar%setData((self%previous_time + current_time_step) * self%time_unit_factor / 2, (/self%counter/))
      case(2)
        call tvar%setData(current_time_step * self%time_unit_factor, (/self%counter/))
    end select
    ! add bounds (with current time at end)
    tvar = self%nc%getVariable("time_bnds")
    call tvar%setData(self%previous_time * self%time_unit_factor, (/1, self%counter/))
    call tvar%setData(current_time_step * self%time_unit_factor, (/2, self%counter/))
    self%previous_time = current_time_step

    do ii = 1, size(self%vars)
      call self%vars(ii)%writeVariableTimestep(self%counter)
    end do

  end subroutine writeTimestep

  !> \brief Close the file
  !> \details Close the file associated with variable of type(OutputDataset)
  !!
  !! Modifications:
  !! - Stephan Thober Oct  2015 - adapted to mRM
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !!
  !> \authors Rohini Kumar & Stephan Thober
  !> \date August 2013
  subroutine close(self)

    use mo_String_utils, only : num2str
    use mo_message, only : message

    implicit none

    class(OutputDataset) :: self

    call self%nc%close()
    call message('  OUTPUT: saved netCDF file for domain', trim(num2str(self%iDomain)))
    call message('    to ', trim(self%nc%fname))

  end subroutine close

  !> \brief Create and initialize output file for X & Y coordinate system
  !> \details Create output file, write all non-dynamic variables
  !!       and global attributes for the given domain for X & Y coordinate system
  !!
  !! Modifications:
  !! - Stephan Thober  Oct 2015 - adapted to mRM
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !! - Pallav Shrestha Mar 2020 - output file lat and lon are 1d or 2d based on coordinate system
  !!
  !> \return type(NcDataset)
  !> \authors David Schaefer
  !> \date June 2015
  function createOutputFile(iDomain, level, file_name, double_precision, outputs_frequence, time_reference) result(nc)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_grid, only : geoCoordinates, mapCoordinates
    use mo_julian, only : dec2date

    implicit none

    ! -> domain id
    integer(i4), intent(in) :: iDomain !< selected domain
    type(Grid), dimension(:), allocatable, intent(in) :: level !< level definitions for all domains
    character(*), intent(in) :: file_name !< long name of the variable
    logical, intent(in) :: double_precision !< mask on desired level
    integer(i4), intent(in) :: outputs_frequence !< write out frequence (-3, -2, -1, 0, >0)
    integer(i4), intent(in) :: time_reference !< time stamp reference (0: begin, 1: center, 2: end of time interval)

    type(NcDataset) :: nc
    type(NcDimension), dimension(4) :: dimids1
    type(NcVariable) :: var
    integer(i4) :: day, month, year
    character(1028) :: fname
    character(128) :: unit, date, time, datetime
    real(dp), allocatable, dimension(:) :: easting, northing
    real(dp), allocatable, dimension(:, :) :: x_bnds, y_bnds
    real(dp) :: half_step
    real(dp), allocatable, dimension(:) :: lat1d, lon1d    ! 1D lat lon vectors. Used if coordinate system is lat & lon
    real(dp), allocatable, dimension(:, :) :: lat2d, lon2d ! temporary storage of mHM's 2D latlon array.
                                                           ! Used as 2d lat lon arrays if coordinate system is X & Y
    character(3) :: dtype

    dtype = data_dtype(double_precision)

    ! half cell step to calculate cell bounds from center
    half_step = level(iDomain)%cellsize / 2.0_dp

    fname = trim(dirOut(iDomain)) // trim(file_name)
    call geoCoordinates(level(iDomain), lat2d, lon2d)

    nc = NcDataset(trim(fname), "w")

    ! set the horizonal dimensions
    if (iFlag_cordinate_sys == 0) then

      ! X & Y coordinate system; 2D lat lon!
      !============================================================
      call mapCoordinates(level(iDomain), northing, easting)
      allocate(x_bnds(2, size(easting)))
      allocate(y_bnds(2, size(northing)))
      x_bnds(1, :) = easting - half_step
      x_bnds(2, :) = easting + half_step
      y_bnds(1, :) = northing - half_step
      y_bnds(2, :) = northing + half_step

      dimids1 = (/&
        nc%setDimension("easting", size(easting)), &
        nc%setDimension("northing", size(northing)), &
        nc%setDimension("time", 0), &
        nc%setDimension("bnds", 2) &
      /)
      ! easting
      var = nc%setVariable("easting", dtype, (/ dimids1(1) /))
      call var%setData(easting)
      call var%setAttribute("axis", "X")
      call var%setAttribute("units", "m")
      call var%setAttribute("long_name", "x-coordinate in the given coordinate system")
      call var%setAttribute("standard_name", "projection_x_coordinate")
      call var%setAttribute("bounds", "easting_bnds")
      var = nc%setVariable("easting_bnds", dtype, (/ dimids1(4), dimids1(1) /))
      call var%setData(x_bnds)
      ! northing
      var = nc%setVariable("northing", dtype, (/ dimids1(2) /))
      call var%setData(northing)
      call var%setAttribute("axis", "Y")
      call var%setAttribute("units", "m")
      call var%setAttribute("long_name", "y-coordinate in the given coordinate system")
      call var%setAttribute("standard_name", "projection_y_coordinate")
      call var%setAttribute("bounds", "northing_bnds")
      var = nc%setVariable("northing_bnds", dtype, (/ dimids1(4), dimids1(2) /))
      call var%setData(y_bnds)
      ! lon
      var = nc%setVariable("lon", dtype, dimids1(1 : 2))
      call var%setFillValue(nodata_dp)
      call var%setData(lon2d)
      call var%setAttribute("units", "degrees_east")
      call var%setAttribute("long_name", "longitude")
      call var%setAttribute("standard_name", "longitude")
      call var%setAttribute("missing_value", nodata_dp)
      ! lat
      var = nc%setVariable("lat", dtype, dimids1(1 : 2))
      call var%setFillValue(nodata_dp)
      call var%setData(lat2d)
      call var%setAttribute("units", "degrees_north")
      call var%setAttribute("long_name", "latitude")
      call var%setAttribute("standard_name", "latitude")
      call var%setAttribute("missing_value", nodata_dp)

    else

      ! lat & lon coordinate system; 1D lat lon!
      !============================================================
      lon1d = lon2d(:, 1) ! first column info is sufficient
      lat1d = lat2d(1, :) ! first row info is sufficient
      allocate(x_bnds(2, size(lon1d)))
      allocate(y_bnds(2, size(lat1d)))
      ! cellsize is given in degree in case of lat-lon coordinates
      x_bnds(1, :) = lon1d - half_step
      x_bnds(2, :) = lon1d + half_step
      y_bnds(1, :) = lat1d - half_step
      y_bnds(2, :) = lat1d + half_step

      dimids1 = (/&
        nc%setDimension("lon", size(lon1d)), &
        nc%setDimension("lat", size(lat1d)), &
        nc%setDimension("time", 0), &
        nc%setDimension("bnds", 2) &
      /)
      ! lon
      var = nc%setVariable("lon", dtype, (/ dimids1(1) /)) ! sufficient to store lon as vector
      call var%setData(lon1d)
      call var%setAttribute("axis", "X")
      call var%setAttribute("units", "degrees_east")
      call var%setAttribute("long_name", "longitude")
      call var%setAttribute("standard_name", "longitude")
      call var%setAttribute("bounds", "lon_bnds")
      var = nc%setVariable("lon_bnds", dtype, (/ dimids1(4), dimids1(1) /))
      call var%setData(x_bnds)
      ! lat
      var = nc%setVariable("lat", dtype, (/ dimids1(2) /)) ! sufficient to store lat as vector
      call var%setData(lat1d)
      call var%setAttribute("axis", "Y")
      call var%setAttribute("units", "degrees_north")
      call var%setAttribute("long_name", "latitude")
      call var%setAttribute("standard_name", "latitude")
      call var%setAttribute("bounds", "lat_bnds")
      var = nc%setVariable("lat_bnds", dtype, (/ dimids1(4), dimids1(2) /))
      call var%setData(y_bnds)

    endif

    ! set record dimension
    ! time units
    call dec2date(real(evalPer(iDomain)%julStart, dp), dd = day, mm = month, yy = year)

    ! check if we need minutes instead of hours as time unit
    if ( (outputs_frequence > 0) &                           ! only for output after n steps
         .and. (mod(timestep * outputs_frequence, 2) == 1) & ! only for uneven hours
         .and. (time_reference == 1) &                       ! only when center of time span is used for time stamp
    ) then
      write(unit, "('minutes since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day
    else
      write(unit, "('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day
    end if

    ! time
    var = nc%setVariable("time", "i32", (/ dimids1(3) /))
    call var%setAttribute("axis", "T")
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name", "time")
    call var%setAttribute("standard_name", "time")
    call var%setAttribute("bounds", "time_bnds")
    var = nc%setVariable("time_bnds", "i32", (/ dimids1(4), dimids1(3) /))

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

  !> \brief Write output variable attributes
  !> \details Modifications:
  !! - Robert Schweppe Jun 2018 - refactoring and reformatting
  !!
  !> \authors David Schaefer
  !> \date June 2015
  subroutine writeVariableAttributes(var, long_name, unit)
    implicit none

    type(OutputVariable), intent(inout) :: var !< NetCDF variable
    character(*), intent(in) :: long_name !< long name of the variable
    character(*), intent(in) :: unit !< unit of the variable

    call var%nc%setFillValue(nodata_dp)
    call var%nc%setAttribute("long_name", long_name)
    call var%nc%setAttribute("units", unit)
    call var%nc%setAttribute("scale_factor", 1.0_dp)
    call var%nc%setAttribute("missing_value", nodata_dp)
    call var%nc%setAttribute("coordinates", "lat lon")

  end subroutine writeVariableAttributes

end module mo_nc_output
