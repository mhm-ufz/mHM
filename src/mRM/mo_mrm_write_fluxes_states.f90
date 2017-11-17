!> \file mo_write_fluxes_states.f90

!> \brief Creates NetCDF output for different fluxes and state variables of mHM.

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!  HISTORY
!>     \authors Matthias Zink
!>     \date Apr 2013
!      Modified:
!          David Schaefer, Aug 2015 - major rewrite
!          Stephan Thober, Oct 2015 - adapted to mRM
!    O. Rakovec, R. Kumar, Nov 2017 - added project description for the netcdf outputs
!

module mo_mrm_write_fluxes_states

  use mo_kind                 , only: i4, dp
  use mo_string_utils         , only: num2str
  use mo_mrm_global_variables , only: gridGeoRef, project_details, setup_description, simulation_type, &
         Conventions, contact, mHM_details, history
  use mo_mrm_constants        , only: nodata_dp
  use mo_netcdf               , only: NcDataset, NcDimension, NcVariable

  implicit none

  type OutputVariable
     type(NcVariable)      :: nc                 !> NcDataset which contains the variable
     logical               :: avg = .false.      !> average data before writing
     logical,  pointer     :: mask(:,:)          !> mask to reconstruct data
     real(dp), allocatable :: data(:)            !> store the data between writes
     integer(i4)           :: counter = 0        !> count the number of updateVariable calls
     
   contains
     procedure, public :: updateVariable
     procedure, public :: writeVariableTimestep

  end type OutputVariable

  ! constructor interface
  interface OutputVariable
     procedure newOutputVariable
  end interface OutputVariable

  type OutputDataset
     integer(i4)                       :: ibasin      !> basin id
     type(NcDataset)                   :: nc          !> NcDataset to write
     type(OutputVariable), allocatable :: vars(:)     !> store all created (dynamic) variables
     integer(i4)                       :: counter = 0 !> count written time steps

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
  !     NAME
  !         newOutputVariable
  !
  !     PURPOSE
  !>        \brief Initialize OutputVariable
  !
  !     CALLING SEQUENCE
  !         var = OutputVariable(nc, ncells, mask, avg)
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
  !
  !     Modified
  !         David Schaefer, Nov 2017, added NcVariable initialization
  function newOutputVariable(nc, name, dtype, dims, ncells, mask, avg) result(out)

    type(NcDataset), intent(in)   :: nc
    character(*), intent(in)      :: name
    character(*), intent(in)      :: dtype
    character(16), intent(in)     :: dims(3)
    integer(i4), intent(in)       :: ncells
    logical, intent(in), target   :: mask(:,:)
    logical, intent(in), optional :: avg
    type(OutputVariable)          :: out

    allocate(out%data(ncells))
    out%nc   =  nc%setVariable(name, dtype, dims, deflate_level=1, shuffle=.true.)
    out%data =  0
    out%mask => mask
    if (present(avg)) out%avg = avg
  end function newOutputVariable

  ! function newOutputVariable(nc, ncells, mask, avg) result(out)
  !   type(NcVariable), intent(in) :: nc
  !   integer(i4), intent(in)      :: ncells
  !   logical, target              :: mask(:,:)
  !   logical, optional            :: avg
  !   type(OutputVariable)         :: out

  !   allocate(out%data(ncells))
  !   out%nc   =  nc
  !   out%mask => mask
  !   out%data =  0
  !   if (present(avg)) out%avg = avg
  ! end function newOutputVariable

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
    real(dp)             , intent(in)    :: data(:)

    self%data    = self%data + data
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
    integer(i4)          , intent(in)    :: timestep

    if (self%avg) then
       self%data = self%data / real(self%counter, dp)
    end if

    call self%nc%setData(unpack(self%data, self%mask, nodata_dp), &
         (/1,1,timestep/))
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
  !             Stephan Thober      , Oct. 2015 - adapted to mRM
  !             David Schaefer      , Nov. 2016 - moved NcVariable initialization to newOutputVariable
  function newOutputDataset(ibasin, mask) result(out)

    use mo_mrm_global_variables,  only : outputFlxState_mrm
    use mo_mrm_tools,         only : get_basin_info_mrm

    integer(i4), intent(in) :: ibasin
    logical, target         :: mask(:,:)
    type(OutputDataset)     :: out 
    ! local
    integer(i4)             :: ii, ncols, nrows, ncells
    character(3)            :: dtype
    character(16)           :: dims1(3)
    type(NcDataset)         :: nc
    type(OutputVariable)    :: tmpvars(size(outputFlxState_mrm))

    call get_basin_info_mrm(ibasin, 11, ncols, nrows, ncells=ncells)

    dtype = "f64"
    dims1 = (/"easting ", "northing","time    "/)
    nc    = createOutputFile(ibasin)

    ii = 0

    if (outputFlxState_mrm(1)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariable(&
            nc, "Qrouted", dtype, dims1, ncells, mask, .true.)
       call writeVariableAttributes(tmpvars(ii), "routed streamflow", "m3 s-1")
    end if

    out = OutputDataset(ibasin, nc, tmpvars(1:ii))

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
  !             L11_qMod )
  !             
  !
  !     INTENT(IN)
  !>             \param[in] "sidx"        -> start index of the basin related data in L1_* arguments
  !>             \param[in] "eidx"        -> end index of the basin related data in L1_* arguments
  !>             \param[in] "L11_qMod"
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
  !             Stephan Thober      , Oct  2015 - adapted to mRM
  subroutine updateDataset(self, sidx, eidx, &
       L11_Qmod )
    
    use mo_mrm_global_variables,  only : outputFlxState_mrm

    class(OutputDataset), intent(inout), target :: self
    integer(i4),          intent(in)            :: sidx, eidx
    ! fluxes,
    real(dp),             intent(in)            :: L11_Qmod(:)

    ! local
    type(OutputVariable), pointer :: vars(:)
    integer(i4)       :: ii

    ii = 0
    vars  => self%vars

    if (outputFlxState_mrm(1)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariable(vars(ii), L11_Qmod(sidx:eidx))
#else
       call vars(ii)%updateVariable(L11_Qmod(sidx:eidx))
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
    integer(i4),          intent(in)    :: timestep
    integer(i4)                         :: ii
    type(NcVariable)                    :: tvar

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
  !             Stephan Thober, Oct  2015 - adapted to mRM
  subroutine close(self)

    use mo_String_utils,         only: num2str
    use mo_message,              only: message
    use mo_mrm_global_variables, only: dirOut 

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
  !>        \param[in] "logical, target :: mask1(:,:)"  -> level11 mask
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
  !             Stephan Thober, Oct  2015 - adapted to mRM
  function createOutputFile(ibasin) result(nc)

    use mo_mrm_global_variables,  only: dirOut, evalPer, level11
    use mo_mrm_file,              only: version, file_mrm_output
    use mo_julian,                only: dec2date

    integer(i4), intent(in) :: ibasin
    type(NcDataset)         :: nc
    type(NcDimension)       :: dimids1(3)
    type(NcVariable)        :: var
    integer(i4)             :: day, month, year
    character(1028)         :: fname
    character(128)          :: unit, date, time, datetime
    real(dp), allocatable   :: northing(:), easting(:), lat(:,:), lon(:,:)

    fname = trim(dirOut(ibasin)) // trim(file_mrm_output)
    call mapCoordinates(ibasin, level11, northing, easting)
    call geoCoordinates(ibasin, level11, lat, lon)

    nc = NcDataset(trim(fname),"w")
    dimids1 = (/&
         nc%setDimension("easting",  size(easting)), &
         nc%setDimension("northing", size(northing)), &
         nc%setDimension("time",     0) &
         /)

    ! time units
    call dec2date(real(evalPer(ibasin)%julStart, dp), dd=day, mm=month, yy=year)
    write(unit,"('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    var = nc%setVariable("time", "i32", (/ dimids1(3) /))
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name","time")

    ! northing
    var = nc%setVariable("northing", "f64", (/ dimids1(2) /))
    call var%setData(northing)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","y-coordinate in cartesian coordinates GK4")

    ! easting
    var = nc%setVariable("easting", "f64", (/ dimids1(1) /))
    call var%setData(easting)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","x-coordinate in cartesian coordinates GK4")

    ! lon
    var = nc%setVariable("lon","f64",dimids1(1:2))
    call var%setData(lon)
    call var%setAttribute("units","degerees_east")
    call var%setAttribute("long_name","longitude")
    call var%setAttribute("missing_value","-9999.")

    ! lat
    var = nc%setVariable("lat","f64",dimids1(1:2))
    call var%setData(lat)
    call var%setAttribute("units","degerees_north")
    call var%setAttribute("long_name","latitude")
    call var%setAttribute("missing_value","-9999.")

    ! global attributes
    call date_and_time(date=date, time=time)
    write(datetime,"(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1:4), &
         date(5:6), date(7:8), time(1:2), time(3:4), time(5:6)

    call nc%setAttribute("project", project_details )
    call nc%setAttribute("setup_description",setup_description)
    call nc%setAttribute("simulation_type",simulation_type)
    call nc%setAttribute("Conventions",Conventions)
    call nc%setAttribute("contact",contact)
    call nc%setAttribute("mHM_details",trim(mHM_details)//", release mRMv"//trim(version))
    call nc%setAttribute("history", trim(datetime)//", "//history)

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
    character(*)        , intent(in) :: long_name, unit

    call var%nc%setAttribute("_FillValue",nodata_dp)
    call var%nc%setAttribute("long_name",long_name)
    call var%nc%setAttribute("unit",unit)
    call var%nc%setAttribute("scale_factor",1.0_dp)
    call var%nc%setAttribute("missing_value", "-9999.")
    call var%nc%setAttribute("coordinates","lat lon")

  end subroutine writeVariableAttributes
  

  !------------------------------------------------------------------
  !     NAME
  !         mapCoordinates
  !
  !     PURPOSE
  !>        \brief Generate map coordinates
  !>        \details Generate map coordinate arrays for given basin and level
  !
  !     CALLING SEQUENCE
  !         call mapCoordinates(ibasin, level, y, x)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)      :: iBasin" -> basin number
  !>        \param[in] "type(geoGridRef) :: level"  -> grid reference
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(:)  :: y(:)"          -> y-coordinates
  !>        \param[out] "real(dp) :: x(:)"          -> x-coorindates
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
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             Stephan Thober, Nov 2013 - removed fproj dependency
  !             David Schaefer, Jun 2015 - refactored the former subroutine CoordSystem
  subroutine mapCoordinates(ibasin, level, y, x)

    implicit none

    integer(i4),      intent(in)               :: iBasin
    type(gridGeoRef), intent(in)               :: level
    real(dp),         intent(out), allocatable :: x(:), y(:)
    integer(i4)                                :: ii, ncols, nrows
    real(dp)                                   :: cellsize

    cellsize = level%cellsize(ibasin)
    nrows    = level%nrows(ibasin)
    ncols    = level%ncols(ibasin)

    allocate(x(nrows), y(ncols))

    x(1) =  level%xllcorner(ibasin) + 0.5_dp * cellsize
    do ii = 2, nrows
       x(ii)   =  x(ii-1) + cellsize
    end do

    ! inverse for Panoply, ncview display
    y(ncols) =  level%yllcorner(ibasin) + 0.5_dp * cellsize
    do ii = ncols-1,1,-1
       y(ii)   =  y(ii+1) + cellsize
    end do

  end subroutine mapCoordinates

  !------------------------------------------------------------------
  !     NAME
  !         geoCoordinates
  !
  !     PURPOSE
  !>        \brief Generate geographic coordinates
  !>        \details Generate geographic coordinate arrays for given basin and level
  !
  !     CALLING SEQUENCE
  !         call mapCoordinates(ibasin, level, y, x)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)      :: iBasin"    -> basin number
  !>        \param[in] "type(gridGeoRef) :: level"     -> grid reference 
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: lat(:,:)"         -> lat-coordinates
  !>        \param[out] "real(dp) :: lon(:,:)"         -> lon-coorindates
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
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             Stephan Thober, Nov 2013 - removed fproj dependency
  !             David Schaefer, Jun 2015 - refactored the former subroutine CoordSystem
  !             Stephan Thober, Sep 2015 - using mask to unpack coordinates
  !             Stephan Thober, Oct 2015 - writing full lat/lon again
  !             Stephan Thober, Oct 2015 - adapted to mRM
  subroutine geoCoordinates(ibasin, level, lat, lon)

    use mo_mrm_global_variables, only : L11_rect_latitude, L11_rect_longitude

    implicit none

    integer(i4),      intent(in)               :: iBasin
    type(gridGeoRef), intent(in)               :: level
    real(dp),         intent(out), allocatable :: lat(:,:), lon(:,:)
    integer(i4)                                :: ncols, nrows
    integer(i4)                                :: ii, pos 
    
    nrows = level%nrows(ibasin) 
    ncols = level%ncols(ibasin)
      
    pos = 1 
    if ( ibasin .gt. 1 ) then 
       do ii = 1, ibasin -1 
          pos = pos + level%ncols(ii) * level%nrows(ii) 
       end do
    end if

    lat = reshape(L11_rect_latitude(pos:pos+nrows*ncols-1),  (/nrows, ncols/)) 
    lon = reshape(L11_rect_longitude(pos:pos+nrows*ncols-1), (/nrows, ncols/)) 
    
  end subroutine geoCoordinates

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

    use mo_mrm_global_variables,  only : simPer, NTSTEPDAY, timestep, timeStep_model_outputs_mrm

    integer(i4), intent(in) :: ibasin
    character(16)           :: fluxesUnit
    real(dp)                :: ntsteps

    if ( timestep*timestep_model_outputs_mrm .eq. 1 ) then
       fluxesUnit = 'mm h-1'
    else if (timestep_model_outputs_mrm > 1) then
       fluxesUnit = 'mm '//trim(adjustl(num2str(timestep)))//'h-1'
    else if (timestep_model_outputs_mrm .eq. 0) then
       ntsteps = ( simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1 ) * NTSTEPDAY
       fluxesUnit = 'mm '//trim(adjustl(num2str(nint(ntsteps))))//'h-1'
    else if (timestep_model_outputs_mrm .eq. -1) then
       fluxesUnit = 'mm d-1'
    else if (timestep_model_outputs_mrm .eq. -2) then
       fluxesUnit = 'mm month-1'
    else if (timestep_model_outputs_mrm .eq. -3) then
       fluxesUnit = 'mm a-1'
    else
       fluxesUnit = ''
    end if

  end function fluxesUnit

end module mo_mrm_write_fluxes_states
