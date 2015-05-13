
module mo_write_fluxes_states

  ! mHM writing module

  use mo_kind, ONLY: i4, dp
  use mo_string_utils, only : num2str
  use mo_global_variables,  only : gridGeoRef
  use mo_mhm_constants, only: nodata_dp

  use mo_netcdf, only : NcDataset, NcDimension, NcVariable

  implicit none

  type OutputVariable
     type(NcVariable)                           :: ncvar
     real(dp), dimension(:), pointer            :: data       => null()       ! point to data
     real(dp), dimension(:), pointer            :: multiplier => null()       ! point to data
     real(dp), dimension(:), allocatable        :: accum                      ! the data accumulator
     integer(i4)                                :: counter = 0
     logical                                    :: avg = .false.

   contains
     
     procedure, private :: updateVariable
     
  end type OutputVariable
    
  interface OutputVariable
     procedure initOutputVariable
  end interface OutputVariable

  type OutputDataset
     type(NcDataset)                               :: nc
     integer(i4)                                   :: counter = 1
     type(OutputVariable),dimension(:),allocatable :: vars
     logical,dimension(:,:),allocatable            :: mask1,mask11

   contains
     
     procedure, public :: updateDataset
     procedure, public :: writeTimestep
     procedure, public :: close
   
  end type OutputDataset

  interface OutputDataset
     procedure initOutputDataset
  end interface OutputDataset
  
  private
  public :: initOutput
  public :: OutputVariable
  public :: OutputDataset
  
CONTAINS
  
  function initOutputVariable(ncvar, data, multip, avg) result(out)
    type(NcVariable), intent(in)                                  :: ncvar
    real(dp), dimension(:), target, intent(in)                    :: data
    real(dp), dimension(size(data)), target, intent(in), optional :: multip
    logical, intent(in),optional                                  :: avg
    type(OutputVariable)                                          :: out
    
    allocate(out%accum(size(data)))
    out%ncvar = ncvar
    out%data  => data
    out%accum = 0
    if (present(multip)) then
       out%multiplier => multip
    end if
    if (present(avg)) then
       out%avg = avg
    end if    
  end function initOutputVariable
  
  type(OutputDataset) function initOutputDataset(nc,vars,mask1)
    type(NcDataset), intent(in)                  :: nc
    type(OutputVariable),dimension(:),intent(in) :: vars
    logical,dimension(:,:),intent(in)            :: mask1 !,mask11

    initOutputDataset%nc = nc
    allocate(initOutputDataset%vars(size(vars)),source=vars)
    allocate(initOutputDataset%mask1(size(mask1,1),size(mask1,2)),&
         source=mask1)
    ! allocate(initOutputDataset%mask11(size(mask11,1),size(mask11,2)),&
    !      source=mask11)
  end function initOutputDataset
  
  subroutine updateVariable(self)
    class(OutputVariable), intent(inout) :: self
    if (associated(self%multiplier)) then
       self%accum = self%accum + self%data * self%multiplier
    else
       self%accum = self%accum + self%data
    end if
    self%counter = self%counter + 1
  end subroutine updateVariable

  subroutine updateDataset(self)

    class(OutputDataset), intent(inout) :: self
    integer(i4)                         :: ii

    do ii=1,size(self%vars)
       call self%vars(ii)%updateVariable()
    end do
    
  end subroutine updateDataset

  subroutine writeTimestep(self,timestep)
    class(OutputDataset),intent(inout) :: self
    integer(i4), intent(in)         :: timestep
    type(OutputVariable) :: var
    type(NcVariable)     :: tvar
    integer(i4)          :: ii
    real(dp), dimension(:,:), allocatable :: tmpdata1

    tvar = self%nc%getVariable("time")
    call tvar%putDataDirect(timestep, (/self%counter/))
    
    do ii = 1,size(self%vars)
       var = self%vars(ii)
       if (var%avg) then
          var%accum = var%accum / real(var%counter, dp)
       end if
       tmpdata1 = unpack(var%accum, self%mask1, nodata_dp)
       call var%ncvar%putDataDirect(tmpdata1,(/1,1,self%counter/))
       self%vars(ii)%accum = 0
       self%vars(ii)%counter = 0
    end do
    
    self%counter = self%counter + 1

  end subroutine writeTimestep
  
  subroutine close(self)
    class(OutputDataset) :: self
    call self%nc%close()
  end subroutine close
   
  function initOutputFile(ibasin)

    use mo_global_variables,  only: dirOut, evalPer, level1
    use mo_julian,            only: dec2date
    
    integer(i4),intent(in)                          :: ibasin
    type(NcDataset)                                 :: initOutputFile, nc
    type(NcDimension),dimension(3)                  :: dimids1
    type(NcVariable)                                :: var
    character(128)                                  :: fname, unit, date, time, datetime
    integer(i4)                                     :: day,month,year
    real(dp),dimension(:),allocatable               :: northing,easting
    real(dp),dimension(:,:),allocatable             :: lat,lon
    
    fname = trim(dirOut(ibasin)) // 'mHM_Fluxes_States.nc'
    call mapCoordinates(ibasin,level1,northing,easting)
    call geoCoordinates(ibasin,level1,lat,lon)
    
    nc = NcDataset(trim(fname),"w")
    
    dimids1 = (/&
         nc%createDimension("time",     0), &
         nc%createDimension("northing", size(northing)), &
         nc%createDimension("easting",  size(easting)) &
         /)

    ! time units
    call dec2date(real(evalPer(ibasin)%julStart, dp), dd=day, mm=month, yy=year)
    write(unit,"('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time 
    var = nc%createVariable("time", "i32", (/ dimids1(1) /))
    call var%createAttribute("units", unit)
    call var%createAttribute("long_name","time")
    
    ! northing
    var = nc%createVariable("northing", "f64", (/ dimids1(2) /))
    call var%putDataDirect(northing)
    call var%createAttribute("units","m")
    call var%createAttribute("long_name","y-coordinate in cartesian coordinates GK4")
    
    ! easting
    var = nc%createVariable("easting", "f64", (/ dimids1(3) /))
    call var%putDataDirect(easting)
    call var%createAttribute("units","m")
    call var%createAttribute("long_name","x-coordinate in cartesian coordinates GK4")

    ! lon
    var = nc%createVariable("lon","f64",dimids1(2:3))
    call var%putDataDirect(lon)
    call var%createAttribute("units","degerees_east")
    call var%createAttribute("long_name","longitude")

    !lat
    var = nc%createVariable("lat","f64",dimids1(2:3))
    call var%putDataDirect(lat)
    call var%createAttribute("units","degerees_north")
    call var%createAttribute("long_name","latitude")
        
    ! global attributes
    call date_and_time(date=date, time=time)
    write(datetime,"(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1:4), &
         date(5:6), date(7:8), time(1:2), time(3:4), time(5:6)

    call nc%createAttribute("title","mHMv5 simulation outputs")
    call nc%createAttribute("creation_date",datetime)
    call nc%createAttribute("institution",&
         "Helmholtz Centre for Environmental Research - UFZ, "// &
         "Department Computational Hydrosystems, Stochastic Hydrology Group")
    
    initOutputFile = nc
    
  end function initOutputFile


  subroutine writeVariableAttributes(var,long_name,unit)
    type(NcVariable), intent(in) :: var
    character(*),intent(in)      :: long_name, unit
    
    call var%createAttribute("_FillValue",nodata_dp)    
    call var%createAttribute("long_name",long_name)
    call var%createAttribute("unit",unit)       
    call var%createAttribute("scale_factor",1.0_dp)
    call var%createAttribute("missing_value", "-9999.")
    call var%createAttribute("coordinates","lat lon")   
    
  end subroutine writeVariableAttributes

  function fluxesUnit(ibasin)

    use mo_global_variables,  only : simPer, NTSTEPDAY, timestep, timeStep_model_outputs
    
    integer(i4),intent(in) :: ibasin
    character(16)           :: fluxesUnit
    real(dp)               :: ntsteps
    

    if ( timestep*timestep_model_outputs .eq. 1 ) then
       fluxesUnit = 'mm h-1'
    else if (timestep_model_outputs > 1) then
       fluxesUnit = 'mm '//trim(adjustl(num2str(timestep)))//'h-1'
    else if (timestep_model_outputs .eq. 0) then
       ntsteps = ( simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1 ) * NTSTEPDAY
       fluxesUnit = 'mm '//trim(adjustl(num2str(nint(ntsteps))))//'h-1'
    else if (timestep_model_outputs .eq. -1) then
       fluxesUnit = 'mm d-1'
    else if (timestep_model_outputs .eq. -2) then
       fluxesUnit = 'mm month-1'
    else if (timestep_model_outputs .eq. -3) then
       fluxesUnit = 'mm a-1'
    else
       fluxesUnit = ''
    endif
    
  end function fluxesUnit
  
  function initOutput(&
       ibasin             , & ! basin number
       ! states
       L1_inter           , & ! Interception
       L1_snowPack        , & ! Snowpack
       L1_soilMoist       , & ! Soil moisture of each horizon
       L1_soilMoistVol    , &
       L1_soilMoistVolAvg , &
       L1_sealSTW         , &
       L1_unsatSTW        , &
       L1_satSTW          , &
       L1_neutrons        , &
       ! fluxes
       L1_pet          , &    ! potential evapotranspiration (PET)
       L1_aet          , &    ! actual ET 
       L1_total_runoff , &    ! Generated runoff
       L1_runoffSeal   , &    ! Direct runoff from impervious areas
       L1_fastRunoff   , &    ! Fast runoff component
       L1_slowRunoff   , &    ! Slow runoff component
       L1_baseflow     , &    ! Baseflow
       L1_percol       , &    ! Percolation
       L1_infilSoil    , &    ! Infiltration
       ! helper
       L1_fSealed      , &
       L1_fNotSealed    &                   
       )
    
    use mo_global_variables,  only : outputFlxState, nSoilHorizons_mHM
    use mo_init_states,       only : get_basin_info
    
    integer(i4),              intent(in)            :: ibasin
    ! states
    real(dp), dimension(:),   intent(in),    target :: L1_inter          ! Interception
    real(dp), dimension(:),   intent(in),    target :: L1_snowPack       ! Snowpack
    real(dp), dimension(:,:), intent(in),    target :: L1_soilMoist      ! Soil moisture of each horizon                   
    real(dp), dimension(:,:), intent(in),    target :: L1_soilMoistVol
    real(dp), dimension(:),   intent(in),    target :: L1_soilMoistVolAvg
    real(dp), dimension(:),   intent(in),    target :: L1_sealSTW      ! Retention storage of impervious areas
    real(dp), dimension(:),   intent(in),    target :: L1_unsatSTW     ! Upper soil storage
    real(dp), dimension(:),   intent(in),    target :: L1_satSTW       ! Groundwater storage
    real(dp), dimension(:),   intent(in),    target :: L1_neutrons     ! ground albedo neutrons
    ! fluxes
    real(dp), dimension(:),   intent(inout), target :: L1_pet          ! potential evapotranspiration (PET)
    real(dp), dimension(:),   intent(inout), target :: L1_aet
    real(dp), dimension(:),   intent(inout), target :: L1_total_runoff ! Generated runoff
    real(dp), dimension(:),   intent(inout), target :: L1_runoffSeal   ! Direct runoff from impervious areas
    real(dp), dimension(:),   intent(inout), target :: L1_fastRunoff   ! Fast runoff component
    real(dp), dimension(:),   intent(inout), target :: L1_slowRunoff   ! Slow runoff component
    real(dp), dimension(:),   intent(inout), target :: L1_baseflow     ! Baseflow
    real(dp), dimension(:),   intent(inout), target :: L1_percol       ! Percolation
    real(dp), dimension(:,:), intent(inout), target :: L1_infilSoil    ! Infiltration
    ! helper
    real(dp), dimension(:),   intent(inout), target :: L1_fSealed      ! 
    real(dp), dimension(:),   intent(inout), target :: L1_fNotSealed      ! 
    ! output
    type(OutputDataset)                             :: initOutput

    ! local
    type(NcDataset)                                                        :: nc
    character(8),dimension(3)                                              :: dims1 
    type(NcVariable)                                                       :: var
    character(3)                                                           :: dtype
    character(16)                                                          :: unit
    type(OutputVariable),dimension(size(outputFlxState)*nSoilHorizons_mHM) :: tmpvars
    integer(i4)                                                            :: ii, nn, ncols, nrows
    logical,dimension(:,:),allocatable                                     :: mask1 

    dtype = "f64"
    dims1 = (/"time    ","northing","easting "/)
    nc = initOutputFile(ibasin)    

    call get_basin_info (ibasin,  1,  ncols, nrows, mask=mask1)
    ! call get_basin_info (ibasin, 11,  ncols, nrows, mask=mask11)

    ii = 0
    ! interception
    if (outputFlxState(1)) then
       ii = ii + 1
       var = nc%createVariable("interception",dtype,dims1)
       call writeVariableAttributes(var,"canopy interception storage", "mm")
       tmpvars(ii) = OutputVariable(var,L1_inter,avg=.true.)
    end if

    if (outputFlxState(2)) then       
       ii = ii + 1
       var = nc%createVariable("snowpack",dtype,dims1)
       call writeVariableAttributes(var,"depth of snowpack", "mm")
       tmpvars(ii) =  OutputVariable(var,L1_snowPack,avg=.true.)
    end if

    if (outputFlxState(3)) then       
       do nn = 1, nSoilHorizons_mHM
          ii = ii + 1
          var = nc%createVariable("SWC_L"//trim(num2str(nn,'(i2.2)')),dtype,dims1)
          call writeVariableAttributes(&
               var,'soil water content of soil layer'//trim(num2str(nn)), "mm")
          tmpvars(ii) = OutputVariable(var, L1_soilMoist(:,nn),avg=.true.)
       end do
    end if
    
    if (outputFlxState(4)) then       
       do nn = 1, nSoilHorizons_mHM
          ii = ii + 1
          var = nc%createVariable("SM_L"//trim(num2str(nn,'(i2.2)')), dtype,dims1)
          call writeVariableAttributes(&
               var,'volumetric soil moisture of soil layer'//trim(num2str(nn)), "mm mm-1")
          tmpvars(ii) = OutputVariable(var, L1_soilMoistVol(:,nn), avg=.true.)
       end do
    end if

    if (outputFlxState(5)) then
       ii = ii + 1
       var = nc%createVariable("SM_Lall", dtype, dims1)
       call writeVariableAttributes(var,"average soil moisture over all layers", "mm mm-1")
       tmpvars(ii) = OutputVariable(var,L1_soilMoistVolAvg, avg=.true.)       
    end if

    if (outputFlxState(6)) then
       ii = ii + 1
       var = nc%createVariable("sealedSTW", dtype, dims1)
       call writeVariableAttributes(var,"reservoir of sealed areas (sealedSTW)", "mm")
       tmpvars(ii) = OutputVariable(var,L1_sealSTW, avg=.true.)       
    end if

    if (outputFlxState(7)) then
       ii = ii + 1
       var = nc%createVariable("unsatSTW", dtype, dims1)
       call writeVariableAttributes(var, "reservoir of unsaturated zone", "mm")
       tmpvars(ii) = OutputVariable(var,L1_unsatSTW, avg=.true.)       
    end if

    if (outputFlxState(8)) then
       ii = ii + 1
       var = nc%createVariable("satSTW", dtype, dims1)
       call writeVariableAttributes(var, "water level in groundwater reservoir", "mm")
       tmpvars(ii) = OutputVariable(var,L1_satSTW, avg=.true.)       
    end if

    if (outputFlxState(18)) then
       ii = ii + 1
       var = nc%createVariable("Neutrons", dtype, dims1)
       call writeVariableAttributes(var, "ground albedo neutrons", "cph")
       tmpvars(ii) = OutputVariable(var,L1_neutrons, avg=.true.)       
    end if

    ! fluxes
    unit = fluxesUnit(ibasin)

    if (outputFlxState(9)) then
       ii = ii + 1
       var = nc%createVariable("PET", dtype, dims1)
       call writeVariableAttributes(var, "potential Evapotranspiration", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_pet)       
    end if

    if (outputFlxState(10)) then
       ii = ii + 1
       var = nc%createVariable("aET", dtype, dims1)
       call writeVariableAttributes(var, "actual Evapotranspiration", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_aet)       
    end if

    if (outputFlxState(11)) then
       ii = ii + 1
       var = nc%createVariable("Q", dtype, dims1)
       call writeVariableAttributes(var, "total runoff generated by every cell", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_total_runoff)       
    end if

    if (outputFlxState(12)) then
       ii = ii + 1
       var = nc%createVariable("QD", dtype, dims1)
       call writeVariableAttributes(&
            var, "direct runoff generated by every cell (runoffSeal)", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_runoffSeal, multip=L1_fSealed )       
    end if

    if (outputFlxState(13)) then
       ii = ii + 1
       var = nc%createVariable("QIf", dtype, dims1)
       call writeVariableAttributes(&
            var, "fast interflow generated by every cell (fastRunoff)", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_fastRunoff, multip=L1_fNotSealed)       
    end if

    if (outputFlxState(14)) then
       ii = ii + 1
       var = nc%createVariable("QIs", dtype, dims1)
       call writeVariableAttributes(&
            var, "slow interflow generated by every cell (slowRunoff)", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_slowRunoff, multip=L1_fNotSealed)       
    end if

    if (outputFlxState(15)) then
       ii = ii + 1
       var = nc%createVariable("QB", dtype, dims1)
       call writeVariableAttributes(&
            var, "baseflow generated by every cell", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_baseflow, multip=L1_fNotSealed)       
    end if

    if (outputFlxState(16)) then
       ii = ii + 1
       var = nc%createVariable("recharge", dtype, dims1)
       call writeVariableAttributes(&
            var, "groundwater recharge", trim(unit))
       tmpvars(ii) = OutputVariable(var, L1_percol, multip=L1_fNotSealed)       
    end if

    if (outputFlxState(17)) then
       do nn = 1, nSoilHorizons_mHM
          ii = ii + 1
          var = nc%createVariable("soil_infil_L"//trim(num2str(nn,'(i2.2)')), dtype,dims1)
          call writeVariableAttributes(&
               var,"infiltration flux from soil layer"//trim(num2str(nn)), unit)
          tmpvars(ii) = OutputVariable(var, L1_infilSoil(:,nn), multip=L1_fNotSealed)
       end do
    end if
    
    initOutput = OutputDataset(nc,tmpvars(1:ii),mask1)
    
  end function initOutput
  
  subroutine mapCoordinates(ibasin, level, y, x)

    implicit none
    
    integer(i4),      intent(in) :: iBasin
    type(gridGeoRef), intent(in) :: level
    real(dp),dimension(:), allocatable, intent(out) :: x, y
    real(dp)                    :: cellsize
    integer(i4)                 :: ii, ncols, nrows
    
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

  subroutine geoCoordinates(ibasin, level, lat, lon)

    use mo_global_variables, only : latitude,longitude
    
    implicit none
    
    integer(i4),      intent(in) :: iBasin
    type(gridGeoRef), intent(in) :: level
    real(dp),dimension(:,:), allocatable, intent(out) :: lat, lon
    integer(i4)                 :: ii, ncols, nrows, pos
    
    nrows    = level%nrows(ibasin)
    ncols    = level%ncols(ibasin)

    allocate(lat(nrows,ncols), lon(nrows,ncols))
    
    pos = 1
    if ( ibasin .gt. 1 ) then
       do ii = 1, ibasin -1
          pos = pos + level%ncols(ii) * level%nrows(ii)
       end do
    end if

    lat = reshape( latitude(pos:pos+nrows*ncols-1),  shape(lat))
    lon = reshape( longitude(pos:pos+nrows*ncols-1), shape(lon))
    
  end subroutine geoCoordinates

END MODULE mo_write_fluxes_states

