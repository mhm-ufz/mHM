
module mo_write_fluxes_states

  ! mHM writing module

  use mo_kind, ONLY: i4, dp
  use mo_string_utils, only : num2str
  use mo_global_variables,  only : outputFlxState !, outputVariable, outputVariableNames
  use mo_mhm_constants, only: nodata_i4, nodata_dp
  ! use netcdf,           only: nf90_create, NF90_NOERR, NF90_NETCDF4,nf90_def_dim,NF90_UNLIMITED, &
  !                            nf90_close,NF90_DOUBLE,nf90_def_var, nf90_put_att, nf90_put_var

  use mo_netcdf, only : NcDataset, NcDimension, NcVariable

  use mo_init_states,         only : get_basin_info

  IMPLICIT NONE

  type OutputVariable
     type(NcVariable)                           :: ncvar
     real(dp), dimension(:), pointer            :: data       => null()       ! point to data
     real(dp), dimension(:), pointer            :: multiplier => null()       ! point to data
     real(dp), dimension(:), allocatable        :: accum                   ! the data accumulator
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
     procedure, public :: writeDataset
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
  
  function initOutputVariable(ncvar, data, multip, avg)!,multip)
    type(NcVariable), intent(in)                                  :: ncvar
    real(dp), dimension(:), target, intent(in)                    :: data
    real(dp), dimension(size(data)), target, intent(in), optional :: multip
    logical, intent(in),optional                                             :: avg
    type(OutputVariable)                                          :: initOutputVariable

    allocate(initOutputVariable%accum(size(data)))
    initOutputVariable%ncvar = ncvar
    initOutputVariable%data  => data
    initOutputVariable%accum = 0
    if (present(multip)) then
       initOutputVariable%multiplier => multip
    end if
    if (present(avg)) then
       initOutputVariable%avg = avg
    end if
  end function initOutputVariable

  function initOutputDataset(nc,vars,mask1)

    type(NcDataset), intent(in)                  :: nc
    type(OutputVariable),dimension(:),intent(in) :: vars
    logical,dimension(:,:),intent(in)            :: mask1 !,mask11
    type(OutputDataset)                          :: initOutputDataset

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

  subroutine writeDataset(self)
    class(OutputDataset) :: self
    type(OutputVariable) :: var
    integer(i4)          :: ii
    real(dp),dimension(:,:),allocatable :: tmpdata1

    do ii = 1,size(self%vars)
       var = self%vars(ii)
       if (var%avg) then
          var%accum = var%accum / real(var%counter, dp)
       end if
       tmpdata1 = unpack(var%accum, self%mask1, nodata_dp)
       call var%ncvar%putData(transpose(tmpdata1),(/self%counter,1,1/))
       self%vars(ii)%accum = 0
       self%vars(ii)%counter = 0
    end do
    
    self%counter = self%counter + 1

  end subroutine writeDataset
  
  subroutine close(self)
    class(OutputDataset) :: self
    call self%nc%close()
  end subroutine close
  
  function initOutputFile(ibasin)

    use mo_global_variables,  only : dirOut !, outputVariable, outputVariableNames

    ! creates the output netcdf file based on the variables passed in outvars
    ! return the netcdf file id
    integer(i4),intent(in)                          :: ibasin
    type(NcDataset)                                 :: initOutputFile, nc
    type(NcDimension),dimension(3)                  :: dimids1
    type(NcVariable)                                :: easting, northing, time
    character(256)                                  :: fname
    integer(i4)                                     :: nrows1,ncols1 !,nrows11,ncols11,varid

    call get_basin_info (ibasin,  1,  ncols1, nrows1)
    
    fname = trim(dirOut(ibasin)) // 'mHM_Fluxes_States.nc'
    nc = NcDataset(trim(fname),"w")    
    
    dimids1 = (/&
         nc%createDimension("time",     0), &
         nc%createDimension("northing", nrows1), &
         nc%createDimension("easting",  ncols1) &
         /)

    ! time -> fill units
    time = nc%createVariable("time", "i32", (/ dimids1(1) /))
    call time%createAttribute("units", "hours since ... ")
    call time%createAttribute("long_name","time")
    
    ! northing
    northing = nc%createVariable("northing", "f64", (/ dimids1(2) /))
    call northing%createAttribute("units","m")
    call northing%createAttribute("long_name","y-coordinate in cartesian coordinates GK4")
    
    ! easting
    easting = nc%createVariable("easting", "f64", (/ dimids1(3) /))
    call easting%createAttribute("units","m")
    call easting%createAttribute("long_name","x-coordinate in cartesian coordinates GK4")
                
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

    integer(i4), intent(in)              :: ibasin
    ! states
    real(dp), dimension(:),   intent(in) :: L1_inter          ! Interception
    real(dp), dimension(:),   intent(in) :: L1_snowPack       ! Snowpack
    real(dp), dimension(:,:), intent(in) :: L1_soilMoist      ! Soil moisture of each horizon                   
    real(dp), dimension(:,:), intent(in) :: L1_soilMoistVol
    real(dp), dimension(:),   intent(in) :: L1_soilMoistVolAvg
    real(dp), dimension(:),   intent(in) :: L1_sealSTW      ! Retention storage of impervious areas
    real(dp), dimension(:),   intent(in) :: L1_unsatSTW     ! Upper soil storage
    real(dp), dimension(:),   intent(in) :: L1_satSTW       ! Groundwater storage
    real(dp), dimension(:),   intent(in) :: L1_neutrons     ! ground albedo neutrons
    ! fluxes
    real(dp), dimension(:),   intent(inout) :: L1_pet          ! potential evapotranspiration (PET)
    real(dp), dimension(:), intent(inout) :: L1_aet
    real(dp), dimension(:),   intent(inout) :: L1_total_runoff ! Generated runoff
    real(dp), dimension(:),   intent(inout) :: L1_runoffSeal   ! Direct runoff from impervious areas
    real(dp), dimension(:),   intent(inout) :: L1_fastRunoff   ! Fast runoff component
    real(dp), dimension(:),   intent(inout) :: L1_slowRunoff   ! Slow runoff component
    real(dp), dimension(:),   intent(inout) :: L1_baseflow     ! Baseflow
    real(dp), dimension(:),   intent(inout) :: L1_percol       ! Percolation
    real(dp), dimension(:,:), intent(inout) :: L1_infilSoil    ! Infiltration
    ! helper
    real(dp), dimension(:),   intent(inout) :: L1_fSealed      ! 
    real(dp), dimension(:),   intent(inout) :: L1_fNotSealed      ! 
    ! output
    type(OutputDataset)                  :: initOutput

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
    ! SWC per horizon

    if (outputFlxState(3)) then       
       do nn = 1, nSoilHorizons_mHM
          ii = ii + 1
          var = nc%createVariable("SWC_L"//trim(num2str(nn,'(i2.2)')),dtype,dims1)
          call writeVariableAttributes(&
               var,'soil water content of soil layer'//trim(num2str(nn)), "mm")
          tmpvars(ii) = OutputVariable(var, L1_soilMoist(:,nn),avg=.true.)
       end do
    end if
    
    ! SM per horizon

    if (outputFlxState(4)) then       
       do nn = 1, nSoilHorizons_mHM
          ii = ii + 1
          var = nc%createVariable("SM_L"//trim(num2str(nn,'(i2.2)')), dtype,dims1)
          call writeVariableAttributes(&
               var,'volumetric soil moisture of soil layer'//trim(num2str(nn)), "mm mm-1")
          tmpvars(ii) = OutputVariable(var, L1_soilMoistVol(:,nn), avg=.true.)
       end do
    end if
    ! SM - average

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
       tmpvars(ii) = OutputVariable(var, L1_pet, avg=.true.)       
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

END MODULE mo_write_fluxes_states

