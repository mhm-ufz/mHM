!> \file mo_write_fluxes_states.f90

!> \brief Creates NetCDF output for different fluxes and state variabels of mHM.

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!> \authors Matthias Zink
!> \date Apr 2013

MODULE mo_write_fluxes_states

  ! This module creates the output for mHM.

  ! Written Matthias Zink, Apr 2012

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: WriteFluxStateInit
  PUBLIC :: WriteFluxState
  PUBLIC :: CloseFluxState_file

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          WriteFluxStateInit

  !>        \brief Initialization of the NetCDF for mHM outputs.

  !>        \details The NetCDF file is set up with its dimensions, variables and variable attributes. Additionally 
  !>                 the arrays for aggregating the output to the output time step are allocated.

  !     INTENT(IN)
  !>         \param[in]  "integer(i4), intent(in) :: iBasin"           ! mumber of subbasin
  !>         \param[in]  "integer(i4), intent(in) :: output_timeStep"  ! timestep (e.g. hour, day,..) of the output

  !     INTENT(INOUT)
  !>         \param[inout] "real(dp), allocatable :: L1_inter_out(:)"        ! Interception
  !>         \param[inout] "real(dp), allocatable :: L1_snowPack_out(:)"     ! Snowpack
  !>         \param[inout] "real(dp), allocatable :: L1_soilMoist_out(:,:)"  ! Soil moisture of each horizon
  !>         \param[inout] "real(dp), allocatable :: L1_sealSTW_out(:)"      ! Retention storage of impervious areas
  !>         \param[inout] "real(dp), allocatable :: L1_unsatSTW_out(:)"     ! Upper soil storage
  !>         \param[inout] "real(dp), allocatable :: L1_satSTW_out(:)"       ! Groundwater storage
  !>         \param[inout] "real(dp), allocatable :: L1_neutrons_out(:)"     ! ground albedo neutrons
  !>         \param[inout] "real(dp), allocatable :: L1_pet_out(:)"          ! potential evapotranspiration (PET) 
  !>         \param[inout] "real(dp), allocatable :: L1_aETSoil_out(:,:)"    ! actual ET of each horizon
  !>         \param[inout] "real(dp), allocatable :: L1_aETCanopy_out(:)"    ! Real evaporation intensity from canopy
  !>         \param[inout] "real(dp), allocatable :: L1_aETSealed_out(:)"    ! Actual ET from free-water surfaces
  !>         \param[inout] "real(dp), allocatable :: L1_total_runoff_out(:)" ! Generated runoff
  !>         \param[inout] "real(dp), allocatable :: L1_runoffSeal_out(:)"   ! Direct runoff from impervious areas
  !>         \param[inout] "real(dp), allocatable :: L1_fastRunoff_out(:)"   ! Fast runoff component
  !>         \param[inout] "real(dp), allocatable :: L1_slowRunoff_out(:)"   ! Slow runoff component
  !>         \param[inout] "real(dp), allocatable :: L1_baseflow_out(:)"     ! Baseflow
  !>         \param[inout] "real(dp), allocatable :: L1_percol_out(:)"       ! Percolation
  !>         \param[inout] "real(dp), allocatable :: L1_infilSoil_out(:)"    ! Infiltration  

  !     INTENT(OUT)
  !>         \param[out] "integer(i4), intent(out) :: ncid" ! ID of NetCDF to be written in

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified: R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                                     during writing of the netcdf file
  !                   Matthias Zink       , Feb. 2014 - added aditional output: pet
  !                   V. Prykhodk, J. Mai,  Nov. 2014 - adding new variable infilSoil - case 16

  Subroutine WriteFluxStateInit(&
       ! Input
       iBasin                 , &
       output_timeStep        , &
       ! Inout: States L1
       L1_inter_out           , & ! Interception
       L1_snowPack_out        , & ! Snowpack
       L1_soilMoist_out       , & ! Soil moisture of each horizon
       L1_sealSTW_out         , & ! Retention storage of impervious areas
       L1_unsatSTW_out        , & ! Upper soil storage
       L1_satSTW_out          , & ! Groundwater storage
       L1_neutrons_out        , & ! Ground albedo neutrons
       ! Inout: Fluxes L1
       L1_pet_out             , & ! potential evapotranspiration (PET)
       L1_aETSoil_out         , & ! actual ET
       L1_aETCanopy_out       , & ! Real evaporation intensity from canopy
       L1_aETSealed_out       , & ! Actual ET from free-water surfaces
       L1_total_runoff_out    , & ! Generated runoff
       L1_runoffSeal_out      , & ! Direct runoff from impervious areas
       L1_fastRunoff_out      , & ! Fast runoff component
       L1_slowRunoff_out      , & ! Slow runoff component
       L1_baseflow_out        , & ! Baseflow
       L1_percol_out          , & ! Percolation 
       L1_infilSoil_out       , & ! Infiltrationf
       ! Output
       ncid)

    use mo_global_variables,  only : & 
                                     dirOut            , & ! output directory
                                     outputFlxState    , & ! definition which output to write
                                     yCoor, xCoor      , & ! kartesian coordinates
                                     lons, lats        , & ! geographic coordinates
                                     level1            , & ! grid information for level 1
                                     evalPer           , & ! information about time period
                                     timestep          , & ! model time step
                                     nSoilHorizons_mHM , & ! Number of horizons to model
                                     simPer            , & ! simulation period
                                     NTSTEPDAY             ! time steps per day
    use mo_mhm_constants,      only: nodata_dp             ! global no data value
    !
    use mo_init_states,        only: get_basin_info
    use mo_julian,             only: dec2date
    use mo_message,            only: message 
    use mo_ncwrite,            only: create_netCDF, write_static_netCDF, V, GAtt
    use mo_set_netcdf_outputs, only: set_netCDF
    use mo_string_utils,       only: num2str
    !   
    implicit none
    
    !
    integer(i4),                           intent(in)    :: iBasin    ! mumber of subbasin
    integer(i4),                           intent(out)   :: ncid      ! ID of NetCDF
    !
    integer(i4),                           intent(in)    :: output_timeStep
    !
    ! States L1
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_inter_out        ! Interception
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_snowPack_out     ! Snowpack
    real(dp), dimension(:,:), allocatable, intent(inout) :: L1_soilMoist_out    ! Soil moisture of each horizon
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_sealSTW_out      ! Retention storage of impervious areas
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_unsatSTW_out     ! Upper soil storage
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_satSTW_out       ! Groundwater storage
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_neutrons_out     ! ground albedo neutrons
    ! Fluxes L1
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_pet_out          ! potential evapotranspiration (PET)
    real(dp), dimension(:,:), allocatable, intent(inout) :: L1_aETSoil_out      ! actual ET of each horizon
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_aETCanopy_out    ! Real evaporation intensity from canopy
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_aETSealed_out    ! Actual ET from free-water surfaces
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_total_runoff_out ! Generated runoff
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_runoffSeal_out   ! Direct runoff from impervious areas
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_fastRunoff_out   ! Fast runoff component
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_slowRunoff_out   ! Slow runoff component
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_baseflow_out     ! Baseflow
    real(dp), dimension(:),   allocatable, intent(inout) :: L1_percol_out       ! Percolation
    real(dp), dimension(:,:), allocatable, intent(inout) :: L1_infilSoil_out    ! Infiltration 
    !
    ! local
    integer(i4)                                      :: i
    integer(i4)                                      :: VarNo      ! counter for variable
    integer(i4)                                      :: totalVarNo ! total number of variables
    character(8)                                     :: date
    character(10)                                    :: time
    character(20)                                    :: unit
    character(256)                                   :: fName, dummy
    !
    integer(i4)                                      :: nrows, ncols, nCells
    real(dp)                                         :: jday_frac
    integer(i4)                                      :: day, month, year
    real(dp)                                         :: nTimeSteps
    !  
    ! first get number of L1-cells for this iBasin
    call get_basin_info(iBasin, 1, nrows, ncols, ncells=nCells) 
    !
    VarNo = 1

    ! Initialize NetCDF dimensions and Variables incl. attributes

    ! ---------
    ! 1. States
    ! ---------

    ! determine total number of variables to print
    totalVarNo = count(outputFlxState)
    
    ! output for soil, every layer increases number of output parameters
    if (outputFlxState(3))  totalVarNo = totalVarNo + nSoilHorizons_mHM - 1
    if (outputFlxState(4))  totalVarNo = totalVarNo + nSoilHorizons_mHM - 1
    if (outputFlxState(17)) totalVarNo = totalVarNo + nSoilHorizons_mHM - 1 
    
    ! Initialize NetCDF dimensions and Variables incl. attributes
    call set_netCDF(totalVarNo, level1%nrows(iBasin), level1%ncols(iBasin))
    !
    if (outputFlxState(1)) then
    
       allocate( L1_inter_out(nCells) )
       L1_inter_out = 0.0_dp
    
       ! name
       V(VarNo+5)%name          = "interception"
       ! unit
       V(VarNo+5)%att(1)%values = "mm"
       ! long_name
       V(VarNo+5)%att(2)%values = "canopy interception storage"
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(2)) then
    
       allocate( L1_snowPack_out(nCells) )
       L1_snowPack_out = 0.0_dp
       
       ! name
       V(VarNo+5)%name          = "snowpack"
       ! unit
       V(VarNo+5)%att(1)%values = "mm"
       ! long_name
       V(VarNo+5)%att(2)%values = "depth of snowpack"
       VarNo = VarNo + 1
    end if


    ! allocate if any of the following condition is true
    if( outputFlxState(3) .OR. outputFlxState(4) .OR. outputFlxState(5) ) then
      allocate( L1_soilMoist_out(nCells, nSoilHorizons_mHM) )
      L1_soilMoist_out = 0.0_dp
    end if

    !
    ! total water amount in soil layers
    if (outputFlxState(3)) then
       do i = 1, nSoilHorizons_mHM
          ! all Vars same attributes
          ! unit
          V(VarNo+5)%att(1)%values = "mm"
          write(V(VarNo+5)%name, "('SWC_L', i2.2)") i
          write(V(VarNo+5)%att(2)%values,"('soil water content of soil layer',i2)") i
          VarNo = VarNo + 1
       end do
    end if
    !
    ! volumetric soil moisture in every soil layer
    if (outputFlxState(4)) then
       do i = 1, nSoilHorizons_mHM
          ! unit
          V(VarNo+5)%att(1)%values = "mm mm-1"
          write(V(VarNo+5)%name, "('SM_L', i2.2)") i
          write(V(VarNo+5)%att(2)%values,"('volumetric soil moisture of soil layer',i2)") i
          VarNo = VarNo + 1
       end do
       ! 
    end if
    !
    ! volumetric soil moisture - layer average
    if (outputFlxState(5)) then
       ! name
       V(VarNo+5)%name          = "SM_Lall"
       ! unit
       V(VarNo+5)%att(1)%values = "mm mm-1"
       ! long_name
       V(VarNo+5)%att(2)%values = "average soil moisture over all layers"
       VarNo = VarNo + 1
    end if
    
    !
    if (outputFlxState(6)) then

       allocate( L1_sealSTW_out(nCells) )
       L1_sealSTW_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "sealedSTW"
       ! unit
       V(VarNo+5)%att(1)%values = "mm"
       ! long_name
       V(VarNo+5)%att(2)%values = "reservoir of sealed areas (sealedSTW)"
       VarNo = VarNo + 1
    end if
    
    !
    if (outputFlxState(7)) then

       allocate( L1_unsatSTW_out(nCells) )
       L1_unsatSTW_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "unsatSTW"
       ! unit
       V(VarNo+5)%att(1)%values = "mm"
       ! long_name
       V(VarNo+5)%att(2)%values = "reservoir of unsaturated zone"
       VarNo = VarNo + 1
    end if
    
    !
    if (outputFlxState(8)) then

       allocate( L1_satSTW_out(nCells) )
       L1_satSTW_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "satSTW"
       ! unit
       V(VarNo+5)%att(1)%values = "mm"
       ! long_name
       V(VarNo+5)%att(2)%values = "water level in groundwater reservoir"
       VarNo = VarNo + 1
    end if
    
    if (outputFlxState(18)) then

       allocate( L1_neutrons_out(nCells) )
       L1_neutrons_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "Neutrons"
       ! unit
       V(VarNo+5)%att(1)%values = "cph"
       ! long_name
       V(VarNo+5)%att(2)%values = "ground albedo neutrons"
       VarNo = VarNo + 1
    end if

    ! ---------
    ! 2. Fluxes
    ! ---------
    nTimeSteps = ( simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1 ) * NTSTEPDAY
    if (output_timeStep > 0) then
       if ( (timestep*output_timeStep)     .EQ. 1 ) then
          unit = 'mm h-1'
       else
          ! write(unit,"('mm ', i2, 'h-1')") timestep
          unit = 'mm '//trim(adjustl(num2str(timestep)))//'h-1'
       end if
    else
       select case(output_timeStep)
       case(0) ! only at last time step
          unit = 'mm '//trim(adjustl(num2str(nint(nTimeSteps))))//'h-1'
       case(-1) ! daily
          unit = 'mm d-1'
       case(-2) ! monthly
          unit = 'mm month-1'
       case(-3) ! yearly
          unit = 'mm a-1'
       case default ! no output at all
          unit = ''
       end select
    endif

    !
    if (outputFlxState(9)) then

       allocate( L1_pet_out(nCells))
       L1_pet_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "PET"
       !unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "potential Evapotranspiration"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(10)) then

       allocate( L1_aETSoil_out(nCells, nSoilHorizons_mHM), L1_aETCanopy_out(nCells), &
                 L1_aETSealed_out(nCells) )
       L1_aETSoil_out   = 0.0_dp
       L1_aETCanopy_out = 0.0_dp
       L1_aETSealed_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "aET"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "actual Evapotranspiration"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(11)) then

       allocate( L1_total_runoff_out(nCells) )
       L1_total_runoff_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "Q"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "total runoff generated by every cell"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(12)) then

       allocate( L1_runoffSeal_out(nCells) )
       L1_runoffSeal_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "QD"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "direct runoff generated by every cell (runoffSeal)"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(13)) then

       allocate( L1_fastRunoff_out (nCells) )
       L1_fastRunoff_out  = 0.0_dp
       
       ! name
       V(VarNo+5)%name          = "QIf"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "fast interflow generated by every cell (fastRunoff)"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(14)) then

       allocate( L1_slowRunoff_out(nCells) )
       L1_slowRunoff_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "QIs"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "slow interflow generated by every cell (slowRunoff)"
       VarNo = VarNo + 1
       ! 
    end if

    !
    if (outputFlxState(15)) then

       allocate( L1_baseflow_out(nCells) )
       L1_baseflow_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "QB"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "baseflow generated by every cell"
       VarNo = VarNo + 1        
    end if

    !
    if (outputFlxState(16)) then

       allocate( L1_percol_out(nCells) )
       L1_percol_out = 0.0_dp

       ! name
       V(VarNo+5)%name          = "recharge"
       ! unit
       V(VarNo+5)%att(1)%values = trim(unit)
       ! long_name
       V(VarNo+5)%att(2)%values = "groundwater recharge"
       VarNo = VarNo + 1
    end if

    ! 
    if (outputFlxState(17)) then 
       allocate( L1_infilSoil_out(nCells, nSoilHorizons_mHM) ) 
       L1_infilSoil_out = 0.0_dp 
       do i = 1, nSoilHorizons_mHM 
          ! all Vars same attributes 
          ! name 
          write(V(VarNo+5)%name, "('soil_infil_L', i2.2)") i 
          ! unit 
          V(VarNo+5)%att(1)%values = trim(unit) 
          ! long_name 
          write(V(VarNo+5)%att(2)%values,"('infiltration flux from soil layer',i2)") i 
          VarNo = VarNo + 1 
       end do
    end if

    !
    !*******************************************************
    ! here come the attributes which are equal for all vars
    ! up to Var 5 are static Vars
    !*******************************************************
    write(dummy,*) nodata_dp
    do i = 6, totalVarNo + 5
       ! scale_factor
       V(i)%att(3)%values = "1."
       ! _FillValue
       V(i)%att(4)%values = trim(adjustl(dummy))
       ! missing_value
       V(i)%att(5)%values = trim(adjustl(dummy))
       ! coordinates
       V(i)%att(6)%values = "lat lon"
    end do

    !
    !*******************************************************
    ! set Attributes for dimensions
    ! set unit for dimension time 
    ! (e.g. hours since 2008-01-01 00:00:00)
    !*******************************************************
    jday_frac = real( (evalPer(iBasin)%julStart), dp )
    call dec2date(jday_frac, dd=day, mm=month, yy=year)
    write(dummy,"('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day
    V(3)%att(1)%values = trim(dummy)

    !
    ! global attributes
    GAtt(1)%name    = "title"
    GAtt(1)%values  = "mHMv5 simulation outputs"
    !
    GAtt(2)%name    = "creating_date"
    call DATE_AND_TIME(DATE=date, TIME=time)
    dummy = trim(date) // trim(time)
    write(dummy,"(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1:4), &
      date(5:6), date(7:8), time(1:2), time(3:4), time(5:6)
    GAtt(2)%values   = trim(dummy)
    !
    GAtt(3)%name    = "institution"
    GAtt(3)%values   = trim("Helmholtz Centre for Environmental Research - UFZ,") // &
                       trim("Department of Computational Hydrosystems, Stochastic Hydrology Group")
    ! to create a new netCDF
    fName = trim(dirOut(iBasin)) // 'mHM_Fluxes_States.nc'
    !
    call message('')
    call message('  OUTPUT: Writing NetCDF file in')
    call message('     to ', trim(fName))
    !
    call create_netCDF(trim(fName), ncid)   

    !
    !*******************************************************
    ! write static variables  
    ! put coordinates sytem to the NetCDF
    !*******************************************************
    call CoordSystem(iBasin)
    !
    V(1)%G1_d        => xCoor
    V(2)%G1_d        => yCoor
    V(4)%G2_d        => lons
    V(5)%G2_d        => lats
    call write_static_netCDF(ncid) 
    !
  end Subroutine WriteFluxStateInit

  ! ------------------------------------------------------------------

  !      NAME
  !          WriteFluxState

  !>        \brief Ouput is written to NetCDF file.

  !>        \details All the outputs defined in the outputs namlist are written to a NetCDF file. 

  !     INTENT(IN)
  !>         \param[in] "integer(i4) :: hours_eval"       ! hours since beginning of evaluation period
  !>         \param[in] "integer(i4) :: curr_wstep"       ! current write out step
  !>         \param[in] "integer(i4) :: ncid"             ! ID of NetCDF to be written in
  !>         \param[in] "integer(i4) :: iBasin"           ! mumber of subbasin
  !>         \param[in] "integer(i4) :: output_timeStep"  ! timestep (e.g. hour, day,..) of the output
  !>         \param[in] "logical     :: mask(:,:)"        ! mask for unpacking vectorized data to array 
  !>         \param[in] "real(dp), allocatable    :: L1_inter_out(:)"        ! Interception
  !>         \param[in] "real(dp), allocatable    :: L1_snowPack_out(:)"     ! Snowpack
  !>         \param[in] "real(dp), allocatable    :: L1_soilMoist_out(:,:)"  ! Soil moisture of each horizon
  !>         \param[in] "real(dp), allocatable    :: L1_sealSTW_out(:)"      ! Retention storage of impervious areas
  !>         \param[in] "real(dp), allocatable    :: L1_unsatSTW_out(:)"     ! Upper soil storage
  !>         \param[in] "real(dp), allocatable    :: L1_satSTW_out(:)"       ! Groundwater storage
  !>         \param[in] "real(dp), allocatable    :: L1_neutrons_out(:)"     ! Ground albedo neutrons
  !>         \param[in] "real(dp), allocatable    :: L1_pet_out(:)"          ! potential evapotranspiration (PET) 
  !>         \param[in] "real(dp), allocatable    :: L1_aETSoil_out(:,:)"    ! actual ET of each horizon
  !>         \param[in] "real(dp), allocatable    :: L1_aETCanopy_out(:)"    ! Real evaporation intensity from canopy
  !>         \param[in] "real(dp), allocatable    :: L1_aETSealed_out(:)"    ! Actual ET from free-water surfaces
  !>         \param[in] "real(dp), allocatable    :: L1_total_runoff_out(:)" ! Generated runoff
  !>         \param[in] "real(dp), allocatable    :: L1_runoffSeal_out(:)"   ! Direct runoff from impervious areas
  !>         \param[in] "real(dp), allocatable    :: L1_fastRunoff_out(:)"   ! Fast runoff component
  !>         \param[in] "real(dp), allocatable    :: L1_slowRunoff_out(:)"   ! Slow runoff component
  !>         \param[in] "real(dp), allocatable    :: L1_baseflow_out(:)"     ! Baseflow
  !>         \param[in] "real(dp), allocatable    :: L1_percol_out(:)"       ! Percolation
  !>         \param[in] "real(dp), allocatable    :: L1_infilSoil_out(:)"    ! Infiltration  

  !     INTENT(INOUT)
  !          None

  !     INTENT(OUT)
  !          None         

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified: R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                                     during writing of the netcdf file
  !                   L. Samaniego et al.,  Dec  2013 - nullify pointer
  !                   Matthias Zink,        Feb. 2014 - added aditional output: pet
  !                   V. Prykhodk, J. Mai,  Nov. 2014 - adding new variable infilSoil - case 16

  Subroutine WriteFluxState( &
       ! input
       hours_eval          , &
       curr_wstep          , &
       ncid                , &
       iBasin              , &
       mask                , &
       ! input: States L1
       L1_inter_out        , & ! Interception
       L1_snowPack_out     , & ! Snowpack
       L1_soilMoist_out    , & ! Soil moisture of each horizon
       L1_sealSTW_out      , & ! Retention storage of impervious areas
       L1_unsatSTW_out     , & ! Upper soil storage
       L1_satSTW_out       , & ! Groundwater storage
       L1_neutrons_out     , & ! Ground albedo neutrons
       ! input: Fluxes L1
       L1_pet_out             , & ! potential evapotranspiration (PET)
       L1_aETSoil_out      , & ! actual ET
       L1_aETCanopy_out    , & ! Real evaporation intensity from canopy
       L1_aETSealed_out    , & ! Actual ET from free-water surfaces
       L1_total_runoff_out , & ! Generated runoff
       L1_runoffSeal_out   , & ! Direct runoff from impervious areas
       L1_fastRunoff_out   , & ! Fast runoff component
       L1_slowRunoff_out   , & ! Slow runoff component
       L1_baseflow_out     , & ! Baseflow
       L1_percol_out       , & ! Percolation
       L1_soilMoistSat_out , & ! Saturation soil moisture for each horizon [mm]
       !                       ! part (s1:e1) from global variable
       L1_infilSoil_out      & ! Infiltration for each horizon
       ! output
       )
                             
    use mo_global_variables, only : & 
         outputFlxState    , & ! definition which output to write
         level1            , & ! grid information for level 1
         nSoilHorizons_mHM     ! Number of horizons to model

    use mo_mhm_constants,    only: & 
         nodata_dp             ! global no data value
    !
    use mo_ncwrite,          only: write_dynamic_netCDF, V

    implicit none
    !
    integer(i4),                           intent(in)    :: hours_eval ! hours since beginning of evaluation period
    integer(i4),                           intent(in)    :: curr_wstep ! current write out step
    integer(i4),                           intent(in)    :: ncid       ! ID of NetCDF
    integer(i4),                           intent(in)    :: iBasin     ! mumber of subbasin
    logical,   dimension(:,:),             intent(in)    :: mask       ! mask for unpacking
    
    ! States L1
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_inter_out        ! Interception
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_snowPack_out     ! Snowpack
    real(dp), dimension(:,:), allocatable, intent(in)    :: L1_soilMoist_out    ! Soil moisture of each horizon
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_sealSTW_out      ! Retention storage of impervious areas
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_unsatSTW_out     ! Upper soil storage
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_satSTW_out       ! Groundwater storage
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_neutrons_out     ! Ground albedo neutrons
    ! Fluxes L1
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_pet_out          ! potential evapotranspiration (PET)
    real(dp), dimension(:,:), allocatable, intent(in)    :: L1_aETSoil_out      ! actual ET of each horizon
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_aETCanopy_out    ! Real evaporation intensity from canopy
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_aETSealed_out    ! Actual ET from free-water surfaces
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_total_runoff_out ! Generated runoff
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_runoffSeal_out   ! Direct runoff from impervious areas
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_fastRunoff_out   ! Fast runoff component
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_slowRunoff_out   ! Slow runoff component
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_baseflow_out     ! Baseflow
    real(dp), dimension(:),   allocatable, intent(in)    :: L1_percol_out       ! Percolation
    real(dp), dimension(:,:),              intent(in)    :: L1_soilMoistSat_out ! Saturation soil moisture for each 
    !                                                                           ! horizon [mm]
    real(dp), dimension(:,:),  allocatable, intent(in)   :: L1_infilSoil_out    ! Infiltration for each horizon

    ! local variables
    integer(i4)                                     :: totalVarNo ! total number of variables
    integer(i4)                                     :: VarNo      ! counter for variable
    integer(i4)                                     :: i
    integer(i4),                             target :: tstep      ! time step of output
    real(dp), dimension(:,:,:), allocatable, target :: OutPut
    !
    VarNo = 1
    !
    totalVarNo = count(outputFlxState)
    if( outputFlxState(3) )  totalVarNo = totalVarNo + nSoilHorizons_mHM - 1
    if( outputFlxState(4) )  totalVarNo = totalVarNo + nSoilHorizons_mHM - 1
    if( outputFlxState(17) ) totalVarNo = totalVarNo + nSoilHorizons_mHM - 1 
    
    ! Write files
    allocate( OutPut(level1%nrows(iBasin), level1%ncols(iBasin), totalVarNo) )

    OutPut = nodata_dp

    ! ---------
    ! 1. States    
    ! ---------
    if (outputFlxState(1)) then

       OutPut(:,:,VarNo)  = unpack(L1_inter_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if

    !
    if (outputFlxState(2)) then

       OutPut(:,:,VarNo)  = unpack(L1_snowPack_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    ! total water amount in soil layers
    if (outputFlxState(3)) then
       do i = 1, nSoilHorizons_mHM

          OutPut(:,:,VarNo) = unpack(L1_soilMoist_out(:,i), mask, nodata_dp)

          V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
          VarNo = VarNo + 1
       end do
    end if
    !
    ! volumetric soil moisture in every soil layer
    if (outputFlxState(4)) then
       do i = 1, nSoilHorizons_mHM
          OutPut(:,:,VarNo) = unpack(L1_soilMoist_out(:,i)    / &
                                     L1_soilMoistSat_out(:,i) , &
                                     mask, nodata_dp )
          !
          V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
          VarNo = VarNo + 1
       end do
       !
    end if
    !
    ! volumetric soil moisture - layer average
    if (outputFlxState(5)) then

       OutPut(:,:,VarNo) = unpack(sum(L1_soilMoist_out(:,:),    dim=2) / &
                                  sum(L1_soilMoistSat_out(:,:), dim=2) , &
                                  mask, nodata_dp)                                 

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1

    end if
    !
    if (outputFlxState(6)) then

       OutPut(:,:,VarNo)  = unpack(L1_sealSTW_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(7)) then

       OutPut(:,:,VarNo)  = unpack(L1_unsatSTW_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(8)) then

       OutPut(:,:,VarNo)  = unpack(L1_satSTW_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(18)) then

       OutPut(:,:,VarNo)  = unpack(L1_neutrons_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if

    ! ---------
    ! 2. Fluxes
    ! ---------
    if (outputFlxState(9)) then
       
       OutPut(:,:,VarNo)  = unpack(L1_pet_out(:), mask, nodata_dp)
       !
       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(10)) then
       
       OutPut(:,:,VarNo)  = unpack(sum(L1_aETSoil_out  (:,:), dim=2) + &
                                       L1_aETCanopy_out(:)           + & 
                                       L1_aETSealed_out(:)           , & 
                                       mask, nodata_dp)
       !
       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(11)) then

       OutPut(:,:,VarNo)  = unpack(L1_total_runoff_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(12)) then

       OutPut(:,:,VarNo)  = unpack(L1_runoffSeal_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(13)) then

       OutPut(:,:,VarNo)  = unpack(L1_fastRunoff_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(14)) then

       OutPut(:,:,VarNo)  = unpack(L1_SlowRunoff_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(15)) then

       OutPut(:,:,VarNo)  = unpack(L1_baseflow_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if
    !
    if (outputFlxState(16)) then

       OutPut(:,:,VarNo)  = unpack(L1_percol_out(:), mask, nodata_dp)

       V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo)
       VarNo = VarNo + 1
    end if

    if (outputFlxState(17)) then 
       do i = 1, nSoilHorizons_mHM 
 
          OutPut(:,:,VarNo) = unpack(L1_infilSoil_out(:,i), mask, nodata_dp) 
          
          V(VarNo+5)%G2_d =>  OutPut(:,:,VarNo) 
          VarNo = VarNo + 1 
       end do
    end if
    
    ! timestep
    ! "-1" to start at timestep 0
    tstep = hours_eval
    V(3)%G0_i => tstep
    call write_dynamic_netCDF(ncid, curr_wstep)
    nullify(V(3)%G0_i)

    !
    VarNo = 0
    do i = 1, size( outputFlxState )
      if (.not. outputFlxState(i)) cycle
      VarNo = VarNo + 1
       nullify ( V(VarNo+5)%G2_d)
    end do

    deallocate(OutPut)
    !
  end subroutine WriteFluxState
  
  
  ! ------------------------------------------------------------------

  !      NAME
  !          CloseFluxState_file

  !>        \brief Close the netcdf file containing flux and states

  !>        \details  Close the netcdf file containing flux and states

  !     INDENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"   - ID number of basin
  !>        \param[in] "integer(i4)    :: ncid_out" - ID of NetCDF file

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !>        \author Rohini Kumar & Stephan Thober
  !>        \date August 2013
  !         Modified, 

 subroutine CloseFluxState_file(iBasin, ncid_out)
    !
    use mo_ncwrite,          only: close_netCDF
    use mo_message,          only: message                        ! For print out
    !
    implicit none
    !
    integer(i4), intent(in)                          :: iBasin    ! mumber of subbasin
    integer(i4), intent(in)                          :: ncid_out  ! ID of NetCDF
    ! local
    character(256)                                   :: dummy
    !
    call message()
    write(dummy,*) iBasin
    call message('  OUTPUT: saved netCDF file for basin', trim(dummy))
    call close_netCDF(ncid_out)
    !
 end subroutine CloseFluxState_file

  ! ------------------------------------------------------------------

  !      NAME
  !          CoordSystem

  !         \brief Setting Gauss Krueger 3 and WGS84 coordinates for output

  !         \details Calculating the coordinates of the projected coordinate system with
  !                  the help of the coordinates of the lower left corner, the cellsize and the number of
  !                  columns and rows. Extracting lats and lons from respective global variables

  !     INTENT(IN)
  !       integer(i4), intent(in) :: iBasin        ! basin index

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !         \author Matthias Zink
  !         \date Apr 2013
  !         Modified, Stephan Thober, Nov 2013 - removed fproj dependency

  subroutine CoordSystem(iBasin) !(xll, yll, cz, nrows, ncols)
    !
    use mo_global_variables, only : & 
         yCoor, xCoor       , & ! kartesian coordinates
         lons, lats         , & ! geographic coordinates
         longitude, latitude, & ! global 1d variables
         level1                 ! level 1 grid information

    !
    implicit none
    ! 
    integer(i4), intent(in) :: iBasin
    
    !
    real(dp)                :: xll, yll      ! coordinates of the lower left corner of
    !  projected coordinate system
    real(dp)                :: cz            ! cellsize of the  projected coordinate system
    integer(i4)             :: nrows, ncols  ! number row and columns of array V accrording
    !
    integer(i4)             :: i, j
    integer(i4)             :: sPos          ! starting position in global array
    !
    xll   = level1%xllcorner(iBasin)
    yll   = level1%yllcorner(iBasin)
    cz    = level1%cellsize(iBasin)
    nrows = level1%nrows(iBasin)
    ncols = level1%ncols(iBasin)
    !
    if (allocated(xCoor) ) deallocate ( xCoor ) ! northing
    if (allocated(yCoor) ) deallocate ( yCoor ) ! easting
    if (allocated(lons))   deallocate ( lons  )
    if (allocated(lats))   deallocate ( lats  )
    allocate ( xCoor(nrows)       )
    allocate ( yCoor(ncols)       ) 
    allocate ( lons(nrows, ncols) )
    allocate ( lats(nrows, ncols) )
    !
    ! def northings and eastings arrays
    xCoor(1)     =  xll + 0.5_dp * cz
    do i = 2, nrows
       xCoor(i)   =  xCoor(i-1) + cz
    end do
    ! inverse for Panoply, ncview display
    yCoor(ncols) =  yll + 0.5_dp * cz
    do j = ncols-1,1,-1 
       yCoor(j)   =  yCoor(j+1) + cz
    end do
    ! calculate position where lats and lons are located in global variables
    sPos = 1
    if ( iBasin .gt. 1 ) then
       do i = 1, iBasin -1 
          sPos = sPos + level1%ncols(i) * level1%nrows(i)
       end do
    end if
    ! extract lats and lons from global field
    lats = reshape( latitude(sPos:sPos+nrows*ncols-1), shape(lats))
    lons = reshape( longitude(sPos:sPos+nrows*ncols-1), shape(lons))
    !
  end subroutine CoordSystem

END MODULE mo_write_fluxes_states



