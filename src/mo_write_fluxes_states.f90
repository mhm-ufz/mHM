!> \file mo_write_fluxes_states.f90

!> \brief Creates NetCDF output for different fluxes and state variabels of mHM.

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!> \authors Matthias Zink
!> \date Apr 2013

MODULE mo_write_fluxes_states

  ! This module creates the output for mHM.

  ! Written Matthias Zink, Apr 2012

  USE mo_kind, ONLY: i4, dp
  USE mo_global_variables, ONLY: outputVariables

  IMPLICIT NONE

  PRIVATE

  public :: writeFluxState
  public :: allocateOutput
  public :: deallocateOutput
  public :: averageStates
  public :: updateOutput
  public :: clearOutput

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
  !>         \param[inout] "real(dp), allocatable :: L11_Qmod(:)"            ! Modelled discharge  

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


  ! function writeTimestep(model_timestep)
  !   integer(i4)  , intent(in) :: model_itimestep
  !   logical :: writeTimestep

  !   if (timeStep_model_outputs .gt. 0) then
  !      if ((mod(model_timestep, timeStep_model_outputs) .eq. 0) .or. (tt .eq. nTimeSteps)) writeout = .true.
  !   else
  !      select case(timeStep_model_outputs)
  !      case(0) ! only at last time step
  !         if (tt .eq. nTimeSteps) writeout = .true.
  !      case(-1) ! daily
  !         if (((model_timestep .gt. 1) .and. (day_counter .ne. day)) .or. (tt .eq. nTimeSteps))     writeout = .true.
  !      case(-2) ! monthly
  !         if (((model_timestep .gt. 1) .and. (month_counter .ne. month)) .or. (tt .eq. nTimeSteps)) writeout = .true.
  !      case(-3) ! yearly
  !         if (((model_timestep .gt. 1) .and. (year_counter .ne. year)) .or. (tt .eq. nTimeSteps))   writeout = .true.
  !      case default ! no output at all
  !         continue
  !      end select
  !   endif


  ! end function writeTimestep


  subroutine averageStates(multiplier,outvars)

    implicit none
    
    real(dp), intent(in) :: multiplier
    type(outputVariables), intent(inout) :: outvars

    if (allocated(outvars%L1_inter))     &
         outvars%L1_inter(:)       = outvars%L1_inter(:)   * multiplier
    if (allocated(outvars%L1_snowPack))   &
         outvars%L1_snowPack       = outvars%L1_snowPack(:)    * multiplier
    if (allocated(outvars%L1_soilMoist)) &
         outvars%L1_soilMoist(:,:) = outvars%L1_soilMoist(:,:) * multiplier
    if (allocated(outvars%L1_sealSTW))   &
         outvars%L1_sealSTW(:)     = outvars%L1_sealSTW(:)     * multiplier
    if (allocated(outvars%L1_unsatSTW))  &
         outvars%L1_unsatSTW(:)    = outvars%L1_unsatSTW(:)    * multiplier
    if (allocated(outvars%L1_satSTW))    &
         outvars%L1_satSTW(:)      = outvars%L1_satSTW(:)      * multiplier
    if (allocated(outvars%L1_neutrons))  &
         outvars%L1_neutrons(:)    = outvars%L1_neutrons(:)    * multiplier
    
  end subroutine averageStates


  subroutine clearOutput(outvars)
    
    implicit none

    type(outputVariables), intent(inout) :: outvars
    !! States L1
    if (allocated(outvars%L1_inter))      outvars%L1_inter(:)        = 0.0_dp       
    if (allocated(outvars%L1_snowPack))   outvars%L1_snowPack(:)     = 0.0_dp      
    if (allocated(outvars%L1_soilMoist))  outvars%L1_soilMoist(:,:)  = 0.0_dp       
    if (allocated(outvars%L1_sealSTW))    outvars%L1_sealSTW(:)      = 0.0_dp      
    if (allocated(outvars%L1_unsatSTW))   outvars%L1_unsatSTW(:)     = 0.0_dp      
    if (allocated(outvars%L1_satSTW))     outvars%L1_satSTW(:)       = 0.0_dp      
    if (allocated(outvars%L1_neutrons))   outvars%L1_neutrons(:)     = 0.0_dp      
    !! Fluxes L1    if (allocated(outvars%L1_pet))        outvars%L1_pet(:)          = 0.0_dp     
    if (allocated(outvars%L1_aETSoil))    outvars%L1_aETSoil(:,:)    = 0.0_dp     
    if (allocated(outvars%L1_aETCanopy))  outvars%L1_aETCanopy(:)    = 0.0_dp    
    if (allocated(outvars%L1_aETSealed))  outvars%L1_aETSealed(:)    = 0.0_dp    
    if (allocated(outvars%L1_total_runoff))      outvars%L1_total_runoff(:) = 0.0_dp   
    if (allocated(outvars%L1_runoffSeal)) outvars%L1_runoffSeal(:)   = 0.0_dp   
    if (allocated(outvars%L1_fastRunoff)) outvars%L1_fastRunoff(:)   = 0.0_dp   
    if (allocated(outvars%L1_slowRunoff)) outvars%L1_slowRunoff(:)   = 0.0_dp   
    if (allocated(outvars%L1_baseflow))   outvars%L1_baseflow(:)     = 0.0_dp   
    if (allocated(outvars%L1_percol))     outvars%L1_percol(:)       = 0.0_dp
    if (allocated(outvars%L1_infilSoil))  outvars%L1_infilSoil(:,:)  = 0.0_dp   

  end subroutine clearOutput

  subroutine allocateOutput(dim1,dim2,outvars)
    
    use mo_global_variables,  only : outputFlxState
    
    implicit none

    integer(i4),           intent(in)    :: dim1,dim2
    type(outputVariables), intent(inout) :: outvars

    !! allocate arrays
    if (outputFlxState(1))  allocate( outvars%L1_inter(dim1) )
    if (outputFlxState(2))  allocate( outvars%L1_snowPack(dim1) )
    if (outputFlxState(3) .or. outputFlxState(4) .or. outputFlxState(5)) allocate( outvars%L1_soilmoist(dim1,dim2) )
    if (outputFlxState(6))  allocate( outvars%L1_sealSTW(dim1) )
    if (outputFlxState(7))  allocate( outvars%L1_unsatSTW(dim1) )
    if (outputFlxState(8))  allocate( outvars%L1_satSTW(dim1) )
    if (outputFlxState(18)) allocate( outvars%L1_neutrons(dim1) )
    if (outputFlxState(9))  allocate( outvars%L1_pet(dim1) )
    if (outputFlxState(10)) then 
       allocate( outvars%L1_aETSoil(dim1, dim2) )
       allocate( outvars%L1_aETCanopy(dim1) )
       allocate( outvars%L1_aETSealed(dim1) )
    end if
    if (outputFlxState(11))  allocate( outvars%L1_total_runoff(dim1) )
    if (outputFlxState(12))  allocate( outvars%L1_runoffSeal(dim1) )
    if (outputFlxState(13))  allocate( outvars%L1_fastRunoff(dim1) )
    if (outputFlxState(14))  allocate( outvars%L1_slowRunoff(dim1) )
    if (outputFlxState(15))  allocate( outvars%L1_baseflow(dim1) )
    if (outputFlxState(16))  allocate( outvars%L1_percol(dim1) )
    if (outputFlxState(17))  allocate( outvars%L1_infilSoil(dim1, dim2))

    !! initialize arrays
    if (outputFlxState(1))  outvars%L1_inter(dim1) = 0.0_dp
    if (outputFlxState(2))  outvars%L1_snowPack(dim1) = 0.0_dp
    if (outputFlxState(3) .or. outputFlxState(4) .or. outputFlxState(5)) outvars%L1_soilmoist(dim1,dim2) = 0.0_dp
    if (outputFlxState(6))  outvars%L1_sealSTW(dim1) = 0.0_dp
    if (outputFlxState(7))  outvars%L1_unsatSTW(dim1) = 0.0_dp
    if (outputFlxState(8))  outvars%L1_satSTW(dim1) = 0.0_dp
    if (outputFlxState(18)) outvars%L1_neutrons(dim1) = 0.0_dp
    if (outputFlxState(9))  outvars%L1_pet(dim1) = 0.0_dp
    if (outputFlxState(10)) then 
       outvars%L1_aETSoil(dim1, dim2) = 0.0_dp
       outvars%L1_aETCanopy(dim1) = 0.0_dp
       outvars%L1_aETSealed(dim1) = 0.0_dp
    end if
    if (outputFlxState(11))  outvars%L1_total_runoff(dim1) = 0.0_dp
    if (outputFlxState(12))  outvars%L1_runoffSeal(dim1) = 0.0_dp
    if (outputFlxState(13))  outvars%L1_fastRunoff(dim1) = 0.0_dp
    if (outputFlxState(14))  outvars%L1_slowRunoff(dim1) = 0.0_dp
    if (outputFlxState(15))  outvars%L1_baseflow(dim1) = 0.0_dp
    if (outputFlxState(16))  outvars%L1_percol(dim1) = 0.0_dp
    if (outputFlxState(17))  outvars%L1_infilSoil(dim1, dim2)= 0.0_dp

  end subroutine allocateOutput


  subroutine deallocateOutput(outvars)
    
    implicit none
    
    type(outputVariables), intent(inout) :: outvars

    if (allocated(outvars%L1_inter))        deallocate(outvars%L1_inter)
    if (allocated(outvars%L1_snowPack))     deallocate(outvars%L1_snowPack)
    if (allocated(outvars%L1_soilmoist))    deallocate(outvars%L1_soilmoist)
    if (allocated(outvars%L1_sealSTW))      deallocate(outvars%L1_sealSTW)
    if (allocated(outvars%L1_unsatSTW))     deallocate(outvars%L1_unsatSTW)
    if (allocated(outvars%L1_satSTW))       deallocate(outvars%L1_satSTW)
    if (allocated(outvars%L1_neutrons))     deallocate(outvars%L1_neutrons)
    if (allocated(outvars%L1_pet))          deallocate(outvars%L1_pet)
    if (allocated(outvars%L1_aETSoil))      deallocate(outvars%L1_aETSoil)
    if (allocated(outvars%L1_aETCanopy))    deallocate(outvars%L1_aETCanopy)
    if (allocated(outvars%L1_aETSealed))    deallocate(outvars%L1_aETSealed)
    if (allocated(outvars%L1_total_runoff)) deallocate(outvars%L1_total_runoff)
    if (allocated(outvars%L1_runoffSeal))   deallocate(outvars%L1_runoffSeal)
    if (allocated(outvars%L1_fastRunoff))   deallocate(outvars%L1_fastRunoff)
    if (allocated(outvars%L1_slowRunoff))   deallocate(outvars%L1_slowRunoff)
    if (allocated(outvars%L1_baseflow))     deallocate(outvars%L1_baseflow)
    if (allocated(outvars%L1_percol))       deallocate(outvars%L1_percol)
    if (allocated(outvars%L1_infilSoil))    deallocate(outvars%L1_infilSoil)                

  end subroutine deallocateOutput

  subroutine updateOutput(&
       L1_inter           , & ! Interception
       L1_snowPack        , & ! Snowpack
       L1_soilMoist       , & ! Soil moisture of each horizon
       L1_sealSTW         , & ! Retention storage of impervious areas
       L1_unsatSTW        , & ! Upper soil storage
       L1_satSTW          , & ! Groundwater storage
       L1_neutrons        , & ! Ground albedo neutrons
       !! Inout: Fluxes L1
       L1_pet             , & ! potential evapotranspiration (PET)
       L1_aETSoil         , & ! actual ET
       L1_aETCanopy       , & ! Real evaporation intensity from canopy
       L1_aETSealed       , & ! Actual ET from free-water surfaces
       L1_total_runoff    , & ! Generated runoff
       L1_runoffSeal      , & ! Direct runoff from impervious areas
       L1_fastRunoff      , & ! Fast runoff component
       L1_slowRunoff      , & ! Slow runoff component
       L1_baseflow        , & ! Baseflow
       L1_percol          , & ! Percolation 
       L1_infilSoil       , & ! Infiltration       
       L1_fSealed         , &
       outvars              &
    )

    use mo_global_variables,  only : nSoilHorizons_mHM    ! Number of horizons to model


    implicit none

    real(dp), dimension(:),   intent(in) :: L1_inter        ! Interception
    real(dp), dimension(:),   intent(in) :: L1_snowPack     ! Snowpack
    real(dp), dimension(:,:), intent(in) :: L1_soilMoist    ! Soil moisture of each horizon
    real(dp), dimension(:),   intent(in) :: L1_sealSTW      ! Retention storage of impervious areas
    real(dp), dimension(:),   intent(in) :: L1_unsatSTW     ! Upper soil storage
    real(dp), dimension(:),   intent(in) :: L1_satSTW       ! Groundwater storage
    real(dp), dimension(:),   intent(in) :: L1_neutrons     ! ground albedo neutrons
    ! Fluxes L1
    real(dp), dimension(:),   intent(in) :: L1_pet          ! potential evapotranspiration (PET)
    real(dp), dimension(:,:), intent(in) :: L1_aETSoil      ! actual ET of each horizon
    real(dp), dimension(:),   intent(in) :: L1_aETCanopy    ! Real evaporation intensity from canopy
    real(dp), dimension(:),   intent(in) :: L1_aETSealed    ! Actual ET from free-water surfaces
    real(dp), dimension(:),   intent(in) :: L1_total_runoff ! Generated runoff
    real(dp), dimension(:),   intent(in) :: L1_runoffSeal   ! Direct runoff from impervious areas
    real(dp), dimension(:),   intent(in) :: L1_fastRunoff   ! Fast runoff component
    real(dp), dimension(:),   intent(in) :: L1_slowRunoff   ! Slow runoff component
    real(dp), dimension(:),   intent(in) :: L1_baseflow     ! Baseflow
    real(dp), dimension(:),   intent(in) :: L1_percol       ! Percolation
    real(dp), dimension(:,:), intent(in) :: L1_infilSoil    ! Infiltration

    real(dp), dimension(:),   intent(in) :: L1_fSealed       
 
    type(outputVariables), intent(inout) :: outvars    

    !! local
    integer(i4) :: hh

    !! States L1 --> AVERAGE
    if (allocated(outvars%L1_inter     )) outvars%L1_inter         = outvars%L1_inter(:)  + L1_inter
    if (allocated(outvars%L1_snowPack  )) outvars%L1_snowPack      = outvars%L1_snowPack  + L1_snowPack
    if (allocated(outvars%L1_soilMoist )) outvars%L1_soilMoist     = outvars%L1_soilMoist + L1_soilMoist
    if (allocated(outvars%L1_sealSTW   )) outvars%L1_sealSTW       = outvars%L1_sealSTW   + L1_sealSTW
    if (allocated(outvars%L1_unsatSTW  )) outvars%L1_unsatSTW      = outvars%L1_unsatSTW  + L1_unsatSTW
    if (allocated(outvars%L1_satSTW    )) outvars%L1_satSTW        = outvars%L1_satSTW    + L1_satSTW
    if (allocated(outvars%L1_neutrons  )) outvars%L1_neutrons      = outvars%L1_neutrons  + L1_neutrons
    
    !! Fluxes L1  --> AGGREGATED
    if (allocated(outvars%L1_pet)) outvars%L1_pet = L1_pet

    if (allocated(outvars%L1_aETSoil)) then
       ! is this loop necessray??
       do hh = 1, nSoilHorizons_mHM
          outvars%L1_aETSoil(:,hh)  = outvars%L1_aETSoil(:,hh) + L1_aETSoil(:,hh)*(1.0_dp - L1_fSealed)
       end do
    end if

    if (allocated(outvars%L1_aETCanopy    ))  outvars%L1_aETCanopy    = outvars%L1_aETCanopy + L1_aETCanopy
    if (allocated(outvars%L1_aETSealed    ))  outvars%L1_aETSealed    = outvars%L1_aETSealed + L1_aETSealed * L1_fSealed
    if (allocated(outvars%L1_total_runoff ))  outvars%L1_total_runoff = outvars%L1_total_runoff + L1_total_runoff
    if (allocated(outvars%L1_runoffSeal   ))  outvars%L1_runoffSeal   = outvars%L1_runoffSeal + L1_runoffSeal*L1_fSealed
    if (allocated(outvars%L1_fastRunoff   ))  outvars%L1_fastRunoff   = outvars%L1_fastRunoff + L1_fastRunoff*(1.0_dp - L1_fSealed)
    if (allocated(outvars%L1_slowRunoff   ))  outvars%L1_slowRunoff   = outvars%L1_slowRunoff + L1_slowRunoff*(1.0_dp - L1_fSealed)
    if (allocated(outvars%L1_baseflow     ))  outvars%L1_baseflow     = outvars%L1_baseflow + L1_baseflow*(1.0_dp - L1_fSealed)
    if (allocated(outvars%L1_percol       ))  outvars%L1_percol       = outvars%L1_percol + L1_percol*(1.0_dp - L1_fSealed)

    if (allocated(outvars%L1_infilSoil)) then
       ! is this loop necessray??
       do hh = 1, nSoilHorizons_mHM
          outvars%L1_infilSoil(:,hh)                 = outvars%L1_infilSoil(:,hh) + L1_infilSoil(:,hh)*(1.0_dp - L1_fSealed) 
       end do
    end if
    
  end subroutine updateOutput

  subroutine writeFluxState(&
       !! Input
       iTimestep              , &
       iBasin                 , &
       mask1                  , &
       L1_soilMoistSat        , &
       outvars                  &
       )
    
    use mo_global_variables,  only : & 
         dirOut            , & ! output directory
         outputFlxState        ! definition which output to write

    use mo_mhm_constants,      only: nodata_dp             ! global no data value
    
    use mo_ncwrite, only: var2nc


    implicit none
    
    !
    integer(i4),                           intent(in)    :: iTimestep    ! mumber of subbasin
    integer(i4),                           intent(in)    :: iBasin       ! mumber of subbasin
    real(dp), dimension(:,:),              intent(in)    :: L1_soilMoistSat
    logical, dimension(:,:),               intent(in)    :: mask1
    type(outputVariables),                 intent(in)    :: outvars
     
    !! local
    character(256)                                       :: fname
    character(256),dimension(3)                          :: dims1, dims11
    character(256),dimension(6,2)                        :: attrs
    logical                                              :: createvar
    integer(i4)                                          :: ii 
    real(dp), dimension(1,size(mask1,1),size(mask1,2))   :: tmp1
    integer(i4), dimension(1461)                         :: tmptime
    !!
    dims1(1) = "time"
    dims1(2) = "northing"
    dims1(3) = "easting"
    ! dims11(1) = "time"
    ! dims11(2) = "y11"
    ! dims11(3) = "x11"

    attrs(1,1) = "name"
    attrs(2,1) = "unit"
    attrs(3,1) = "long_name"
    attrs(4,1) = "scale_factor"
    attrs(4,2) = "1."
    attrs(5,1) = "missing_value"
    attrs(5,2) = "-9999."
    attrs(6,1) = "coordinates"
    attrs(6,2) = "lat lon"
    
    fname = trim(dirOut(iBasin)) // 'mHM_Fluxes_States.nc'
    createvar = iTimestep .eq. 1

    if (createvar) then
       tmptime=0
       call var2nc(fname, tmptime, (/ "time" /), "time", create=createvar)
    end if

    !! ---------
    !! 2. States
    !! ---------
    if (outputFlxState(1)) then
       attrs(1,2) = "interception"
       attrs(2,2) = "2"
       attrs(3,2) = "canopy interception storage"
       tmp1(1,:,:) = unpack(outvars%L1_inter, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value = nodata_dp)
    end if

    if (outputFlxState(2)) then
       attrs(1,2) = "snowpack"
       attrs(2,2) = "mm"
       attrs(3,2) = "depth of snowpack"      
       tmp1(1,:,:) = unpack(outvars%L1_snowPack, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value = nodata_dp)
    end if

    if (outputFlxState(3)) then
       attrs(2,2) = "mm"
       do ii = 1, size(outvars%L1_soilMoist,2)
          write(attrs(1,2), "('SWC_L', i2.2)") ii
          write(attrs(3,2), "('soil water content of soil layer',i2)") ii
          tmp1(1,:,:) = unpack(outvars%L1_soilMoist(:,ii), mask1, nodata_dp)
          call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)
       end do
    end if

    if (outputFlxState(4)) then
       attrs(2,2) = "mm"
       do ii = 1, size(outvars%L1_soilMoist,2)
          write(attrs(1,2), "('SM_L', i2.2)") ii
          write(attrs(3,2), "('volumetric soil moisture of soil layer',i2)") ii
          tmp1(1,:,:) = unpack(outvars%L1_soilMoist(:,ii)/L1_soilMoistSat(:,ii), mask1, nodata_dp)
          call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
       end do
    end if

    if (outputFlxState(5)) then
       attrs(1,2) = "SM_Lall"
       attrs(2,2) = "mm mm-1"
       attrs(2,2) = "average soil moisture over all layers"
       tmp1(1,:,:) = unpack(sum(outvars%L1_soilMoist,dim=2)/sum(L1_soilMoistSat,dim=2), mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(6)) then
       attrs(1,2) = "sealedSTW"
       attrs(2,2) = "mm"
       attrs(2,2) = "reservoir of sealed areas (sealedSTW)"
       tmp1(1,:,:) = unpack(outvars%L1_sealSTW, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)
    end if

    if (outputFlxState(7)) then
       attrs(1,2) = "unsatSTW"
       attrs(2,2) = "mm"
       attrs(2,2) = "reservoir of unsaturated zone"
       tmp1(1,:,:) = unpack(outvars%L1_unsatSTW, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(8)) then
       attrs(1,2) = "satSTW"
       attrs(2,2) = "mm"
       attrs(2,2) = "water level in groudwater reservoir"
       tmp1(1,:,:) = unpack(outvars%L1_satSTW, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(18)) then
       attrs(1,2) = "Neutrons"
       attrs(2,2) = "cph"
       attrs(2,2) = "ground albedo neutrons"
       tmp1(1,:,:) = unpack(outvars%L1_neutrons, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    ! ! ---------
    ! ! 2. Fluxes
    ! ! ---------
    attrs(2,2) = trim(fluxesUnit(iBasin,iTimestep))

    if (outputFlxState(9)) then
       attrs(1,2) = "PET"
       attrs(2,2) = "potential Evapotranspiration"
       tmp1(1,:,:) = unpack(outvars%L1_pet, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)
    end if

    if (outputFlxState(10)) then
       attrs(1,2) = "aET"
       attrs(2,2) = "actual Evapotranspiration"
       tmp1(1,:,:) =  unpack(sum(outvars%L1_aetSoil, dim=2) + outvars%L1_aETCanopy + outvars%L1_aETSealed, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)
    end if

    if (outputFlxState(11)) then
       attrs(1,2) = "Q"
       attrs(2,2) = "total runoff generated by every cell"
       tmp1(1,:,:) = unpack(outvars%L1_total_runoff, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if
    
    if (outputFlxState(12)) then
       attrs(1,2) = "QD"
       attrs(2,2) = "direct runoff generated by every cell (runoffSeal)"
       tmp1(1,:,:) = unpack(outvars%L1_runoffSeal, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(13)) then
       attrs(1,2) = "QIf"
       attrs(2,2) = "fast interflow generated by every cell (fastRunoff)"
       tmp1(1,:,:) = unpack(outvars%L1_fastRunoff, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(14)) then
       attrs(1,2) = "QIs"
       attrs(2,2) = "slow interflow generated by every cell (slowRunoff)"
       tmp1(1,:,:) = unpack(outvars%L1_SlowRunoff, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if
    
    if (outputFlxState(15)) then
       attrs(1,2) = "QB"
       attrs(2,2) = "baseflow generated by every cell"
       tmp1(1,:,:) = unpack(outvars%L1_baseflow, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(16)) then
       attrs(1,2) = "recharge"
       attrs(2,2) = "groundwater recharge"
       tmp1(1,:,:) = unpack(outvars%L1_percol, mask1, nodata_dp)
       call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)       
    end if

    if (outputFlxState(17)) then
       do ii = 1, size(outvars%L1_infilSoil,2)
          write(attrs(1,2), "('soil_infil_L', i2.2)") ii 
          write(attrs(3,2), "('infiltration flux from soil layer',i2)") ii
          tmp1(1,:,:) = unpack(outvars%L1_infilSoil(:,ii), mask1, nodata_dp)
          call var2nc(fname, tmp1, dims1, attrs(1,2), attributes=attrs, create=createvar, dim_unlimited=1, missing_value=nodata_dp)
       end do
    end if

  end subroutine writeFluxState

  function fluxesUnit(ibasin,itimestep)

    use mo_global_variables,  only : timestep, simPer, NTSTEPDAY
    use mo_string_utils,      only: num2str

    integer(i4), intent(in) :: ibasin    
    integer(i4), intent(in) :: itimestep
    ! local
    character(256)          :: fluxesUnit
    real(dp)             :: nTimeSteps


    nTimeSteps = ( simPer(ibasin)%julEnd - simPer(ibasin)%julStart + 1 ) * NTSTEPDAY

    if (itimestep > 0) then
       if ((timestep*itimestep) .eq. 1 ) then
          fluxesUnit = 'mm h-1'
       else
          fluxesUnit = 'mm '//trim(adjustl(num2str(timestep)))//'h-1'
       end if
    else
       select case(itimestep)
       case(0) ! only at last time step
          fluxesUnit = 'mm '//trim(adjustl(num2str(nint(nTimeSteps))))//'h-1'
       case(-1) ! daily
          fluxesUnit = 'mm d-1'
       case(-2) ! monthly
          fluxesUnit = 'mm month-1'
       case(-3) ! yearly
          fluxesUnit = 'mm a-1'
       case default ! no output at all
          fluxesUnit = ''
       end select
    endif

  end function fluxesUnit

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

