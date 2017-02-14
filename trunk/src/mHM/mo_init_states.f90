!> \file mo_init_states.f90

!> \brief Initialization of all state variables of mHM.

!> \details This module initializes all state variables required to run mHM.
!> Two options are provided:
!>     - (1) default values
!>     - (2) from nc file

!> \author Luis Samaniego & Rohini Kumar
!> \date Dec 2012

MODULE mo_init_states

  ! This module provides the startup routines for mHM.

  ! Written Luis Samaniego & Rohini Kumar, Dec 2012

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: variables_alloc                 ! allocation of space for state variables/fluxes/effective parameters
  PUBLIC :: variables_default_init          ! initialization for state variables/fluxes/effective parameters

  PUBLIC :: get_basin_info                  ! get basin info
  PUBLIC :: calculate_grid_properties       ! generate Lx (x = 1, 11, 2) grid properties using L0 grid

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          variables_alloc

  !>        \brief Allocation of space for mHM related L1 and L11 variables.

  !>        \details Allocation of space for mHM related L1 and L11 variables (e.g., states,
  !>                 fluxes, and parameters) for a given basin. Variables allocated here is
  !>                 defined in them mo_global_variables.f90 file. After allocating any variable
  !>                 in this routine, initalize them in the following variables_default_init
  !>                 subroutine:
  !>
  !
  !     CALLING SEQUENCE
  !         call variables_alloc(iBasin)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: iBasin"        - basin id

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL

  !     RETURN

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date Jan 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template
  !                   S. Thober, Aug 2015  - removed routing related variables

  subroutine variables_alloc(iBasin)

    use mo_global_variables, only: nSoilHorizons_mHM,            &
         L1_fSealed, L1_fForest, L1_fPerm, L1_inter, L1_snowPack, L1_sealSTW,   &
         L1_soilMoist, L1_unsatSTW, L1_satSTW,                                  &
         L1_pet_calc, L1_aETSoil, L1_aETCanopy, L1_aETSealed,                   &
         L1_baseflow, L1_infilSoil, L1_fastRunoff, L1_melt,                     &
         L1_percol, L1_preEffect, L1_rain, L1_runoffSeal, L1_slowRunoff,        &
         L1_snow, L1_Throughfall, L1_total_runoff, L1_alpha, L1_degDayInc,      &
         L1_degDayMax, L1_degDayNoPre, L1_degDay, L1_karstLoss, L1_fAsp,        &
         L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_surfResist,         &
         L1_fRoots, L1_maxInter, L1_kfastFlow, L1_kSlowFlow, L1_kBaseFlow,      &
         L1_kPerco, L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp,           &
         L1_tempThresh, L1_unsatThresh, L1_sealedThresh, L1_wiltingPoint,       &
         L1_neutrons

    use mo_mhm_constants,    only: YearMonths_i4
    use mo_append,           only: append                      ! append vector

    implicit none

    integer(i4), intent(in)                   :: iBasin

    ! local variables
    integer(i4)                               :: nrows1, ncols1
    integer(i4)                               :: ncells1
    real(dp), dimension(:),   allocatable     :: dummy_Vector
    real(dp), dimension(:,:), allocatable     :: dummy_Matrix
    real(dp), dimension(:,:), allocatable     :: dummy_Matrix_months

    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1 )

    ! for appending and intialization
    allocate( dummy_Vector        (nCells1)                        )
    allocate( dummy_Matrix        (nCells1, nSoilHorizons_mHM)     )
    allocate( dummy_Matrix_months (nCells1,     YearMonths_i4)     )

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------
    dummy_Vector(:) = 0.0_dp
    call append( L1_fSealed,  dummy_Vector )
    call append( L1_fForest,  dummy_Vector )
    call append( L1_fPerm,    dummy_Vector )

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    dummy_Vector(:) = 0.0_dp
    call append( L1_inter,  dummy_Vector )

    !Snowpack
    dummy_Vector(:) = 0.0_dp
    call append( L1_snowPack,  dummy_Vector )

    !Retention storage of impervious areas
    dummy_Vector(:) = 0.0_dp
    call append( L1_sealSTW,  dummy_Vector )

    ! Soil moisture of each horizon
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_soilMoist, dummy_Matrix )

    ! upper soil storage
    dummy_Vector(:) = 0.0_dp
    call append( L1_unsatSTW,  dummy_Vector )

    ! groundwater storage
    dummy_Vector(:) = 0.0_dp
    call append( L1_satSTW,  dummy_Vector )

    ! ground albedo neutrons
    dummy_Vector(:) = 0.0_dp
    call append( L1_neutrons,  dummy_Vector )

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------

    ! calculated / corrected potential evapotranspiration
    dummy_Vector(:) = 0.0_dp
    call append( L1_pet_calc,  dummy_Vector )

    !  soil actual ET
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_aETSoil,  dummy_Matrix )

    ! canopy actual ET
    dummy_Vector(:) = 0.0_dp
    call append( L1_aETCanopy,  dummy_Vector )

    ! sealed area actual ET
    dummy_Vector(:) = 0.0_dp
    call append( L1_aETSealed,  dummy_Vector )

    ! baseflow
    dummy_Vector(:) = 0.0_dp
    call append( L1_baseflow,  dummy_Vector )

    !  soil in-exfiltration
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_infilSoil,  dummy_Matrix )

    ! fast runoff
    dummy_Vector(:) = 0.0_dp
    call append( L1_fastRunoff,  dummy_Vector )

    ! snow melt
    dummy_Vector(:) = 0.0_dp
    call append( L1_melt,  dummy_Vector )

    ! percolation
    dummy_Vector(:) = 0.0_dp
    call append( L1_percol,  dummy_Vector )

    ! effective precip. depth (snow melt + rain)
    dummy_Vector(:) = 0.0_dp
    call append( L1_preEffect,  dummy_Vector )

    ! rain (liquid water)
    dummy_Vector(:) = 0.0_dp
    call append( L1_rain,  dummy_Vector )

    ! runoff from impervious area
    dummy_Vector(:) = 0.0_dp
    call append( L1_runoffSeal,  dummy_Vector )

    ! slow runoff
    dummy_Vector(:) = 0.0_dp
    call append( L1_slowRunoff,  dummy_Vector )

    ! snow (solid water)
    dummy_Vector(:) = 0.0_dp
    call append( L1_snow,  dummy_Vector )

    ! throughfall
    dummy_Vector(:) = 0.0_dp
    call append( L1_Throughfall,  dummy_Vector )

    ! throughfall
    dummy_Vector(:) = 0.0_dp
    call append( L1_total_runoff,  dummy_Vector )

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    dummy_Vector(:) = 0.0_dp
    call append( L1_alpha,  dummy_Vector )

    ! increase of the Degree-day factor per mm of increase in precipitation
    dummy_Vector(:) = 0.0_dp
    call append( L1_degDayInc,  dummy_Vector )

    ! maximum degree-day factor
    dummy_Vector(:) = 0.0_dp
    call append( L1_degDayMax,  dummy_Vector )

    ! degree-day factor with no precipitation
    dummy_Vector(:) = 0.0_dp
    call append( L1_degDayNoPre,  dummy_Vector )

    ! degree-day factor
    dummy_Vector(:) = 0.0_dp
    call append( L1_degDay,  dummy_Vector )

    ! Karstic percolation loss
    dummy_Vector(:) = 0.0_dp
    call append( L1_karstLoss,  dummy_Vector )

    ! PET correction factor due to terrain aspect
    dummy_Vector(:) = 0.0_dp
    call append( L1_fAsp,  dummy_Vector )

    ! PET Hargreaves Samani coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_HarSamCoeff,   dummy_Vector )

    ! PET Prietley Taylor coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_PrieTayAlpha,  dummy_Matrix_months )

    ! PET aerodynamical resistance
    dummy_Matrix_months = 0.0_dp
    call append( L1_aeroResist,   dummy_Matrix_months )

    ! PET bulk surface resistance
    dummy_Matrix_months = 0.0_dp
    call append( L1_surfResist,   dummy_Matrix_months )

    ! Fraction of roots in soil horizons
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_fRoots,  dummy_Matrix )

    ! Maximum interception
    dummy_Vector(:) = 0.0_dp
    call append( L1_maxInter,  dummy_Vector(:) )

    ! fast interflow recession coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_kfastFlow,  dummy_Vector )

    ! slow interflow recession coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_kSlowFlow,  dummy_Vector )

    ! baseflow recession coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_kBaseFlow,  dummy_Vector )

    ! percolation coefficient
    dummy_Vector(:) = 0.0_dp
    call append( L1_kPerco,  dummy_Vector )

    ! Soil moisture below which actual ET is reduced linearly till PWP
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_soilMoistFC,  dummy_Matrix )

    ! Saturation soil moisture for each horizon [mm]
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_soilMoistSat,  dummy_Matrix )

    ! Exponential parameter to how non-linear is the soil water retention
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_soilMoistExp,  dummy_Matrix )

    ! Threshold temperature for snow/rain
    dummy_Vector(:) = 0.0_dp
    call append( L1_tempThresh,  dummy_Vector )

    ! Threshhold water depth controlling fast interflow
    dummy_Vector(:) = 0.0_dp
    call append( L1_unsatThresh,  dummy_Vector )

    ! Threshhold water depth for surface runoff in sealed surfaces
    dummy_Vector(:) = 0.0_dp
    call append( L1_sealedThresh,  dummy_Vector )

    ! Permanent wilting point
    dummy_Matrix(:,:) = 0.0_dp
    call append( L1_wiltingPoint,  dummy_Matrix )

    ! free space
    if ( allocated( dummy_Vector          ) ) deallocate( dummy_Vector          )
    if ( allocated( dummy_Matrix          ) ) deallocate( dummy_Matrix          )
    if ( allocated( dummy_Matrix_months   ) ) deallocate( dummy_Matrix_months   )

  end subroutine variables_alloc

  ! ------------------------------------------------------------------

  !      NAME
  !          variables_default_init

  !>        \brief Default initalization mHM related L1 variables

  !>        \details Default initalization of mHM related L1 variables (e.g., states,
  !>                 fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !>                 Variables initalized here is defined in the mo_global_variables.f90 file.
  !>                 Only Variables that are defined in the variables_alloc subroutine are
  !>                 intialized here.
  !
  !>                 If a variable is added or removed here, then it also has to be added or removed
  !>                 in the subroutine state_variables_set in the module mo_restart and in the
  !>                 subroutine set_state in the module mo_set_netcdf_restart.
  !
  !     CALLING SEQUENCE
  !         call variables_default_init()

  !     INTENT(IN)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL

  !     RETURN

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         \authors  R. Kumar & J. Mai
  !         \date    Sep. 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template
  !              Stephan Thober Aug 2015   - moved routing variables to mRM

  subroutine variables_default_init()

    use mo_global_variables, only:                                              &
         nSoilHorizons_mHM, HorizonDepth_mHM,                                   &
         L1_fSealed, L1_fForest, L1_fPerm, L1_inter, L1_snowPack, L1_sealSTW,   &
         L1_soilMoist, L1_unsatSTW, L1_satSTW,                                  &
         L1_pet_calc, L1_aETSoil, L1_aETCanopy, L1_aETSealed,                   &
         L1_baseflow, L1_infilSoil, L1_fastRunoff, L1_melt,                     &
         L1_percol, L1_preEffect, L1_rain, L1_runoffSeal, L1_slowRunoff,        &
         L1_snow, L1_Throughfall, L1_total_runoff, L1_alpha, L1_degDayInc,      &
         L1_degDayMax, L1_degDayNoPre, L1_degDay, L1_karstLoss, L1_fAsp,        &
         L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_surfResist,         &
         L1_fRoots, L1_maxInter, L1_kfastFlow, L1_kSlowFlow, L1_kBaseFlow,      &
         L1_kPerco, L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp,           &
         L1_tempThresh, L1_unsatThresh, L1_sealedThresh, L1_wiltingPoint,       &
         L1_neutrons

    use mo_mhm_constants,    only:               &
         P1_InitStateFluxes, P2_InitStateFluxes, &
         P3_InitStateFluxes, P4_InitStateFluxes, &
         P5_InitStateFluxes, C1_InitStateSM

    implicit none

    ! inital values
    integer(i4)                               :: i

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------
    L1_fSealed = P1_InitStateFluxes
    L1_fForest = P1_InitStateFluxes
    L1_fPerm   = P1_InitStateFluxes

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    L1_inter = P1_InitStateFluxes

    !Snowpack
    L1_snowPack = P2_InitStateFluxes

    !Retention storage of impervious areas
    L1_sealSTW = P1_InitStateFluxes

    ! Soil moisture of each horizon
    do i = 1, nSoilHorizons_mHM-1
       if (i .eq. 1) then
          L1_soilMoist(:,i) = HorizonDepth_mHM(i)*C1_InitStateSM
       else
          L1_soilMoist(:,i)  = ( HorizonDepth_mHM(i) - HorizonDepth_mHM(i-1) )*C1_InitStateSM
       end if
    end do
    L1_soilMoist(:,nSoilHorizons_mHM) = ( P5_InitStateFluxes - &
         HorizonDepth_mHM(nSoilHorizons_mHM-1) )*C1_InitStateSM

    ! upper soil storage
    L1_unsatSTW = P3_InitStateFluxes

    ! groundwater storage
    L1_satSTW = P4_InitStateFluxes

    ! ground albedo neutrons, initially zero
    L1_neutrons = P1_InitStateFluxes

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------

    ! corrected / calculated potential ET
    L1_pet_calc = P1_InitStateFluxes

    !  soil actual ET
    L1_aETSoil = P1_InitStateFluxes

    ! canopy actual ET
    L1_aETCanopy = P1_InitStateFluxes

    ! sealed area actual ET
    L1_aETSealed = P1_InitStateFluxes

    ! baseflow
    L1_baseflow = P1_InitStateFluxes

    !  soil in-exfiltration
    L1_infilSoil = P1_InitStateFluxes

    ! fast runoff
    L1_fastRunoff = P1_InitStateFluxes

    ! snow melt
    L1_melt = P1_InitStateFluxes

    ! percolation
    L1_percol = P1_InitStateFluxes

    ! effective precip. depth (snow melt + rain)
    L1_preEffect = P1_InitStateFluxes

    ! rain (liquid water)
    L1_rain = P1_InitStateFluxes

    ! runoff from impervious area
    L1_runoffSeal = P1_InitStateFluxes

    ! slow runoff
    L1_slowRunoff = P1_InitStateFluxes

    ! snow (solid water)
    L1_snow = P1_InitStateFluxes

    ! throughfall
    L1_Throughfall = P1_InitStateFluxes

    ! throughfall
    L1_total_runoff = P1_InitStateFluxes

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    L1_alpha = P1_InitStateFluxes

    ! increase of the Degree-day factor per mm of increase in precipitation
    L1_degDayInc = P1_InitStateFluxes

    ! maximum degree-day factor
    L1_degDayMax = P1_InitStateFluxes

    ! degree-day factor with no precipitation
    L1_degDayNoPre = P1_InitStateFluxes

    ! degree-day factor
    L1_degDay = P1_InitStateFluxes

    ! Karstic percolation loss
    L1_karstLoss = P1_InitStateFluxes

    ! PET correction factor due to terrain aspect
    L1_fAsp = P1_InitStateFluxes

    ! PET Hargreaves Samani Coefficient
    L1_HarSamCoeff = P1_InitStateFluxes

    ! PET Priestley Taylor coefficient
    L1_PrieTayAlpha = P1_InitStateFluxes

    ! PET aerodynamical resistance
    L1_aeroResist = P1_InitStateFluxes

    ! PET bulk surface resistance
    L1_surfResist = P1_InitStateFluxes


    ! Fraction of roots in soil horizons
    L1_fRoots = P1_InitStateFluxes

    ! Maximum interception
    L1_maxInter = P1_InitStateFluxes

    ! fast interflow recession coefficient
    L1_kfastFlow = P1_InitStateFluxes

    ! slow interflow recession coefficient
    L1_kSlowFlow = P1_InitStateFluxes

    ! baseflow recession coefficient
    L1_kBaseFlow = P1_InitStateFluxes

    ! percolation coefficient
    L1_kPerco = P1_InitStateFluxes

    ! Soil moisture below which actual ET is reduced linearly till PWP
    L1_soilMoistFC = P1_InitStateFluxes

    ! Saturation soil moisture for each horizon [mm]
    L1_soilMoistSat = P1_InitStateFluxes

    ! Exponential parameter to how non-linear is the soil water retention
    L1_soilMoistExp = P1_InitStateFluxes

    ! Threshold temperature for snow/rain
    L1_tempThresh = P1_InitStateFluxes

    ! Threshhold water depth controlling fast interflow
    L1_unsatThresh = P1_InitStateFluxes

    ! Threshhold water depth for surface runoff in sealed surfaces
    L1_sealedThresh = P1_InitStateFluxes

    ! Permanent wilting point
    L1_wiltingPoint = P1_InitStateFluxes

  end subroutine variables_default_init
  ! ------------------------------------------------------------------

  !      NAME
  !          get_basin_info

  !>        \brief Get basic basin information (e.g., nrows, ncols, indices, mask)

  !>        \details Get basic basin information (e.g., nrows, ncols, indices, mask) for
  !>                 different levels (L0, L1, and L2).
  !
  !     CALLING SEQUENCE
  !         call get_basin_info(iBasin, iLevel,nrows,ncols, ncells, iStart, iEnd, &
  !                             iStartMask, iEndMask, mask, xllcorner, yllcorner, cellsize)

  !     INTENT(IN)
  !>        \param[in] "integer(i4)             :: iBasin"    basin id
  !>        \param[in] "integer(i4)             :: iLevel"    level id (e.g., 0, 1, 11, 2)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)             :: nRows"      no. of rows
  !>        \param[out] "integer(i4)             :: nCols"      no. of coloums
  !>        \param[out] "integer(i4)             :: ncells"     no. of cells
  !>        \param[out] "integer(i4)             :: iStart"     start cell index of a given basin at a given level
  !>        \param[out] "integer(i4)             :: iEnd"       end cell index of a given basin at a given level
  !>        \param[out] "integer(i4)             :: iStartMask" start cell index of mask a given basin at a given level
  !>        \param[out] "integer(i4)             :: iEndMask"   end cell index of mask a given basin at a given level

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "logical, optional       :: mask"       Mask at a given level
  !>        \param[out] "real(dp), optional      :: xllcorner"  x coordinate of the lowerleft corner at a given level
  !>        \param[out] "real(dp), optional      :: yllcorner"  y coordinate of the lowerleft corner at a given level
  !>        \param[out] "real(dp), optional      :: cellsize"   cell size at a given level

  !     RETURN

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         \authors  Rohini Kumar, Luis Samaniego
  !         \date     Jan 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template
  !                   Stephan Thober, Aug 2015 - moved L11 and L110 to mRM

  subroutine get_basin_info(iBasin, iLevel, nrows, ncols, ncells, iStart, iEnd, &
                            iStartMask, iEndMask, mask, xllcorner, yllcorner, cellsize)

    use mo_message, only: message
    use mo_global_variables, only: basin, level0, level1, level2
    implicit none

    integer(i4), intent(in)                                      :: iBasin
    integer(i4), intent(in)                                      :: iLevel
    integer(i4), intent(out)                                     :: nrows, ncols
    integer(i4), optional, intent(out)                           :: ncells
    integer(i4), optional, intent(out)                           :: iStart, iEnd
    integer(i4), optional, intent(out)                           :: iStartMask, iEndMask
    logical, optional, dimension(:,:), allocatable,  intent(out) :: mask
    real(dp), optional, intent(out)                              :: xllcorner, yllcorner, cellsize

    ! level information
    select case (iLevel)
    case (0)
       nrows = level0%nrows(iBasin)
       ncols = level0%ncols(iBasin)
       if (present(ncells)) ncells = basin%L0_iEnd(iBasin) - basin%L0_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin%L0_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin%L0_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin%L0_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin%L0_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin%L0_Mask( basin%L0_iStartMask(iBasin): basin%L0_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level0%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level0%yllcorner(iBasin)
       if (present(cellsize))  cellsize  = level0%cellsize(iBasin)

    case (1)
       nrows = level1%nrows(iBasin)
       ncols = level1%ncols(iBasin)
       if (present(ncells)) ncells = basin%L1_iEnd(iBasin) - basin%L1_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin%L1_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin%L1_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin%L1_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin%L1_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin%L1_Mask( basin%L1_iStartMask(iBasin): basin%L1_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level1%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level1%yllcorner(iBasin)
       if (present(cellsize))  cellsize = level1%cellsize(iBasin)

    case (11)
       call message('***ERROR: get_basin_info has been called for level 11 that does not exist within mHM')
       call message('***ERROR: use get_basin_info_mrm from mRM instead')
       stop

    case (110)
       call message('***ERROR: get_basin_info has been called for level 110 that does not exist within mHM')
       call message('***ERROR: use get_basin_info_mrm from mRM instead')
       stop

    case (2)
       nrows = level2%nrows(iBasin)
       ncols = level2%ncols(iBasin)
       if (present(ncells)) ncells = basin%L2_iEnd(iBasin) - basin%L2_iStart(iBasin) + 1
       if (present(iStart)) iStart = basin%L2_iStart(iBasin)
       if (present(iEnd))   iEnd   = basin%L2_iEnd(iBasin)
       if (present(iStartMask)) iStartMask = basin%L2_iStartMask(iBasin)
       if (present(iEndMask))   iEndMask   = basin%L2_iEndMask(iBasin)
       if (present(Mask)) then
          allocate ( mask(nrows, ncols) )
          mask(:,:) = .FALSE.
          mask(:,:) = RESHAPE( basin%L2_Mask( basin%L2_iStartMask(iBasin): basin%L2_iEndMask(iBasin)),&
               (/nrows,ncols/) )
       end if
       if (present(xllcorner)) xllcorner = level2%xllcorner(iBasin)
       if (present(yllcorner)) yllcorner = level2%yllcorner(iBasin)
       if (present(cellsize)) cellsize   = level2%cellsize(iBasin)

    end select

  end subroutine get_basin_info

  ! ------------------------------------------------------------------

  !      NAME
  !         calculate_grid_properties

  !     PURPOSE
  !>        \brief Calculates basic grid properties at a required coarser level using
  !>              information of a given finer level.

  !>        \brief Calculates basic grid properties at a required coarser level (e.g., L11) using
  !>              information of a given finer level (e.g., L0). Basic grid properties such as
  !>              nrows, ncols, xllcorner, yllcorner cellsize are estimated in this
  !>              routine.

  !     CALLING SEQUENCE
  !         call calculate_grid_properties( nrowsIn, ncolsIn,  xllcornerIn,                     &
  !                                         yllcornerIn,  cellsizeIn, nodata_valueIn,           &
  !                                         aimingResolution, nrowsOut, ncolsOut, xllcornerOut, &
  !                                         yllcornerOut, cellsizeOut, nodata_valueOut )
  !     INTENT(IN)
  !>        \param[in] "integer(i4)             :: nrowsIn"           no. of rows at an input level
  !>        \param[in] "integer(i4)             :: ncolsIn"           no. of cols at an input level
  !>        \param[in] "real(dp)                :: xllcornerIn"       xllcorner at an input level
  !>        \param[in] "real(dp)                :: yllcornerIn"       yllcorner at an input level
  !>        \param[in] "real(dp)                :: cellsizeIn"        cell size at an input level
  !>        \param[in] "real(dp)                :: nodata_valueIn"    nodata value at an input level
  !>        \param[in] "real(dp)                :: aimingResolution"  resolution of an output level

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)             :: nrowsOut"         no. of rows at an output level
  !>        \param[out] "integer(i4)             :: ncolsOut"         no. of cols at an output level
  !>        \param[out] "real(dp)                :: xllcornerOut"      xllcorner at an output level
  !>        \param[out] "real(dp)                :: yllcornerOut"      yllcorner at an output level
  !>        \param[out] "real(dp)                :: cellsizeOut"       cell size at an output level
  !>        \param[out] "real(dp)                :: nodata_valueOut"   nodata value at an output level

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note resolutions of input and output levels should confirm each other.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Zink & Rohini Kumar
  !>        \date Feb 2013
  !         Modified, R. Kumar, Sep 2013   - documentation added according to the template

  subroutine calculate_grid_properties( &
       nrowsIn, ncolsIn,  xllcornerIn,  yllcornerIn,  cellsizeIn, nodata_valueIn,  &
       aimingResolution,                                                           &
       nrowsOut, ncolsOut, xllcornerOut, yllcornerOut, cellsizeOut, nodata_valueOut )

    use mo_message,      only: message       ! for print out
    use mo_string_utils, only: num2str

    implicit none

    integer(i4), intent(in) :: nrowsIn
    integer(i4), intent(in) :: ncolsIn
    real(dp), intent(in)    :: xllcornerIn
    real(dp), intent(in)    :: yllcornerIn
    real(dp), intent(in)    :: cellsizeIn
    real(dp), intent(in)    :: nodata_valueIn
    real(dp), intent(in)    :: aimingResolution

    integer(i4), intent(out) :: nrowsOut
    integer(i4), intent(out) :: ncolsOut
    real(dp), intent(out)    :: xllcornerOut
    real(dp), intent(out)    :: yllcornerOut
    real(dp), intent(out)    :: cellsizeOut
    real(dp), intent(out)    :: nodata_valueOut

    ! local variables
    real(dp)                 :: cellfactor

    cellFactor = aimingResolution / cellsizeIn

    if ( nint(mod(aimingResolution, cellsizeIn)) .ne. 0) then
       call message()
       call message('***ERROR: Two resolutions size do not confirm: ',   &
            trim(adjustl(num2str(nint(AimingResolution)))),              &
            trim(adjustl(num2str(nint(cellsizeIn))))         )
       stop
    end if

    cellsizeOut     = cellsizeIn * cellFactor
    ncolsOut        = ceiling( real(ncolsIn, dp) / cellFactor)
    nrowsOut        = ceiling( real(nrowsIn, dp) / cellFactor)
    xllcornerOut    = xllcornerIn + real(ncolsIn,dp) * cellsizeIn - real(ncolsOut,dp) * cellsizeOut
    yllcornerOut    = yllcornerIn + real(nrowsIn,dp) * cellsizeIn - real(nrowsOut,dp) * cellsizeOut
    nodata_valueOut  = nodata_valueIn

  end subroutine calculate_grid_properties

END MODULE mo_init_states
