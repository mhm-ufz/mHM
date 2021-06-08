!>       \file mo_startup.f90

!>       \brief Startup procedures for mHM.

!>       \details This module initializes all variables required to run mHM. This
!>       module needs to be run only one time at the beginning of a simulation if
!>       re-starting files do not exist.

!>       \authors Luis Samaniego, Rohini Kumar

!>       \date Dec 2012

! Modifications:

MODULE mo_startup

  ! This module provides the startup routines for mHM.

  ! Written Luis Samaniego, Rohini Kumar, Dec 2012

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_initialize        ! initialization sequence

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mhm_initialize

  !    PURPOSE
  !>       \brief Initialize main mHM variables

  !>       \details Initialize main mHM variables for a given domain.
  !>       Calls the following procedures in this order:
  !>       - Constant initialization.
  !>       - Generate soil database.
  !>       - Checking inconsistencies input fields.
  !>       - Variable initialization at level-0.
  !>       - Variable initialization at level-1.
  !>       - Variable initialization at level-11.
  !>       - Space allocation of remaining variable/parameters.
  !>       Global variables will be used at this stage.

  !    HISTORY
  !>       \authors Luis Samaniego, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Luis Samaniego Mar 2008 - fully distributed multilayer
  ! Rohini Kumar   Oct 2010 - matrix to vector version 
  !                         - openmp parallelization 
  !                         - routing level 11
  ! Luis Samaniego Jul 2012 - removal of IMSL dependencies
  ! Luis Samaniego Dec 2012 - modular version
  ! Rohini Kumar   May 2013 - code cleaned and error checks
  ! Rohini Kumar   Nov 2013 - updated documentation
  ! Stephan Thober Jun 2014 - copied L2 initialization from mo_meteo_forcings
  ! Stephan Thober Jun 2014 - updated flag for read_restart
  ! Stephan Thober Aug 2015 - removed initialisation of routing
  ! Rohini Kumar   Mar 2016 - changes for handling multiple soil database options
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mhm_initialize(parameterValues, parameterNames, opti_domain_indices)

    use mo_grid, only : read_grid_info, set_domain_indices, calculate_grid_properties,  infer_grid_info, Grid
    use mo_common_variables, only : level0, level1, domainMeta, mhmFileRestartIn, read_restart, nuniqueL0Domains
    use mo_global_variables, only : level2
    use mo_init_states, only : variables_alloc
    use mo_global_variables, only : dirPrecipitation
    use mo_string_utils, only : num2str
    use mo_message, only : message
    use mo_restart, only: read_restart_states
    use mo_file, only: file_namelist_mhm, file_namelist_mhm_param, unamelist_mhm, unamelist_mhm_param
    use mo_mhm_mpr_interface, only: call_mpr

    implicit none

    real(dp), dimension(:), intent(in) :: parameterValues
    character(64), dimension(:), intent(in) :: parameterNames
    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices
    integer(i4) :: iDomain, domainID
    type(Grid) :: dummy
    type(Grid), pointer :: level1_iDomain

    ! constants initialization
    allocate(level2(domainMeta%nDomains))
    allocate(level1(domainMeta%nDomains))

    if (read_restart) then
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)

        ! this reads only the domain properties
        ! domainID, inputFile, level_name, new_grid
        call read_grid_info(domainID, mhmFileRestartIn(iDomain), "1", level1(iDomain))
        ! Parameter fields have to be allocated in any case
        call init_eff_params(level1(iDomain)%nCells)

        call read_restart_states(iDomain, mhmFileRestartIn(iDomain), .false.)

      end do
    else
      call call_mpr(parameterValues, parameterNames, level1, .true., opti_domain_indices)
    end if

    call constants_init()

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      ! State variables and fluxes
      ! have to be allocated and initialised in any case
      call variables_alloc(level1(iDomain)%nCells)

      ! L2 inialization
      call infer_grid_info(trim(dirPrecipitation(iDomain)) // 'pre.nc', 'x', 'y', 'pre', level2(iDomain))

      level1_iDomain => level1(iDomain)
      call calculate_grid_properties(level1_iDomain%nrows, level1_iDomain%ncols, &
        level1_iDomain%xllcorner, level1_iDomain%yllcorner, level1_iDomain%cellsize, &
        level2(iDomain)%cellsize, &
        dummy%nrows, dummy%ncols, &
        dummy%xllcorner, dummy%yllcorner, dummy%cellsize)

      ! check
      if ((dummy%ncols     /=  level2(iDomain)%ncols)         .or. &
              (dummy%nrows /= level2(iDomain)%nrows)         .or. &
              (abs(dummy%xllcorner - level2(iDomain)%xllcorner) > tiny(1.0_dp))     .or. &
              (abs(dummy%yllcorner - level2(iDomain)%yllcorner) > tiny(1.0_dp))     .or. &
              (abs(dummy%cellsize - level2(iDomain)%cellsize)  > tiny(1.0_dp))) then
        call message('   ***ERROR: size mismatch in grid file for meteo input in domain ', &
                trim(adjustl(num2str(iDomain))), '!')
        call message('  Provided (in precipitation file):')
        call message('... rows:     ', trim(adjustl(num2str(level2(iDomain)%nrows))), ', ')
        call message('... cols:     ', trim(adjustl(num2str(level2(iDomain)%ncols))), ', ')
        call message('... cellsize: ', trim(adjustl(num2str(level2(iDomain)%cellsize))), ', ')
        call message('... xllcorner:', trim(adjustl(num2str(level2(iDomain)%xllcorner))), ', ')
        call message('... yllcorner:', trim(adjustl(num2str(level2(iDomain)%yllcorner))), ', ')
        call message('  Expected to have following properties (based on level 1):')
        call message('... rows:     ', trim(adjustl(num2str(dummy%nrows))), ', ')
        call message('... cols:     ', trim(adjustl(num2str(dummy%ncols))), ', ')
        call message('... cellsize: ', trim(adjustl(num2str(dummy%cellsize))), ', ')
        call message('... xllcorner:', trim(adjustl(num2str(dummy%xllcorner))), ', ')
        call message('... yllcorner:', trim(adjustl(num2str(dummy%yllcorner))), ', ')
        stop 1
      end if


    end do

    call set_domain_indices(level2)

  end subroutine mhm_initialize

  ! ------------------------------------------------------------------

  !    NAME
  !        constants_init

  !    PURPOSE
  !>       \brief Initialize mHM constants

  !>       \details transformation of time units & initialize constants

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2012

  ! Modifications:
  ! Rohini Kumar                 Jan 2013 - 
  ! Juliane Mai & Matthias Cuntz Nov 2013 - check timeStep
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine constants_init

    use mo_common_variables, only : processMatrix, c2TSTu
    use mo_common_datetime_type, only: timestep
    use mo_file, only : file_namelist_mhm_param
    use mo_global_variables, only : neutron_integral_AFast
    use mo_message, only : message
    use mo_mpr_file, only : file_hydrogeoclass
    use mo_mpr_global_variables, only : GeoUnitList
    use mo_neutrons, only : TabularIntegralAFast
    use mo_string_utils, only : num2str

    implicit none


    !Fill Tabular for neutron flux integral
    if (processMatrix(10, 1) .eq. 2) then
      allocate(neutron_integral_AFast(10000 + 2))
      call TabularIntegralAFast(neutron_integral_AFast, 20.0_dp)
    else
      allocate(neutron_integral_AFast(1))
      neutron_integral_AFast(:) = 0.0_dp
    endif

    c2TSTu = real(timeStep, dp) / 24.0_dp   ! from per timeStep to per day

  end subroutine constants_init

    subroutine init_eff_params(ncells1)

    use mo_append, only : append
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_common_variables, only: nLandCoverPeriods
    use mo_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_alpha, &
                                        L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, &
                                        L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, &
                                        L1_kFastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                        L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, &
                                        L1_wiltingPoint, L1_latitude, nLAIs, nSoilHorizons

    implicit none

    integer(i4), intent(in) :: ncells1

    real(dp), dimension(:, :, :), allocatable :: dummy_3D
    real(dp), dimension(:, :), allocatable :: dummy_2D
    real(dp), dimension(:), allocatable :: dummy_1D

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    ! for appending and intialization
    allocate(dummy_3D(nCells1, nSoilHorizons, nLandCoverPeriods))
    dummy_3D = P1_InitStateFluxes
    ! Fraction of roots in soil horizons
    call append(L1_fRoots, dummy_3D)
    ! Soil moisture below which actual ET is reduced linearly till PWP
    call append(L1_soilMoistFC, dummy_3D)
    ! Saturation soil moisture for each horizon [mm]
    call append(L1_soilMoistSat, dummy_3D)
    ! Exponential parameter to how non-linear is the soil water retention
    call append(L1_soilMoistExp, dummy_3D)
    ! Permanent wilting point
    call append(L1_wiltingPoint, dummy_3D)
    deallocate(dummy_3D)

    allocate(dummy_3D(nCells1, nLAIs, nLandCoverPeriods))
    dummy_3D = P1_InitStateFluxes
    ! PET correction factor due to LAI
    call append(L1_petLAIcorFactor, dummy_3D)
    ! PET aerodynamical resistance
    call append(L1_aeroResist, dummy_3D)
    deallocate(dummy_3D)

    allocate(dummy_2D(nCells1, nLandCoverPeriods))
    dummy_2D = P1_InitStateFluxes
    call append(L1_fSealed, dummy_2D)
    ! increase of the Degree-day factor per mm of increase in precipitation
    call append(L1_degDayInc, dummy_2D)
    ! maximum degree-day factor
    call append(L1_degDayMax, dummy_2D)
    ! degree-day factor with no precipitation
    call append(L1_degDayNoPre, dummy_2D)
    ! fast interflow recession coefficient
    call append(L1_kFastFlow, dummy_2D)
    ! Threshold temperature for snow/rain
    call append(L1_tempThresh, dummy_2D)
    ! Threshold water depth controlling fast interflow
    call append(L1_unsatThresh, dummy_2D)
    ! slow interflow recession coefficient
    call append(L1_kSlowFlow, dummy_2D)
    ! percolation coefficient
    call append(L1_kPerco, dummy_2D)
    ! exponent for the upper reservoir
    call append(L1_alpha, dummy_2D)
    ! baseflow recession coefficient
    call append(L1_kBaseFlow, dummy_2D)
    deallocate(dummy_2D)

    allocate(dummy_2D(nCells1, nLAIs))
    dummy_2D = P1_InitStateFluxes
    ! PET Prietley Taylor coefficient
    call append(L1_PrieTayAlpha, dummy_2D)
    ! PET bulk surface resistance
    call append(L1_surfResist, dummy_2D)
    ! Maximum interception
    call append(L1_maxInter, dummy_2D)
    deallocate(dummy_2D)

    allocate(dummy_1D(nCells1))
    dummy_1D = P1_InitStateFluxes
    ! Karstic percolation loss
    call append(L1_karstLoss, dummy_1D)
    ! PET correction factor due to terrain aspect
    call append(L1_fAsp, dummy_1D)
    ! latitude
    call append(L1_latitude, dummy_1D)
    ! PET Hargreaves Samani coefficient
    call append(L1_HarSamCoeff, dummy_1D)
    ! jarvis critical value for normalized soil water content
    call append(L1_jarvis_thresh_c1, dummy_1D)
    ! Threshhold water depth for surface runoff in sealed surfaces
    call append(L1_sealedThresh, dummy_1D)
    deallocate(dummy_1D)


  end subroutine init_eff_params

END MODULE mo_startup
