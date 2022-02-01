!>       \file mo_init_states.f90

!>       \brief Initialization of all state variables of mHM.

!>       \details This module initializes all state variables required to run mHM.
!>       Two options are provided:
!>       - (1) default values
!>       - (2) from nc file

!>       \authors Luis Samaniego & Rohini Kumar

!>       \date Dec 2012

! Modifications:

MODULE mo_init_states

  ! This module provides the startup routines for mHM.

  ! Written Luis Samaniego & Rohini Kumar, Dec 2012

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: variables_alloc                 ! allocation of space for state variables/fluxes/effective parameters
  PUBLIC :: variables_default_init          ! initialization for state variables/fluxes/effective parameters

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        variables_alloc

  !    PURPOSE
  !>       \brief Allocation of space for mHM related L1 and L11 variables.

  !>       \details Allocation of space for mHM related L1 and L11 variables (e.g., states,
  !>       fluxes, and parameters) for a given domain. Variables allocated here is
  !>       defined in them mo_global_variables.f90 file. After allocating any variable
  !>       in this routine, initalize them in the following variables_default_init
  !>       subroutine:

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: ncells1"

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Jan 2013

  ! Modifications:
  ! R. Kumar           Sep 2013 - documentation added according to the template
  ! S. Thober          Aug 2015 - removed routing related variables
  ! Zink M. Demirel C. Mar 2017 - Init Jarvis soil water stress variable at SM process(3)
  ! Robert Schweppe    Dec 2017 - restructured allocation in variables_alloc, expanded dimensions of effective parameters
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine variables_alloc(ncells1)

    use mo_append, only : append
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_global_variables, only : L1, &
                                    soilHorizonBoundaries, nSoilHorizons
    use mo_mhm_constants, only : C1_InitStateSM, P2_InitStateFluxes, P3_InitStateFluxes, &
                                 P4_InitStateFluxes, P5_InitStateFluxes

    implicit none

    integer(i4), intent(in) :: ncells1

    integer(i4) :: i

    real(dp), dimension(:), allocatable :: dummy_1D

    real(dp), dimension(:, :), allocatable :: dummy_2D


    ! for appending and intialization
    allocate(dummy_1D(nCells1))
    allocate(dummy_2D(nCells1, nSoilHorizons))

    dummy_1D = P1_InitStateFluxes
    dummy_2D = P1_InitStateFluxes

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------
    ! calculated / corrected potential evapotranspiration
    call append(L1%pet_calc, dummy_1D)
    !  soil actual ET
    call append(L1%aETSoil, dummy_2D)
    ! canopy actual ET
    call append(L1%aETCanopy, dummy_1D)
    ! sealed area actual ET
    call append(L1%aETSealed, dummy_1D)
    ! baseflow
    call append(L1%baseflow, dummy_1D)
    !  soil in-exfiltration
    call append(L1%infilSoil, dummy_2D)
    ! fast runoff
    call append(L1%fastRunoff, dummy_1D)
    ! snow melt
    call append(L1%melt, dummy_1D)
    ! percolation
    call append(L1%percol, dummy_1D)
    ! effective precip. depth (snow melt + rain)
    call append(L1%preEffect, dummy_1D)
    ! rain (liquid water)
    call append(L1%rain, dummy_1D)
    ! runoff from impervious area
    call append(L1%runoffSeal, dummy_1D)
    ! slow runoff
    call append(L1%slowRunoff, dummy_1D)
    ! snow (solid water)
    call append(L1%snow, dummy_1D)
    ! throughfall
    call append(L1%Throughfall, dummy_1D)
    ! throughfall
    call append(L1%total_runoff, dummy_1D)

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    call append(L1%inter, dummy_1D)
    ! Degree-day factor
    call append(L1%degDay, dummy_1D)
    !Retention storage of impervious areas
    call append(L1%sealSTW, dummy_1D)
    ! ground albedo neutrons
    call append(L1%neutrons, dummy_1D)
    !Snowpack
    dummy_1D = P2_InitStateFluxes
    call append(L1%snowPack, dummy_1D)
    ! upper soil storage
    dummy_1D = P3_InitStateFluxes
    call append(L1%unsatSTW, dummy_1D)
    ! groundwater storage
    dummy_1D = P4_InitStateFluxes
    call append(L1%satSTW, dummy_1D)
    ! Soil moisture of each horizon
    do i = 1, nSoilHorizons - 1
      dummy_2D(:, i) = (soilHorizonBoundaries(i+1) - soilHorizonBoundaries(i)) * C1_InitStateSM * 1000_dp
    end do
    dummy_2D(:, nSoilHorizons) = (P5_InitStateFluxes - &
            soilHorizonBoundaries(nSoilHorizons+1)) * C1_InitStateSM * 1000_dp
    call append(L1%soilMoist, dummy_2D)

    ! free space
    if (allocated(dummy_1D)) deallocate(dummy_1D)
    if (allocated(dummy_2D)) deallocate(dummy_2D)

  end subroutine variables_alloc


  ! ------------------------------------------------------------------

  !    NAME
  !        variables_default_init

  !    PURPOSE
  !>       \brief Default initalization mHM related L1 variables

  !>       \details Default initalization of mHM related L1 variables (e.g., states,
  !>       fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !>       Variables initalized here is defined in the mo_global_variables.f90 file.
  !>       Only Variables that are defined in the variables_alloc subroutine are
  !>       intialized here.
  !>       If a variable is added or removed here, then it also has to be added or removed
  !>       in the subroutine state_variables_set in the module mo_restart and in the
  !>       subroutine set_state in the module mo_set_netcdf_restart.

  !    HISTORY
  !>       \authors R. Kumar & J. Mai

  !>       \date Sep 2013

  ! Modifications:
  ! R. Kumar       Sep 2013 - documentation added according to the template
  ! Stephan Thober Aug 2015 - moved routing variables to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine variables_default_init

    use mo_common_constants, only : P1_InitStateFluxes
    use mo_global_variables, only : L1, &
                                    nSoilHorizons, soilHorizonBoundaries, &
                                    are_parameter_initialized
    use mo_mhm_constants, only : C1_InitStateSM, P2_InitStateFluxes, P3_InitStateFluxes, &
                                 P4_InitStateFluxes, P5_InitStateFluxes

    implicit none

    integer(i4) :: i

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    L1%inter = P1_InitStateFluxes

    !Snowpack
    L1%snowPack = P2_InitStateFluxes

    !Retention storage of impervious areas
    L1%sealSTW = P1_InitStateFluxes

    ! Soil moisture of each horizon
    do i = 1, nSoilHorizons - 1
      L1%soilMoist(:, i) = (soilHorizonBoundaries(i+1) - soilHorizonBoundaries(i)) * C1_InitStateSM  * 1000_dp
    end do
    L1%soilMoist(:, nSoilHorizons) = (P5_InitStateFluxes - &
            soilHorizonBoundaries(nSoilHorizons) * 1000_dp) * C1_InitStateSM

    ! upper soil storage
    L1%unsatSTW = P3_InitStateFluxes

    ! groundwater storage
    L1%satSTW = P4_InitStateFluxes

    ! ground albedo neutrons, initially zero
    L1%neutrons = P1_InitStateFluxes

    ! degree-day factor
    L1%degDay = P1_InitStateFluxes

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------

    ! corrected / calculated potential ET
    L1%pet_calc = P1_InitStateFluxes

    !  soil actual ET
    L1%aETSoil = P1_InitStateFluxes

    ! canopy actual ET
    L1%aETCanopy = P1_InitStateFluxes

    ! sealed area actual ET
    L1%aETSealed = P1_InitStateFluxes

    ! baseflow
    L1%baseflow = P1_InitStateFluxes

    !  soil in-exfiltration
    L1%infilSoil = P1_InitStateFluxes

    ! fast runoff
    L1%fastRunoff = P1_InitStateFluxes

    ! snow melt
    L1%melt = P1_InitStateFluxes

    ! percolation
    L1%percol = P1_InitStateFluxes

    ! effective precip. depth (snow melt + rain)
    L1%preEffect = P1_InitStateFluxes

    ! rain (liquid water)
    L1%rain = P1_InitStateFluxes

    ! runoff from impervious area
    L1%runoffSeal = P1_InitStateFluxes

    ! slow runoff
    L1%slowRunoff = P1_InitStateFluxes

    ! snow (solid water)
    L1%snow = P1_InitStateFluxes

    ! throughfall
    L1%Throughfall = P1_InitStateFluxes

    ! throughfall
    L1%total_runoff = P1_InitStateFluxes

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    if (.not. are_parameter_initialized) then
      ! sealed fraction of LCover
      L1%fSealed = P1_InitStateFluxes

      ! exponent for the upper reservoir
      L1%alpha = P1_InitStateFluxes

      ! increase of the Degree-day factor per mm of increase in precipitation
      L1%degDayInc = P1_InitStateFluxes

      ! maximum degree-day factor
      L1%degDayMax = P1_InitStateFluxes

      ! degree-day factor with no precipitation
      L1%degDayNoPre = P1_InitStateFluxes

      ! Karstic percolation loss
      L1%karstLoss = P1_InitStateFluxes

      ! PET correction factor due to LAI
      L1%petLAIcorFactor = P1_InitStateFluxes

      ! PET correction factor due to terrain aspect
      L1%fAsp = P1_InitStateFluxes

      ! PET Hargreaves Samani Coefficient
      L1%HarSamCoeff = P1_InitStateFluxes

      ! PET Priestley Taylor coefficient
      L1%PrieTayAlpha = P1_InitStateFluxes

      ! PET aerodynamical resistance
      L1%aeroResist = P1_InitStateFluxes

      ! PET bulk surface resistance
      L1%surfResist = P1_InitStateFluxes

      ! Fraction of roots in soil horizons
      L1%fRoots = P1_InitStateFluxes

      ! Maximum interception
      L1%maxInter = P1_InitStateFluxes

      ! fast interflow recession coefficient
      L1%kFastFlow = P1_InitStateFluxes

      ! slow interflow recession coefficient
      L1%kSlowFlow = P1_InitStateFluxes

      ! baseflow recession coefficient
      L1%kBaseFlow = P1_InitStateFluxes

      ! percolation coefficient
      L1%kPerco = P1_InitStateFluxes

      ! Soil moisture below which actual ET is reduced linearly till PWP
      L1%soilMoistFC = P1_InitStateFluxes

      ! Saturation soil moisture for each horizon [mm]
      L1%soilMoistSat = P1_InitStateFluxes

      ! Exponential parameter to how non-linear is the soil water retention
      L1%soilMoistExp = P1_InitStateFluxes

      ! jarvis critical value for normalized soil water content
      L1%jarvis_thresh_c1 = P1_InitStateFluxes

      ! Threshold temperature for snow/rain
      L1%tempThresh = P1_InitStateFluxes

      ! Threshhold water depth controlling fast interflow
      L1%unsatThresh = P1_InitStateFluxes

      ! Threshhold water depth for surface runoff in sealed surfaces
      L1%sealedThresh = P1_InitStateFluxes

      ! Permanent wilting point
      L1%wiltingPoint = P1_InitStateFluxes
    end if

  end subroutine variables_default_init

END MODULE mo_init_states
