!> \file mo_init_states.f90
!> \brief \copybrief mo_init_states
!> \details \copydetails mo_init_states

!> \brief Initialization of all state variables of mHM.
!> \details This module initializes all state variables required to run mHM.
!!
!!       Two options are provided:
!!       - (1) default values
!!       - (2) from nc file
!> \authors Luis Samaniego & Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_init_states

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: variables_alloc                 ! allocation of space for state variables/fluxes/effective parameters
  PUBLIC :: variables_default_init          ! initialization for state variables/fluxes/effective parameters
  PUBLIC :: fluxes_states_default_init      ! initialization for state/fluxes variables

CONTAINS


  !> \brief Allocation of space for mHM related L1 and L11 variables.
  !> \details Allocation of space for mHM related L1 and L11 variables (e.g., states,
  !! fluxes, and parameters) for a given domain. Variables allocated here is
  !! defined in them mo_global_variables.f90 file. After allocating any variable
  !! in this routine, initalize them in the following variables_default_init subroutine.
  !> \changelog
  !! - R. Kumar           Sep 2013
  !!   - documentation added according to the template
  !! - S. Thober          Aug 2015
  !!   - removed routing related variables
  !! - Zink M. Demirel C. Mar 2017
  !!   - Init Jarvis soil water stress variable at SM process(3)
  !! - Robert Schweppe    Dec 2017
  !!   - restructured allocation in variables_alloc, expanded dimensions of effective parameters
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine variables_alloc(ncells1)

    use mo_append, only : append
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_global_variables, only : L1_Throughfall, L1_aETCanopy, L1_aETSealed, L1_aETSoil, L1_baseflow, &
                                    L1_fastRunoff, L1_infilSoil, L1_inter, L1_melt, L1_neutrons, L1_percol, &
                                    L1_pet_calc, L1_temp_calc, L1_prec_calc, &
                                    L1_preEffect, L1_rain, L1_runoffSeal, L1_satSTW, L1_sealSTW, &
                                    L1_slowRunoff, L1_snow, L1_snowPack, L1_soilMoist, L1_total_runoff, L1_unsatSTW
    use mo_mpr_constants, only : C1_InitStateSM, P2_InitStateFluxes, P3_InitStateFluxes, &
                                 P4_InitStateFluxes, P5_InitStateFluxes
    use mo_mpr_global_variables, only : HorizonDepth_mHM, nSoilHorizons_mHM

    implicit none

    integer(i4), intent(in) :: ncells1 !< number of level-1 cells

    integer(i4) :: i

    real(dp), dimension(:), allocatable :: dummy_1D

    real(dp), dimension(:, :), allocatable :: dummy_2D


    ! for appending and intialization
    allocate(dummy_1D(nCells1))
    allocate(dummy_2D(nCells1, nSoilHorizons_mHM))

    dummy_1D = P1_InitStateFluxes
    dummy_2D = P1_InitStateFluxes

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------
    ! calculated / corrected potential evapotranspiration
    call append(L1_pet_calc, dummy_1D)
    ! temperature for current time step
    call append(L1_temp_calc, dummy_1D)
    ! precipitation for current time step
    call append(L1_prec_calc, dummy_1D)
    !  soil actual ET
    call append(L1_aETSoil, dummy_2D)
    ! canopy actual ET
    call append(L1_aETCanopy, dummy_1D)
    ! sealed area actual ET
    call append(L1_aETSealed, dummy_1D)
    ! baseflow
    call append(L1_baseflow, dummy_1D)
    !  soil in-exfiltration
    call append(L1_infilSoil, dummy_2D)
    ! fast runoff
    call append(L1_fastRunoff, dummy_1D)
    ! snow melt
    call append(L1_melt, dummy_1D)
    ! percolation
    call append(L1_percol, dummy_1D)
    ! effective precip. depth (snow melt + rain)
    call append(L1_preEffect, dummy_1D)
    ! rain (liquid water)
    call append(L1_rain, dummy_1D)
    ! runoff from impervious area
    call append(L1_runoffSeal, dummy_1D)
    ! slow runoff
    call append(L1_slowRunoff, dummy_1D)
    ! snow (solid water)
    call append(L1_snow, dummy_1D)
    ! throughfall
    call append(L1_Throughfall, dummy_1D)
    ! throughfall
    call append(L1_total_runoff, dummy_1D)

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------
    ! Interception
    call append(L1_inter, dummy_1D)
    !Retention storage of impervious areas
    call append(L1_sealSTW, dummy_1D)
    ! ground albedo neutrons
    call append(L1_neutrons, dummy_1D)
    !Snowpack
    dummy_1D = P2_InitStateFluxes
    call append(L1_snowPack, dummy_1D)
    ! upper soil storage
    dummy_1D = P3_InitStateFluxes
    call append(L1_unsatSTW, dummy_1D)
    ! groundwater storage
    dummy_1D = P4_InitStateFluxes
    call append(L1_satSTW, dummy_1D)
    ! Soil moisture of each horizon
    do i = 1, nSoilHorizons_mHM - 1
      if (i == 1) then
        dummy_2D(:, i) = HorizonDepth_mHM(i) * C1_InitStateSM
      else
        dummy_2D(:, i) = (HorizonDepth_mHM(i) - HorizonDepth_mHM(i - 1)) * C1_InitStateSM
      end if
    end do
    dummy_2D(:, nSoilHorizons_mHM) = (P5_InitStateFluxes - &
            HorizonDepth_mHM(nSoilHorizons_mHM - 1)) * C1_InitStateSM
    call append(L1_soilMoist, dummy_2D)

    ! free space
    if (allocated(dummy_1D)) deallocate(dummy_1D)
    if (allocated(dummy_2D)) deallocate(dummy_2D)

  end subroutine variables_alloc


  !> \brief Default initalization mHM related L1 variables
  !> \details Default initalization of mHM related L1 variables (e.g., states,
  !! fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !! Variables initalized here is defined in the mo_global_variables.f90 file.
  !! Only Variables that are defined in the variables_alloc subroutine are
  !! intialized here.
  !! If a variable is added or removed here, then it also has to be added or removed
  !! in the subroutine state_variables_set in the module mo_restart and in the
  !! subroutine set_state in the module mo_set_netcdf_restart.
  !> \changelog
  !! - R. Kumar       Sep 2013
  !!   - documentation added according to the template
  !! - Stephan Thober Aug 2015
  !!   - moved routing variables to mRM
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - added separate fluxes_states_default_init
  !> \authors R. Kumar & J. Mai
  !> \date Sep 2013
  subroutine variables_default_init

    use mo_common_constants, only : P1_InitStateFluxes
    use mo_mpr_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, &
                                        L1_aeroResist, L1_alpha, L1_degDay, L1_degDayInc, L1_degDayMax,&
                                        L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, L1_jarvis_thresh_c1, &
                                        L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, L1_kfastFlow, &
                                        L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                        L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, &
                                        L1_unsatThresh, L1_wiltingPoint

    implicit none

    ! init fluxes and states
    call fluxes_states_default_init()

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! sealed fraction of LCover
    L1_fSealed = P1_InitStateFluxes
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
    ! PET correction factor due to LAI
    L1_petLAIcorFactor = P1_InitStateFluxes
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
    ! jarvis critical value for normalized soil water content
    L1_jarvis_thresh_c1 = P1_InitStateFluxes
    ! Threshold temperature for snow/rain
    L1_tempThresh = P1_InitStateFluxes
    ! Threshhold water depth controlling fast interflow
    L1_unsatThresh = P1_InitStateFluxes
    ! Threshhold water depth for surface runoff in sealed surfaces
    L1_sealedThresh = P1_InitStateFluxes
    ! Permanent wilting point
    L1_wiltingPoint = P1_InitStateFluxes

  end subroutine variables_default_init


  !> \brief initialize fluxes and states with default values
  !> \authors Sebastian Müller
  !> \date Mar 2023
  subroutine fluxes_states_default_init

    use mo_common_constants, only : P1_InitStateFluxes
    use mo_mpr_constants, only : C1_InitStateSM, P2_InitStateFluxes, P3_InitStateFluxes, &
                                 P4_InitStateFluxes, P5_InitStateFluxes
    use mo_global_variables, only : L1_Throughfall, L1_aETCanopy, L1_aETSealed, L1_aETSoil, L1_baseflow, &
                                    L1_fastRunoff, L1_infilSoil, &
                                    L1_inter, L1_melt, L1_neutrons, L1_percol, L1_pet_calc, L1_temp_calc, L1_prec_calc, &
                                    L1_preEffect, L1_rain, &
                                    L1_runoffSeal, L1_satSTW, L1_sealSTW, L1_slowRunoff, L1_snow, L1_snowPack, &
                                    L1_soilMoist, L1_total_runoff, L1_unsatSTW
    use mo_mpr_global_variables, only : HorizonDepth_mHM, nSoilHorizons_mHM

    implicit none

    integer(i4) :: i

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
    do i = 1, nSoilHorizons_mHM - 1
      if (i .eq. 1) then
        L1_soilMoist(:, i) = HorizonDepth_mHM(i) * C1_InitStateSM
      else
        L1_soilMoist(:, i) = (HorizonDepth_mHM(i) - HorizonDepth_mHM(i - 1)) * C1_InitStateSM
      end if
    end do
    L1_soilMoist(:, nSoilHorizons_mHM) = (P5_InitStateFluxes - &
            HorizonDepth_mHM(nSoilHorizons_mHM - 1)) * C1_InitStateSM

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
    ! temperature for current time step
    L1_temp_calc = P1_InitStateFluxes
    ! precipitation for current time step
    L1_prec_calc = P1_InitStateFluxes
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
    ! total runoff
    L1_total_runoff = P1_InitStateFluxes

  end subroutine fluxes_states_default_init

END MODULE mo_init_states
