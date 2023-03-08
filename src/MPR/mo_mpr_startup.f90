!> \file mo_mpr_startup.f90
!> \brief \copybrief mo_mpr_startup
!> \details \copydetails mo_mpr_startup

!> \brief Startup procedures for mHM.
!> \details This module initializes all variables required to run mHM. This
!! module needs to be run only one time at the beginning of a simulation if re-starting files do not exist.
!> \authors Luis Samaniego, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_mpr_startup

  USE mo_kind, ONLY : i4, dp
  use mo_common_constants, only : nodata_i4, nodata_dp   ! global nodata values (i4, dp)

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mpr_initialize, init_eff_params                      ! initialization sequence

CONTAINS

  !> \brief Initialize main mHM variables
  !> \details Initialize main mHM variables for a given domain.
  !! Calls the following procedures in this order:
  !! - Constant initialization.
  !! - Generate soil database.
  !! - Checking inconsistencies input fields.
  !! - Variable initialization at level-0.
  !! - Variable initialization at level-1.
  !! - Variable initialization at level-11.
  !! - Space allocation of remaining variable/parameters.
  !! Global variables will be used at this stage.
  !> \changelog
  !! - Luis Samaniego Mar 2008
  !!   - fully distributed multilayer
  !! - Rohini Kumar   Oct 2010
  !!   - matrix to vector version
  !!   - openmp parallelization
  !!   - routing level 11
  !! - Luis Samaniego Jul 2012
  !!   - removal of IMSL dependencies
  !! - Luis Samaniego Dec 2012
  !!   - modular version
  !! - Rohini Kumar   May 2013
  !!   - code cleaned and error checks
  !! - Rohini Kumar   Nov 2013
  !!   - updated documentation
  !! - Stephan Thober Jun 2014
  !!   - copied L2 initialization from mo_meteo_forcings
  !! - Stephan Thober Jun 2014
  !!   - updated flag for read_restart
  !! - Stephan Thober Aug 2015
  !!   - removed initialisation of routing
  !! - Rohini Kumar   Mar 2016
  !!   - changes for handling multiple soil database options
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Luis Samaniego, Rohini Kumar
  !> \date Dec 2012
  subroutine mpr_initialize

    use mo_common_variables, only : l0_l1_remap, level0, level1, domainMeta, resolutionHydrology
    use mo_grid, only : init_lowres_level, set_domain_indices
    use mo_kind, only : i4
    use mo_read_latlon, only : read_latlon
    use mo_soil_database, only : generate_soil_database

    implicit none

    integer(i4) :: iDomain


    ! soilDB common for all domains
    call generate_soil_database()

    allocate(level1(domainMeta%nDomains))
    allocate(l0_l1_remap(domainMeta%nDomains))

    ! L0 and L1 initialization
    do iDomain = 1, domainMeta%nDomains
      if (iDomain .eq. 1) then
        call L0_check_input(domainMeta%L0DataFrom(iDomain))
        call L0_variable_init(domainMeta%L0DataFrom(iDomain))
      ! ToDo: adopt to parallel
      ! ToDo: check change
      else if (domainMeta%L0DataFrom(iDomain) == iDomain) then
        ! this needs only be done if there is new input
        call L0_check_input(domainMeta%L0DataFrom(iDomain))
        call L0_variable_init(domainMeta%L0DataFrom(iDomain))
      end if

      call init_lowres_level(level0(domainMeta%L0DataFrom(iDomain)), resolutionHydrology(iDomain), &
      level1(iDomain), l0_l1_remap(iDomain))
      ! read lat lon coordinates for level 1
      call read_latlon(iDomain, "lon", "lat", "level1", level1(iDomain))

      ! Parameter fields have to be allocated in any case
      call init_eff_params(level1(iDomain)%nCells)

    end do

    call set_domain_indices(level1)

  end subroutine mpr_initialize


  !> \brief Check for errors in L0 input data
  !> \details Check for possible errors in input data (morphological and land cover) at level-0
  !> \changelog
  !! - Rohini Kumar   Aug  2013
  !!   - added iFlag_LAI_data_format to handle LAI options, and changed within the code made accordingly
  !! - Rohini  Kumar  Sep 2013
  !!   - read input data for routing processes according & Stephan Thober, to process_matrix flag
  !! - Stephan Thober Aug 2015
  !!   - moved check of L0 routing variables to mRM
  !! - Rohini Kumar   Mar 2016
  !!   - changes for handling multiple soil database options
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine L0_check_input(iDomain)

    use mo_common_constants, only : eps_dp
    use mo_common_variables, only : L0_LCover, L0_elev, level0, nLCoverScene
    use mo_message, only : error_message
    use mo_mpr_global_variables, only : L0_asp, L0_geoUnit, L0_gridded_LAI, &
                                        L0_slope, L0_soilId, iFlag_soilDB, nSoilHorizons_mHM, timeStep_LAI_input
    use mo_string_utils, only : num2str
    use mo_utils, only : eq

    implicit none

    !> domain id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: k, n, nH

    CHARACTER(len=1024) :: message_text = ''

    ! START CHECKING VARIABLES
    do k = level0(iDomain)%iStart, level0(iDomain)%iEnd

      ! elevation [m]
      if (abs(L0_elev(k) - nodata_dp) .lt. eps_dp) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)'))
        call error_message(' Error: elevation has missing value within the valid masked area at cell in domain ', &
                trim(message_text))
      end if

      ! slope [%]
      if (abs(L0_slope(k) - nodata_dp) .lt. eps_dp) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)'))
        call error_message(' Error: slope has missing value within the valid masked area at cell in domain ', &
                trim(message_text))
      end if

      ! aspect [degree]
      if (abs(L0_asp(k) - nodata_dp) .lt. eps_dp) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)'))
        call error_message(' Error: aspect has missing values within the valid masked area at cell in domain ', &
                trim(message_text))
      end if

      ! soil-Id [-]
      nH = 1 !> by default; when iFlag_soilDB = 0
      if (iFlag_soilDB .eq. 1) nH = nSoilHorizons_mHM
      ! another option to handle multiple soil horizons properties
      do n = 1, nH
        if (L0_soilId(k, n) .eq. nodata_i4) then
          message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)')) // ',' // trim(num2str(n, '(I5)'))
          call error_message(' Error: soil id has missing values within the valid masked area at cell in domain and horizon ', &
                  trim(message_text))
        end if
      end do

      ! geological-Id [-]
      if (L0_geoUnit(k) .eq. nodata_i4) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)'))
        call error_message(' Error: geological formation id has missing values within the valid masked area at cell in domain ', &
                trim(message_text))
      end if

      ! landcover scenes
      do  n = 1, nLCoverScene
        if (L0_LCover(k, n) .eq. nodata_i4) then
          message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(iDomain, '(I5)')) // ',' // trim(num2str(n, '(I5)'))
          call error_message(' Error: land cover id has missing values within the valid masked area at cell in domain and scene ', &
                  trim(message_text))
        end if
      end do

      ! land cover scenes related to LAI
      if(timeStep_LAI_input .EQ. 0) then
        if (eq(L0_gridded_LAI(k, 1), nodata_dp)) then
          message_text = trim(num2str(k, '(G5.3)')) // ',' // trim(num2str(iDomain, '(I5)'))
          call error_message(' Error: gridded LAI has missing values within the valid masked area at cell in domain ', &
                  trim(message_text))
        end if
      end if

    end do

  end subroutine L0_check_input


  !> \brief level 0 variable initialization
  !> \details following tasks are performed for L0 data sets
  !! -  cell id & numbering
  !! -  storage of cell cordinates (row and coloum id)
  !! -  empirical dist. of terrain slope
  !! -  flag to determine the presence of a particular soil id
  !! in this configuration of the model run
  !! If a variable is added or removed here, then it also has to
  !! be added or removed in the subroutine config_variables_set in
  !! module mo_restart and in the subroutine set_config in module
  !! mo_set_netcdf_restart
  !> \changelog
  !! - Rohini Kumar & Matthias Cuntz  May 2014
  !!   - cell area calulation based on a regular lat-lon grid or on a regular X-Y coordinate system
  !! - Matthias Cuntz           May 2014
  !!   - changed empirical distribution function so that doubles get the same value
  !! - Matthias Zink & Matthias Cuntz Feb 2016
  !!   - code speed up due to reformulation of CDF calculation
  !! - Rohini Kumar             Mar 2016
  !!   - changes for handling multiple soil database options
  !! - Maren Kaluza             Feb 2018
  !!   - removed slope_val, temp, only sort the index to speed up finding the empirical distribution slope_emp
  !! - Robert Schweppe          Jun 2018
  !!   - refactoring and reformatting
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine L0_variable_init(iDomain)

    use mo_append, only : append
    use mo_common_variables, only : level0
    use mo_grid, only : L0_grid_setup
    use mo_mpr_global_variables, only : L0_slope, L0_slope_emp, L0_soilId, iFlag_soilDB, nSoilHorizons_mHM, &
            nSoilTypes, soilDB
    use mo_orderpack, only : sort_index
    use mo_utils, only : eq

    implicit none

    !> domain id
    integer(i4), intent(in) :: iDomain

    real(dp), dimension(:), allocatable :: slope_emp

    integer(i4), dimension(:), allocatable :: slope_sorted_index

    integer(i4) :: i, j, k, nH, i_sort, i_sortpost

    ! STEPS ::


    !--------------------------------------------------------
    ! 1) Estimate each variable locally for a given domain
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    !------------------------------------------------------
    ! Assign whether a given soil type is present or not
    !------------------------------------------------------
    if (iDomain .eq. 1) then
      allocate(soilDB%is_present(nSoilTypes))
      soilDB%is_present(:) = 0_i4
    end if

    call L0_grid_setup(level0(iDomain))

    !---------------------------------------------------
    ! Estimate empirical distribution of slope
    !---------------------------------------------------
    allocate(slope_emp(level0(iDomain)%nCells), slope_sorted_index(level0(iDomain)%nCells))

    ! get sorted data and sorted indexes to remap later
    slope_sorted_index = sort_index(L0_slope(level0(iDomain)%iStart : level0(iDomain)%iEnd))

    ! empirical distribution of slopes = cumulated number points with slopes that are <= the slope at this point
    !
    !       sorted data                     emp. CDF
    ! 9 |             x x       7/8 |             x x
    !   |                           |
    ! 8 |           x           5/8 |           x
    !   |                           |
    ! 5 |     x x x             4/8 |     x x x
    !   |                           |
    ! 2 |  x                    1/8 |  x
    !   |__________________         |__________________
    !
    ! highest slope value = highest rank or No. of data points / (data points + 1)
    slope_emp(slope_sorted_index(level0(iDomain)%nCells)) = real(level0(iDomain)%nCells, dp) / &
            real(level0(iDomain)%nCells + 1_i4, dp)

    ! backward loop to check if the preceding data point has the same slope value
    do i = level0(iDomain)%nCells - 1, 1, -1
      i_sort=slope_sorted_index(i)
      i_sortpost=slope_sorted_index(i+1)
      if (eq(L0_slope(level0(iDomain)%iStart-1_i4+i_sort), L0_slope(level0(iDomain)%iStart-1_i4+i_sortpost))) then
        ! if yes: assign the same probabitity
        slope_emp(i_sort) = slope_emp(i_sortpost)
      else
        ! if not: assign rank / (data points + 1)
        slope_emp(i_sort) = real(i, dp) / real(level0(iDomain)%nCells + 1_i4, dp)
      end if
    end do

    ! EXAMPLE
    ! in      = [  7, 20, 31, 31, 12, 31, 42 ]
    ! sorted  = [  7, 12, 20, 31, 31, 31, 42 ]
    ! index   = [  1,  5,  2,  3,  4,  6,  7 ]
    ! temp    = [  1,  2,  3,  6,  6,  6,  7 ]
    ! out     = [  1,  3,  6,  6,  2,  6,  7 ] / (len(out) + 1 )

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    call append(L0_slope_emp, slope_emp)

    nH = 1_i4 !> by default; when iFlag_soilDB = 0
    if (iFlag_soilDB .eq. 1) nH = nSoilHorizons_mHM
    do i = 1, nH
      do k = level0(iDomain)%iStart, level0(iDomain)%iEnd
        j = L0_soilId(k, i)
        soilDB%is_present(j) = 1_i4
      end do
    end do

    ! free space
    deallocate(slope_emp, slope_sorted_index)


  end subroutine L0_variable_init


  !> \brief Allocation of space for mHM related L1 and L11 variables.
  !> \details Allocation of space for mHM related L1 and L11 variables (e.g., states,
  !! fluxes, and parameters) for a given domain. Variables allocated here is
  !! defined in them mo_global_variables.f90 file. After allocating any variable
  !! in this routine, initalize them in the following variables_default_init
  !! subroutine.
  !> \changelog
  !! - R. Kumar           Sep 2013
  !!   - documentation added according to the template
  !! - S. Thober          Aug 2015
  !!   - removed routing related variables
  !! - Zink M. Demirel C. Mar 2017
  !!   - Init Jarvis soil water stress variable at SM process(3)
  !! - Robert Schweppe    Dec 2017
  !!   - restructured allocation in variables_alloc, expanded dimensions of effective parameters
  !! - Robert Schweppe    Jun 2018
  !!   - refactoring and reformatting
  !! - Rohini Kumar       Oct 2021
  !!   - Added Neutron count module to mHM integrate into develop branch (5.11.2)
  !! - Sebastian Müller Mar 2023
  !!   - made L1_alpha, L1_kSlowFlow, L1_kBaseFlow and L1_kPerco land cover dependent
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine init_eff_params(ncells1)

    use mo_append, only : append
    use mo_constants, only : YearMonths
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_common_variables, only : nLCoverScene
    use mo_mpr_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_alpha, L1_degDay, &
                                        L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, &
                                        L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, &
                                        L1_kfastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                        L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, &
                                        L1_wiltingPoint, nLAI, nSoilHorizons_mHM, &
                                        L1_No_Count, L1_bulkDens, L1_latticeWater, L1_COSMICL3

    implicit none

    !> number of L1 cells
    integer(i4), intent(in) :: ncells1

    real(dp), dimension(:, :, :), allocatable :: dummy_3D

    integer(i4) :: max_extent


    ! get maximum extent of one dimension 2 or 3
    max_extent = max(nSoilHorizons_mHM, int(YearMonths, i4), nLCoverScene, nLAI)

    ! for appending and intialization
    allocate(dummy_3D(nCells1, max_extent, nLCoverScene))

    dummy_3D = P1_InitStateFluxes

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    call append(L1_fSealed, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! exponent for the upper reservoir
    call append(L1_alpha, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! increase of the Degree-day factor per mm of increase in precipitation
    call append(L1_degDayInc, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! maximum degree-day factor
    call append(L1_degDayMax, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! degree-day factor with no precipitation
    call append(L1_degDayNoPre, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! degree-day factor
    call append(L1_degDay, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! Karstic percolation loss
    call append(L1_karstLoss, dummy_3D(:, 1 : 1, 1 : 1))
    ! PET correction factor due to terrain aspect
    call append(L1_fAsp, dummy_3D(:, 1 : 1, 1 : 1))
    ! PET correction factor due to LAI
    call append(L1_petLAIcorFactor, dummy_3D(:, 1 : nLAI, 1 : nLCoverScene))
    ! PET Hargreaves Samani coefficient
    call append(L1_HarSamCoeff, dummy_3D(:, 1 : 1, 1 : 1))
    ! PET Prietley Taylor coefficient
    call append(L1_PrieTayAlpha, dummy_3D(:, 1 : nLAI, 1 : 1))
    ! PET aerodynamical resistance
    call append(L1_aeroResist, dummy_3D(:, 1 : nLAI, 1 : nLCoverScene))
    ! PET bulk surface resistance
    call append(L1_surfResist, dummy_3D(:, 1 : nLAI, 1 : 1))
    ! Fraction of roots in soil horizons
    call append(L1_fRoots, dummy_3D(:, 1 : nSoilHorizons_mHM, 1 : nLCoverScene))
    ! Maximum interception
    call append(L1_maxInter, dummy_3D(:, 1 : nLAI, 1 : 1))
    ! fast interflow recession coefficient
    call append(L1_kfastFlow, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! slow interflow recession coefficient
    call append(L1_kSlowFlow, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! baseflow recession coefficient
    call append(L1_kBaseFlow, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! percolation coefficient
    call append(L1_kPerco, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! Soil moisture below which actual ET is reduced linearly till PWP
    call append(L1_soilMoistFC, dummy_3D(:, 1 : nSoilHorizons_mHM, 1 : nLCoverScene))
    ! Saturation soil moisture for each horizon [mm]
    call append(L1_soilMoistSat, dummy_3D(:, 1 : nSoilHorizons_mHM, 1 : nLCoverScene))
    ! jarvis critical value for normalized soil water content
    call append(L1_jarvis_thresh_c1, dummy_3D(:, 1 : 1, 1 : 1))
    ! Exponential parameter to how non-linear is the soil water retention
    call append(L1_soilMoistExp, dummy_3D(:, 1 : nSoilHorizons_mHM, 1 : nLCoverScene))
    ! Threshold temperature for snow/rain
    call append(L1_tempThresh, dummy_3D(:, 1 : 1, 1 : nLCoverScene))
    ! Threshold water depth controlling fast interflow
    call append(L1_unsatThresh, dummy_3D(:, 1 : 1, 1 : 1))
    ! Threshold water depth for surface runoff in sealed surfaces
    call append(L1_sealedThresh, dummy_3D(:, 1 : 1, 1 : 1))
    ! Permanent wilting point
    call append(L1_wiltingPoint, dummy_3D(:, 1 : nSoilHorizons_mHM, 1 : nLCoverScene))
    ! neutron count related parameters
    call append(L1_No_Count,     dummy_3D(:, 1:1,                 1:1))            ! N0 count
    call append(L1_bulkDens,     dummy_3D(:, 1:nSoilHorizons_mHM, 1:nLCoverScene)) ! bulk density
    call append(L1_latticeWater, dummy_3D(:, 1:nSoilHorizons_mHM, 1:nLCoverScene)) ! lattice water
    call append(L1_COSMICL3,     dummy_3D(:, 1:nSoilHorizons_mHM, 1:nLCoverScene)) ! cosmic L3 parameter

    ! free space
    if (allocated(dummy_3D)) deallocate(dummy_3D)

  end subroutine init_eff_params

END MODULE mo_mpr_startup
