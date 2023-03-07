!> \file mo_startup.f90
!> \brief \copybrief mo_startup
!> \details \copydetails mo_startup

!> \brief Startup procedures for mHM.
!> \details This module initializes all variables required to run mHM. This
!> module needs to be run only one time at the beginning of a simulation if
!> re-starting files do not exist.
!> \authors Luis Samaniego, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_startup

  USE mo_kind, ONLY : i4, dp
  use mo_message, only: message, error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_initialize        ! initialization sequence

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
  !! - Sebastian Müller Mar 2023
  !!   - added separate read_nLAI_and_check_dims to correctly read nLAI from restart
  !> \authors Luis Samaniego, Rohini Kumar
  !> \date Dec 2012
  subroutine mhm_initialize

    use mo_common_mHM_mRM_variables, only : mhmFileRestartIn, read_restart
    use mo_common_restart, only : read_grid_info, read_nLAI_and_check_dims
    use mo_common_variables, only : level0, level1, domainMeta
    use mo_global_variables, only : level2
    use mo_grid, only : set_domain_indices
    use mo_init_states, only : variables_alloc
    use mo_kind, only : i4
    use mo_mpr_startup, only : init_eff_params, mpr_initialize
    use mo_mpr_global_variables, only: nLAI

    implicit none

    integer(i4) :: iDomain, domainID

    ! constants initialization
    call constants_init()
    allocate(level2(domainMeta%nDomains))

    if (read_restart) then
      allocate(level1(domainMeta%nDomains))
      allocate(level0(domainMeta%nDomains))
      ! read nLAI from restart files (-1 indicates first reading)
      nLAI = -1_i4
    else
      call mpr_initialize()
    end if

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      if (read_restart) then
        ! this reads only the domain properties
        if (domainMeta%L0DataFrom(iDomain) == iDomain) then
          ! only read level0 data if it is new
          ! similar to mo_common_read_data::read_dem
          call read_grid_info(domainID, mhmFileRestartIn(iDomain), "0", level0(iDomain))
        endif
        call read_grid_info(domainID, mhmFileRestartIn(iDomain), "1", level1(iDomain))
        ! read nLAI from restart
        call read_nLAI_and_check_dims(iDomain, mhmFileRestartIn(iDomain))
        ! Parameter fields have to be allocated in any case
        call init_eff_params(level1(iDomain)%nCells)
      end if

      ! State variables and fluxes
      ! have to be allocated and initialised in any case
      call variables_alloc(level1(iDomain)%nCells)

      ! L2 inialization
      call L2_variable_init(iDomain, level0(domainMeta%L0DataFrom(iDomain)), level2(iDomain))

    end do

    ! if no restart, this is done already in MPR
    if (read_restart) then
      call set_domain_indices(level0, indices=domainMeta%L0DataFrom)
      call set_domain_indices(level1)
    end if

    call set_domain_indices(level2)

  end subroutine mhm_initialize


  !> \brief Initialize mHM constants
  !> \details transformation of time units & initialize constants
  !> \changelog
  !! - Rohini Kumar                 Jan 2013
  !! - Juliane Mai & Matthias Cuntz Nov 2013
  !!   - check timeStep
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Luis Samaniego
  !> \date Dec 2012
  subroutine constants_init
    use mo_common_mHM_mRM_variables, only : timestep, c2TSTu, read_restart
    use mo_common_variables, only : processMatrix
    use mo_file, only : file_namelist_mhm_param
    use mo_global_variables, only : neutron_integral_AFast
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

    ! if reading restart, we don't need GeoUnitList
    if (.not. read_restart) then
      ! check if enough geoparameter are defined in mhm_parameter.nml
      ! this was formerly done after reading of data, but mHM and MPR are now seperate processes
      if ((processMatrix(9, 2)) .NE.  size(GeoUnitList, 1)) then
        call error_message('***ERROR: Mismatch: Number of geological units in ', trim(adjustl(file_hydrogeoclass)), &
                ' is ', trim(adjustl(num2str(size(GeoUnitList, 1)))), raise=.false.)
        call error_message('          while it is ', trim(num2str(processMatrix(9, 2))), &
                ' in ', trim(file_namelist_mhm_param), '!')
      end if
    end if

    c2TSTu = real(timeStep, dp) / 24.0_dp   ! from per timeStep to per day

  end subroutine constants_init


  !> \brief Initalize Level-2 meteorological forcings data
  !> \details following tasks are performed
  !! 1)  cell id & numbering
  !! 2)  mask creation
  !! 3)  append variable of intrest to global ones
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !> \authors Rohini Kumar
  !> \date Feb 2013
  subroutine L2_variable_init(iDomain, level0_iDomain, level2_iDomain)

    use mo_common_variables, only : Grid
    use mo_global_variables, only : dirPrecipitation
    use mo_grid, only : init_lowres_level
    use mo_mpr_file, only : file_meteo_header, umeteo_header
    use mo_read_spatial_data, only : read_header_ascii
    use mo_string_utils, only : num2str

    implicit none

    !> domain Id
    integer(i4), intent(in) :: iDomain
    !> associated level-0 grid
    type(Grid), intent(in) :: level0_iDomain
    !> associated level-2 grid to populate
    type(Grid), intent(inout) :: level2_iDomain

    integer(i4) :: nrows2, ncols2
    real(dp) :: xllcorner2, yllcorner2
    real(dp) :: cellsize2, nodata_dummy
    character(256) :: fName

    !--------------------------------------------------------
    ! 1) Estimate each variable locally for a given domain
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! read header file
    ! NOTE: assuming the header file for all meteo variables are same as that of precip.
    fName = trim(adjustl(dirPrecipitation(iDomain))) // trim(adjustl(file_meteo_header))
    call read_header_ascii(trim(fName), umeteo_header, &
            nrows2, ncols2, xllcorner2, &
            yllcorner2, cellsize2, nodata_dummy)

    call init_lowres_level(level0_iDomain, cellsize2, level2_iDomain)

    ! check
    if ((ncols2     .ne.  level2_iDomain%ncols)         .or. &
            (nrows2     .ne.  level2_iDomain%nrows)         .or. &
            (abs(xllcorner2 - level2_iDomain%xllcorner) .gt. tiny(1.0_dp))     .or. &
            (abs(yllcorner2 - level2_iDomain%yllcorner) .gt. tiny(1.0_dp))     .or. &
            (abs(cellsize2 - level2_iDomain%cellsize)  .gt. tiny(1.0_dp))) then
      call error_message('   ***ERROR: subroutine L2_variable_init: size mismatch in grid file for level2 in domain ', &
              trim(adjustl(num2str(iDomain))), '!', raise=.false.)
      call error_message('  Expected to have following properties (based on L0):', raise=.false.)
      call error_message('... rows:     ', trim(adjustl(num2str(level2_iDomain%nrows))), ', ', raise=.false.)
      call error_message('... cols:     ', trim(adjustl(num2str(level2_iDomain%ncols))), ', ', raise=.false.)
      call error_message('... cellsize: ', trim(adjustl(num2str(level2_iDomain%cellsize))), ', ', raise=.false.)
      call error_message('... xllcorner:', trim(adjustl(num2str(level2_iDomain%xllcorner))), ', ', raise=.false.)
      call error_message('... yllcorner:', trim(adjustl(num2str(level2_iDomain%yllcorner))), ', ', raise=.false.)
      call error_message('  Provided (in precipitation file):', raise=.false.)
      call error_message('... rows:     ', trim(adjustl(num2str(nrows2))), ', ', raise=.false.)
      call error_message('... cols:     ', trim(adjustl(num2str(ncols2))), ', ', raise=.false.)
      call error_message('... cellsize: ', trim(adjustl(num2str(cellsize2))), ', ', raise=.false.)
      call error_message('... xllcorner:', trim(adjustl(num2str(xllcorner2))), ', ', raise=.false.)
      call error_message('... yllcorner:', trim(adjustl(num2str(yllcorner2))), ', ')
    end if

  end subroutine L2_variable_init

END MODULE mo_startup
