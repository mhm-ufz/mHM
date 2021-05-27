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

  subroutine mhm_initialize

    use mo_common_mHM_mRM_variables, only : mhmFileRestartIn, read_restart
    use mo_grid, only : read_grid_info, Grid, set_domain_indices
    use mo_common_variables, only : level0, level1, domainMeta
    use mo_global_variables, only : level2
    use mo_init_states, only : variables_alloc
    use mo_mpr_startup, only : init_eff_params, mpr_initialize
    use mo_global_variables, only : dirPrecipitation

    implicit none

    integer(i4) :: iDomain, domainID
    type(Grid) :: dummy
    type(Grid), pointer :: level0_iDomain


    ! constants initialization
    call constants_init()
    allocate(level2(domainMeta%nDomains))

    if (read_restart) then
      allocate(level1(domainMeta%nDomains))
    else
      call mpr_initialize()
    end if

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      if (read_restart) then
        ! this reads only the domain properties
        call read_grid_info(domainID, mhmFileRestartIn(iDomain), "1", level1(iDomain))
        ! Parameter fields have to be allocated in any case
        call init_eff_params(level1(iDomain)%nCells)
      end if

      ! State variables and fluxes
      ! have to be allocated and initialised in any case
      call variables_alloc(level1(iDomain)%nCells)

      ! L2 inialization
      call infer_grid_info(trim(dirPrecipitation(iDomain)) // 'pre/pre.nc', 'xc', 'yc', 'pre', level2(iDomain))

      level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
      call calculate_grid_properties(level0_iDomain%nrows, level0_iDomain%ncols, &
        level0_iDomain%xllcorner, level0_iDomain%yllcorner, level0_iDomain%cellsize, &
        level2(iDomain)%cellsize, &
        dummy%nrows, dummy%ncols, &
        dummy%xllcorner, dummy%yllcorner, dummy%cellsize)

      ! check
      if ((dummy%ncols     /=  level2(iDomain)%ncols)         .or. &
              (dummy%nrows /= level2(iDomain)%nrows)         .or. &
              (abs(dummy%xllcorner - level2(iDomain)%xllcorner) > tiny(1.0_dp))     .or. &
              (abs(dummy%yllcorner - level2(iDomain)%yllcorner) > tiny(1.0_dp))     .or. &
              (abs(dummy%cellsize - level2(iDomain)%cellsize)  > tiny(1.0_dp))) then
        call message('   ***ERROR: subroutine L2_variable_init: size mismatch in grid file for level2 in domain ', &
                trim(adjustl(num2str(iDomain))), '!')
        call message('  Provided (in precipitation file):')
        call message('... rows:     ', trim(adjustl(num2str(level2(iDomain)%nrows))), ', ')
        call message('... cols:     ', trim(adjustl(num2str(level2(iDomain)%ncols))), ', ')
        call message('... cellsize: ', trim(adjustl(num2str(level2(iDomain)%cellsize))), ', ')
        call message('... xllcorner:', trim(adjustl(num2str(level2(iDomain)%xllcorner))), ', ')
        call message('... yllcorner:', trim(adjustl(num2str(level2(iDomain)%yllcorner))), ', ')
        call message('  Expected to have following properties (based on L0):')
        call message('... rows:     ', trim(adjustl(num2str(dummy%nrows))), ', ')
        call message('... cols:     ', trim(adjustl(num2str(dummy%ncols))), ', ')
        call message('... cellsize: ', trim(adjustl(num2str(dummy%cellsize))), ', ')
        call message('... xllcorner:', trim(adjustl(num2str(dummy%xllcorner))), ', ')
        call message('... yllcorner:', trim(adjustl(num2str(dummy%yllcorner))), ', ')
        stop 1
      end if

    end do

    ! if no restart, this is done already in MPR
    if (read_restart) then
      call set_domain_indices(level1)
    end if

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

    use mo_common_mHM_mRM_variables, only : timestep, c2TSTu
    use mo_common_variables, only : processMatrix
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

    ! check if enough geoparameter are defined in mhm_parameter.nml
    ! this was formerly done after reading of data, but mHM and MPR are now seperate processes
    if ((processMatrix(9, 2)) .NE.  size(GeoUnitList, 1)) then
      call message()
      call message('***ERROR: Mismatch: Number of geological units in ', trim(adjustl(file_hydrogeoclass)), &
              ' is ', trim(adjustl(num2str(size(GeoUnitList, 1)))))
      call message('          while it is ', trim(num2str(processMatrix(9, 2))), &
              ' in ', trim(file_namelist_mhm_param), '!')
      stop 1
    end if

    c2TSTu = real(timeStep, dp) / 24.0_dp   ! from per timeStep to per day

  end subroutine constants_init

END MODULE mo_startup
