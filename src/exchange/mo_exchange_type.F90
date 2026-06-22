!> \dir exchange
!> \copydoc f_exchange

!> \defgroup   f_exchange exchange - Fortran modules
!> \brief      Modules to deal with data exchange between mHM components.
!> \details    This module provides different types to enable the data exchange between components.

!> \file    mo_exchange_type.f90
!> \copydoc mo_exchange_type

!> \brief   Module to provide the exchange type.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Mar 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_exchange_type
  use mo_logging
  use mo_grid, only: grid_t
  use mo_geology_classdefinition, only: geology_classdefinition_t
  use mo_datetime, only: datetime, timedelta
  use mo_kind, only: dp, i4, i8
  use mo_string_utils, only: n2s=>num2str
  use mo_main_config, only: parameters_t
  use mo_utils, only: optval
  use nml_config_time, only: nml_config_time_t
  use nml_config_project, only: nml_config_project_t
  use nml_config_domain, only: nml_config_domain_t
  use nml_config_resolution, only: nml_config_resolution_t
  use nml_helper, only: NML_OK

  implicit none
  private

  public :: standard_path, get_n_domains

  !> \class   exchange_config_t
  !> \brief   Exchange-owned run and domain metadata namelists.
  type, public :: exchange_config_t
    type(nml_config_project_t) :: project     !< top-level/run project configuration
    type(nml_config_domain_t) :: domain       !< top-level domain directory configuration
    type(nml_config_time_t) :: time           !< domain-local time configuration
    type(nml_config_resolution_t) :: resolution !< domain-local target grid resolutions
  end type exchange_config_t

  !> \name Level Selectors
  !> \brief Constants to specify the grid for levels in mHM: L0, L1, L2 and L3.
  !!@{
  integer(i4), public, parameter :: nogrid = -1_i4 !< no grid (yet) defined
  integer(i4), public, parameter :: l0 = 0_i4      !< level0 - morphology
  integer(i4), public, parameter :: l1 = 1_i4      !< level1 - hydrology
  integer(i4), public, parameter :: l2 = 2_i4      !< level2 - meteorology
  integer(i4), public, parameter :: l3 = 3_i4      !< level3 - routing
  !!@}

  !> \class   variable_abc
  !> \brief   Abstract base class for a variable in the exchange type.
  type, abstract :: variable_abc
    character(:), allocatable :: name          !< variable name
    character(:), allocatable :: units         !< variable unit
    character(:), allocatable :: long_name     !< long name of the variable
    character(:), allocatable :: standard_name !< standard name of the variable
    integer(i4) :: grid = nogrid               !< ID of the grid the data is defined on
    integer(i4) :: stepping = 0_i4             !< time-step size of this variable in hours (0 - static)
    logical :: static = .false.                !< flag to indicated static data (.false. by default)
    logical :: provided = .false.              !< flag to indicate that data is provided by a component (.false. by default)
    logical :: required = .false.              !< flag to indicate that data is required by a component (.false. by default)
  contains
    procedure(variable_has_data_i), public, deferred :: has_data
    procedure(variable_data_shape_i), public, deferred :: data_shape
    procedure(variable_clear_data_i), public, deferred :: clear_data
    procedure, public :: available => variable_available
    procedure, public :: require => variable_require
    procedure, public :: expect_handoff => variable_expect_handoff
    procedure, public :: clear => variable_clear
  end type variable_abc

  !> \class   var_dp
  !> \brief   Class for a double precision variable in the exchange type.
  type, public, extends(variable_abc) :: var_dp
    real(dp), dimension(:), pointer :: data => null() !< 1D real pointer (n-cells)
  contains
    procedure, public :: has_data => var_dp_has_data
    procedure, public :: data_shape => var_dp_data_shape
    procedure, public :: clear_data => var_dp_clear_data
    procedure, public :: publish_local => var_dp_publish_local
    procedure, public :: publish_alias => var_dp_publish_alias
  end type var_dp

  !> \class   var_i4
  !> \brief   Class for a 32bit integer variable in the exchange type.
  type, public, extends(variable_abc) :: var_i4
    integer(i4), dimension(:), pointer :: data => null() !< 1D integer pointer (n-cells)
  contains
    procedure, public :: has_data => var_i4_has_data
    procedure, public :: data_shape => var_i4_data_shape
    procedure, public :: clear_data => var_i4_clear_data
    procedure, public :: publish_local => var_i4_publish_local
    procedure, public :: publish_alias => var_i4_publish_alias
  end type var_i4

  !> \class   var_lg
  !> \brief   Class for a logical variable in the exchange type.
  type, public, extends(variable_abc) :: var_lg
    logical, dimension(:), pointer :: data => null() !< 1D logical pointer (n-cells)
  contains
    procedure, public :: has_data => var_lg_has_data
    procedure, public :: data_shape => var_lg_data_shape
    procedure, public :: clear_data => var_lg_clear_data
    procedure, public :: publish_local => var_lg_publish_local
    procedure, public :: publish_alias => var_lg_publish_alias
  end type var_lg

  !> \class   var2d_dp
  !> \brief   Class for a double precision variable for each horizon in the exchange type.
  type, public, extends(variable_abc) :: var2d_dp
    real(dp), dimension(:,:), pointer :: data => null() !< 2D real pointer (n-cells, horizons)
  contains
    procedure, public :: has_data => var2d_dp_has_data
    procedure, public :: data_shape => var2d_dp_data_shape
    procedure, public :: clear_data => var2d_dp_clear_data
    procedure, public :: publish_local => var2d_dp_publish_local
    procedure, public :: publish_alias => var2d_dp_publish_alias
  end type var2d_dp

  !> \class   var2d_i4
  !> \brief   Class for a 32bit integer variable for each horizon in the exchange type.
  type, public, extends(variable_abc) :: var2d_i4
    integer(i4), dimension(:,:), pointer :: data => null() !< 2D integer pointer (n-cells, horizons)
  contains
    procedure, public :: has_data => var2d_i4_has_data
    procedure, public :: data_shape => var2d_i4_data_shape
    procedure, public :: clear_data => var2d_i4_clear_data
    procedure, public :: publish_local => var2d_i4_publish_local
    procedure, public :: publish_alias => var2d_i4_publish_alias
  end type var2d_i4

  !> \class   var2d_lg
  !> \brief   Class for a logical variable for each horizon in the exchange type.
  type, public, extends(variable_abc) :: var2d_lg
    logical, dimension(:,:), pointer :: data => null() !< 2D logical pointer (n-cells, horizons)
  contains
    procedure, public :: has_data => var2d_lg_has_data
    procedure, public :: data_shape => var2d_lg_data_shape
    procedure, public :: clear_data => var2d_lg_clear_data
    procedure, public :: publish_local => var2d_lg_publish_local
    procedure, public :: publish_alias => var2d_lg_publish_alias
  end type var2d_lg

  abstract interface
    !> \brief Return whether an exchange variable currently has an associated data pointer.
    logical function variable_has_data_i(self)
      import :: variable_abc
      class(variable_abc), intent(in) :: self
    end function variable_has_data_i

    !> \brief Return the shape of an associated data pointer, or a zero-length shape if unassociated.
    function variable_data_shape_i(self) result(shape)
      import :: variable_abc, i8
      class(variable_abc), intent(in) :: self
      integer(i8), allocatable :: shape(:)
    end function variable_data_shape_i

    !> \brief Clear the data pointer of an exchange variable.
    subroutine variable_clear_data_i(self)
      import :: variable_abc
      class(variable_abc), intent(inout) :: self
    end subroutine variable_clear_data_i
  end interface

  !> \class   exchange_t
  !> \brief   Class for dynamically exchanging variables in mHM.
  type, public :: exchange_t
    integer(i4) :: step_count               !< current time step
    type(datetime) :: time                  !< time-stamp for the current time step
    type(datetime) :: start_time            !< start time of simulation
    type(datetime) :: eval_start_time       !< start time of evaluation
    type(datetime) :: end_time              !< end time of simulation
    type(timedelta) :: step                 !< time step of the simulation
    type(exchange_config_t) :: config        !< exchange-owned namelist configuration
    type(parameters_t) :: parameters        !< parameters container
    integer(i4) :: domain_id                !< stable domain ID in the run
    integer(i4) :: n_domains = 1_i4         !< number of domains in the run
    integer(i4) :: nml_n_domains = 1_i4     !< domain dimension used when reading namelists
    integer(i4) :: nml_domain_id = 1_i4     !< domain index used when reading namelist arrays
    integer(i4) :: n_geo_units = 25_i4      !< number of geological units in the global parameter set
    integer(i4) :: max_layers = 10_i4       !< maximum number of soil-layer entries in the run
    integer(i4) :: n_layers = 1_i4          !< number of soil layers for this domain
    character(:), allocatable :: cwd        !< current working directory to set relative paths
    character(:), allocatable :: main_file  !< effective main namelist file for this domain

    ! grids
    type(grid_t), pointer :: level0 => null() !< level0 grid of the morphology
    type(grid_t), pointer :: level1 => null() !< level1 grid of the hydrology
    type(grid_t), pointer :: level2 => null() !< level2 grid of the meteorology
    type(grid_t), pointer :: level3 => null() !< level3 grid of the river network
    real(dp), dimension(:), pointer :: soil_horizon_bounds => null() !< soil-horizon boundary depths [mm] for mHM metadata

    ! grid resolutions (for deriving grids after configuration)
    real(dp) :: level0_resolution = -1.0_dp !< level0 resolution of the morphology
    real(dp) :: level1_resolution = -1.0_dp !< level1 resolution of the hydrology
    real(dp) :: level2_resolution = -1.0_dp !< level2 resolution of the meteorology
    real(dp) :: level3_resolution = -1.0_dp !< level3 resolution of the river network

    ! variables
    ! raw meteorology (level2)
    type(var_dp) :: raw_pre             !< raw precipitation [mm] on level l2
    type(var_dp) :: raw_temp            !< raw air temperature [degC] on level l2
    type(var_dp) :: raw_ssrd            !< raw solar short wave radiation downward [W m-2] on level l2
    type(var_dp) :: raw_strd            !< raw surface thermal radiation downward [W m-2] on level l2
    type(var_dp) :: raw_tann            !< raw annual mean air temperature [degC] on level l2
    type(var_dp) :: raw_tmin            !< raw minimum daily temperature [degC] on level l2
    type(var_dp) :: raw_tmax            !< raw maximum daily temperature [degC] on level l2
    type(var_dp) :: raw_netrad          !< raw net radiation [W m-2] on level l2
    type(var_dp) :: raw_eabs            !< raw vapor pressure [Pa] on level l2
    type(var_dp) :: raw_wind            !< raw wind speed [m s-1] on level l2
    type(var_dp) :: raw_pet             !< input potential evapotranspiration [mm] on level l2

    ! processed meteorology (level1)
    type(var_dp) :: pre                 !< precipitation [mm] on level l1
    type(var_dp) :: temp                !< air temperature [degC] on level l1
    type(var_dp) :: pet                 !< potential evapotranspiration [mm] on level l1
    type(var_dp) :: ssrd                !< solar short wave radiation downward [W m-2] on level l1
    type(var_dp) :: strd                !< surface thermal radiation downward [W m-2] on level l1
    type(var_dp) :: tann                !< annual mean air temperature [degC] on level l1

    ! morphology (level0)
    type(var_dp) :: dem                 !< elevation [m] on level l0 (static)
    type(var_dp) :: slope               !< slope [%] on level l0 (static)
    type(var_dp) :: aspect              !< aspect [degree] on level l0 (static)
    type(var_i4) :: fdir                !< flow direction [1] on level l0 (static)
    type(var_i4) :: facc                !< flow accumulation [1] on level l0 (static)
    type(var2d_i4) :: soil_id           !< soil class ID on level l0 (static)
    type(var_i4) :: geo_unit            !< geological unit ID on level l0 (static)
    type(var_i4) :: lai_class           !< LAI class ID on level l0 (static)
    type(var_dp) :: slope_emp           !< empirical slope distribution on level l0 (static)
    type(geology_classdefinition_t), pointer :: geo_class_def => null() !< geology class lookup

    ! hydrology (level1)
    ! canopy
    type(var_dp) :: interception        !< canopy interception storage [mm] on level l1
    type(var_dp) :: throughfall         !< throughfall amount [mm] on level l1
    ! storage and SM
    type(var2d_dp) :: soil_moisture     !< soil water content of soil layer [mm] on level l1
    type(var_dp) :: sealed_storage      !< reservoir of sealed areas [mm] on level l1
    type(var_dp) :: unsat_storage       !< reservoir of unsaturated zone [mm] on level l1
    type(var_dp) :: sat_storage         !< water level in groundwater reservoir [mm] on level l1
    ! type(var_dp) :: water_table_depth   !< depth to water table in groundwater reservoir [m] on level l1
    ! AET
    type(var_dp) :: aet_canopy          !< actual evapotranspiration from canopy [mm] on level l1
    type(var_dp) :: aet_sealed          !< actual evapotranspiration from free water surfaces [mm] on level l1
    type(var2d_dp) :: aet_soil          !< actual evapotranspiration from soil layer [mm] on level l1
    ! rain/snow
    type(var_dp) :: snowpack            !< depth of snowpack [mm] on level l1
    type(var_dp) :: rain                !< rain precipitation [mm] on level l1
    type(var_dp) :: snow                !< snow precipitation [mm] on level l1
    type(var_dp) :: melt                !< melting snow [mm] on level l1
    type(var_dp) :: pre_eff             !< effective precipitation [mm] on level l1 (rain + melt)
    ! vertical soil water movement
    type(var2d_dp) :: infiltration      !< infiltration intensity in soil layer [mm] on level l1
    type(var_dp) :: percolation         !< percolation [mm] on level l1
    ! type(var_dp) :: loss                !< gain/loss flux in a leaking linear reservoir [mm] on level l1
    ! lateral water movement
    type(var_dp) :: runoff_total        !< total runoff [mm] on level l1
    type(var_dp) :: runoff_sealed       !< direct runoff from impervious areas [mm] on level l1
    type(var_dp) :: interflow_fast      !< fast runoff component [mm] on level l1
    type(var_dp) :: interflow_slow      !< slow runoff component [mm] on level l1
    type(var_dp) :: baseflow            !< baseflow [mm] on level l1
    ! neutrons
    type(var_dp) :: neutrons            !< ground albedo neutrons [count h-1] on level l1
    ! degday calculated by mHM from MPR degday_X variables
    type(var_dp) :: degday              !< Degree-day factor [mm d-1 degC-1] on level l1

    ! MPR results (level1)
    ! PET
    type(var_dp) :: pet_coeff_pt        !< PET calculation coefficient for Priestley Taylor (alpha) [1] on level l1
    type(var_dp) :: pet_coeff_hs        !< PET calculation coefficient for Hargreaves Samani [1] on level l1
    type(var_dp) :: pet_fac_aspect      !< PET correction based on aspect [1] on level l1
    type(var_dp) :: pet_fac_lai         !< PET correction based on LAI [1] on level l1
    type(var_dp) :: resist_aero         !< aerodynamical resistance [s m-1] on level l1
    type(var_dp) :: resist_surf         !< bulk surface resistance [s m-1] on level l1
    ! canopy
    type(var_dp) :: max_interception    !< Maximum interception [mm] on level l1
    ! snow
    type(var_dp) :: degday_inc          !< Increase of the degree-day factor per precipitation [d-1 degC-1] on level l1
    type(var_dp) :: degday_max          !< Maximum degree-day factor [mm d-1 degC-1] on level l1
    type(var_dp) :: degday_dry          !< Degree-day factor for no precipitation [mm d-1 degC-1] on level l1
    type(var_dp) :: thresh_temp         !< Threshold temperature for phase transition snow and rain [degC] on level l1
    ! soil moisture
    type(var_dp) :: f_sealed            !< Fraction of sealed area [1] on level l1
    type(var2d_dp) :: f_roots           !< Fraction of roots in soil horizons [1] on level l1
    type(var2d_dp) :: sm_saturation     !< Saturation soil moisture [mm] on level l1
    type(var2d_dp) :: sm_exponent       !< Exponential parameter controlling non-linearity of soil water retention [1] on level l1
    type(var2d_dp) :: sm_field_capacity !< Field capacity - soil moisture below which actual ET is reduced [mm] on level l1
    type(var2d_dp) :: wilting_point     !< permanent wilting point [mm] on level l1
    type(var_dp) :: thresh_jarvis       !< Jarvis critical value (C1) for normalized soil water content [1] on level l1
    ! runoff
    type(var_dp) :: alpha               !< Exponent for the upper reservoir [1] on level l1
    type(var_dp) :: k_fastflow          !< Fast interflow recession coefficient [d-1] on level l1
    type(var_dp) :: k_slowflow          !< Slow interflow recession coefficient [d-1] on level l1
    type(var_dp) :: k_baseflow          !< Baseflow recession coefficient [d-1] on level l1
    type(var_dp) :: k_percolation       !< Percolation coefficient [d-1] on level l1
    type(var_dp) :: f_karst_loss        !< Fraction of karstic percolation loss [1] on level l1
    type(var_dp) :: thresh_unsat        !< Threshold water depth for fast interflow [mm] on level l1
    type(var_dp) :: thresh_sealed       !< Threshold water depth for runoff on sealed surfaces [mm] on level l1
    ! neutrons
    type(var_dp) :: desilets_n0         !< neutron count rate under dry reference conditions (N_0 in Desilets eq.) [count h-1] on level l1
    type(var2d_dp) :: bulk_density      !< bulk density [g cm-3] on level l1
    type(var2d_dp) :: lattice_water     !< Ratio of structurally bound water [g g-1] on level l1
    type(var2d_dp) :: cosmic_l3         !< cosmic L3 parameter [g cm-2] on level l1

    ! routing (level3)
    ! type(var_dp) :: q_out               !< accumulated runoff [m3 s-1] on level L3
    ! type(var_dp) :: e_out               !< accumulated source energy [W] on level L3
    type(var_dp) :: discharge           !< modelled discharge [m3 s-1] on level L3
    type(var_dp) :: energy_flux         !< modelled routed energy [W] on level L3
    type(var_dp) :: river_temp          !< simulated river temperature [degC] on level L3

    ! groundwater (level0)
    type(var_dp) :: riverhead           !< simulated riverhead [m] on level l0

  contains
    procedure, public  :: init => exchange_init
    procedure, public  :: set_dims => exchange_set_dims
    procedure, public  :: configure => exchange_configure
    procedure, public  :: get_grid => exchange_get_grid
    procedure, public  :: has_grid => exchange_has_grid
    procedure, public :: get_meta => exchange_get_var_meta
    procedure, public :: get_path => exchange_get_path
    procedure, private :: get_var_class => exchange_get_var_class
    procedure, private  :: get_data_1d_dp => exchange_get_data_1d_dp
    procedure, private  :: get_data_1d_i4 => exchange_get_data_1d_i4
    procedure, private  :: get_data_1d_lg => exchange_get_data_1d_lg
    procedure, private  :: get_data_2d_dp => exchange_get_data_2d_dp
    procedure, private  :: get_data_2d_i4 => exchange_get_data_2d_i4
    procedure, private  :: get_data_2d_lg => exchange_get_data_2d_lg
    generic, public :: get_data => get_data_1d_dp, get_data_1d_i4, get_data_1d_lg, get_data_2d_dp, get_data_2d_i4, get_data_2d_lg
    procedure, private  :: set_data_1d => exchange_set_data_1d
    procedure, private  :: set_data_2d => exchange_set_data_2d
    generic, public :: set_data => set_data_1d, set_data_2d
  end type exchange_t

contains

  !> \brief Format a path by prepending the current working directory and appending a file name if given.
  function standard_path(cwd, path, file) result(std_path)
    use mo_os, only: path_join, path_normpath, get_cwd
    character(len=*), optional, intent(in) :: cwd !< current working directory (current directory by default)
    character(len=*), optional, intent(in) :: path !< path to be formatted
    character(len=*), optional, intent(in) :: file !< file to be appended (omitted by default)
    character(:), allocatable :: std_path !< formatted absolute path
    character(:), allocatable :: work
    if (present(cwd)) then
      work = cwd
    else
      call get_cwd(work)
    end if
    std_path = path_normpath(path_join(work, path, file))
  end function standard_path

  !> \brief Return the number of domains configured in a meta namelist file.
  function get_n_domains(main_file) result(n_domains)
    character(*), intent(in) :: main_file !< file containing the main namelist
    integer(i4) :: n_domains !< number of configured domains
    type(nml_config_project_t) :: project
    character(1024) :: errmsg
    integer :: status

    status = project%from_file(file=main_file, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Error reading project config: ", trim(errmsg)
      error stop 1
    end if
    status = project%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Project config not valid: ", trim(errmsg)
      error stop 1
    end if
    n_domains = project%n_domains
  end function get_n_domains

  !> \brief Initialize the exchange type
  subroutine exchange_init(self, meta_file, domain, cwd)
    use mo_os, only: path_abspath, check_path_isdir
    class(exchange_t), intent(inout) :: self
    character(*), intent(in), optional :: meta_file !< file containing the run metadata namelists
    integer(i4), intent(in), optional :: domain !< domain ID of the current domain in the configuration arrays (1 by default)
    character(len=*), intent(in), optional :: cwd !< current working directory to set relative paths
    logical :: from_dirs
    character(1024) :: errmsg
    type(nml_config_project_t) :: local_project
    integer(i4) :: id(1)
    integer :: status
    log_info(*) "Configure exchange."

    self%domain_id = optval(domain, 1_i4) ! 1 by default for single domain initialization
    self%cwd = path_abspath(optval(cwd, "."))
    call check_path_isdir(self%cwd, raise=.true.)

    if (present(meta_file)) then
      log_info(*) "Read project attributes: ", meta_file
      status = self%config%project%from_file(file=meta_file, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading project config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%project%is_configured) then
      log_fatal(*) "Project config not set."
      error stop 1
    end if
    status = self%config%project%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Project config not valid: ", trim(errmsg)
      error stop 1
    end if

    self%n_domains = self%config%project%n_domains
    self%n_geo_units = self%config%project%n_geo_units
    self%max_layers = self%config%project%max_layers
    if (self%domain_id < 1_i4 .or. self%domain_id > self%n_domains) then
      log_fatal(*) "Invalid domain id ", n2s(self%domain_id), " for run n_domains=", n2s(self%n_domains), "."
      error stop 1
    end if

    from_dirs = self%config%project%read_domains_from_dirs
    if (from_dirs) then
      if (present(meta_file)) then
        status = self%config%domain%set_dims(n_domains=self%n_domains, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error setting domain config dimensions: ", trim(errmsg)
          error stop 1
        end if
        log_info(*) "Read domain config: ", meta_file
        status = self%config%domain%from_file(file=meta_file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading domain config: ", trim(errmsg)
          error stop 1
        end if
      end if
      if (.not.self%config%domain%is_configured) then
        log_fatal(*) "Domain config not set for read_domains_from_dirs."
        error stop 1
      end if
      status = self%config%domain%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Domain config not valid: ", trim(errmsg)
        error stop 1
      end if

      id(1) = self%domain_id
      status = self%config%domain%is_set("domain_dirs", idx=id, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Directory not specified for domain ", n2s(self%domain_id), ": ", trim(errmsg)
        error stop 1
      end if
      self%cwd = standard_path(cwd=self%cwd, path=trim(self%config%domain%domain_dirs(self%domain_id)))
      call check_path_isdir(self%cwd, raise=.true.)
      self%main_file = standard_path(cwd=self%cwd, file=trim(self%config%domain%domain_nmls(self%domain_id)))

      log_info(*) "Read local project attributes: ", self%main_file
      status = local_project%from_file(file=self%main_file, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading local project config: ", trim(errmsg)
        error stop 1
      end if
      status = local_project%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Local project config not valid: ", trim(errmsg)
        error stop 1
      end if
      if (local_project%n_domains /= 1_i4) then
        log_fatal(*) "Domain-local project config must set n_domains = 1 when read_domains_from_dirs is enabled."
        error stop 1
      end if
      if (local_project%read_domains_from_dirs) then
        log_fatal(*) "Domain-local project config must set read_domains_from_dirs = .false."
        error stop 1
      end if
      if (local_project%n_geo_units > self%n_geo_units) then
        log_fatal(*) "Domain-local n_geo_units=", n2s(local_project%n_geo_units), &
          " exceeds top-level n_geo_units=", n2s(self%n_geo_units), "."
        error stop 1
      end if
      if (local_project%max_layers > self%max_layers) then
        log_fatal(*) "Domain-local max_layers=", n2s(local_project%max_layers), &
          " exceeds top-level max_layers=", n2s(self%max_layers), "."
        error stop 1
      end if
    end if

    if (from_dirs) then
      self%nml_n_domains = 1_i4
      self%nml_domain_id = 1_i4
    else
      self%nml_n_domains = self%n_domains
      self%nml_domain_id = self%domain_id
    end if
    if (self%nml_domain_id < 1_i4 .or. self%nml_domain_id > self%nml_n_domains) then
      log_fatal(*) "Invalid namelist domain index ", n2s(self%nml_domain_id), " for n_domains=", n2s(self%nml_n_domains), "."
      error stop 1
    end if

    ! variables
    ! raw meteorology (level2)
    self%raw_pre    = var_dp(grid=l2, name="pre",       units="mm",    long_name="precipitation", standard_name="precipitation_amount")
    self%raw_pet    = var_dp(grid=l2, name="pet",       units="mm",    long_name="potential evapotranspiration", standard_name="water_potential_evapotranspiration_amount")
    self%raw_temp   = var_dp(grid=l2, name="temp",      units="degC",  long_name="air temperature", standard_name="air_temperature")
    self%raw_ssrd   = var_dp(grid=l2, name="ssrd",      units="W m-2", long_name="solar short wave radiation downward", standard_name="surface_downwelling_shortwave_flux")
    self%raw_strd   = var_dp(grid=l2, name="strd",      units="W m-2", long_name="surface thermal radiation downward", standard_name="surface_downwelling_longwave_flux")
    self%raw_tann   = var_dp(grid=l2, name="tann",      units="degC",  long_name="annual mean air temperature", standard_name="air_temperature")
    self%raw_tmin   = var_dp(grid=l2, name="tmin",      units="degC",  long_name="minimum daily temperature", standard_name="air_temperature")
    self%raw_tmax   = var_dp(grid=l2, name="tmax",      units="degC",  long_name="maximum daily temperature", standard_name="air_temperature")
    self%raw_netrad = var_dp(grid=l2, name="netrad",    units="W m-2", long_name="net radiation", standard_name="surface_net_downward_radiative_flux")
    self%raw_eabs   = var_dp(grid=l2, name="eabs",      units="Pa",    long_name="vapor pressure", standard_name="water_vapor_pressure")
    self%raw_wind   = var_dp(grid=l2, name="windspeed", units="m s-1", long_name="wind speed", standard_name="wind_speed")

    ! processed meteorology (level1)
    self%pre  = var_dp(grid=l1, name="pre",  units="mm",    long_name="precipitation", standard_name="precipitation_amount")
    self%temp = var_dp(grid=l1, name="temp", units="degC",  long_name="air temperature", standard_name="air_temperature")
    self%pet  = var_dp(grid=l1, name="pet",  units="mm",    long_name="potential evapotranspiration", standard_name="water_potential_evapotranspiration_amount")
    self%ssrd = var_dp(grid=l1, name="ssrd", units="W m-2", long_name="solar short wave radiation downward", standard_name="surface_downwelling_shortwave_flux")
    self%strd = var_dp(grid=l1, name="strd", units="W m-2", long_name="surface thermal radiation downward", standard_name="surface_downwelling_longwave_flux")
    self%tann = var_dp(grid=l1, name="tann", units="degC",  long_name="annual mean air temperature", standard_name="air_temperature")

    ! morphology (level0)
    self%dem    = var_dp(static=.true., grid=l0, name="dem",    units="m",      long_name="elevation", standard_name="height_above_mean_sea_level")
    self%slope  = var_dp(static=.true., grid=l0, name="slope",  units="%",      long_name="slope", standard_name="ground_slope_angle")
    self%aspect = var_dp(static=.true., grid=l0, name="aspect", units="degree", long_name="aspect", standard_name="ground_slope_direction")
    self%fdir   = var_i4(static=.true., grid=l0, name="fdir",   units="1",      long_name="flow direction")
    self%facc   = var_i4(static=.true., grid=l0, name="facc",   units="1",      long_name="flow accumulation")
    self%soil_id = var2d_i4(static=.true., grid=l0, name="soil_id", units="1", long_name="soil class ID")
    self%geo_unit = var_i4(static=.true., grid=l0, name="geo_unit", units="1", long_name="geological unit ID")
    self%lai_class = var_i4(static=.true., grid=l0, name="lai_class", units="1", long_name="LAI class ID")
    self%slope_emp = var_dp(static=.true., grid=l0, name="slope_emp", units="1", long_name="empirical slope distribution")

    ! hydrology (level1)
    ! canopy
    self%interception      =   var_dp(grid=l1, name="interception",      units="mm",  long_name="canopy interception storage")
    self%throughfall       =   var_dp(grid=l1, name="throughfall",       units="mm",  long_name="throughfall amount")
    ! storage and SM
    self%soil_moisture     = var2d_dp(grid=l1, name="soil_moisture",     units="mm",  long_name="soil water content of soil layer")
    self%sealed_storage    =   var_dp(grid=l1, name="sealedSTW",         units="mm",  long_name="reservoir of sealed areas")
    self%unsat_storage     =   var_dp(grid=l1, name="unsatSTW",          units="mm",  long_name="reservoir of unsaturated zone")
    self%sat_storage       =   var_dp(grid=l1, name="satSTW",            units="mm",  long_name="water level in groundwater reservoir")
    ! AET
    self%aet_canopy        =   var_dp(grid=l1, name="aet_canopy",        units="mm",  long_name="actual evapotranspiration from canopy")
    self%aet_sealed        =   var_dp(grid=l1, name="aet_sealed",        units="mm",  long_name="actual evapotranspiration from free water surfaces")
    self%aet_soil          = var2d_dp(grid=l1, name="aet_soil",          units="mm",  long_name="actual evapotranspiration from soil layer")
    ! rain/snow
    self%snowpack          =   var_dp(grid=l1, name="snowpack",          units="mm",  long_name="depth of snowpack", standard_name="surface_snow_amount")
    self%rain              =   var_dp(grid=l1, name="rain",              units="mm",  long_name="rain precipitation", standard_name="rainfall_amount")
    self%snow              =   var_dp(grid=l1, name="snow",              units="mm",  long_name="snow precipitation", standard_name="snowfall_amount")
    self%melt              =   var_dp(grid=l1, name="melt",              units="mm",  long_name="melting snow", standard_name="surface_snow_melt_amount")
    self%pre_eff           =   var_dp(grid=l1, name="pre_eff",           units="mm",  long_name="effective precipitation") ! rain + melt
    ! vertical soil water movement
    self%infiltration      = var2d_dp(grid=l1, name="infiltration",      units="mm",  long_name="infiltration intensity in soil layer")
    self%percolation       =   var_dp(grid=l1, name="percolation",       units="mm",  long_name="percolation")
    ! lateral water movement
    self%runoff_total      =   var_dp(grid=l1, name="Q",                 units="mm",  long_name="total runoff", standard_name="runoff_amount")
    self%runoff_sealed     =   var_dp(grid=l1, name="QD",                units="mm",  long_name="direct runoff from impervious areas", standard_name="surface_runoff_amount")
    self%interflow_fast    =   var_dp(grid=l1, name="QIf",               units="mm",  long_name="fast runoff component", standard_name="subsurface_runoff_amount")
    self%interflow_slow    =   var_dp(grid=l1, name="QIs",               units="mm",  long_name="slow runoff component", standard_name="subsurface_runoff_amount")
    self%baseflow          =   var_dp(grid=l1, name="QB",                units="mm",  long_name="baseflow", standard_name="baseflow_amount")
    ! neutrons
    self%neutrons          =   var_dp(grid=l1, name="neutrons",          units="cph", long_name="ground albedo neutrons")

    ! MPR results (level1)
    ! PET
    self%pet_coeff_pt      =   var_dp(grid=l1, name="pet_coeff_pt",      units="1",                 long_name="PET calculation coefficient for Priestley Taylor (alpha)")
    self%pet_coeff_hs      =   var_dp(grid=l1, name="pet_coeff_hs",      units="1", static=.true.,  long_name="PET calculation coefficient for Hargreaves Samani")
    self%pet_fac_aspect    =   var_dp(grid=l1, name="pet_fac_aspect",    units="1", static=.true.,  long_name="PET correction factor based on aspect")
    self%pet_fac_lai       =   var_dp(grid=l1, name="pet_fac_lai",       units="1",                 long_name="PET correction factor based on LAI")
    self%resist_aero       =   var_dp(grid=l1, name="resist_aero",       units="s m-1",             long_name="aerodynamical resistance")
    self%resist_surf       =   var_dp(grid=l1, name="resist_surf",       units="s m-1",             long_name="bulk surface resistance")
    ! canopy
    self%max_interception  =   var_dp(grid=l1, name="max_interception",  units="mm",                long_name="Maximum interception")
    ! snow
    self%degday_inc        =   var_dp(grid=l1, name="degday_inc",        units="d-1 degC-1",       long_name="Increase of the degree-day factor per precipitation")
    self%degday_max        =   var_dp(grid=l1, name="degday_max",        units="mm d-1 degC-1",    long_name="Maximum degree-day factor")
    self%degday_dry        =   var_dp(grid=l1, name="degday_dry",        units="mm d-1 degC-1",    long_name="Degree-day factor for no precipitation")
    self%thresh_temp       =   var_dp(grid=l1, name="thresh_temp",       units="degC",              long_name="Threshold temperature for phase transition snow and rain")
    ! soil moisture
    self%f_sealed          =   var_dp(grid=l1, name="f_sealed",          units="1",                 long_name="Fraction of sealed area")
    self%f_roots           = var2d_dp(grid=l1, name="f_roots",           units="1",                 long_name="Fraction of roots in soil horizons")
    self%sm_saturation     = var2d_dp(grid=l1, name="sm_saturation",     units="mm",                long_name="Saturation soil moisture")
    self%sm_exponent       = var2d_dp(grid=l1, name="sm_exponent",       units="1",                 long_name="Exponential parameter controlling non-linearity of soil water retention")
    self%sm_field_capacity = var2d_dp(grid=l1, name="sm_field_capacity", units="mm",                long_name="Field capacity - soil moisture below which actual ET is reduced")
    self%wilting_point     = var2d_dp(grid=l1, name="wilting_point",     units="mm",                long_name="permanent wilting point")
    self%thresh_jarvis     =   var_dp(grid=l1, name="thresh_jarvis",     units="1",  static=.true., long_name="Jarvis critical value (C1) for normalized soil water content")
    ! runoff
    self%alpha             =   var_dp(grid=l1, name="alpha",             units="1",                 long_name="Exponent for the upper reservoir")
    self%k_fastflow        =   var_dp(grid=l1, name="k_fastflow",        units="d-1",              long_name="Fast interflow recession coefficient")
    self%k_slowflow        =   var_dp(grid=l1, name="k_slowflow",        units="d-1",              long_name="Slow interflow recession coefficient")
    self%k_baseflow        =   var_dp(grid=l1, name="k_baseflow",        units="d-1",              long_name="Baseflow recession coefficient")
    self%k_percolation     =   var_dp(grid=l1, name="k_percolation",     units="d-1",              long_name="Percolation coefficient")
    self%f_karst_loss      =   var_dp(grid=l1, name="f_karst_loss",      units="1",  static=.true., long_name="Fraction of karstic percolation loss")
    self%thresh_unsat      =   var_dp(grid=l1, name="thresh_unsat",      units="mm", static=.true., long_name="Threshold water depth for fast interflow")
    self%thresh_sealed     =   var_dp(grid=l1, name="thresh_sealed",     units="mm", static=.true., long_name="Threshold water depth for runoff on sealed surfaces")
    ! neutrons
    self%desilets_n0       =   var_dp(grid=l1, name="desilets_n0",       units="count h-1", static=.true., long_name="neutron count rate under dry reference conditions (N_0 in Desilets eq.)")
    self%bulk_density      = var2d_dp(grid=l1, name="bulk_density",      units="g cm-3",            long_name="bulk density")
    self%lattice_water     = var2d_dp(grid=l1, name="lattice_water",     units="g g-1",             long_name="Ratio of structurally bound water")
    self%cosmic_l3         = var2d_dp(grid=l1, name="cosmic_l3",         units="g cm-2",            long_name="cosmic L3 parameter")

    ! routing (level3)
    ! self%q_out             =   var_dp(grid=l3, name="q_out",            units="m3 s-1",            long_name="accumulated runoff")
    ! self%e_out             =   var_dp(grid=l3, name="e_out",            units="W",                 long_name="accumulated source energy")
    self%discharge         =   var_dp(grid=l3, name="discharge",        units="m3 s-1",            long_name="modelled discharge", standard_name="outgoing_water_volume_transport_along_river_channel")
    self%energy_flux       =   var_dp(grid=l3, name="e_mod",            units="W",                 long_name="modelled routed energy")
    self%river_temp        =   var_dp(grid=l3, name="river_temp",       units="degC",              long_name="simulated river temperature")

    ! groundwater (level0)
    self%riverhead         =   var_dp(grid=l0,  name="riverhead",        units="m",                 long_name="simulated riverhead")
  end subroutine exchange_init

  !> \brief Set runtime dimensions for generated exchange-owned namelists.
  subroutine exchange_set_dims(self)
    class(exchange_t), intent(inout) :: self
    character(1024) :: errmsg
    integer :: status

    if (.not.self%config%time%is_configured) then
      status = self%config%time%set_dims(n_domains=self%nml_n_domains, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error setting time config dimensions: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%resolution%is_configured) then
      status = self%config%resolution%set_dims(n_domains=self%nml_n_domains, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error setting resolution config dimensions: ", trim(errmsg)
        error stop 1
      end if
    end if
    call self%parameters%set_dims(self%n_geo_units)
  end subroutine exchange_set_dims

  !> \brief Read exchange-owned namelists after runtime dimensions have been set.
  subroutine exchange_configure(self, para_file)
    class(exchange_t), intent(inout) :: self
    character(*), intent(in), optional :: para_file !< file containing the parameter namelists
    integer(i4) :: id(1)
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer :: status

    if (allocated(self%main_file)) then
      path = self%main_file
    end if
    if (allocated(path)) then
      log_info(*) "Read time config: ", path
      status = self%config%time%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading time config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%time%is_configured) then
      log_fatal(*) "Time config not set."
      error stop 1
    end if
    status = self%config%time%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Time config not valid: ", trim(errmsg)
      error stop 1
    end if

    ! parameters are created redundantly for each exchange instance
    ! but this simplifies the code structure
    if (allocated(path)) then
      call self%parameters%configure(main_file=path, para_file=para_file)
    else
      call self%parameters%configure(para_file=para_file)
    end if

    if (allocated(self%main_file)) then
      path = self%main_file
    else if (allocated(path)) then
      deallocate(path)
    end if
    if (allocated(path)) then
      log_info(*) "Read resolution config: ", path
      status = self%config%resolution%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading resolution config: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%resolution%is_configured) then
      log_fatal(*) "Resolution configuration not set."
      error stop 1
    end if
    status = self%config%resolution%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Resolution config not valid: ", trim(errmsg)
      error stop 1
    end if

    id(1) = self%nml_domain_id
    status = self%config%resolution%is_set("hydro", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      self%level1_resolution = self%config%resolution%hydro(id(1))
      log_info(*) "Set hydro resolution for domain ", n2s(id(1)), ": ", n2s(self%level1_resolution)
    else if (self%parameters%meteo_active() .or. self%parameters%mhm_active()) then
      log_fatal(*) "Hydro resolution not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg)
      error stop 1
    end if

    status = self%config%resolution%is_set("route", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      self%level3_resolution = self%config%resolution%route(id(1))
      log_info(*) "Set route resolution for domain ", n2s(id(1)), ": ", n2s(self%level3_resolution)
    else if (self%parameters%mrm_active()) then
      log_fatal(*) "Route resolution not set for domain ", n2s(id(1)), ". Error: ", trim(errmsg)
      error stop 1
    end if

    ! time settings
    id(1) = self%nml_domain_id
    if (self%config%time%share_time_period) id(1) = 1_i4
    status = self%config%time%is_set("sim_start", idx=id, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Simulation start time input error: ", trim(errmsg)
      error stop 1
    end if
    self%start_time = datetime(self%config%time%sim_start(id(1)))

    status = self%config%time%is_set("sim_end", idx=id, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Simulation end time input error: ", trim(errmsg)
      error stop 1
    end if
    self%end_time = datetime(self%config%time%sim_end(id(1)))

    status = self%config%time%is_set("eval_start", idx=id, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Evaluation start time input error: ", trim(errmsg)
      error stop 1
    end if
    self%eval_start_time = datetime(self%config%time%eval_start(id(1)))

    id(1) = self%nml_domain_id
    if (self%config%time%share_time_step) id(1) = 1_i4
    self%step = timedelta(hours=self%config%time%time_step(id(1)))

    self%step_count = 0_i4
    self%time = self%start_time
  end subroutine exchange_configure

  !> \brief get the grid specifications for the selected level
  subroutine exchange_get_grid(self, selector, grid)
    use mo_message, only: error_message
    class(exchange_t), intent(in) :: self
    integer(i4), intent(in) :: selector !< level selector (0: L0, 1: L1, 2: L2, 3: L3, -1: nogrid)
    type(grid_t), pointer, intent(out) :: grid !< resulting pointer to the selected grid
    select case(selector)
      case(nogrid)
        grid => null() ! exchangable
      case(l0)
        grid => self%level0
      case(l1)
        grid => self%level1
      case(l2)
        grid => self%level2
      case(l3)
        grid => self%level3
      case default
        log_fatal(*) "exchange%get_grid: unknown grid selector '", n2s(selector), "'."
        error stop 1
    end select
  end subroutine exchange_get_grid

  !> \brief get the grid specifications for the selected level
  logical function exchange_has_grid(self, selector)
    use mo_message, only: error_message
    class(exchange_t), intent(in) :: self
    integer(i4), intent(in) :: selector !< level selector (0: l0, 1: l1, 2: l2, 3: L3, -1: nogrid)
    select case(selector)
      case(l0)
        exchange_has_grid = associated(self%level0)
      case(l1)
        exchange_has_grid = associated(self%level1)
      case(l2)
        exchange_has_grid = associated(self%level2)
      case(l3)
        exchange_has_grid = associated(self%level3)
      case default
        exchange_has_grid = .false.
    end select
  end function exchange_has_grid

  !> \brief get class pointer to a variable
  subroutine exchange_get_var_class(self, var, var_pnt)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    class(*), pointer, intent(out) :: var_pnt !< resulting pointer to the selected variable
    select case(var)
      case("raw_pre")
        var_pnt => self%raw_pre
      case("raw_temp")
        var_pnt => self%raw_temp
      case("raw_ssrd")
        var_pnt => self%raw_ssrd
      case("raw_strd")
        var_pnt => self%raw_strd
      case("raw_tann")
        var_pnt => self%raw_tann
      case("raw_tmin")
        var_pnt => self%raw_tmin
      case("raw_tmax")
        var_pnt => self%raw_tmax
      case("raw_netrad")
        var_pnt => self%raw_netrad
      case("raw_eabs")
        var_pnt => self%raw_eabs
      case("raw_wind")
        var_pnt => self%raw_wind
      ! processed meteorology (level1)
      case("pre")
        var_pnt => self%pre
      case("temp")
        var_pnt => self%temp
      case("pet")
        var_pnt => self%pet
      case("ssrd")
        var_pnt => self%ssrd
      case("strd")
        var_pnt => self%strd
      case("tann")
        var_pnt => self%tann
      ! morphology (level0)
      case("dem")
        var_pnt => self%dem
      case("slope")
        var_pnt => self%slope
      case("aspect")
        var_pnt => self%aspect
      case("fdir")
        var_pnt => self%fdir
      case("facc")
        var_pnt => self%facc
      case("soil_id")
        var_pnt => self%soil_id
      case("geo_unit")
        var_pnt => self%geo_unit
      case("lai_class")
        var_pnt => self%lai_class
      case("slope_emp")
        var_pnt => self%slope_emp
      ! hydrology (level1)
      ! canopy
      case("interception")
        var_pnt => self%interception
      case("throughfall")
        var_pnt => self%throughfall
      ! storage and SM
      case("soil_moisture")
        var_pnt => self%soil_moisture
      case("sealed_storage")
        var_pnt => self%sealed_storage
      case("unsat_storage")
        var_pnt => self%unsat_storage
      case("sat_storage")
        var_pnt => self%sat_storage
      ! case("water_table_depth")
      !   var => self%water_table_depth
      ! AET
      case("aet_canopy")
        var_pnt => self%aet_canopy
      case("aet_sealed")
        var_pnt => self%aet_sealed
      case("aet_soil")
        var_pnt => self%aet_soil
      ! rain/snow
      case("snowpack")
        var_pnt => self%snowpack
      case("rain")
        var_pnt => self%rain
      case("snow")
        var_pnt => self%snow
      case("melt")
        var_pnt => self%melt
      case("pre_eff")
        var_pnt => self%pre_eff
      ! vertical soil water movement
      case("infiltration")
        var_pnt => self%infiltration
      case("percolation")
        var_pnt => self%percolation
      ! case("loss")
      !   var => self%loss
      ! lateral water movement
      case("runoff_total")
        var_pnt => self%runoff_total
      case("runoff_sealed")
        var_pnt => self%runoff_sealed
      case("interflow_fast")
        var_pnt => self%interflow_fast
      case("interflow_slow")
        var_pnt => self%interflow_slow
      case("baseflow")
        var_pnt => self%baseflow
      ! neutrons
      case("neutrons")
        var_pnt => self%neutrons
      ! degday calculated by mHM from MPR degday_X variables
      case("degday")
        var_pnt => self%degday
      ! MPR results (level1)
      ! PET
      case("pet_coeff_pt")
        var_pnt => self%pet_coeff_pt
      case("pet_coeff_hs")
        var_pnt => self%pet_coeff_hs
      case("pet_fac_aspect")
        var_pnt => self%pet_fac_aspect
      case("pet_fac_lai")
        var_pnt => self%pet_fac_lai
      case("resist_aero")
        var_pnt => self%resist_aero
      case("resist_surf")
        var_pnt => self%resist_surf
      ! canopy
      case("max_interception")
        var_pnt => self%max_interception
      ! snow
      case("degday_inc")
        var_pnt => self%degday_inc
      case("degday_max")
        var_pnt => self%degday_max
      case("degday_dry")
        var_pnt => self%degday_dry
      case("thresh_temp")
        var_pnt => self%thresh_temp
      ! soil moisture
      case("f_sealed")
        var_pnt => self%f_sealed
      case("f_roots")
        var_pnt => self%f_roots
      case("sm_saturation")
        var_pnt => self%sm_saturation
      case("sm_exponent")
        var_pnt => self%sm_exponent
      case("sm_field_capacity")
        var_pnt => self%sm_field_capacity
      case("wilting_point")
        var_pnt => self%wilting_point
      case("thresh_jarvis")
        var_pnt => self%thresh_jarvis
      ! runoff
      case("alpha")
        var_pnt => self%alpha
      case("k_fastflow")
        var_pnt => self%k_fastflow
      case("k_slowflow")
        var_pnt => self%k_slowflow
      case("k_baseflow")
        var_pnt => self%k_baseflow
      case("k_percolation")
        var_pnt => self%k_percolation
      case("f_karst_loss")
        var_pnt => self%f_karst_loss
      case("thresh_unsat")
        var_pnt => self%thresh_unsat
      case("thresh_sealed")
        var_pnt => self%thresh_sealed
      ! neutrons
      case("desilets_n0")
        var_pnt => self%desilets_n0
      case("bulk_density")
        var_pnt => self%bulk_density
      case("lattice_water")
        var_pnt => self%lattice_water
      case("cosmic_l3")
        var_pnt => self%cosmic_l3
      ! routing (level3)
      ! case("q_out")
      !   var_pnt => self%q_out
      ! case("e_out")
      !   var_pnt => self%e_out
      case("discharge")
        var_pnt => self%discharge
      case("energy_flux")
        var_pnt => self%energy_flux
      case("river_temp")
        var_pnt => self%river_temp
      ! groundwater (level0)
      case("riverhead")
        var_pnt => self%riverhead
      case default
        log_fatal(*) "exchange%get_var: variable '", var, "' not available."
        error stop 1
    end select
  end subroutine exchange_get_var_class

  !> \brief get var_dp pointer to a variable
  subroutine exchange_get_var_meta(self, var, name, units, long_name, standard_name, grid, static, provided, required)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var                                   !< name of the variable (attribute name)
    character(:), allocatable, intent(out), optional :: name          !< variable name
    character(:), allocatable, intent(out), optional :: units         !< variable unit
    character(:), allocatable, intent(out), optional :: long_name     !< long name of the variable
    character(:), allocatable, intent(out), optional :: standard_name !< standard name of the variable
    integer(i4), intent(out), optional :: grid                        !< ID of the grid the data is defined on
    logical, intent(out), optional :: static                          !< flag to indicated static data
    logical, intent(out), optional :: provided                        !< flag to indicate that data is provided by a component
    logical, intent(out), optional :: required                        !< flag to indicate that data is required by a component
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (variable_abc)
        if (present(name)          .and. allocated(tmp%name))          name          = tmp%name
        if (present(units)         .and. allocated(tmp%units))         units         = tmp%units
        if (present(long_name)     .and. allocated(tmp%long_name))     long_name     = tmp%long_name
        if (present(standard_name) .and. allocated(tmp%standard_name)) standard_name = tmp%standard_name
        if (present(grid))     grid     = tmp%grid
        if (present(static))   static   = tmp%static
        if (present(provided)) provided = tmp%provided
        if (present(required)) required = tmp%required
    end select
  end subroutine exchange_get_var_meta

  !> \brief get pointer to the 1D variable data
  subroutine exchange_get_data_1d_dp(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    real(dp), pointer, intent(out) :: data(:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var_dp)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 1D real(dp)."
        error stop 1
    end select
  end subroutine exchange_get_data_1d_dp

  !> \brief get pointer to the 1D variable data
  subroutine exchange_get_data_1d_i4(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    integer(i4), pointer, intent(out) :: data(:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var_i4)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 1D integer(i4)."
        error stop 1
    end select
  end subroutine exchange_get_data_1d_i4

  !> \brief get pointer to the 1D variable data
  subroutine exchange_get_data_1d_lg(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    logical, pointer, intent(out) :: data(:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var_lg)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 1D logical."
        error stop 1
    end select
  end subroutine exchange_get_data_1d_lg

  !> \brief get pointer to the 2D variable data
  subroutine exchange_get_data_2d_dp(self, var, data)
  use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    real(dp), pointer, intent(out) :: data(:,:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var2d_dp)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 2D real(dp)."
        error stop 1
    end select
  end subroutine exchange_get_data_2d_dp

  !> \brief get pointer to the 2D variable data
  subroutine exchange_get_data_2d_i4(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    integer(i4), pointer, intent(out) :: data(:,:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var2d_i4)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 2D integer(i4)."
        error stop 1
    end select
  end subroutine exchange_get_data_2d_i4

  !> \brief get pointer to the 2D variable data
  subroutine exchange_get_data_2d_lg(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(in) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    logical, pointer, intent(out) :: data(:,:) !< resulting pointer to the selected variable data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var2d_lg)
        data => tmp%data
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not 2D logical."
        error stop 1
    end select
  end subroutine exchange_get_data_2d_lg

  !> \brief set pointer to the 1D variable data
  subroutine exchange_set_data_1d(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(inout) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    class(*), target, intent(in) :: data(:) !< target data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var_dp)
        select type (data)
          type is (real(dp))
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type real(dp)."
            error stop 1
        end select
      class is (var_i4)
        select type (data)
          type is (integer(i4))
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type integer(i4)."
            error stop 1
        end select
      class is (var_lg)
        select type (data)
          type is (logical)
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type logical."
            error stop 1
        end select
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not one dimensional."
        error stop 1
    end select
  end subroutine exchange_set_data_1d

  !> \brief set target for variable data 2D
  subroutine exchange_set_data_2d(self, var, data)
    use mo_message, only: error_message
    class(exchange_t), target, intent(inout) :: self ! target attribute valid here since Fortran 2003
    character(*), intent(in) :: var !< name of the variable (attribute name)
    class(*), target, intent(in) :: data(:,:) !< target data
    class(*), pointer :: tmp
    call self%get_var_class(var, tmp)
    select type (tmp)
      class is (var2d_dp)
        select type (data)
          type is (real(dp))
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type real(dp)."
            error stop 1
        end select
      class is (var2d_i4)
        select type (data)
          type is (integer(i4))
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type integer(i4)."
            error stop 1
        end select
      class is (var2d_lg)
        select type (data)
          type is (logical)
            tmp%data => data
          class default
            log_fatal(*) "exchange%get_var: variable data of '", var, "' is of type logical."
            error stop 1
        end select
      class default
        log_fatal(*) "exchange%get_var: variable data of '", var, "' not two dimensional."
        error stop 1
    end select
  end subroutine exchange_set_data_2d

  !> \brief Format a path by prepending the current working directory and appending a file name if given.
  function exchange_get_path(self, path, file) result(norm_path)
    class(exchange_t), intent(in) :: self
    character(len=*), intent(in) :: path !< path to be formatted
    character(len=*), optional, intent(in) :: file !< file to be appended
    character(:), allocatable :: norm_path !< formatted path
    norm_path = standard_path(self%cwd, path, file)
  end function exchange_get_path

  !> \brief Return whether an exchange variable is already available or will be owned by the caller.
  logical function variable_available(self, owned)
    class(variable_abc), intent(in) :: self
    logical, intent(in), optional :: owned !< caller owns the field and may bypass provided/data availability checks

    variable_available = optval(owned, .false.) .or. (self%provided .and. self%has_data())
  end function variable_available

  !> \brief Mark an exchange variable as required and validate source/data availability when needed.
  subroutine variable_require(self, component, required, expected_shape, check_data)
    class(variable_abc), intent(inout) :: self
    character(*), intent(in) :: component           !< calling component name for diagnostics
    logical, intent(in) :: required                 !< whether this variable is required by the caller
    integer(i8), intent(in), optional :: expected_shape(:) !< expected data shape for contract validation
    logical, intent(in), optional :: check_data     !< enforce data-pointer association check when .true.

    self%required = self%required .or. required
    if (.not.required) return
    call variable_validate(self, component, .false., expected_shape, optval(check_data, .true.))
  end subroutine variable_require

  !> \brief Mark an exchange variable as a required handoff and validate data availability when needed.
  subroutine variable_expect_handoff(self, component, required, expected_shape)
    class(variable_abc), intent(inout) :: self
    character(*), intent(in) :: component           !< calling component name for diagnostics
    logical, intent(in) :: required                 !< whether a handoff is required by the caller
    integer(i8), intent(in), optional :: expected_shape(:) !< expected data shape for handoff validation

    self%required = self%required .or. required
    if (.not.required) return
    call variable_validate(self, component, .true., expected_shape, .true.)
  end subroutine variable_expect_handoff

  !> \brief Clear a variable publication when the caller owns the exchange field.
  subroutine variable_clear(self, owned)
    class(variable_abc), intent(inout) :: self
    logical, intent(in), optional :: owned !< caller ownership flag; only owned publications are cleared

    if (.not.optval(owned, .false.)) return
    call self%clear_data()
    self%provided = .false.
  end subroutine variable_clear

  !> \brief Validate the provided/data/shape contract of an exchange variable.
  subroutine variable_validate(self, component, handoff, expected_shape, check_data)
    class(variable_abc), intent(in) :: self
    character(*), intent(in) :: component
    logical, intent(in) :: handoff
    integer(i8), intent(in), optional :: expected_shape(:)
    logical, intent(in) :: check_data
    integer(i8), allocatable :: actual_shape(:)
    character(:), allocatable :: var_name
    logical :: shape_ok

    var_name = variable_name(self)
    if (.not.self%provided) then
      if (handoff) then
        log_fatal(*) trim(component), ": required handoff not provided: ", var_name, "."
      else
        log_fatal(*) trim(component), ": ", var_name, " not provided."
      end if
      error stop 1
    end if
    if (.not.check_data) return
    if (.not.self%has_data()) then
      if (handoff) then
        log_fatal(*) trim(component), ": required handoff data not connected: ", var_name, "."
      else
        log_fatal(*) trim(component), ": ", var_name, " data not connected."
      end if
      error stop 1
    end if
    if (present(expected_shape)) then
      actual_shape = self%data_shape()
      shape_ok = size(actual_shape) == size(expected_shape)
      if (shape_ok) shape_ok = all(actual_shape == expected_shape)
      if (.not.shape_ok) then
        if (handoff) then
          log_fatal(*) trim(component), ": handoff ", var_name, " has unexpected shape. Expected ", &
            variable_shape_string(expected_shape), ", got ", variable_shape_string(actual_shape), "."
        else
          log_fatal(*) trim(component), ": ", var_name, " has unexpected shape. Expected ", &
            variable_shape_string(expected_shape), ", got ", variable_shape_string(actual_shape), "."
        end if
        error stop 1
      end if
    end if
  end subroutine variable_validate

  !> \brief Validate that a publication target is not already occupied.
  subroutine variable_validate_publish_target(self, component)
    class(variable_abc), intent(in) :: self
    character(*), intent(in) :: component

    if (self%provided .or. self%has_data()) then
      log_fatal(*) trim(component), ": exchange field already provided before publication: ", variable_name(self), "."
      error stop 1
    end if
  end subroutine variable_validate_publish_target

  !> \brief Validate that an alias source is already connected.
  subroutine variable_validate_alias_source(source, component, target)
    class(variable_abc), intent(in) :: source
    character(*), intent(in) :: component
    class(variable_abc), intent(in) :: target

    if (.not.source%provided .or. .not.source%has_data()) then
      log_fatal(*) trim(component), ": pass-through source not provided for ", variable_name(target), "."
      error stop 1
    end if
  end subroutine variable_validate_alias_source

  !> \brief Return the configured variable name or a fallback for diagnostics.
  function variable_name(self) result(name)
    class(variable_abc), intent(in) :: self
    character(:), allocatable :: name
    name = "<unnamed>"
    if (allocated(self%name)) then
      if (len_trim(self%name) > 0) name = trim(self%name)
    end if
  end function variable_name

  !> \brief Format a shape vector for diagnostics.
  function variable_shape_string(shape) result(text)
    integer(i8), intent(in) :: shape(:)
    character(:), allocatable :: text
    integer :: i

    text = "["
    do i = 1, size(shape)
      if (i > 1) text = text // ", "
      text = text // trim(n2s(shape(i)))
    end do
    text = text // "]"
  end function variable_shape_string

  !> \brief Return whether a 1D real exchange variable has data connected.
  logical function var_dp_has_data(self)
    class(var_dp), intent(in) :: self

    var_dp_has_data = associated(self%data)
  end function var_dp_has_data

  !> \brief Return the data shape of a 1D real exchange variable.
  function var_dp_data_shape(self) result(shape)
    class(var_dp), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var_dp_data_shape

  !> \brief Clear the data pointer of a 1D real exchange variable.
  subroutine var_dp_clear_data(self)
    class(var_dp), intent(inout) :: self

    nullify(self%data)
  end subroutine var_dp_clear_data

  !> \brief Publish a local 1D real field through the exchange variable.
  subroutine var_dp_publish_local(self, component, local)
    class(var_dp), intent(inout) :: self
    character(*), intent(in) :: component    !< publishing component name for diagnostics
    real(dp), intent(inout), target :: local(:) !< local 1D real field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var_dp_publish_local

  !> \brief Publish a 1D real alias through the exchange variable.
  subroutine var_dp_publish_alias(self, component, source)
    class(var_dp), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var_dp), intent(in) :: source    !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var_dp_publish_alias

  !> \brief Return whether a 1D integer exchange variable has data connected.
  logical function var_i4_has_data(self)
    class(var_i4), intent(in) :: self

    var_i4_has_data = associated(self%data)
  end function var_i4_has_data

  !> \brief Return the data shape of a 1D integer exchange variable.
  function var_i4_data_shape(self) result(shape)
    class(var_i4), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var_i4_data_shape

  !> \brief Clear the data pointer of a 1D integer exchange variable.
  subroutine var_i4_clear_data(self)
    class(var_i4), intent(inout) :: self

    nullify(self%data)
  end subroutine var_i4_clear_data

  !> \brief Publish a local 1D integer field through the exchange variable.
  subroutine var_i4_publish_local(self, component, local)
    class(var_i4), intent(inout) :: self
    character(*), intent(in) :: component       !< publishing component name for diagnostics
    integer(i4), intent(inout), target :: local(:) !< local 1D integer field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var_i4_publish_local

  !> \brief Publish a 1D integer alias through the exchange variable.
  subroutine var_i4_publish_alias(self, component, source)
    class(var_i4), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var_i4), intent(in) :: source    !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var_i4_publish_alias

  !> \brief Return whether a 1D logical exchange variable has data connected.
  logical function var_lg_has_data(self)
    class(var_lg), intent(in) :: self

    var_lg_has_data = associated(self%data)
  end function var_lg_has_data

  !> \brief Return the data shape of a 1D logical exchange variable.
  function var_lg_data_shape(self) result(shape)
    class(var_lg), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var_lg_data_shape

  !> \brief Clear the data pointer of a 1D logical exchange variable.
  subroutine var_lg_clear_data(self)
    class(var_lg), intent(inout) :: self

    nullify(self%data)
  end subroutine var_lg_clear_data

  !> \brief Publish a local 1D logical field through the exchange variable.
  subroutine var_lg_publish_local(self, component, local)
    class(var_lg), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    logical, intent(inout), target :: local(:) !< local 1D logical field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var_lg_publish_local

  !> \brief Publish a 1D logical alias through the exchange variable.
  subroutine var_lg_publish_alias(self, component, source)
    class(var_lg), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var_lg), intent(in) :: source    !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var_lg_publish_alias

  !> \brief Return whether a 2D real exchange variable has data connected.
  logical function var2d_dp_has_data(self)
    class(var2d_dp), intent(in) :: self

    var2d_dp_has_data = associated(self%data)
  end function var2d_dp_has_data

  !> \brief Return the data shape of a 2D real exchange variable.
  function var2d_dp_data_shape(self) result(shape)
    class(var2d_dp), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8), size(self%data, 2, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var2d_dp_data_shape

  !> \brief Clear the data pointer of a 2D real exchange variable.
  subroutine var2d_dp_clear_data(self)
    class(var2d_dp), intent(inout) :: self

    nullify(self%data)
  end subroutine var2d_dp_clear_data

  !> \brief Publish a local 2D real field through the exchange variable.
  subroutine var2d_dp_publish_local(self, component, local)
    class(var2d_dp), intent(inout) :: self
    character(*), intent(in) :: component   !< publishing component name for diagnostics
    real(dp), intent(inout), target :: local(:, :) !< local 2D real field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var2d_dp_publish_local

  !> \brief Publish a 2D real alias through the exchange variable.
  subroutine var2d_dp_publish_alias(self, component, source)
    class(var2d_dp), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var2d_dp), intent(in) :: source  !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var2d_dp_publish_alias

  !> \brief Return whether a 2D integer exchange variable has data connected.
  logical function var2d_i4_has_data(self)
    class(var2d_i4), intent(in) :: self

    var2d_i4_has_data = associated(self%data)
  end function var2d_i4_has_data

  !> \brief Return the data shape of a 2D integer exchange variable.
  function var2d_i4_data_shape(self) result(shape)
    class(var2d_i4), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8), size(self%data, 2, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var2d_i4_data_shape

  !> \brief Clear the data pointer of a 2D integer exchange variable.
  subroutine var2d_i4_clear_data(self)
    class(var2d_i4), intent(inout) :: self

    nullify(self%data)
  end subroutine var2d_i4_clear_data

  !> \brief Publish a local 2D integer field through the exchange variable.
  subroutine var2d_i4_publish_local(self, component, local)
    class(var2d_i4), intent(inout) :: self
    character(*), intent(in) :: component      !< publishing component name for diagnostics
    integer(i4), intent(inout), target :: local(:, :) !< local 2D integer field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var2d_i4_publish_local

  !> \brief Publish a 2D integer alias through the exchange variable.
  subroutine var2d_i4_publish_alias(self, component, source)
    class(var2d_i4), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var2d_i4), intent(in) :: source  !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var2d_i4_publish_alias

  !> \brief Return whether a 2D logical exchange variable has data connected.
  logical function var2d_lg_has_data(self)
    class(var2d_lg), intent(in) :: self

    var2d_lg_has_data = associated(self%data)
  end function var2d_lg_has_data

  !> \brief Return the data shape of a 2D logical exchange variable.
  function var2d_lg_data_shape(self) result(shape)
    class(var2d_lg), intent(in) :: self
    integer(i8), allocatable :: shape(:)

    if (associated(self%data)) then
      shape = [size(self%data, 1, kind=i8), size(self%data, 2, kind=i8)]
    else
      allocate(shape(0))
    end if
  end function var2d_lg_data_shape

  !> \brief Clear the data pointer of a 2D logical exchange variable.
  subroutine var2d_lg_clear_data(self)
    class(var2d_lg), intent(inout) :: self

    nullify(self%data)
  end subroutine var2d_lg_clear_data

  !> \brief Publish a local 2D logical field through the exchange variable.
  subroutine var2d_lg_publish_local(self, component, local)
    class(var2d_lg), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    logical, intent(inout), target :: local(:, :) !< local 2D logical field to publish

    call variable_validate_publish_target(self, component)
    self%data => local
    self%provided = .true.
  end subroutine var2d_lg_publish_local

  !> \brief Publish a 2D logical alias through the exchange variable.
  subroutine var2d_lg_publish_alias(self, component, source)
    class(var2d_lg), intent(inout) :: self
    character(*), intent(in) :: component !< publishing component name for diagnostics
    type(var2d_lg), intent(in) :: source  !< already-published source variable to alias

    call variable_validate_alias_source(source, component, self)
    call variable_validate_publish_target(self, component)
    self%data => source%data
    self%provided = .true.
  end subroutine var2d_lg_publish_alias

end module mo_exchange_type
