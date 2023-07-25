!> \dir coupling
!> \brief \copybrief f_coupling
!> \details \copydetails f_coupling

!> \defgroup   f_coupling coupling - Fortran modules
!> \brief      Core modules to couple mHM.
!> \details    These modules provide the core components to couple mHM.

!> \file    mo_coupling_type.f90
!> \brief   \copybrief mo_coupling_type
!> \details \copydetails mo_coupling_type

!> \brief   Types to specify the coupling configuration of mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Apr 2023
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_coupling
module mo_coupling_type

  use mo_kind, only: i4, dp

  implicit none

  private

  !> \class   couple_cfg_type
  !> \brief   This is a container to hold all coupling configurations for mHM.
  type, public :: couple_cfg_type
    character(256) :: nml_name = 'coupling' !< namelist name in mhm namelist
    logical :: read_nml = .true.            !< whether to read the namelist (or set via routine)
    integer(i4) :: case                     !< coupling case (0: no coupling, 1: single domain coupling)
    integer(i4) :: meteo_timestep           !< timestep for meteo-data from coupling
    logical :: meteo_time_ref_endpoint      !< expect meteo has time reference point at end of associated time interval
    logical :: meteo_expect_pre             !< expect meteo from coupling: [mm]      Precipitation
    logical :: meteo_expect_temp            !< expect meteo from coupling: [degC]    Air temperature
    logical :: meteo_expect_pet             !< expect meteo from coupling: [mm TS-1] Potential evapotranspiration
    logical :: meteo_expect_tmin            !< expect meteo from coupling: [degC]    minimum daily air temperature
    logical :: meteo_expect_tmax            !< expect meteo from coupling: [degC]    maximum daily air temperature
    logical :: meteo_expect_netrad          !< expect meteo from coupling: [W m2]    net radiation
    logical :: meteo_expect_absvappress     !< expect meteo from coupling: [Pa]      absolute vapour pressure
    logical :: meteo_expect_windspeed       !< expect meteo from coupling: [m s-1]   windspeed
    ! riv-temp related
    logical :: meteo_expect_ssrd            !< expect meteo from coupling: [W m2]    short wave radiation
    logical :: meteo_expect_strd            !< expect meteo from coupling: [W m2]    long wave radiation
    logical :: meteo_expect_tann            !< expect meteo from coupling: [degC]    annual mean air temperature

  contains
    !> \copydoc mo_coupling_type::read_config
    procedure :: read_config !< \see mo_coupling_type::read_config
    !> \copydoc mo_coupling_type::set_config
    procedure :: set_config !< \see mo_coupling_type::set_config
    !> \copydoc mo_coupling_type::any_meteo_expected
    procedure :: any_meteo_expected !< \see mo_coupling_type::any_meteo_expected
    !> \copydoc mo_coupling_type::clean_up
    procedure :: clean_up !< \see mo_coupling_type::clean_up
    !> \copydoc mo_coupling_type::check
    procedure :: check !< \see mo_coupling_type::check
    !> \copydoc mo_coupling_type::active
    procedure :: active !< \see mo_coupling_type::active
  end type couple_cfg_type

contains

  !> \brief clean up
  subroutine clean_up(self)
    implicit none

    class(couple_cfg_type), intent(inout) :: self

    ! restore defaults
    call self%set_config(read_nml=.true.)

  end subroutine clean_up

  !> \brief read configuration for the \ref couple_cfg_type class from the mhm namelist
  subroutine read_config(self, file_namelist, unamelist)
    use mo_nml, only : close_nml, open_nml, position_nml
    implicit none

    class(couple_cfg_type), intent(inout) :: self
    character(*), intent(in) :: file_namelist !< mhm namelist file
    integer, intent(in) :: unamelist !< unit to open namelist file

    integer(i4) :: case                 ! coupling case
    integer(i4) :: meteo_timestep       ! timestep for meteo-data from coupling
    logical :: meteo_time_ref_endpoint  ! expect meteo has time reference point at end of associated time interval
    logical :: meteo_expect_pre         ! expect meteo from coupling: [mm]      Precipitation
    logical :: meteo_expect_temp        ! expect meteo from coupling: [degC]    Air temperature
    logical :: meteo_expect_pet         ! expect meteo from coupling: [mm TS-1] Potential evapotranspiration
    logical :: meteo_expect_tmin        ! expect meteo from coupling: [degC]    minimum daily air temperature
    logical :: meteo_expect_tmax        ! expect meteo from coupling: [degC]    maximum daily air temperature
    logical :: meteo_expect_netrad      ! expect meteo from coupling: [W m2]    net radiation
    logical :: meteo_expect_absvappress ! expect meteo from coupling: [Pa]      absolute vapour pressure
    logical :: meteo_expect_windspeed   ! expect meteo from coupling: [m s-1]   windspeed
    logical :: meteo_expect_ssrd        ! expect meteo from coupling: [W m2]    short wave radiation
    logical :: meteo_expect_strd        ! expect meteo from coupling: [W m2]    long wave radiation
    logical :: meteo_expect_tann        ! expect meteo from coupling: [degC]    annual mean air temperature

    integer(i4) :: status ! nml status

    ! namelist coupling
    namelist /coupling/ &
      case, &
      meteo_timestep, &
      meteo_time_ref_endpoint, &
      meteo_expect_pre, &
      meteo_expect_temp, &
      meteo_expect_pet, &
      meteo_expect_tmin, &
      meteo_expect_tmax, &
      meteo_expect_netrad, &
      meteo_expect_absvappress, &
      meteo_expect_windspeed, &
      meteo_expect_ssrd, &
      meteo_expect_strd, &
      meteo_expect_tann

    ! set defaults of namelist should be read
    if (self%read_nml) then
      call self%set_config(read_nml=.true.)
    else
      return ! config already set via set_config
    end if
    ! get defaults for local variables
    case = self%case
    meteo_timestep = self%meteo_timestep
    meteo_time_ref_endpoint = self%meteo_time_ref_endpoint
    meteo_expect_pre = self%meteo_expect_pre
    meteo_expect_temp = self%meteo_expect_temp
    meteo_expect_pet = self%meteo_expect_pet
    meteo_expect_tmin = self%meteo_expect_tmin
    meteo_expect_tmax = self%meteo_expect_tmax
    meteo_expect_netrad = self%meteo_expect_netrad
    meteo_expect_absvappress = self%meteo_expect_absvappress
    meteo_expect_windspeed = self%meteo_expect_windspeed
    meteo_expect_ssrd = self%meteo_expect_ssrd
    meteo_expect_strd = self%meteo_expect_strd
    meteo_expect_tann = self%meteo_expect_tann

    ! open the namelist file
    call open_nml(file_namelist, unamelist, quiet=.true.)
    call position_nml(self%nml_name, unamelist, status=status)

    ! coupling namelist can be missing
    if (status == 0_i4) then
      ! read namelist if present
      read(unamelist, nml=coupling)
    end if

    ! closing the namelist file
    call close_nml(unamelist)

    self%case = case
    self%meteo_timestep = meteo_timestep
    self%meteo_time_ref_endpoint = meteo_time_ref_endpoint
    self%meteo_expect_pre = meteo_expect_pre
    self%meteo_expect_temp = meteo_expect_temp
    self%meteo_expect_pet = meteo_expect_pet
    self%meteo_expect_tmin = meteo_expect_tmin
    self%meteo_expect_tmax = meteo_expect_tmax
    self%meteo_expect_netrad = meteo_expect_netrad
    self%meteo_expect_absvappress = meteo_expect_absvappress
    self%meteo_expect_windspeed = meteo_expect_windspeed
    self%meteo_expect_ssrd = meteo_expect_ssrd
    self%meteo_expect_strd = meteo_expect_strd
    self%meteo_expect_tann = meteo_expect_tann

  end subroutine read_config

  !> \brief set configuration for the \ref couple_cfg_type class
  subroutine set_config( &
    self, &
    case, &
    meteo_timestep, &
    meteo_time_ref_endpoint, &
    meteo_expect_pre, &
    meteo_expect_temp, &
    meteo_expect_pet, &
    meteo_expect_tmin, &
    meteo_expect_tmax, &
    meteo_expect_netrad, &
    meteo_expect_absvappress, &
    meteo_expect_windspeed, &
    meteo_expect_ssrd, &
    meteo_expect_strd, &
    meteo_expect_tann, &
    read_nml &
  )
    implicit none

    class(couple_cfg_type), intent(inout) :: self
    integer(i4), intent(in), optional :: case                 !< coupling case
    integer(i4), intent(in), optional :: meteo_timestep       !< timestep for meteo-data from coupling
    logical, intent(in), optional :: meteo_time_ref_endpoint  !< expect meteo has time reference point at end of time interval
    logical, intent(in), optional :: meteo_expect_pre         !< expect meteo from coupling: [mm]      Precipitation
    logical, intent(in), optional :: meteo_expect_temp        !< expect meteo from coupling: [degC]    Air temperature
    logical, intent(in), optional :: meteo_expect_pet         !< expect meteo from coupling: [mm TS-1] Potential evapotranspiration
    logical, intent(in), optional :: meteo_expect_tmin        !< expect meteo from coupling: [degC]    minimum daily air temperature
    logical, intent(in), optional :: meteo_expect_tmax        !< expect meteo from coupling: [degC]    maximum daily air temperature
    logical, intent(in), optional :: meteo_expect_netrad      !< expect meteo from coupling: [W m2]    net radiation
    logical, intent(in), optional :: meteo_expect_absvappress !< expect meteo from coupling: [Pa]      absolute vapour pressure
    logical, intent(in), optional :: meteo_expect_windspeed   !< expect meteo from coupling: [m s-1]   windspeed
    logical, intent(in), optional :: meteo_expect_ssrd        !< expect meteo from coupling: [W m2]    short wave radiation
    logical, intent(in), optional :: meteo_expect_strd        !< expect meteo from coupling: [W m2]    long wave radiation
    logical, intent(in), optional :: meteo_expect_tann        !< expect meteo from coupling: [degC]    annual mean air temperature
    logical, intent(in), optional :: read_nml                 !< whether to read the namelist (or set via routine)

    ! defaults
    self%case = 0_i4 ! no coupling by default
    self%meteo_timestep = 0_i4 ! only valid if no meteo expected
    self%meteo_time_ref_endpoint = .false. ! meteo data usually given at begin of time interval (i.e. 00:00 for current day)
    self%meteo_expect_pre = .false.
    self%meteo_expect_temp = .false.
    self%meteo_expect_pet = .false.
    self%meteo_expect_tmin = .false.
    self%meteo_expect_tmax = .false.
    self%meteo_expect_netrad = .false.
    self%meteo_expect_absvappress = .false.
    self%meteo_expect_windspeed = .false.
    self%meteo_expect_ssrd = .false.
    self%meteo_expect_strd = .false.
    self%meteo_expect_tann = .false.
    ! don't read namelist if configured via this routine (by default)
    ! if set_config is only used to set defaults, set this to true
    self%read_nml = .false.

    if (present(case)) self%case = case
    if (present(meteo_timestep)) self%meteo_timestep = meteo_timestep
    if (present(meteo_time_ref_endpoint)) self%meteo_time_ref_endpoint = meteo_time_ref_endpoint
    if (present(meteo_expect_pre)) self%meteo_expect_pre = meteo_expect_pre
    if (present(meteo_expect_temp)) self%meteo_expect_temp = meteo_expect_temp
    if (present(meteo_expect_pet)) self%meteo_expect_pet = meteo_expect_pet
    if (present(meteo_expect_tmin)) self%meteo_expect_tmin = meteo_expect_tmin
    if (present(meteo_expect_tmax)) self%meteo_expect_tmax = meteo_expect_tmax
    if (present(meteo_expect_netrad)) self%meteo_expect_netrad = meteo_expect_netrad
    if (present(meteo_expect_absvappress)) self%meteo_expect_absvappress = meteo_expect_absvappress
    if (present(meteo_expect_windspeed)) self%meteo_expect_windspeed = meteo_expect_windspeed
    if (present(meteo_expect_ssrd)) self%meteo_expect_ssrd = meteo_expect_ssrd
    if (present(meteo_expect_strd)) self%meteo_expect_strd = meteo_expect_strd
    if (present(meteo_expect_tann)) self%meteo_expect_tann = meteo_expect_tann
    if (present(read_nml)) self%read_nml = read_nml

  end subroutine set_config

  !> \brief whether any meteo data is expected
  !> \return True if any meteo data expected, else False
  logical function any_meteo_expected(self)
    implicit none

    class(couple_cfg_type), intent(in) :: self

    any_meteo_expected = &
      self%meteo_expect_pre .or. &
      self%meteo_expect_temp .or. &
      self%meteo_expect_pet .or. &
      self%meteo_expect_tmin .or. &
      self%meteo_expect_tmax .or. &
      self%meteo_expect_netrad .or. &
      self%meteo_expect_absvappress .or. &
      self%meteo_expect_windspeed .or. &
      self%meteo_expect_ssrd .or. &
      self%meteo_expect_strd .or. &
      self%meteo_expect_tann

  end function any_meteo_expected

  !> \brief check configuration
  subroutine check(self, domainMeta, optimize)
    use mo_message, only : error_message
    use mo_string_utils, only : num2str
    use mo_common_types, only : domain_meta
    implicit none

    class(couple_cfg_type), intent(inout) :: self
    type(domain_meta), intent(in) :: domainMeta !< domain general description
    logical, intent(in) :: optimize             !< Optimization (.true. ) or Evaluation run (.false.)

    if (.not. any(self%case == [0, 1])) &
      call error_message("Coupling: case needs to be 0 or 1. Got: ", num2str(self%case))

    ! coupling case 1 for single domain coupling
    if (self%case == 1 .and. domainMeta%nDomains > 1) &
      call error_message("Coupling: Only one domain allowed when coupling.")

    ! coupling only without internal optimization
    if (self%case /= 0 .and. optimize) &
      call error_message("Coupling: no internal optimization allowed when coupling.")

    ! if meteo data is expected, an associated time-step (1 or 24) needs to be given
    if (self%case /= 0 .and. self%any_meteo_expected() .and. .not. any(self%meteo_timestep == [1, 24])) &
      call error_message("Coupling: meteo data expected but no valid time-step (1 or 24) given for it.")

  end subroutine check

  !> \brief whether coupling is actived
  !> \return True if any case > 0, else False
  logical function active(self)
    implicit none
    class(couple_cfg_type), intent(in) :: self
    active = self%case > 0_i4
  end function active

end module mo_coupling_type
