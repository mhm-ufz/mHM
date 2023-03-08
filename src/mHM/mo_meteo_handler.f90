!> \file    mo_meteo_handler.f90
!> \brief   \copybrief mo_meteo_handler
!> \details \copydetails mo_meteo_handler

!> \brief   Class for the meteo handler
!> \details River temperature routing on top of mRM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Mar 2023
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
module mo_meteo_handler

  USE mo_kind, ONLY : i4, dp
  USE mo_constants, ONLY : YearMonths
  use mo_common_types, only: Grid

  implicit none

  private

  !> \class   meteo_handler_type
  !> \brief   This is a handler for the meteorological forcings
  !> \details This class provides all procedures to handle meteorological forcings from file or interfaces.
  type, public :: meteo_handler_type
    ! config settings
    character(256) :: dir_nml_name = 'directories_mHM'                     !< namelist name in mhm namelist
    character(256) :: weight_nml_name = 'nightDayRatio'                    !< namelist name in mhm namelist
    ! -------------------------------------------------------------------
    ! GRID description
    ! -------------------------------------------------------------------
    type(Grid), dimension(:), allocatable, public :: level2                !< Reference of the metereological variables
    ! -------------------------------------------------------------------
    ! INPUT variables for configuration of mHM
    ! -------------------------------------------------------------------
    integer(i4), dimension(:), allocatable, public :: timeStep_model_inputs !< frequency for reading meteo input
    logical, public :: read_meteo_weights                                   !< read weights for tavg and pet
    character(256), public :: inputFormat_meteo_forcings                    !< format of meteo input data (nc)
    ! ------------------------------------------------------------------
    ! DIRECTORIES
    ! ------------------------------------------------------------------
    ! has the dimension of nDomains
    character(256), dimension(:), allocatable, public :: dirPrecipitation  !< Directory where precipitation files are located
    character(256), dimension(:), allocatable, public :: dirTemperature    !< Directory where temperature files are located
    character(256), dimension(:), allocatable, public :: dirMinTemperature !< Directory where minimum temp. files are located
    character(256), dimension(:), allocatable, public :: dirMaxTemperature !< Directory where maximum temp. files are located
    character(256), dimension(:), allocatable, public :: dirNetRadiation   !< Directory where abs. vap. pressure files are located
    character(256), dimension(:), allocatable, public :: dirabsVapPressure !< Directory where abs. vap. pressure files are located
    character(256), dimension(:), allocatable, public :: dirwindspeed      !< Directory where windspeed files are located
    character(256), dimension(:), allocatable, public :: dirReferenceET    !< Directory where reference-ET files are located
    ! riv-temp releated
    character(256), dimension(:), allocatable, public :: dirRadiation      !< Directory where short/long-wave rad. files are located
    ! -------------------------------------------------------------------
    ! L1 DOMAIN description
    ! -------------------------------------------------------------------
    ! Forcings
    ! dim1 = number grid cells L1
    ! dim2 = number of meteorological time steps
    ! dim2 = month of year
    real(dp), public, dimension(:, :, :), allocatable :: L1_temp_weights   !< hourly temperature weights for daily values
    real(dp), public, dimension(:, :, :), allocatable :: L1_pet_weights    !< hourly pet weights for daily values
    real(dp), public, dimension(:, :, :), allocatable :: L1_pre_weights    !< hourly pre weights for daily values
    real(dp), public, dimension(:, :), allocatable :: L1_pre               !< [mm]    Precipitation
    real(dp), public, dimension(:, :), allocatable :: L1_temp              !< [degC]  Air temperature
    real(dp), public, dimension(:, :), allocatable :: L1_pet               !< [mm TS-1] Potential evapotranspiration
    real(dp), public, dimension(:, :), allocatable :: L1_tmin              !< [degC]  minimum daily air temperature
    real(dp), public, dimension(:, :), allocatable :: L1_tmax              !< [degC]  maximum daily air temperature
    real(dp), public, dimension(:, :), allocatable :: L1_netrad            !< [W m2]  net radiation
    real(dp), public, dimension(:, :), allocatable :: L1_absvappress       !< [Pa]    absolute vapour pressure
    real(dp), public, dimension(:, :), allocatable :: L1_windspeed         !< [m s-1] windspeed
    ! riv-temp related
    real(dp), public, dimension(:, :), allocatable :: L1_ssrd              !< [W m2]  short wave radiation
    real(dp), public, dimension(:, :), allocatable :: L1_strd              !< [W m2]  long wave radiation
    real(dp), public, dimension(:, :), allocatable :: L1_tann              !< [degC]  annual mean air temperature
    ! -------------------------------------------------------------------
    ! Monthly day/night variation of Meteorological variables
    ! for temporal disaggregation
    ! -------------------------------------------------------------------
    ! dim1 = number of months in a year
    real(dp), public, dimension(int(YearMonths, i4)) :: fday_prec          !< [-] Day ratio precipitation < 1
    real(dp), public, dimension(int(YearMonths, i4)) :: fnight_prec        !< [-] Night ratio precipitation < 1
    real(dp), public, dimension(int(YearMonths, i4)) :: fday_pet           !< [-] Day ratio PET  < 1
    real(dp), public, dimension(int(YearMonths, i4)) :: fnight_pet         !< [-] Night ratio PET  < 1
    real(dp), public, dimension(int(YearMonths, i4)) :: fday_temp          !< [-] Day factor mean temp
    real(dp), public, dimension(int(YearMonths, i4)) :: fnight_temp        !< [-] Night factor mean temp
    real(dp), public, dimension(int(YearMonths, i4)) :: fday_ssrd          !< [-] Day factor short-wave rad.
    real(dp), public, dimension(int(YearMonths, i4)) :: fnight_ssrd        !< [-] Night factor short-wave rad.
    real(dp), public, dimension(int(YearMonths, i4)) :: fday_strd          !< [-] Day factor long-wave rad.
    real(dp), public, dimension(int(YearMonths, i4)) :: fnight_strd        !< [-] Night factor long-wave rad.

  contains
    !> \copydoc mo_meteo_handler::config
    procedure :: config !< \see mo_meteo_handler::config
    !> \copydoc mo_meteo_handler::clean_up
    procedure :: clean_up !< \see mo_meteo_handler::clean_up
  end type meteo_handler_type

contains

  !> \brief clean up
  subroutine clean_up(self)
    implicit none

    class(meteo_handler_type), intent(inout) :: self

    if ( allocated(self%timeStep_model_inputs) ) deallocate(self%timeStep_model_inputs)
    if ( allocated(self%dirPrecipitation) ) deallocate(self%dirPrecipitation)
    if ( allocated(self%dirTemperature) ) deallocate(self%dirTemperature)
    if ( allocated(self%dirMinTemperature) ) deallocate(self%dirMinTemperature)
    if ( allocated(self%dirMaxTemperature) ) deallocate(self%dirMaxTemperature)
    if ( allocated(self%dirNetRadiation) ) deallocate(self%dirNetRadiation)
    if ( allocated(self%dirabsVapPressure) ) deallocate(self%dirabsVapPressure)
    if ( allocated(self%dirwindspeed) ) deallocate(self%dirwindspeed)
    if ( allocated(self%dirReferenceET) ) deallocate(self%dirReferenceET)
    if ( allocated(self%dirRadiation) ) deallocate(self%dirRadiation)
    if ( allocated(self%level2) ) deallocate(self%level2)
    if ( allocated(self%L1_temp_weights) ) deallocate(self%L1_temp_weights)
    if ( allocated(self%L1_pet_weights) ) deallocate(self%L1_pet_weights)
    if ( allocated(self%L1_pre_weights) ) deallocate(self%L1_pre_weights)
    if ( allocated(self%L1_pre) ) deallocate(self%L1_pre)
    if ( allocated(self%L1_temp) ) deallocate(self%L1_temp)
    if ( allocated(self%L1_pet) ) deallocate(self%L1_pet)
    if ( allocated(self%L1_tmin) ) deallocate(self%L1_tmin)
    if ( allocated(self%L1_tmax) ) deallocate(self%L1_tmax)
    if ( allocated(self%L1_netrad) ) deallocate(self%L1_netrad)
    if ( allocated(self%L1_absvappress) ) deallocate(self%L1_absvappress)
    if ( allocated(self%L1_windspeed) ) deallocate(self%L1_windspeed)
    if ( allocated(self%L1_ssrd) ) deallocate(self%L1_ssrd)
    if ( allocated(self%L1_strd) ) deallocate(self%L1_strd)
    if ( allocated(self%L1_tann) ) deallocate(self%L1_tann)

  end subroutine clean_up

  !> \brief configure the \ref meteo_handler_type class from the mhm namelist
  subroutine config(self, file_namelist, unamelist, optimize, domainMeta)

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_types, only : domain_meta
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_check, only : check_dir
    USE mo_string_utils, ONLY : num2str

    implicit none

    class(meteo_handler_type), intent(inout) :: self
    character(*), intent(in) :: file_namelist !< mhm namelist file
    integer, intent(in) :: unamelist !< unit to open namelist file
    logical, intent(in) :: optimize !< Optimization flag
    type(domain_meta), intent(in) :: domainMeta !< domain general description

    integer(i4), dimension(maxNoDomains) :: time_step_model_inputs
    character(256), dimension(maxNoDomains) :: dir_Precipitation
    character(256), dimension(maxNoDomains) :: dir_Temperature
    character(256), dimension(maxNoDomains) :: dir_MinTemperature
    character(256), dimension(maxNoDomains) :: dir_MaxTemperature
    character(256), dimension(maxNoDomains) :: dir_NetRadiation
    character(256), dimension(maxNoDomains) :: dir_windspeed
    character(256), dimension(maxNoDomains) :: dir_absVapPressure
    character(256), dimension(maxNoDomains) :: dir_ReferenceET
    character(256), dimension(maxNoDomains) :: dir_Radiation ! riv-temp related
    character(256) :: inputFormat_meteo_forcings

    logical :: read_meteo_weights
    real(dp), dimension(int(YearMonths, i4)) :: fday_prec
    real(dp), dimension(int(YearMonths, i4)) :: fnight_prec
    real(dp), dimension(int(YearMonths, i4)) :: fday_pet
    real(dp), dimension(int(YearMonths, i4)) :: fnight_pet
    real(dp), dimension(int(YearMonths, i4)) :: fday_temp
    real(dp), dimension(int(YearMonths, i4)) :: fnight_temp
    real(dp), dimension(int(YearMonths, i4)) :: fday_ssrd
    real(dp), dimension(int(YearMonths, i4)) :: fnight_ssrd
    real(dp), dimension(int(YearMonths, i4)) :: fday_strd
    real(dp), dimension(int(YearMonths, i4)) :: fnight_strd

    ! namelist directories
    namelist /directories_mHM/ &
      inputFormat_meteo_forcings, &
      dir_Precipitation, &
      dir_Temperature, &
      dir_ReferenceET, &
      dir_MinTemperature, &
      dir_MaxTemperature, &
      dir_absVapPressure, &
      dir_windspeed, &
      dir_NetRadiation, &
      dir_Radiation, &
      time_step_model_inputs

    ! namelist for night-day ratio of precipitation, referenceET and temperature
    namelist /nightDayRatio/ &
      read_meteo_weights, &
      fnight_prec, &
      fnight_pet, &
      fnight_temp, &
      fnight_ssrd, &
      fnight_strd

    ! allocate variables
    allocate(self%dirPrecipitation(domainMeta%nDomains))
    allocate(self%dirTemperature(domainMeta%nDomains))
    allocate(self%dirwindspeed(domainMeta%nDomains))
    allocate(self%dirabsVapPressure(domainMeta%nDomains))
    allocate(self%dirReferenceET(domainMeta%nDomains))
    allocate(self%dirMinTemperature(domainMeta%nDomains))
    allocate(self%dirMaxTemperature(domainMeta%nDomains))
    allocate(self%dirNetRadiation(domainMeta%nDomains))
    allocate(self%dirRadiation(domainMeta%nDomains))
    ! allocate time periods
    allocate(timestep_model_inputs(domainMeta%nDomains))

    ! open the namelist file
    call open_nml(file_namelist, unamelist, quiet=.true.)

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call position_nml(self%dir_nml_name, unamelist)
    read(unamelist, nml = directories_mHM)

    self%inputFormat_meteo_forcings = inputFormat_meteo_forcings
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      self%timestep_model_inputs(iDomain) = time_step_model_inputs(domainID)
      self%dirPrecipitation(iDomain) = dir_Precipitation(domainID)
      self%dirTemperature(iDomain) = dir_Temperature(domainID)
      self%dirReferenceET(iDomain) = dir_ReferenceET(domainID)
      self%dirMinTemperature(iDomain) = dir_MinTemperature(domainID)
      self%dirMaxTemperature(iDomain) = dir_MaxTemperature(domainID)
      self%dirNetRadiation(iDomain) = dir_NetRadiation(domainID)
      self%dirwindspeed(iDomain) = dir_windspeed(domainID)
      self%dirabsVapPressure(iDomain) = dir_absVapPressure(domainID)
      ! riv-temp related
      self%dirRadiation(iDomain) = dir_Radiation(domainID)
    end do

    ! consistency check for timestep_model_inputs
    if (any(timestep_model_inputs .ne. 0) .and. .not. all(timestep_model_inputs .ne. 0)) then
      call error_message('***ERROR: timestep_model_inputs either have to be all zero or all non-zero')
    end if
    ! check for optimzation and timestep_model_inputs options
    if (optimize .and. (any(timestep_model_inputs .ne. 0))) then
      call error_message('***ERROR: optimize and chunk read is switched on! (set timestep_model_inputs to zero)')
    end if

    !===============================================================
    ! Read night-day ratios and pan evaporation
    !===============================================================
    ! default values for long/shortwave rad.
    fnight_ssrd = 0.0_dp
    fnight_strd = 0.45_dp

    ! namelist for night-day ratio of precipitation, referenceET and temperature
    call position_nml(self%weight_nml_name, unamelist)
    read(unamelist, nml = nightDayRatio)
    self%read_meteo_weights = read_meteo_weights
    self%fnight_prec = fnight_prec
    self%fnight_pet = fnight_pet
    self%fnight_temp = fnight_temp
    self%fnight_ssrd = fnight_ssrd
    self%fnight_strd = fnight_strd
    self%fday_prec = 1.0_dp - self%fnight_prec
    self%fday_pet = 1.0_dp - self%fnight_pet
    self%fday_temp = -1.0_dp * self%fnight_temp
    self%fday_ssrd = 1.0_dp - self%fnight_ssrd
    self%fday_strd = 1.0_dp - self%fnight_strd

    ! TODO-RIV-TEMP:
    ! - add short- and long-wave raidiation weights (nc files)

    ! closing the namelist file
    call close_nml(unamelist)

  end subroutine config

end module mo_meteo_handler
