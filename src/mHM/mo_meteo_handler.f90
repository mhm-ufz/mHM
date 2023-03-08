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
  use mo_common_types, only: Grid, period

  implicit none

  private

  !> \class   meteo_handler_type
  !> \brief   This is a handler for the meteorological forcings
  !> \details This class provides all procedures to handle meteorological forcings from file or interfaces.
  type, public :: meteo_handler_type
    ! config settings
    character(256) :: dir_nml_name = 'directories_mHM'                     !< namelist name in mhm namelist
    character(256) :: weight_nml_name = 'nightDayRatio'                    !< namelist name in mhm namelist
    integer(i4) :: pet_case                                                !< process case for PET (processCase(5))
    integer(i4) :: riv_temp_case                                           !< process case for river temperature (processCase(11))
    type(period), public :: readPer                                        !< start and end dates of read period
    integer(i4), public :: nTstepForcingDay                                !< Number of forcing intervals per day
    ! -------------------------------------------------------------------
    ! level 2 description
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

    ! -------------------------------------------------------------------
    ! current mHM array indizes
    ! -------------------------------------------------------------------
    !> meteorological time step for process 5 (PET)
    integer(i4), dimension(6) :: iMeteo_p5
    !> process 5: start index of vectors
    !! index 1: pet
    !! index 2: tmin
    !! index 3: tmax
    !! index 4: netrad
    !! index 5: absolute vapour pressure
    !! index 6: windspeed
    integer(i4), dimension(6) :: s_p5
    !> process 5: end index of vectors
    !! index 1: pet
    !! index 2: tmin
    !! index 3: tmax
    !! index 4: netrad
    !! index 5: absolute vapour pressure
    !! index 6: windspeed
    integer(i4), dimension(6) :: e_p5
    !> start and end index of meteo variables
    integer(i4) :: s_meteo, e_meteo
    !> index of meteo time-step
    integer(i4) :: iMeteoTS

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
  subroutine config(self, file_namelist, unamelist, optimize, domainMeta, processMatrix)

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
    integer(i4), dimension(nProcesses, 3), intent(in) :: processMatrix !< Info about which process runs in which option

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

    ! # init of number of forcing timesteps, will be set when reading forcings
    self%nTStepForcingDay = nodata_i4

    ! store important process cases
    self%pet_case = processMatrix(5,1)
    self%riv_temp_case = processMatrix(11,1)

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

  !> \brief update the current time-step of the \ref meteo_handler_type class
  subroutine update_timestep(self, tt, iDomain, domainID, domainMeta, level1, nTstepDay, simPer, timestep)
    use mo_common_types, only : domain_meta
    use mo_common_variables, only : nProcesses

    implicit none

    class(meteo_handler_type), intent(inout) :: self
    integer(i4), intent(in) :: tt !< current time step
    integer(i4), intent(in) :: iDomain !< current domain
    integer(i4), intent(in) :: domainID !< current domain ID
    type(domain_meta), intent(in) :: domainMeta !< domain general description
    !> grid information at hydrologic level
    type(Grid), dimension(:), intent(in) :: level1
    !> Number of time intervals per day
    integer(i4), intent(in) :: nTstepDay
    !> warmPer + evalPer
    type(period), dimension(:), intent(in) :: simPer
    !> [h] simulation time step (= TS) in [h]
    integer(i4), intent(in) :: timeStep

    ! time increment is done right after call to mrm (and initially before looping)
    if (self%timeStep_model_inputs(iDomain) .eq. 0_i4) then
      ! whole meteorology is already read

      ! set start and end of meteo position
      self%s_meteo = level1(iDomain)%iStart
      self%e_meteo = level1(iDomain)%iEnd

      ! time step for meteorological variable (daily values)
      ! iMeteoTS = ceiling(real(tt, dp) / real(nTstepDay, dp))
      self%iMeteoTS = ceiling(real(tt, dp) / real(nint( 24._dp / real(self%nTstepForcingDay, dp)), dp))
    else
      ! read chunk of meteorological forcings data (reading, upscaling/downscaling)
      call self%prepare_data(iDomain, domainID, tt, domainMeta, level1, nTstepDay, simPer, timestep)
      ! set start and end of meteo position
      self%s_meteo = 1
      self%e_meteo = e1 - s1 + 1
      ! time step for meteorological variable (daily values)
      self%iMeteoTS = &
        ceiling(real(tt, dp) / real(nint( 24._dp / real(self%nTstepForcingDay, dp)), dp)) &
        - (self%readPer%julStart - simPer(iDomain)%julStart)
    end if

    ! preapare vector length specifications depending on the process case
    ! process 5 - PET
    ! customize iMeteoTS for process 5 - PET
    select case (self%pet_case)
      ! s_p5/e_p5:  [pet,        tmax,    tmin,  netrad, absVapP,windspeed]
      case(-1 : 0) ! PET is input
        self%s_p5 = [self%s_meteo, 1, 1, 1, 1, 1]
        self%e_p5 = [self%e_meteo, 1, 1, 1, 1, 1]
        self%iMeteo_p5 = [self%iMeteoTS, 1, 1, 1, 1, 1 ]
      case(1) ! Hargreaves-Samani
        self%s_p5 = [self%s_meteo, self%s_meteo, self%s_meteo, 1, 1, 1]
        self%e_p5 = [self%e_meteo, self%e_meteo, self%e_meteo, 1, 1, 1]
        self%iMeteo_p5 = [self%iMeteoTS, self%iMeteoTS, self%iMeteoTS, 1, 1, 1 ]
      case(2) ! Priestely-Taylor
        self%s_p5 = [self%s_meteo, 1, 1, self%s_meteo, 1, 1]
        self%e_p5 = [self%e_meteo, 1, 1, self%e_meteo, 1, 1]
        self%iMeteo_p5 = [self%iMeteoTS, 1, 1, self%iMeteoTS, 1, 1 ]
      case(3) ! Penman-Monteith
        self%s_p5 = [self%s_meteo, 1, 1, self%s_meteo, self%s_meteo, self%s_meteo]
        self%e_p5 = [self%e_meteo, 1, 1, self%e_meteo, self%e_meteo, self%e_meteo]
        self%iMeteo_p5 = [self%iMeteoTS, 1, 1, self%iMeteoTS, self%iMeteoTS, self%iMeteoTS ]
    end select

  end subroutine update_timestep

  !> \brief Prepare meteorological forcings data for a given variable
  !> \details Prepare meteorological forcings data for a given variable.
  !! Internally this subroutine calls another routine meteo_wrapper
  !! for different meterological variables
  !> \changelog
  !! - Matthias Zink,   Jun 2013
  !!   - addded NetCDf reader
  !! - Rohini Kumar,    Aug 2013
  !!   - name changed "inputFormat" to inputFormat_meteo_forcings
  !! - Matthias Zink,   Feb 2014
  !!   - added read in for different PET processes (process 5)
  !! - Stephan Thober,  Jun 2014
  !!   - add chunk_config for chunk read, copied L2 initialization to mo_startup
  !! - Stephan Thober,  Nov 2016
  !!   - moved processMatrix to common variables
  !! - Stephan Thober,  Jan 2017
  !!   - added subroutine for meteo_weights
  !! - Robert Schweppe  Jun 2018
  !!   - refactoring and reformatting
  !! - Sebastian Müller Mar 2023
  !!   - converted routine to meteo-handler method
  !> \authors Rohini Kumar
  !> \date Jan 2013
  subroutine prepare_data(self, iDomain, domainID, tt, &
    domainMeta, level1, nTstepDay, simPer, timestep)

    use mo_common_types, only : domain_meta
    use mo_common_variables, only : nProcesses
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, timer_stop
    use mo_meteo_forcings, only : meteo_forcings_wrapper, meteo_weights_wrapper, chunk_config

    implicit none

    class(meteo_handler_type), intent(inout) :: self
    integer(i4), intent(in) :: iDomain !< Domain number
    integer(i4), intent(in) :: domainID !< Domain ID
    integer(i4), intent(in) :: tt !< current timestep
    type(domain_meta), intent(in) :: domainMeta !< domain general description
    !> grid information at hydrologic level
    type(Grid), dimension(:), intent(in) :: level1
    !> Number of time intervals per day
    integer(i4), intent(in) :: nTstepDay
    !> warmPer + evalPer
    type(period), dimension(:), intent(in) :: simPer
    !> [h] simulation time step (= TS) in [h]
    integer(i4), intent(in) :: timeStep

    ! indicate whether data should be read
    logical :: read_flag

    ! configuration of chunk_read
    call chunk_config(iDomain, tt, nTstepDay, simPer, timestep, self%timeStep_model_inputs, read_flag, self%readPer)

    ! only read, if read_flag is true
    if (read_flag) then

      ! read weights for hourly disaggregation of temperature
      if (tt .eq. 1) then
        ! TODO-RIV-TEMP: No NC files for weights for radiation at the moment
        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for tavg     ...')
        call meteo_weights_wrapper(iDomain, self%read_meteo_weights, self%dirTemperature(iDomain), &
          self%L1_temp_weights, level1=level1, level2=self%level2, ncvarName = 'tavg_weight')

        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for pet     ...')
        call meteo_weights_wrapper(iDomain, self%read_meteo_weights, self%dirReferenceET(iDomain), &
          self%L1_pet_weights, level1=level1, level2=self%level2, ncvarName = 'pet_weight')

        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read meteo weights for pre     ...')
        call meteo_weights_wrapper(iDomain, self%read_meteo_weights, self%dirPrecipitation(iDomain), &
          self%L1_pre_weights, level1=level1, level2=self%level2, ncvarName = 'pre_weight')
      end if

      ! free L1 variables if chunk read is activated
      if (self%timeStep_model_inputs(iDomain) .ne. 0) then
        if (allocated(self%L1_pre)) deallocate(self%L1_pre)
        if (allocated(self%L1_temp)) deallocate(self%L1_temp)
        if (allocated(self%L1_pet)) deallocate(self%L1_pet)
        if (allocated(self%L1_tmin)) deallocate(self%L1_tmin)
        if (allocated(self%L1_tmax)) deallocate(self%L1_tmax)
        if (allocated(self%L1_netrad)) deallocate(self%L1_netrad)
        if (allocated(self%L1_absvappress)) deallocate(self%L1_absvappress)
        if (allocated(self%L1_windspeed)) deallocate(self%L1_windspeed)
        ! TODO-RIV-TEMP: deallocate riv-temp related vars also
      end if

      !  Domain characteristics and read meteo header
      if (self%timeStep_model_inputs(iDomain) .eq. 0) then
        call message('  Reading meteorological forcings for Domain: ', trim(adjustl(num2str(domainID))), ' ...')
        call timer_start(1)
      end if

      ! precipitation
      if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read precipitation        ...')
      call meteo_forcings_wrapper(iDomain, self%dirPrecipitation(iDomain), self%inputFormat_meteo_forcings, &
        dataOut1=self%L1_pre, &
        readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
        lower = 0.0_dp, upper = 1000._dp, ncvarName = 'pre')

      ! temperature
      if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read temperature          ...')
      call meteo_forcings_wrapper(iDomain, self%dirTemperature(iDomain), self%inputFormat_meteo_forcings, &
        dataOut1=self%L1_temp, &
        readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
        lower = -100._dp, upper = 100._dp, ncvarName = 'tavg')

      ! read input for PET (process 5) depending on specified option
      ! 0 - input, 1 - Hargreaves-Samani, 2 - Priestley-Taylor, 3 - Penman-Monteith
      select case (self%pet_case)
        case(-1 : 0) ! pet is input
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read pet                  ...')
          call meteo_forcings_wrapper(iDomain, self%dirReferenceET(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_pet, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = 0.0_dp, upper = 1000._dp, ncvarName = 'pet')
          ! allocate PET and dummies for mhm_call
          if ((iDomain.eq.domainMeta%nDomains) .OR. (self%timeStep_model_inputs(iDomain) .NE. 0)) then
            allocate(self%L1_tmin(1, 1))
            allocate(self%L1_tmax(1, 1))
            allocate(self%L1_netrad(1, 1))
            allocate(self%L1_absvappress(1, 1))
            allocate(self%L1_windspeed(1, 1))
          end if

        case(1) ! Hargreaves-Samani formulation (input: minimum and maximum Temperature)
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read min. temperature     ...')
          call meteo_forcings_wrapper(iDomain, self%dirMinTemperature(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_tmin, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = -100.0_dp, upper = 100._dp, ncvarName = 'tmin')
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read max. temperature     ...')
          call meteo_forcings_wrapper(iDomain, self%dirMaxTemperature(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_tmax, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = -100.0_dp, upper = 100._dp, ncvarName = 'tmax')
          ! allocate PET and dummies for mhm_call
          if ((iDomain .eq. domainMeta%nDomains) .OR. (self%timeStep_model_inputs(iDomain) .NE. 0)) then
            allocate(self%L1_pet    (size(self%L1_tmax, dim = 1), size(self%L1_tmax, dim = 2)))
            allocate(self%L1_netrad(1, 1))
            allocate(self%L1_absvappress(1, 1))
            allocate(self%L1_windspeed(1, 1))
          end if

        case(2) ! Priestley-Taylor formulation (input: net radiation)
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read net radiation        ...')
          call meteo_forcings_wrapper(iDomain, self%dirNetRadiation(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_netrad, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = -500.0_dp, upper = 1500._dp, ncvarName = 'net_rad')
          ! allocate PET and dummies for mhm_call
          if ((iDomain .eq. domainMeta%nDomains) .OR. (self%timeStep_model_inputs(iDomain) .NE. 0)) then
            allocate(self%L1_pet    (size(self%L1_netrad, dim = 1), size(self%L1_netrad, dim = 2)))
            allocate(self%L1_tmin(1, 1))
            allocate(self%L1_tmax(1, 1))
            allocate(self%L1_absvappress(1, 1))
            allocate(self%L1_windspeed(1, 1))
          end if

        case(3) ! Penman-Monteith formulation (input: net radiationm absulute vapour pressure, windspeed)
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read net radiation        ...')
          call meteo_forcings_wrapper(iDomain, self%dirNetRadiation(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_netrad, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = -500.0_dp, upper = 1500._dp, ncvarName = 'net_rad')
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read absolute vapour pressure  ...')
          call meteo_forcings_wrapper(iDomain, self%dirabsVapPressure(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_absvappress, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = 0.0_dp, upper = 15000.0_dp, ncvarName = 'eabs')
          if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read windspeed            ...')
          call meteo_forcings_wrapper(iDomain, self%dirwindspeed(iDomain), self%inputFormat_meteo_forcings, &
            dataOut1=self%L1_windspeed, &
            readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
            lower = 0.0_dp, upper = 250.0_dp, ncvarName = 'windspeed')
          ! allocate PET and dummies for mhm_call
          if ((iDomain.eq.domainMeta%nDomains) .OR. (self%timeStep_model_inputs(iDomain) .NE. 0)) then
            allocate(self%L1_pet    (size(self%L1_absvappress, dim = 1), size(self%L1_absvappress, dim = 2)))
            allocate(self%L1_tmin(1, 1))
            allocate(self%L1_tmax(1, 1))
          end if
      end select

      ! long/short-wave radiation and annual mean temperature for river-temperature routing
      if ( self%riv_temp_case .ne. 0 ) then
        ! free L1 variables if chunk read is activated
        if (self%timeStep_model_inputs(iDomain) .ne. 0) then
          if (allocated(self%L1_ssrd)) deallocate(self%L1_ssrd)
          if (allocated(self%L1_strd)) deallocate(self%L1_strd)
          if (allocated(self%L1_tann)) deallocate(self%L1_tann)
        end if
        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read short-wave radiation ...')
        call meteo_forcings_wrapper( &
          iDomain, self%dirRadiation(iDomain), self%inputFormat_meteo_forcings, &
          dataOut1=self%L1_ssrd, &
          readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
          lower = 0.0_dp, upper = 1500._dp, ncvarName = 'ssrd')
        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read long-wave radiation ...')
        call meteo_forcings_wrapper( &
          iDomain, self%dirRadiation(iDomain), self%inputFormat_meteo_forcings, &
          dataOut1=self%L1_strd, &
          readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
          lower = 0.0_dp, upper = 1500._dp, ncvarName = 'strd')
        if (self%timeStep_model_inputs(iDomain) .eq. 0) call message('    read annual mean temperature ...')
        call meteo_forcings_wrapper( &
          iDomain, self%dirTemperature(iDomain), self%inputFormat_meteo_forcings, &
          dataOut1=self%L1_tann, &
          readPer=self%readPer, nTstepForcingDay=self%nTstepForcingDay, level1=level1, level2=self%level2, &
          lower = -100.0_dp, upper = 100._dp, ncvarName = 'tann')
      end if

      if (self%timeStep_model_inputs(iDomain) .eq. 0) then
        call timer_stop(1)
        call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')
      end if
    end if

  end subroutine prepare_data

end module mo_meteo_handler
