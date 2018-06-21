!>       \file mhm_driver.f90

!>       \brief Distributed precipitation-runoff model mHM

!>       \details This is the main driver of mHM, which calls
!>       one instance of mHM for a multiple basins and a given period.
!>       \image html  mhm5-logo.png "Typical mHM cell"
!>       \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>       \authors Luis Samaniego & Rohini Kumar (UFZ)

!>       \date Jun 2018

!>       \version 5.9

!>       \copyright (c)2005-2018, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!>       All rights reserved.

!>       This code is a property of:

!>       ----------------------------------------------------------

!>       Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ
!>       Registered Office: Leipzig
!>       Registration Office: Amtsgericht Leipzig
!>       Trade Register: Nr. B 4703
!>       Chairman of the Supervisory Board: MinDirig Wilfried Kraus
!>       Scientific Director: Prof. Dr. Georg Teutsch
!>       Administrative Director: Dr. Heike Grassmann

!>       ----------------------------------------------------------

!>       NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!>       EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!>       OF THIS SOFTWARE. If software is modified to produce
!>       derivative works, such modified software should be
!>       clearly marked, so as not to confuse it with the version
!>       available from UFZ.  This code can be used for research
!>       purposes ONLY provided that the following sources are
!>       acknowledged:

!>       Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!>       parameter regionalization of a grid-based hydrologic
!>       model at the mesoscale.  Water Resour. Res., 46,
!>       W05523, doi:10.1029/2008WR007327.

!>       Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!>       of distributed hydrologic model parameterization on water
!>       fluxes at multiple scales and locations, Water Resour. Res.,
!>       49, doi:10.1029/2012WR012195.

!>       For commercial applications you have to consult the
!>       authorities of the UFZ.

! Modifications:
! Stephan Thober                Nov 2013 - added read in of latitude longitude fields
! Matthias Zink                 Mar 2013 - edited screen output for gauges added inflow gauges
! Matthias Cuntz & Juliane Mai  Mar 2014 - Likelihood Kavetski uses 2 more parameters for the
!                                          error model global_parameters -> local_parameters
! Rohini Kumar                  Apr 2014 - implementation of the mHM run on a single cell configuration
!                                          also for the routing mode.
!                                        - run mHM at the input data level i.e. L0 grid
! Rohini Kumar                  May 2014 - model run on a regular lat-lon grid or on a regular X-Y coordinate system
! Stephan Thober                May 2014 - moved read meteo forcings to mo_mhm_eval
! Matthias Cuntz & Juliane Mai  Nov 2014 - LAI input from daily, monthly or yearly files
! Matthias Zink                 Mar 2015 - added optional soil mositure read in for calibration
! Luis Samaniego                Jul 2015 - added temporal directories for optimization
! Stephan Thober                Aug 2015 - removed routing related variables
! Stephan Thober                Oct 2015 - reorganized optimization (now compatible with mRM)
! Oldrich Rakovec, Rohini Kumar Oct 2015 - added reading of basin averaged TWS and objective function 15
!                                          for simultaneous calibration based on runoff and TWS
! Rohini Kumar                  Mar 2016 - options to handle different soil databases modified MPR to included
!                                          soil horizon specific properties/parameters
! Stephan Thober                Nov 2016 - implemented adaptive timestep for routing
! Rohini Kumar                  Dec 2016 - options to read (monthly mean) LAI fields
! Robert Schweppe               Jun 2018 - refactoring and reformatting

PROGRAM mhm_driver

  USE mo_file, ONLY : &
          version, version_date, file_main, &      ! main info
          file_namelist_mhm, unamelist_mhm, &      ! filename of namelist: main setup
          file_namelist_mhm_param, unamelist_mhm_param, &      ! filename of namelist: mhm model parameter
          file_defOutput                                               ! filename of namelist: output setup
  USE mo_finish, ONLY : finish                         ! Finish with style
  USE mo_global_variables, ONLY : &
          dirPrecipitation, &      ! directories
          dirTemperature, &      ! directories
          dirReferenceET, &      ! PET input path  if process 5 is 'PET is input' (case 0)
          dirMinTemperature, dirMaxTemperature, &      ! PET input paths if process 5 is Hargreaves-Samani  (case 1)
          dirNetRadiation, &      ! PET input paths if process 5 is Priestley-Taylor (case 2)
          dirabsVapPressure, dirwindspeed, &      ! PET input paths if process 5 is Penman-Monteith  (case 3)
          timestep_model_inputs !frequency of input read
  USE mo_common_mHM_mRM_variables, ONLY : &
          nTstepDay, &      ! number of timesteps per day (former: NAGG)
          simPer, &      ! simulation period
          optimize, opti_function, &                                   ! optimization on/off and optimization method
          mrm_coupling_mode
  USE mo_common_variables, ONLY : &
          write_restart, &      ! restart writing flags
          dirRestartOut, &
          dirConfigOut, &
          dirMorpho, dirLCover, &                                         ! directories
          dirOut, &      ! directories
          nbasins, &      ! number of basins
          processMatrix, &      ! basin information,  processMatrix
          global_parameters, global_parameters_name      ! mhm parameters (gamma) and their clear names
  USE mo_kind, ONLY : i4, dp                         ! number precision
  USE mo_message, ONLY : message, message_text          ! For print out
  USE mo_meteo_forcings, ONLY : prepare_meteo_forcings_data
  USE mo_mhm_eval, ONLY : mhm_eval
  USE mo_read_optional_data, ONLY : read_soil_moisture, &      ! optional soil moisture reader, basin_avg_TWS reader
          read_basin_avg_TWS, &
          read_neutrons, &
          read_evapotranspiration
  USE mo_common_read_config, ONLY : common_read_config                    ! Read main configuration files
  USE mo_common_mHM_mRM_read_config, ONLY : &
          common_mHM_mRM_read_config, check_optimization_settings ! Read main configuration files
  USE mo_mpr_read_config, ONLY : mpr_read_config                    ! Read main configuration files
  USE mo_mhm_read_config, ONLY : mhm_read_config                    ! Read main configuration files
  USE mo_read_wrapper, ONLY : read_data                      ! Read all input data
  USE mo_restart, ONLY : write_restart_files
  USE mo_startup, ONLY : mhm_initialize
  USE mo_string_utils, ONLY : num2str, separator             ! String magic
  USE mo_timer, ONLY : &
          timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  USE mo_write_ascii, ONLY : &
          write_configfile, &      ! Writing Configuration file
          write_optifile, &      ! Writing optimized parameter set and objective
          write_optinamelist     ! Writing optimized parameter set to a namelist
  USE mo_objective_function, ONLY : objective                 ! objective functions and likelihoods
  USE mo_optimization, ONLY : optimization
#ifdef MRM2MHM
  USE mo_mrm_objective_function_runoff, only : single_objective_runoff
  USE mo_mrm_init, ONLY : mrm_init
  USE mo_mrm_write, only : mrm_write

#endif
  !$ USE omp_lib, ONLY : OMP_GET_NUM_THREADS           ! OpenMP routines

  IMPLICIT NONE

  ! local
  integer(i4), dimension(8) :: datetime         ! Date and time
  !$ integer(i4)                        :: n_threads        ! OpenMP number of parallel threads
  integer(i4) :: iBasin               ! Counters
  integer(i4) :: iTimer           ! Current timer number
  integer(i4) :: nTimeSteps
  real(dp) :: funcbest         ! best objective function achivied during optimization
  logical, dimension(:), allocatable :: maskpara ! true  = parameter will be optimized, = parameter(i,4) = 1
  !                                              ! false = parameter will not be optimized = parameter(i,4) = 0
  procedure(mhm_eval), pointer :: eval
  procedure(objective), pointer :: obj_func

  ! --------------------------------------------------------------------------
  ! START
  ! --------------------------------------------------------------------------
  call message(separator)
  call message('              mHM-UFZ')
  call message()
  call message('    MULTISCALE HYDROLOGIC MODEL')
  call message('           Version ', trim(version))
  call message('           ', trim(version_date))
  call message()
  call message('Originally by L. Samaniego & R. Kumar')

  call message(separator)

  call message()
  call date_and_time(values = datetime)
  message_text = trim(num2str(datetime(3), '(I2.2)')) // "." // trim(num2str(datetime(2), '(I2.2)')) &
          // "." // trim(num2str(datetime(1), '(I4.4)')) // " " // trim(num2str(datetime(5), '(I2.2)')) &
          // ":" // trim(num2str(datetime(6), '(I2.2)')) // ":" // trim(num2str(datetime(7), '(I2.2)'))
  call message('Start at ', trim(message_text), '.')
  call message('Using main file ', trim(file_main), ' and namelists: ')
  call message('     ', trim(file_namelist_mhm))
  call message('     ', trim(file_namelist_mhm_param))
  call message('     ', trim(file_defOutput))
  call message()

  !$OMP PARALLEL
  !$ n_threads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ call message('Run with OpenMP with ', trim(num2str(n_threads)), ' threads.')

  call message()
  call message('Read namelist file: ', trim(file_namelist_mhm))
  call message('Read namelist file: ', trim(file_namelist_mhm_param))
  call message('Read namelist file: ', trim(file_defOutput))
  call common_read_config(file_namelist_mhm, unamelist_mhm)
  call mpr_read_config(file_namelist_mhm, unamelist_mhm, file_namelist_mhm_param, unamelist_mhm_param)
  call common_mHM_mRM_read_config(file_namelist_mhm, unamelist_mhm)
  call mhm_read_config(file_namelist_mhm, unamelist_mhm)
  call check_optimization_settings()

  call message()
  call message('# of basins:         ', trim(num2str(nbasins)))
  call message()
  call message('  Input data directories:')
  do iBasin = 1, nbasins
    call message('  --------------')
    call message('      BASIN                   ', num2str(iBasin, '(I3)'))
    call message('  --------------')
    call message('    Morphological directory:    ', trim(dirMorpho(iBasin)))
    call message('    Land cover directory:       ', trim(dirLCover(iBasin)))
    call message('    Precipitation directory:    ', trim(dirPrecipitation(iBasin)))
    call message('    Temperature directory:      ', trim(dirTemperature(iBasin)))
    select case (processMatrix(5, 1))
    case(-1 : 0) ! PET is input
      call message('    PET directory:              ', trim(dirReferenceET(iBasin)))
    case(1) ! Hargreaves-Samani
      call message('    Min. temperature directory: ', trim(dirMinTemperature(iBasin)))
      call message('    Max. temperature directory: ', trim(dirMaxTemperature(iBasin)))
    case(2) ! Priestely-Taylor
      call message('    Net radiation directory:    ', trim(dirNetRadiation(iBasin)))
    case(3) ! Penman-Monteith
      call message('    Net radiation directory:    ', trim(dirNetRadiation(iBasin)))
      call message('    Abs. vap. press. directory: ', trim(dirabsVapPressure(iBasin)))
      call message('    Windspeed directory:        ', trim(dirwindspeed(iBasin)))
    end select
    call message('    Output directory:           ', trim(dirOut(iBasin)))

    call message('')
  end do

  ! Start timings
  call timers_init

  ! --------------------------------------------------------------------------
  ! READ AND INITIALIZE
  ! --------------------------------------------------------------------------
  itimer = 1
  call message()

  call message('  Read data ...')
  call timer_start(itimer)
  ! for DEM, slope, ... define nGvar local
  ! read_data has a basin loop inside
  call read_data(simPer)
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

  ! read data for every basin
  itimer = itimer + 1
  call message('  Initialize basins ...')
  call timer_start(itimer)
  call mhm_initialize()
  call timer_stop(itimer)
  call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

  itimer = itimer + 1
  call message('  Read forcing and optional data ...')
  call timer_start(itimer)

  do iBasin = 1, nbasins
    ! read meteorology now, if optimization is switched on
    ! meteorological forcings (reading, upscaling or downscaling)
    if (timestep_model_inputs(iBasin) .eq. 0_i4) then
      call prepare_meteo_forcings_data(iBasin, 1)
    end if

    ! read optional optional data if necessary
    if (optimize) then
      select case (opti_function)
      case(10 : 13, 28)
        ! read optional spatio-temporal soil mositure data
        call read_soil_moisture(iBasin)
      case(17)
        ! read optional spatio-temporal neutrons data
        call read_neutrons(iBasin)
      case(27, 29, 30)
        ! read optional spatio-temporal evapotranspiration data
        call read_evapotranspiration(iBasin)
      case(15)
        ! read optional basin average TWS data at once, therefore only read it
        ! the last iteration of the basin loop to ensure same time for all basins
        ! note: this is similar to how the runoff is read using mrm below
        if (iBasin == nbasins) then
          call read_basin_avg_TWS()
        end if
      end select
    end if

  end do
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')

#ifdef MRM2MHM
  ! --------------------------------------------------------------------------
  ! READ and INITIALISE mRM ROUTING
  ! --------------------------------------------------------------------------
  mrm_coupling_mode = 2_i4
  if (processMatrix(8, 1) .ne. 0_i4) call mrm_init(file_namelist_mhm, unamelist_mhm, &
          file_namelist_mhm_param, unamelist_mhm_param)
#else
  mrm_coupling_mode = -1_i4
#endif

  !this call may be moved to another position as it writes the master config out file for all basins
  call write_configfile()

  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  iTimer = iTimer + 1
  call message()
  if (optimize) then
    eval => mhm_eval

    select case(opti_function)
#ifdef MRM2MHM
     case(1 : 9, 14, 31)
      ! call optimization against only runoff (no other variables)
      obj_func => single_objective_runoff
      call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
#endif
     case(10 : 13, 15, 17, 27, 28, 29, 30)
      ! call optimization for other variables
      obj_func => objective
      call optimization(eval, obj_func, dirConfigOut, funcBest, maskpara)
    case default
      call message('***ERROR: mhm_driver: The given objective function number ', &
              trim(adjustl(num2str(opti_function))), ' in mhm.nml is not valid!')
      stop 1
    end select

    ! write a file with final objective function and the best parameter set
    call write_optifile(funcbest, global_parameters(:, 3), global_parameters_name(:))
    ! write a file with final best parameter set in a namlist format
    call write_optinamelist(processMatrix, global_parameters, maskpara, global_parameters_name(:))
    deallocate(maskpara)
  else

    ! --------------------------------------------------------------------------
    ! call mHM
    ! get runoff timeseries if possible (i.e. when processMatrix(8,1) > 0)
    ! get other model outputs  (i.e. gridded fields of model output)
    ! --------------------------------------------------------------------------
    call message('  Run mHM')
    call timer_start(iTimer)
    call mhm_eval(global_parameters(:, 3))
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F12.3)')), ' seconds.')

  end if

  ! --------------------------------------------------------------------------
  ! WRITE RESTART files
  ! --------------------------------------------------------------------------
  if (write_restart  .AND. (.NOT. optimize)) then
    itimer = itimer + 1
    call message()
    call message('  Write restart file')
    call timer_start(itimer)
    call write_restart_files(dirRestartOut)
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
  end if

#ifdef MRM2MHM
  ! --------------------------------------------------------------------------
  ! WRITE RUNOFF (INCLUDING RESTART FILES, has to be called after mHM restart
  ! files are written)
  ! --------------------------------------------------------------------------
  if (processMatrix(8, 1) .ne. 0) call mrm_write()
#endif


  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  itimer = itimer + 1
  ! call message()
  ! call message('  Write ouput data')
  ! call timer_start(itimer)
  ! ! call write_data()
  ! call timer_stop(itimer)
  ! call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

  nTimeSteps = maxval(simPer(1 : nBasins)%julEnd - simPer(1 : nBasins)%julStart + 1) * nTstepDay
  call date_and_time(values = datetime)
  call message()
  message_text = 'Done ' // trim(num2str(nTimeSteps, '(I10)')) // " time steps."
  call message(trim(message_text))
  message_text = trim(num2str(datetime(3), '(I2.2)')) // "." // trim(num2str(datetime(2), '(I2.2)')) &
          // "." // trim(num2str(datetime(1), '(I4.4)')) // " " // trim(num2str(datetime(5), '(I2.2)')) &
          // ":" // trim(num2str(datetime(6), '(I2.2)')) // ":" // trim(num2str(datetime(7), '(I2.2)'))
  call message('Finished at ', trim(message_text), '.')
  call message()
  call finish('mHM', 'Finished!')

END PROGRAM mhm_driver
