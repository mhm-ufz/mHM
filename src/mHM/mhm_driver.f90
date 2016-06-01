!> \file mhm_driver.f90
! --------------------------------------------------------------------------
!> \authors   Luis Samaniego & Rohini Kumar (UFZ)
!  CONTACT    luis.samaniego@ufz.de / rohini.kumar@ufz.de
!> \version   5.4
!> \date      Dec 2015

!  PURPOSE
!>            \brief Distributed precipitation-runoff model mHM

!>            \details This is the main driver of mHM, which calls
!>             one instance of mHM for a multiple basins and a given period.

!>              \image html  mhm5-logo.png "Typical mHM cell"
!>              \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>  \copyright (c)2005-2016, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!>             All rights reserved.
!>
!>             This code is a property of:
!>
!>             ----------------------------------------------------------
!>
!>             Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ\n
!>             Registered Office: Leipzig\n
!>             Registration Office: Amtsgericht Leipzig\n
!>             Trade Register: Nr. B 4703\n
!>             Chairman of the Supervisory Board: MinDirig Wilfried Kraus\n
!>             Scientific Director: Prof. Dr. Georg Teutsch\n
!>             Administrative Director: Dr. Heike Grassmann\n
!>
!>             ----------------------------------------------------------
!>
!>             NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!>             EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!>             OF THIS SOFTWARE. If software is modified to produce
!>             derivative works, such modified software should be
!>             clearly marked, so as not to confuse it with the version
!>             available from UFZ.  This code can be used for research
!>             purposes ONLY provided that the following sources are
!>             acknowledged:
!>
!>                Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!>                parameter regionalization of a grid-based hydrologic
!>                model at the mesoscale.  Water Resour. Res., 46,
!>                W05523, doi:10.1029/2008WR007327.
!>
!>                Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!>                of distributed hydrologic model parameterization on water
!>                fluxes at multiple scales and locations, Water Resour. Res.,
!>                49, doi:10.1029/2012WR012195.
!>
!>             For commercial applications you have to consult the
!>             authorities of the UFZ.


! REDISTRIBUTION
!             Redistribution and use in source and binary forms,
!             with or without modification, are permitted provided that
!             the following conditions are met: * Redistributions of
!             source code must retain the above copyright notice, this
!             list of conditions and the following disclaimer.  *
!             Redistributions in binary form must reproduce the above
!             copyright notice, this list of conditions and the
!             following disclaimer in the documentation and/or other
!             materials provided with the distribution.  * Neither the
!             name of Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ,
!             nor the names of its contributors may be used to endorse
!             or promote products derived from this software without
!             specific prior written permission.

! DISCLAIM
!             THIS SOFTWARE IS PROVIDED BY HELMHOLTZ-ZENTRUM FUER
!             UMWELTFORSCHUNG GMBH - UFZ AND CONTRIBUTORS "AS IS" AND
!             ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!             LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
!             FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
!             EVENT SHALL THE HELMHOLTZ-ZENTRUM FUER UMWELTFORSCHUNG
!             GMBH - UFZ AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!             INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!             CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!             PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!             DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!             CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!             CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
!             OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!             SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
!             DAMAGE.

!             This program is distributed in the hope that it will be
!             useful,but WITHOUT ANY WARRANTY; without even the implied
!             warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!             PURPOSE. UFZ and the DEVELOPERS of this code do not take
!             any liabilities on the aplication of this code.

! --------------------------------------------------------------------------
!  HISTORY
!         Written          Sa 24.11.2005    v1.0 main structure
!         Modified         Sa 02.02.2007    v2.0 distributed hydrological
!                                                processes
!                          Sa 06.05.2007    v3.0 Muskingun routing scheme /
!                                                file structure
!                       Sa Ku 14.05.2007    v3.1 QMeasures, parameter limits
!                          Sa 10.09.2007    v3.1 distributed PET/TEMP/PRE
!                                                day/night fraction
!                          Ku 19.03.2008    v3.2 fully distributed
!                       Sa Ku 10.01.2009    v3.2 fully distributed
!                                                with MPR scheme
!                          Ku 01.10.2010    v4.0 matrix into vector version
!                          Ku 15.10.2010    v4.0 paralleized with OPENMP
!                          Ku 07.11.2010    v4.1 inclusion of the routing level
!                          Ku 10.12.2010    v4.2 append trans. funct. param
!                                                with geological parameter set
!                                                (nFixPar + nGeoUnits)
!                          Ku 25.01.2011    v4.3 n-soil Horizons /
!                                                vertical root distribution
!                          Ku 19.04.2011    v4.3 regionalisation of muskingum
!                                                routing parameters
!                          Ku 30.05.2011    v4.4 changes in regionalisation
!                                                functions of unsaturated
!                                                soil layers
!                          Sa 03.07.2012    v4.4 removed IMSL dependencies
!                 Sa Ku Ma Cu 12.12.2012    v5.0 modularization
!                       Sa Cu 12.12.2012    v5.0 automatic documentation
!                          Ku 02.05.2013    v5.0 error/compatability checks
!                     Stephan Thober, Nov 2013 - added read in of latitude longitude fields
!                     Matthias Zink,  Mar 2013 - edited screen output for gauges
!                                                added inflow gauges
!       Matthias Cuntz & Juliane Mai, Mar 2014 - Likelihood Kavetski uses 2 more parameters for the error model
!                                                global_parameters -> local_parameters
!                       Rohini Kumar, Apr 2014 - implementation of the mHM run on a single cell
!                                                configuration also for the routing mode.
!                                              - run mHM at the input data level i.e. L0 grid
!                       Rohini Kumar, May 2014 - model run on a regular lat-lon grid or
!                                                on a regular X-Y coordinate system
!                      Stephan Thober May 2014 - moved read meteo forcings to mo_mhm_eval
!       Matthias Cuntz & Juliane Mai, Nov 2014 - LAI input from daily, monthly or yearly files
!                      Matthias Zink, Mar 2015 - added optional soil mositure read in for calibration
!                     Luis Samaniego, Jul 2015 - added temporal directories for optimization
!                     Stephan Thober, Aug 2015 - removed routing related variables
!                     Stephan Thober, Oct 2015 - reorganized optimization (now compatible with mRM)
!      Oldrich Rakovec, Rohini Kumar, Oct 2015 - added reading of basin averaged TWS and objective function 15
!                                                for simultaneous calibration based on runoff and TWS
!                       Rohini Kumar, Mar 2016 - options to handle different soil databases
!                                                modified MPR to included soil horizon specific properties/parameters
!
! --------------------------------------------------------------------------

PROGRAM mhm_driver

  USE mo_file,                ONLY :                         &
       version, version_date, file_main,                     &      ! main info
       file_namelist,                                        &      ! filename of namelist: main setup
       file_namelist_param,                                  &      ! filename of namelist: mhm model parameter
       file_defOutput                                               ! filename of namelist: output setup
  USE mo_finish,              ONLY : finish                         ! Finish with style
  USE mo_global_variables,    ONLY :                         &
       nbasins, timestep_model_inputs,                       &      ! number of basins, frequency of input read
       write_restart,                                        &      ! restart writing flags
       dirRestartOut,                                        &      ! directories
       dirMorpho, dirLCover,  dirPrecipitation,              &      ! directories
       dirTemperature, dirOut,                               &      ! directories
       dirReferenceET,                                       &      ! PET input path  if process 5 is 'PET is input' (case 0)
       dirMinTemperature, dirMaxTemperature,                 &      ! PET input paths if process 5 is Hargreaves-Samani  (case 1)
       dirNetRadiation,                                      &      ! PET input paths if process 5 is Priestley-Taylor (case 2)
       dirabsVapPressure, dirwindspeed,                      &      ! PET input paths if process 5 is Penman-Monteith  (case 3)
       dirgridded_LAI,                                       &      ! directories
       simPer,                                               &      ! simulation period
       NTSTEPDAY,                                            &      ! number of timesteps per day (former: NAGG)
       timeStep_LAI_input,                                   &      ! LAI option for reading gridded LAI field
       processMatrix,                                        &      ! basin information,  processMatrix
       dirConfigOut,                                         &
       basin,                                                & ! L0_mask for mrm_init call
       L0_elev,                                              & ! L0_elev for mrm_init call
       L0_LCover                                               ! L0_LCover for mrm_init call
  USE mo_common_variables,    ONLY : &
       optimize, opti_function,                              &      ! optimization on/off and optimization method
       global_parameters, global_parameters_name                    ! mhm parameters (gamma) and their clear names
  USE mo_kind,                ONLY : i4, dp                         ! number precision
  USE mo_message,             ONLY : message, message_text          ! For print out
  USE mo_meteo_forcings,      ONLY : prepare_meteo_forcings_data
  USE mo_mhm_eval,            ONLY : mhm_eval
  USE mo_prepare_gridded_LAI, ONLY : prepare_gridded_daily_LAI_data ! prepare daily LAI gridded fields
  USE mo_read_optional_data,  ONLY : read_soil_moisture,     &      ! optional soil moisture reader, basin_avg_TWS reader
                                     read_basin_avg_TWS,     &
                                     read_neutrons
  USE mo_read_config,         ONLY : read_config                    ! Read main configuration files
  USE mo_read_wrapper,        ONLY : read_data                      ! Read all input data
  USE mo_read_latlon,         ONLY : read_latlon
  USE mo_restart,             ONLY : write_restart_files
  USE mo_startup,             ONLY : initialise
  USE mo_string_utils,        ONLY : num2str, separator             ! String magic
  USE mo_timer,               ONLY :                         &
       timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  USE mo_write_ascii,         ONLY :                         &
       write_configfile,                                     &      ! Writing Configuration file
       write_optifile,                                       &      ! Writing optimized parameter set and objective
       write_optinamelist                                           ! Writing optimized parameter set to a namelist
  USE mo_objective_function,  ONLY : objective                 ! objective functions and likelihoods
  USE mo_optimization,        ONLY : optimization
#ifdef mrm2mhm
  USE mo_mrm_objective_function_runoff, only: single_objective_runoff
  USE mo_mrm_init,            ONLY : mrm_init
  USE mo_mrm_write,           only : mrm_write
#endif
  !$ USE omp_lib,             ONLY : OMP_GET_NUM_THREADS           ! OpenMP routines

  IMPLICIT NONE

  ! local
  integer(i4), dimension(8)             :: datetime         ! Date and time
  !$ integer(i4)                        :: n_threads        ! OpenMP number of parallel threads
  integer(i4)                           :: ii               ! Counters
  integer(i4)                           :: iTimer           ! Current timer number
  integer(i4)                           :: nTimeSteps
  real(dp)                              :: funcbest         ! best objective function achivied during optimization
  logical,  dimension(:),   allocatable :: maskpara         ! true  = parameter will be optimized     = parameter(i,4) = 1
  !                                                         ! false = parameter will not be optimized = parameter(i,4) = 0

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
  call date_and_time(values=datetime)
  message_text = trim(num2str(datetime(3),'(I2.2)'))//"."//trim(num2str(datetime(2),'(I2.2)')) &
       //"."//trim(num2str(datetime(1),'(I4.4)'))//" "//trim(num2str(datetime(5),'(I2.2)')) &
       //":"//trim(num2str(datetime(6),'(I2.2)'))//":"//trim(num2str(datetime(7),'(I2.2)'))
  call message('Start at ', trim(message_text), '.')
  call message('Using main file ', trim(file_main), ' and namelists: ')
  call message('     ',trim(file_namelist))
  call message('     ',trim(file_namelist_param))
  call message('     ',trim(file_defOutput))
  call message()

  !$OMP PARALLEL
  !$ n_threads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ call message('Run with OpenMP with ', trim(num2str(n_threads)), ' threads.')

  call message()
  call message('Read namelist file: ', trim(file_namelist))
  call message('Read namelist file: ', trim(file_namelist_param))
  call message('Read namelist file: ', trim(file_defOutput))
  call read_config()

  call message()
  call message('  # of basins:         ', trim(num2str(nbasins)))
  call message()
  call message('  Input data directories:')
  do ii = 1, nbasins
     call message( '  --------------' )
     call message( '      BASIN                   ', num2str(ii,'(I3)') )
     call message( '  --------------' )
     call message('    Morphological directory:    ',   trim(dirMorpho(ii) ))
     call message('    Land cover directory:       ',   trim(dirLCover(ii) ))
     call message('    Precipitation directory:    ',   trim(dirPrecipitation(ii)  ))
     call message('    Temperature directory:      ',   trim(dirTemperature(ii)  ))
     select case (processMatrix(5,1))
     case(0) ! PET is input
        call message('    PET directory:              ', trim(dirReferenceET(ii)  ))
     case(1) ! Hargreaves-Samani
        call message('    Min. temperature directory: ', trim(dirMinTemperature(ii)  ))
        call message('    Max. temperature directory: ', trim(dirMaxTemperature(ii)  ))
     case(2) ! Priestely-Taylor
        call message('    Net radiation directory:    ', trim(dirNetRadiation(ii) ))
     case(3) ! Penman-Monteith
        call message('    Net radiation directory:    ', trim(dirNetRadiation(ii) ))
        call message('    Abs. vap. press. directory: ', trim(dirabsVapPressure(ii)  ))
        call message('    Windspeed directory:        ', trim(dirwindspeed(ii)  ))
     end select
     call message('    Output directory:           ',   trim(dirOut(ii) ))
     if (timeStep_LAI_input < 0) then
        call message('    LAI directory:             ', trim(dirgridded_LAI(ii)) )
     end if

     call message('')
  end do

  ! Start timings
  call timers_init

  ! --------------------------------------------------------------------------
  ! READ AND (RE-)INITIALISE
  ! --------------------------------------------------------------------------
  itimer = 1
  call message()

  call message('  Read data ' )
  call timer_start(itimer)
  ! for DEM, slope, ... define nGvar local
  ! read_data has a basin loop inside
  call read_data()
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

  ! read data for every basin
  do ii=1, nbasins
     itimer = itimer + 1
     call message('  Initialise basin: ', trim(adjustl(num2str(ii))),' ...')
     call timer_start(itimer)
     call initialise(ii)
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

     ! read meteorology now, if optimization is switched on
     ! meteorological forcings (reading, upscaling or downscaling)
     if ( timestep_model_inputs(ii) .eq. 0_i4 ) then
        call prepare_meteo_forcings_data(ii, 1)
     end if

     ! read lat lon coordinates of each basin
     call message('  Reading lat-lon for basin: ', trim(adjustl(num2str(ii))),' ...')
     call timer_start(itimer)
     call read_latlon(ii)
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

     ! daily gridded LAI values
     if (timeStep_LAI_input < 0) then
        call message('  Reading LAI for basin: ', trim(adjustl(num2str(ii))),' ...')
        call timer_start(itimer)
        call prepare_gridded_daily_LAI_data(ii)
        call timer_stop(itimer)
        call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
     endif

     ! read optional optional data
     ! e.g. for optimization against soil mopisture, soil moisture is read
     if ((opti_function .GE. 10) .AND. (opti_function .LE. 13) .AND. optimize) then
        call read_soil_moisture(ii)
     endif

     ! read optional spatio-temporal neutrons data
    if ( (opti_function .EQ. 17) .AND. optimize ) then
        call read_neutrons(ii)
        call message('  neutrons data read')
    endif

  end do

  ! read optional basin average TWS data at once, therefore outside of the basin loop to ensure same time for all basins
  ! note this is similar to how the runoff is read using mrm below
     if ( (opti_function .EQ. 15) .AND. optimize ) then
        call read_basin_avg_TWS()
        call message('  basin_avg TWS data read')
     endif

  ! The following block is for testing of the restart <<<<<<<<<<<<<<<<<<<<<<<<<<
  ! if ( write_restart ) then
  !    itimer = itimer + 1
  !    call message()
  !    call message( '  Write restart config file')
  !    call timer_start(itimer)
  !    call write_restart_config( dirRestartOut )
  !    call timer_stop(itimer)
  !    call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  ! end if
  ! stop 'Test restart' ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef mrm2mhm
  ! --------------------------------------------------------------------------
  ! READ and INITIALISE mRM ROUTING
  ! --------------------------------------------------------------------------
  if (processMatrix(8, 1) .eq. 1) call mrm_init(basin%L0_mask, L0_elev, L0_LCover)
#endif

  !this call may be moved to another position as it writes the master config out file for all basins
  call write_configfile()

  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  iTimer = iTimer + 1
  call message()
  if ( optimize ) then

     select case(opti_function) 
#ifdef mrm2mhm
     case(1:9,14) 
        ! call optimization against only runoff (no other variables)
        call optimization(single_objective_runoff, dirConfigOut, funcBest, maskpara)
#endif
     case(10:13,15,17) 
        ! call optimization for other variables
        call optimization(objective, dirConfigOut, funcBest, maskpara)
     case default 
        call message('mhm_driver: 1: The kind (SO or MO) is not specified for the given objective function!') 
        stop 
     end select

     ! write a file with final objective function and the best parameter set
     call write_optifile(funcbest, global_parameters(:,3), global_parameters_name(:))
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
     call mhm_eval(global_parameters(:,3))
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F12.3)')), ' seconds.')

  end if

  ! --------------------------------------------------------------------------
  ! WRITE RESTART files
  ! --------------------------------------------------------------------------
  if ( write_restart  .AND. (.NOT. optimize)) then
     itimer = itimer + 1
     call message()
     call message( '  Write restart file')
     call timer_start(itimer)
     call write_restart_files( dirRestartOut )
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  end if

#ifdef mrm2mhm
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

  nTimeSteps = maxval( simPer(1:nBasins)%julEnd - simPer(1:nBasins)%julStart + 1 ) * NTSTEPDAY
  call date_and_time(values=datetime)
  call message()
  message_text = 'Done '//trim(num2str(nTimeSteps,'(I10)'))//" time steps."
  call message(trim(message_text))
  message_text = trim(num2str(datetime(3),'(I2.2)'))//"."//trim(num2str(datetime(2),'(I2.2)')) &
       //"."//trim(num2str(datetime(1),'(I4.4)'))//" "//trim(num2str(datetime(5),'(I2.2)')) &
       //":"//trim(num2str(datetime(6),'(I2.2)'))//":"//trim(num2str(datetime(7),'(I2.2)'))
  call message('Finished at ', trim(message_text), '.')
  call message()
  call finish('mHM','Finished!')

END PROGRAM mhm_driver
