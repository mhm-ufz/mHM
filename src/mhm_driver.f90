!> \file mhm_driver.f90
! --------------------------------------------------------------------------
!> \authors   Luis Samaniego & Rohini Kumar (UFZ)
!  CONTACT    luis.samaniego@ufz.de / rohini.kumar@ufz.de
!> \version   5.0
!> \date      Dec 2012

!  PURPOSE
!>            \brief Distributed precipitation-runoff model mHM

!>            \details This is the main driver of mHM, which calls
!>             one instance of mHM for a single basin and a given period.

!>              \image html  mhm5-logo.png "Typical mHM cell"
!>              \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>  \copyright (c)2005-2014, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
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
!                Stephan Thober Nov 2013         added read in of latitude longitude fields
!                Matthias Zink  Mar 2013         edited screen output for gauges
!                                                added inflow gauges
!               Rohini Kumar, Apr 2014         - implementation of the mHM run on a single cell
!                                                configuration that too in the routing mode. 
!                                              - run mHM at the input data level i.e. L0 grid
!               Rohini Kumar, May 2014         - model run on a regular lat-lon grid or 
!                                                on a regular X-Y coordinate system
!
! --------------------------------------------------------------------------

PROGRAM mhm_driver 

  USE mo_anneal,              ONLY : anneal                         ! Optimise with Simulated Annealing SA
  USE mo_dds,                 ONLY : dds                            ! Optimise with Dynam. Dimens. Search DDS
  USE mo_file,                ONLY :                         &
       version, file_main,                                   &      ! main info
       file_namelist,                                        &      ! filename of namelist: main setup
       file_namelist_param,                                  &      ! filename of namelist: mhm model parameter
       file_defOutput                                               ! filename of namelist: output setup
  USE mo_finish,              ONLY : finish                         ! Finish with style
  USE mo_global_variables,    ONLY :                         &
       nbasins,                                              &      ! number of basins
       restart_flag_states_write, restart_flag_config_write, &      ! restart writing flags
       optimize, opti_method,                                &      ! optimization on/off and optimization method
       global_parameters, global_parameters_name,            &      ! mhm parameters (gamma) and their clear names
       dirRestartOut,                                        &      ! directories
       dirMorpho, dirLCover, dirGauges, dirPrecipitation,    &      ! directories
       dirTemperature, dirReferenceET, dirOut,               &      ! directories
       simPer,                                               &      ! simulation period
       NTSTEPDAY,                                            &      ! number of timesteps per day (former: NAGG)
       nIterations, seed,                                    &      ! settings for optimization algorithms
       dds_r, sa_temp, sce_ngs, sce_npg, sce_nps,            &      ! settings for optimization algorithms
       iFlag_LAI_data_format,                                &      ! LAI option for reading gridded LAI field 
       basin, processMatrix                                         ! basin information,  processMatrix
  USE mo_kind,                ONLY : i4, i8, dp                     ! number precision
  USE mo_mcmc,                ONLY : mcmc                           ! Monte Carlo Markov Chain method
  USE mo_message,             ONLY : message, message_text          ! For print out
  USE mo_meteo_forcings,      ONLY : prepare_meteo_forcings_data 
  USE mo_mhm_eval,            ONLY : mhm_eval
  USE mo_objective_function,  ONLY : objective, loglikelihood       ! objective functions and likelihoods
  USE mo_prepare_gridded_LAI, ONLY : prepare_gridded_daily_LAI_data ! prepare daily LAI gridded fields
  USE mo_read_config,         ONLY : read_config                    ! Read main configuration files
  USE mo_read_wrapper,        ONLY : read_data                      ! Read all input data
  USE mo_read_latlon,         ONLY : read_latlon
  USE mo_restart,             ONLY :                         &
       write_restart_states,                                 &      ! Writing states for restart
       write_restart_config                                         ! Writing configuration for restart
  USE mo_sce,                 ONLY : sce                            ! Optimize with Shuffled Complex Evolution SCE
  USE mo_startup,             ONLY : initialise
  USE mo_string_utils,        ONLY : num2str, separator             ! String magic
  USE mo_timer,               ONLY :                         &
       timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  USE mo_write_ascii,         ONLY :                         &
       write_configfile,                                     &      ! Writing Configuration file
       write_optifile,                                       &      ! Writing optimized parameter set and objective
       write_optinamelist                                           ! Writing optimized parameter set to a namelist
  !$ USE omp_lib,             ONLY : OMP_GET_NUM_THREADS           ! OpenMP routines

  IMPLICIT NONE

  ! local
  integer, dimension(8)                 :: datetime        ! Date and time
  !$ integer(i4)                        :: n_threads       ! OpenMP number of parallel threads
  integer(i4)                           :: ii, jj           ! Counters
  integer(i4)                           :: iTimer          ! Current timer number
  integer(i4)                           :: nTimeSteps
  real(dp)                              :: funcbest        ! best objective function achivied during optimization 
  ! model output
  real(dp), allocatable, dimension(:,:) :: riverrun        ! simulated river runoff at all gauges, timepoints
  ! mcmc
  real(dp), dimension(:,:), allocatable :: burnin_paras    ! parameter sets sampled during burnin
  real(dp), dimension(:,:), allocatable :: mcmc_paras      ! parameter sets sampled during proper mcmc
  logical,  dimension(:),   allocatable :: maskpara        ! true  = parameter will be optimized     = parameter(i,4) = 1
  !                                                        ! false = parameter will not be optimized = parameter(i,4) = 0

  ! --------------------------------------------------------------------------
  ! START
  ! --------------------------------------------------------------------------
  call message(separator)
  call message('              mHM-UFZ')
  call message()
  call message('    MULTISCALE HYDROLOGIC MODEL')
  call message('           Revision ', trim(version))
  call message('Originally by L. Samaniego & R. Kumar')
  call message('          February 2013')
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
  call message('  Input data directories:')
  call message('    # of basins  :          ', trim(num2str(nbasins)))
  do ii = 1, nbasins
     call message( '  --------------' )
     call message( '      BASIN ', num2str(ii,'(I3)') )
     call message( '  --------------' )
     call message('    Morphological directory:  ', trim(dirMorpho(ii) ))         
     call message('    Land cover directory:     ', trim(dirLCover(ii) ))          
     call message('    Discharge directory:      ', trim(dirGauges(ii)  ))          
     call message('    Precipitation directory:  ', trim(dirPrecipitation(ii)  ))
     call message('    Temperature directory:    ', trim(dirTemperature(ii)  ))   
     call message('    PET directory:            ', trim(dirReferenceET(ii)  )) 
     call message('    Output directory:         ', trim(dirOut(ii) ))        
     if (processMatrix(8,1) .GT. 0) then
        call message('    Evaluation gauge          ', 'ID')
        do jj = 1 , basin%nGauges(ii)
           call message('    ',trim(adjustl(num2str(jj))),'                         ', &
                trim(adjustl(num2str(basin%gaugeIdList(ii,jj)))))
        end do
     end if
     if (basin%nInflowGauges(ii) .GT. 0) then
        call message('    Inflow gauge              ', 'ID')
        do jj = 1 , basin%nInflowGauges(ii)
           call message('    ',trim(adjustl(num2str(jj))),'                         ', &
                trim(adjustl(num2str(basin%InflowGaugeIdList(ii,jj)))))
        end do
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

  !if (.not. restart_flag_read ) then
     call message('  Read data ' )
     call timer_start(itimer)
     ! for DEM, slope, ... define nGvar local
     ! read_data has a basin loop inside
     call read_data()
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  !end if

  ! read data for every basin
  do ii=1, nbasins
     itimer = itimer + 1
     call message('  Initialise basin: ', trim(adjustl(num2str(ii))),' ...')
     call timer_start(itimer)
     call initialise(ii)
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
     
     ! meteorological forcings (reading, upscaling or downscaling)
     call prepare_meteo_forcings_data(ii)

     ! read lat lon coordinates of each basin
     call message('  Reading lat-lon for basin: ', trim(adjustl(num2str(ii))),' ...')
     call timer_start(itimer)
     call read_latlon(ii)
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
     
     ! daily gridded LAI values
     if(iFlag_LAI_data_format .EQ. 1) call prepare_gridded_daily_LAI_data(ii)
  
  end do

  ! The following block is for testing of the restart <<<<<<<<<<<<<<<<<<<<<<<<<<
  ! if ( restart_flag_config_write ) then
  !    itimer = itimer + 1
  !    call message()
  !    call message( '  Write restart config file')
  !    call timer_start(itimer)
  !    call write_restart_config( dirRestartOut )
  !    call timer_stop(itimer)
  !    call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  ! end if
  ! ! write restart for states
  ! if ( restart_flag_states_write ) then
  !    itimer = itimer + 1
  !    call message()
  !    call message( '  Write restart states file')
  !    call timer_start(itimer)
  !    call write_restart_states( dirRestartOut )
  !    call timer_stop(itimer)
  !    call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  ! end if
  ! stop 'Test restart' ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  !this call may be moved to another position as it writes the master config out file for all basins
  call write_configfile()

  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  iTimer = iTimer + 1
  call message()
  if ( optimize ) then
     call message('  Start optimization')
     call timer_start(iTimer)

     ! mask parameter which have a FLAG=0 in mhm_parameter.nml
     ! maskpara = true : parameter will be optimized
     ! maskpara = false : parameter is discarded during optimization
     allocate(maskpara(size(global_parameters,1)))
     maskpara = .true.
     do ii=1, size(global_parameters,1)
        if ( nint(global_parameters(ii,4),i4) .eq. 0_i4 ) then
           maskpara(ii) = .false.
        end if
     end do

     select case (opti_method)
     case (0)
        call message('    Use MCMC')

        if (seed .gt. 0_i8) then
           ! use fixed user-defined seed 
           call mcmc(loglikelihood, global_parameters(:,3), global_parameters(:,1:2), mcmc_paras, burnin_paras, &
                ParaSelectMode_in=2_i4,tmp_file='mcmc_tmp_parasets.nc',                                         &
                maskpara_in=maskpara,                                                                           &
                seed_in=seed, loglike_in=.true., printflag_in=.true.)
        else
           ! use flexible clock-time seed
           call mcmc(loglikelihood, global_parameters(:,3), global_parameters(:,1:2), mcmc_paras, burnin_paras, &
                ParaSelectMode_in=2_i4,tmp_file='mcmc_tmp_parasets.nc',                                         &
                maskpara_in=maskpara,                                                                           &
                loglike_in=.true., printflag_in=.true.)
        end if
     case (1)
        call message('    Use DDS')
        if (seed .gt. 0_i8) then
           ! use fixed user-defined seed
           global_parameters(:,3) = dds(objective, global_parameters(:,3), global_parameters(:,1:2),            &
                maxiter=int(nIterations,i8), r=dds_r, seed=seed,                                                &
                tmp_file='dds_results.out', mask=maskpara,                                                      &
                funcbest=funcbest)
        else
           ! use flexible clock-time seed
           global_parameters(:,3) = dds(objective, global_parameters(:,3), global_parameters(:,1:2),            &
                maxiter=int(nIterations,i8), r=dds_r,                                                           &
                tmp_file='dds_results.out', mask=maskpara,                                                      &
                funcbest=funcbest)
        end if
     case (2)
        call message('    Use Simulated Annealing')
        if (seed .gt. 0_i8 .and. sa_temp .gt. 0.0_dp) then
           ! use fixed user-defined seed and user-defined initial temperature
           global_parameters(:,3) = anneal(objective, global_parameters(:,3), global_parameters(:,1:2),         &
                temp=sa_temp, seeds=(/seed, seed+1000_i8, seed+2000_i8/), nITERmax=nIterations,                 &
                tmp_file='anneal_results.out', maskpara=maskpara,                                               &
                funcbest=funcbest)
        else if (seed .gt. 0_i8) then
           ! use fixed user-defined seed and adaptive initial temperature
           global_parameters(:,3) = anneal(objective, global_parameters(:,3), global_parameters(:,1:2),         &
                seeds=(/seed, seed+1000_i8, seed+2000_i8/), nITERmax=nIterations,                               &
                tmp_file='anneal_results.out', maskpara=maskpara,                                               &
                funcbest=funcbest)
        else if (sa_temp .gt. 0.0_dp) then
           ! use flexible clock-time seed and user-defined initial temperature
           global_parameters(:,3) = anneal(objective, global_parameters(:,3), global_parameters(:,1:2),         &
                temp=sa_temp, nITERmax=nIterations,                                                             &
                tmp_file='anneal_results.out', maskpara=maskpara,                                               &
                funcbest=funcbest)
        else
           ! use flexible clock-time seed and adaptive initial temperature
           global_parameters(:,3) = anneal(objective, global_parameters(:,3), global_parameters(:,1:2),         &
                nITERmax=nIterations,                                                                           &
                tmp_file='anneal_results.out', maskpara=maskpara,                                               &
                funcbest=funcbest)
        end if
     case (3)
        call message('    Use SCE')
        if (seed .gt. 0_i8) then
           ! use fixed user-defined seed
           global_parameters(:,3) = sce(objective, global_parameters(:,3), global_parameters(:,1:2),            &
                mymaxn=int(nIterations,i8), myseed=seed, myngs=sce_ngs, mynpg=sce_npg, mynps=sce_nps,           &
                parallel=.false., mymask=maskpara,                                                              &
                tmp_file='sce_results.out', popul_file='sce_population.out',                                    &
                bestf=funcbest)
        else
           ! use flexible clock-time seed
           global_parameters(:,3) = sce(objective, global_parameters(:,3), global_parameters(:,1:2),            &
                mymaxn=int(nIterations,i8), myngs=sce_ngs, mynpg=sce_npg, mynps=sce_nps,                        &
                parallel=.false., mymask=maskpara,                                                              &
                tmp_file='sce_results.out', popul_file='sce_population.out',                                    &
                bestf=funcbest)
        end if
     case default
        call finish('mHM','This optimization method is not implemented.')
     end select
     call timer_stop(iTimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
     
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
    if ( processMatrix(8,1) .eq. 0 ) then
       ! call mhm without routing
       call mhm_eval(global_parameters(:,3))
    else
       ! call mhm with routing
       call mhm_eval(global_parameters(:,3), runoff=riverrun)
    end if
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer),'(F12.3)')), ' seconds.')
    !
  end if

  ! --------------------------------------------------------------------------
  ! WRITE RESTART files
  ! --------------------------------------------------------------------------
  ! write restart for configuration
  if ( restart_flag_config_write ) then
     itimer = itimer + 1
     call message()
     call message( '  Write restart config file')
     call timer_start(itimer)
     call write_restart_config( dirRestartOut )
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  end if
  ! write restart for states
  if ( restart_flag_states_write ) then
     itimer = itimer + 1
     call message()
     call message( '  Write restart states file')
     call timer_start(itimer)
     call write_restart_states( dirRestartOut )
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  end if
  
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

  nTimeSteps = ( simPer%julEnd - simPer%julStart + 1 ) * NTSTEPDAY
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
