!>       \file mhm_driver.f90

!>       \brief Distributed precipitation-runoff model mHM

!>       \details This is the main driver of mHM, which calls
!>       one instance of mHM for a multiple domains and a given period.
!>       \image html  mhm5-logo.png "Typical mHM cell"
!>       \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>       \authors Luis Samaniego & Rohini Kumar (UFZ)

!>       \date Jun 2018

!>       \version \htmlinclude version.txt \latexinclude version.txt

!>       \copyright (c) \f$2005 - \the\year{}\f$, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!!       All rights reserved.
!!
!!       This code is a property of:
!!
!!       ----------------------------------------------------------
!!
!!       Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ
!!       Registered Office: Leipzig
!!       Registration Office: Amtsgericht Leipzig
!!       Trade Register: Nr. B 4703
!!       Chairman of the Supervisory Board: MinDirig Wilfried Kraus
!!       Scientific Director: Prof. Dr. Georg Teutsch
!!       Administrative Director: Dr. Heike Grassmann
!!
!!       ----------------------------------------------------------
!!
!!       NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!!       EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!!       OF THIS SOFTWARE. If software is modified to produce
!!       derivative works, such modified software should be
!!       clearly marked, so as not to confuse it with the version
!!       available from UFZ.  This code can be used for research
!!       purposes ONLY provided that the following sources are
!!       acknowledged:
!!
!!       Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!!       parameter regionalization of a grid-based hydrologic
!!       model at the mesoscale.  Water Resour. Res., 46,
!!       W05523, doi:10.1029/2008WR007327.
!!
!!       Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!!       of distributed hydrologic model parameterization on water
!!       fluxes at multiple scales and locations, Water Resour. Res.,
!!       49, doi:10.1029/2012WR012195.
!!
!!       For commercial applications you have to consult the
!!       authorities of the UFZ.

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
! Oldrich Rakovec, Rohini Kumar Oct 2015 - added reading of domain averaged TWS and objective function 15
!                                          for simultaneous calibration based on runoff and TWS
! Rohini Kumar                  Mar 2016 - options to handle different soil databases modified MPR to included
!                                          soil horizon specific properties/parameters
! Stephan Thober                Nov 2016 - implemented adaptive timestep for routing
! Rohini Kumar                  Dec 2016 - options to read (monthly mean) LAI fields
! Robert Schweppe               Jun 2018 - refactoring and reformatting
! Maren Kaluza                  Oct 2019 - TWS to data structure
! M.C. Demirel, Simon Stisen    Jun 2020 - New Soil Moisture Process: Feddes and FC dependency on root fraction coefficient processCase(3) = 4
PROGRAM mhm_driver

  use mo_optimization_types, only: &
          optidata ! type for opti data
  use mo_common_variables, only: &
#ifdef MPI
          comm, &
#endif
          optimize                                   ! optimization on/off and optimization method
  use mo_kind, only: i4                         ! number precision
  use mo_meteo_forcings, only: prepare_meteo_forcings_data
  use mo_mhm_eval, only: mhm_eval
  use mo_read_optional_data, only: readOptidataObs ! read optional observed data
  use mo_common_read_config, only: common_read_config, &       ! Read main configuration files
                                    check_optimization_settings ! Read main configuration files
  use mo_mhm_read_config, only: mhm_read_config                    ! Read main configuration files
  use mo_restart, only: write_restart_files
  use mo_startup, only: mhm_initialize
  use mo_string_utils, only: num2str             ! String magic
  use mo_timer, only: &
          timers_init, timer_start, timer_stop, timer_get              ! Timing of processes
  use mo_write_ascii, only: &
          write_configfile, &      ! Writing Configuration file
          write_optifile, &      ! Writing optimized parameter set and objective
          write_optinamelist     ! Writing optimized parameter set to a namelist
  use mo_objective_function, only: &
#ifdef MPI
          objective_subprocess, &
          objective_master, &
#endif
          objective                 ! objective functions and likelihoods
  use mo_optimization, only: optimization
  use mo_mrm_objective_function_runoff, only: &
#ifdef MPI
          single_objective_runoff_master, &
          single_objective_runoff_subprocess, &
#endif
          single_objective_runoff
  use mo_mrm_init, only: mrm_init, mrm_configuration
  use mo_mrm_write, only : mrm_write

#ifdef MPI
  use mpi_f08
#endif

  use mo_mhm_cli, only: parse_command_line
  use mo_mhm_messages, only: startup_message, domain_dir_check_message, finish_message
  use mo_mhm_interface, only: &
    mhm_interface_init, &
    mhm_interface_run, &
    mhm_interface_run_optimization, &
    mhm_interface_finalize

  IMPLICIT NONE

  procedure(mhm_eval), pointer :: eval
  procedure(objective), pointer :: obj_func

  ! MPI variables
  integer             :: ierror
  integer(i4)         :: nproc, rank, oldrank

#ifdef MPI
  ! Initialize MPI
  call MPI_Init(ierror)
  call MPI_Comm_dup(MPI_COMM_WORLD, comm, ierror)
  ! find number of processes nproc
  call MPI_Comm_size(comm, nproc, ierror)
  ! find the number the process is referred to, called rank
  call MPI_Comm_rank(comm, rank, ierror)
  oldrank = rank
  write(*,*) 'MPI!, comm', rank, nproc
#endif

  ! parse command line arguments
  call parse_command_line()

  ! initialize mhm
  call mhm_interface_init()

  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  if (optimize) then
    call mhm_interface_run_optimization()
  else
    ! single mhm run with current settings
    call mhm_interface_run()
  end if

  ! WRITE RESTART files and RUNOFF
  call mhm_interface_finalize()

  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call finish_message()

#ifdef MPI
  ! find number of processes nproc
  call MPI_Comm_size(comm, nproc, ierror)
  call MPI_Comm_rank(comm, rank, ierror)
  write(*,*) 'MPI finished', rank, nproc
  call MPI_Finalize(ierror)
#endif

END PROGRAM mhm_driver
