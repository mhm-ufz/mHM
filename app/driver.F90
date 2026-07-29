!> \file    app/driver.F90
!> \copydoc driver

!> \brief mHM v6 driver.
!> \details This is the main driver program for mHM v6. It initializes the
!>          model, reads configurations, sets up domains, and runs the time loop.
!> \authors Sebastian Mueller
!> \date    Jan 2026
#include "logging.h"
program driver
  use mo_logging
  use mo_cli, only: cli_parser
  use mo_os, only: path_abspath, check_path_isdir
  use mo_message, only: use_log
  use mo_mhm_cli, only: set_verbosity_level
  use mo_domain, only: domain_t
  use mo_kind, only: i4
  use mo_exchange_type, only: get_n_domains
  !$ use omp_lib, only: omp_get_num_threads
  !$ integer(i4) :: n_threads
  logical :: openmp_enabled = .false.
  integer :: unit
  integer(i4) :: n_domains, i
  character(len=*), parameter :: separator = repeat("-", 72)
  type(domain_t), allocatable, target :: domains(:)
  ! command line interface parser
  type(cli_parser) :: parser

  character(:), allocatable :: cwd, main_file, para_file, out_file

  parser = cli_parser(prog="mhm", description="The mesoscale hydrological model - mHM v6", &
    add_logger_options=.true., add_version_option=.true., version="6.0")

  call parser%add_option(name='log-file', has_value=.true., required=.false., &
    value_name="path", help='Write to log file.')

  call parser%add_option(name="nml", s_name="n", has_value=.true., &
    value_name="path", default="mhm.nml", help="The mHM configuration namelist.")

  call parser%add_option(name="parameters", s_name="p", has_value=.true., &
    value_name="path", default="parameters.nml", help="The mHM parameter namelist.")

  call parser%add_option(name="output", s_name="o", has_value=.true., &
    value_name="path", default="outputs.nml", help="The mHM output namelist.")

  call parser%add_option(name="cwd", blank=.true., help="The desired working directory (optional).")

  ! parse given command line arguments
  call parser%parse()

  use_log = .true. ! enable logging in message module (after parsing for clean help-messages)

  if (parser%option_was_read('log-file')) then
    open(newunit=unit, file=parser%option_value('log-file'), status='replace', action='write')
    ! redirect log units to the opened file
    log_unit = unit
    log_unit_error = unit
  end if

  call set_verbosity_level(3_i4 - parser%option_read_count("quiet")) ! TODO: when switching to logging, not needed anymore

  call startup_message()

  !$omp parallel
  !$ n_threads = omp_get_num_threads()
  !$omp end parallel
  !$ openmp_enabled = .true.
  !$ log_info(*) "OpenMP enabled with ", n_threads, " threads."

  if (.not. openmp_enabled) then
    log_info(*) "OpenMP not enabled."
  end if

  ! get current working directory
  ! we don't change the process working directory, but use cwd for all relative paths
  ! we store the CWD in the exchange type later and use the "get_path" method to get paths
  cwd = "."
  if (parser%option_was_read("cwd")) cwd = parser%option_value("cwd")
  cwd = path_abspath(cwd)
  call check_path_isdir(cwd, raise=.true.)

  ! global config file names; components resolve them against the run root
  main_file = parser%option_value("nml")
  para_file = parser%option_value("parameters")
  out_file  = parser%option_value("output")

  ! get number of domains
  n_domains = get_n_domains(main_file=main_file, cwd=cwd)
  log_info(*) "ALLOCATE DOMAINS: ", n_domains
  allocate(domains(n_domains))

  ! read configs
  do i = 1_i4, size(domains)
    log_text(*) separator
    log_info(*) "DOMAIN: ", i
    ! create new domain and its exchange
    call domains(i)%create(main_file=main_file, domain=i, cwd=cwd)
    ! configure domain components
    log_text(*) separator
    call domains(i)%configure(main_file=main_file, para_file=para_file, out_file=out_file)
    ! check for connections and dependencies
    log_text(*) separator
    call domains(i)%connect()
  end do

  ! simple run
  do i = 1_i4, size(domains)
    log_text(*) separator
    log_info(*) "RUN DOMAIN: ", i
    call domains(i)%initialize()
    log_text(*) separator
    log_info(*) "RUN TIME LOOP"
    do while(domains(i)%exchange%time < domains(i)%exchange%end_time)
      call domains(i)%update()
    end do
    log_text(*) separator
    call domains(i)%finalize()
  end do

contains
  subroutine startup_message()
    use mo_file, only: version, version_date
    log_text(*) ".----------------------------------------------------------------------."
    log_text(*) "|                                                                      |"
    log_text(*) "|                  /MM   /MM /MM      /MM                   /MMMMMM    |"
    log_text(*) "|                 | MM  | MM| MMM    /MMM                  /MM__  MM   |"
    log_text(*) "|    /MMMMMM/MMMM | MM  | MM| MMMM  /MMMM       /MM    /MM| MM  \__/   |"
    log_text(*) "|   | MM_  MM_  MM| MMMMMMMM| MM MM/MM MM      |  MM  /MM/| MMMMMMM    |"
    log_text(*) "|   | MM \ MM \ MM| MM__  MM| MM  MMM| MM       \  MM/MM/ | MM__  MM   |"
    log_text(*) "|   | MM | MM | MM| MM  | MM| MM\  M | MM        \  MMM/  | MM  \ MM   |"
    log_text(*) "|   | MM | MM | MM| MM  | MM| MM \/  | MM         \  M/   |  MMMMMM/   |"
    log_text(*) "|   |__/ |__/ |__/|__/  |__/|__/     |__/          \_/     \______/    |"
    log_text(*) "|                                                                      |"
    log_text(*) "|                               mHM-UFZ                                |"
    log_text(*) "|                     MULTISCALE HYDROLOGIC MODEL                      |"
    log_text(*) "|                                                                      |"
    log_text(*) "|                      Version: " // trim(version)      // repeat(" ", 39-len_trim(version))      // "|"
    log_text(*) "|                      Date:    " // trim(version_date) // repeat(" ", 39-len_trim(version_date)) // "|"
    log_text(*) "|                                                                      |"
    log_text(*) "|                Originally by L. Samaniego & R. Kumar                 |"
    log_text(*) "'----------------------------------------------------------------------'"
  end subroutine startup_message
end program driver
