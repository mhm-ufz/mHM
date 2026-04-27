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
  use mo_os, only: path_abspath, path_join, check_path_isdir
  use mo_message, only: use_log
  use mo_mhm_cli, only: set_verbosity_level
  use mo_string_utils, only: n2s => num2str
  use mo_domain, only: domains, selected_domains, domain_t
  use mo_kind, only: i4
  use mo_exchange_type, only: standard_path
  use nml_config_project, only: nml_config_project_t, NML_OK
  use mo_string_utils, only: n2s => num2str
  !$ use omp_lib, only: omp_get_num_threads
  !$ integer(i4) :: n_threads
  logical :: openmp_enabled = .false.
  integer :: unit
  logical :: from_dirs
  integer(i4) :: n_domains, i, id
  character(len=*), parameter :: separator = repeat("-", 72)
  type(domain_t), pointer :: domain
  ! global configs
  type(nml_config_project_t) :: project
  ! command line interface parser
  type(cli_parser) :: parser

  character(:), allocatable :: cwd, domain_dir, meta_file, main_file, para_file, out_file
  character(1024) :: errmsg
  integer :: status

  parser = cli_parser(prog="mhm", description="The mesoscale hydrological model - mHM v6", &
    add_logger_options=.true., add_version_option=.true., version="6.0")

  call parser%add_option(name='log-file', has_value=.true., required=.false., &
    value_name="path", help='Write to log file.')

  call parser%add_option(name="nml", s_name="n", has_value=.true., &
    value_name="path", default="mhm.nml", help="The mHM configuration namelist.")

  call parser%add_option(name="parameter", s_name="p", has_value=.true., &
    value_name="path", default="parameter.nml", help="The mHM parameter namelist.")

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

  ! global configs
  meta_file = standard_path(cwd=cwd, file=parser%option_value("nml"))
  para_file = standard_path(cwd=cwd, file=parser%option_value("parameter"))
  out_file  = standard_path(cwd=cwd, file=parser%option_value("output"))
  log_info(*) "READ MAIN CONFIG: ", meta_file
  status = project%from_file(file=meta_file, errmsg=errmsg)
  if (status /= NML_OK) then
    log_fatal(*) "Error reading config_project: ", trim(errmsg)
    error stop 1
  end if
  status = project%is_valid(errmsg=errmsg)
  if (status /= NML_OK) then
    log_fatal(*) "Invalid config_project: ", trim(errmsg)
    error stop 1
  end if

  ! determine number of domains
  n_domains = project%n_domains
  from_dirs = project%read_domains_from_dirs
  allocate(selected_domains(n_domains))

  log_info(*) "CREATE DOMAINS: ", n_domains
  if (from_dirs) then
    log_info(*) "Reading domains from separate directories."
  end if

  ! create domain-list
  ! we use a linked list to be able to dynamically add domains
  ! These domains are stored as allocated pointers, which have implicitly the "target" attribute
  ! so components can safely point to "exchange" of their domain
  do i = 1_i4, n_domains
    ! returns the list key (domain id) which can be used to get the domain later
    ! prepared for MPI runs with multiple domains
    selected_domains(i) = domains%add_domain() ! returns domain id
  end do

  log_debug(*) "Selected domains", selected_domains

  ! read configs
  domain_dir = cwd
  main_file = meta_file
  do i = 1_i4, size(selected_domains)
    id = selected_domains(i)
    log_text(*) separator
    log_info(*) "CONFIGURE DOMAIN: ", id
    if (from_dirs) then
      status = project%is_set("domain_dirs", idx=[id], errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Directory not specified for domain ", n2s(id), ": ", trim(errmsg)
        error stop 1
      end if
      domain_dir = standard_path(cwd=cwd, path=project%domain_dirs(id))
      main_file = standard_path(cwd=domain_dir, file=project%domain_nmls(id))
    end if
    ! get domain
    call domains%get_domain(id, domain)
    ! id either from list or 1 if from dirs (always take domain 1 in each sub-dir)
    if (from_dirs) id = 1_i4
    ! create new domain and its exchange
    call domain%init(meta_file, main_file, para_file, id, domain_dir)
    ! configure domain components
    log_text(*) separator
    call domain%configure(main_file, out_file)
    ! check for connections and dependencies
    log_text(*) separator
    call domain%connect()
  end do

  ! simple run
  do i = 1_i4, size(selected_domains)
    id = selected_domains(i)
    call domains%get_domain(id, domain)
    log_text(*) separator
    log_info(*) "RUN DOMAIN: ", id
    call domain%initialize()
    log_text(*) separator
    log_info(*) "RUN TIME LOOP"
    do while(domain%exchange%time < domain%exchange%end_time)
      call domain%update()
    end do
    log_text(*) separator
    call domain%finalize()
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
